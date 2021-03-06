#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TMatrixDSymEigen.h"
#include "Fit/Fitter.h"
#include "TCanvas.h"

#include <vector>
#include <string>

#include "command_line_tools.h"
#include "datasets.h"
#include "stat.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

int n_ff_parameters;

const unsigned int n_dsdt_samples = 200;

// TODO: rename
const unsigned int n_repetitions = 1000;

string fitModel = "";

//----------------------------------------------------------------------------------------------------

struct Constraint
{
	double sqrt_s;
	TVectorD par;
	TMatrixDSym par_V;
	TMatrixDSym par_V_inv;

	Constraint() {}

	Constraint(double _s, const TVectorD &_p, const TMatrixDSym &_V) : sqrt_s(_s), par(_p), par_V(_V)
	{
		par_V_inv.ResizeTo(par_V);
		par_V_inv = par_V;
		par_V_inv.Invert();
	}
};

vector<Constraint> constraints;

//----------------------------------------------------------------------------------------------------

class S2_FCN
{
	public:
		S2_FCN() {}

  		double operator() (const double *par) const;

		string extrapolationModel = "";
		bool useCorrelationMatrix = false;
};

//----------------------------------------------------------------------------------------------------

TVectorD EvaluateParameters(const TVectorD &a, const TVectorD &b, double sqrt_s, const string &extrapolationModel)
{
	double f = 0.;
	if (extrapolationModel == "sqrt_s") f = sqrt_s;
	if (extrapolationModel == "log_sqrt_s") f = log(sqrt_s);

	return a*f + b;
}

//----------------------------------------------------------------------------------------------------

double S2_FCN::operator() (const double *par) const
{
	double S2 = 0.;

	TVectorD a(n_ff_parameters), b(n_ff_parameters);
	for (int i = 0; i < n_ff_parameters; ++i)
	{
		a(i) = par[i];
		b(i) = par[n_ff_parameters + i];
	}

	for (const auto &cnt : constraints)
	{
		TVectorD p = EvaluateParameters(a, b, cnt.sqrt_s, extrapolationModel);

		TVectorD diff = p - cnt.par;

		for (int i = 0; i < n_ff_parameters; i++)
		{
			for (int j = 0; j < n_ff_parameters; j++)
			{
				if (useCorrelationMatrix == false && (i != j))
					continue;

				S2 += diff(i) * cnt.par_V_inv(i, j) * diff(j);
			}
		}
	}

	return S2;
}

//----------------------------------------------------------------------------------------------------

void SaveFitPlots(const ROOT::Fit::FitResult &result, const string &extrapolationModel)
{
	// calculate extrapolated parameters
	const double *par = result.GetParams();

	TVectorD a(n_ff_parameters), b(n_ff_parameters);
	for (int i = 0; i < n_ff_parameters; ++i)
	{
		a(i) = par[i];
		b(i) = par[n_ff_parameters + i];
	}

	// prepare graph vectors
	vector<TGraphErrors *> vg_cnt, vg_cnt_orig;
	vector<TGraph *> vg_fit, vg_fit_pl_unc, vg_fit_mi_unc;
	for (int i = 0; i < n_ff_parameters; ++i)
	{
		TGraphErrors *ge;
		ge = new TGraphErrors(); ge->SetName("g_cnt"); vg_cnt.push_back(ge);
		ge = new TGraphErrors(); ge->SetName("g_cnt_orig"); vg_cnt_orig.push_back(ge);

		TGraph *g;
		g = new TGraph(); g->SetName("g_fit"); vg_fit.push_back(g);
		g = new TGraph(); g->SetName("g_fit_pl_unc"); vg_fit_pl_unc.push_back(g);
		g = new TGraph(); g->SetName("g_fit_mi_unc"); vg_fit_mi_unc.push_back(g);
	}

	// make graph of data points (constraints)
	for (int i = 0; i < n_ff_parameters; ++i)
	{
		for (const auto &cnt : constraints)
		{
			int idx = vg_cnt[i]->GetN();
			vg_cnt[i]->SetPoint(idx, cnt.sqrt_s, cnt.par(i));
			vg_cnt[i]->SetPointError(idx, 0., sqrt(cnt.par_V(i, i)));
		}
	}

	// make graph of original (dip-bump) fit constraints
	for (const auto &ds : datasets)
	{
		for (const auto &c : ds.constraints)
		{
			TGraphErrors *g = vg_cnt_orig[c.first];
			int idx = g->GetN();
			g->SetPoint(idx, ds.sqrt_s, c.second.mean);
			g->SetPointError(idx, 0., c.second.sigma);
		}
	}

	// prepare pertubation generator matrix
	int dim = result.NPar();
	TMatrixDSym V(dim);
	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j)
			V(i, j) = result.CovMatrix(i, j);

	TMatrixDSymEigen eig_decomp(V);
	TVectorD eig_values(eig_decomp.GetEigenValues());
	TMatrixDSym S(V.GetNrows());
	for (int i = 0; i < V.GetNrows(); i++)
		S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
	TMatrixD m_gen = eig_decomp.GetEigenVectors() * S;

	// sample fit (central values)
	for (double W = 1.; W < 14.; W += 0.5)
	{
		// evaluate parameter central values
		TVectorD p = EvaluateParameters(a, b, W, extrapolationModel);

		// add point to graphs
		for (int i = 0; i < n_ff_parameters; ++i)
		{
			int idx = vg_fit[i]->GetN();
			vg_fit[i]->SetPoint(idx, W, p(i));
		}
	}

	// sample fit uncertainty
	vector<Stat> st(vg_fit[0]->GetN(), Stat(n_ff_parameters));

	for (unsigned int ir = 0; ir < n_repetitions; ++ir)
	{
		// generate biased extrapolation
		TVectorD delta(V.GetNrows()), rdm(V.GetNrows());
		for (int i = 0; i < V.GetNrows(); i++)
			rdm(i) = gRandom->Gaus();
		delta = m_gen * rdm;

		TVectorD a_mod(n_ff_parameters), b_mod(n_ff_parameters);
		for (int i = 0; i < n_ff_parameters; ++i)
		{
			a_mod(i) = par[i] + delta(i);
			b_mod(i) = par[n_ff_parameters + i] + delta(n_ff_parameters + i);
		}

		// sample biased parameters at different energy values
		for (int iw = 0; iw < vg_fit[0]->GetN(); ++iw)
		{
			const double W = vg_fit[0]->GetX()[iw];

			TVectorD p_mod = EvaluateParameters(a_mod, b_mod, W, extrapolationModel);

			st[iw].Fill(p_mod);
		}
	}

	// build plots
	for (int i = 0; i < n_ff_parameters; ++i)
	{
		for (int iw = 0; iw < vg_fit[i]->GetN(); ++iw)
		{
			const double W = vg_fit[i]->GetX()[iw];
			const double f = vg_fit[i]->GetY()[iw];

			vg_fit_pl_unc[i]->SetPoint(iw, W, f + st[iw].GetStdDev(i));
			vg_fit_mi_unc[i]->SetPoint(iw, W, f - st[iw].GetStdDev(i));
		}
	}

	// save plots
	for (int i = 0; i < n_ff_parameters; ++i)
	{
		char buf[100];
		sprintf(buf, "c_par%i", i);
		TCanvas *c = new TCanvas();
		c->SetName(buf);

		vg_fit_pl_unc[i]->SetLineStyle(2);
		vg_fit_mi_unc[i]->SetLineStyle(2);

		vg_cnt[i]->Draw("ap");
		vg_cnt_orig[i]->Draw("p");
		vg_fit[i]->Draw("l");
		vg_fit_pl_unc[i]->Draw("l");
		vg_fit_mi_unc[i]->Draw("l");

		c->Write();
	}
}

//----------------------------------------------------------------------------------------------------

void SaveParameterUncertaintyBand(const ROOT::Fit::FitResult &result, const string &extrapolationModel, const Dataset &ds,
		const double *par, const double ep, const TGraph *g_ref, TGraph *g_dsdt_ext_unc, TGraph *g_B_ext_unc)
{
	// prepare pertubation generator matrix
	int dim = result.NPar();
	TMatrixDSym V(dim);
	for (int i = 0; i < dim; ++i)
		for (int j = 0; j < dim; ++j)
			V(i, j) = result.CovMatrix(i, j);

	TMatrixDSymEigen eig_decomp(V);
	TVectorD eig_values(eig_decomp.GetEigenValues());
	TMatrixDSym S(V.GetNrows());
	for (int i = 0; i < V.GetNrows(); i++)
		S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
	TMatrixD m_gen = eig_decomp.GetEigenVectors() * S;

	// uncertainty band
	vector<Stat> stat(n_dsdt_samples, Stat(1)), stat_B(n_dsdt_samples, Stat(1));

	TH1D *h_dsdt_mod_at_t_min = new TH1D("", "dsdt_mod(t_min)", 100., 0., 0.);
	TH1D *h_dsdt_mod_at_t_dip = new TH1D("", "dsdt_mod(t_dip)", 100., 0., 0.);
	TH1D *h_dsdt_mod_at_t_bmp = new TH1D("", "dsdt_mod(t_bmp)", 100., 0., 0.);
	TH1D *h_dsdt_mod_at_t_max = new TH1D("", "dsdt_mod(t_max)", 100., 0., 0.);

	const double dsdt_at_t_min = g_ref->Eval(ds.t_min);
	const double dsdt_at_t_dip = g_ref->Eval(ds.t_dip);
	const double dsdt_at_t_bmp = g_ref->Eval(ds.t_bmp);
	const double dsdt_at_t_max = g_ref->Eval(ds.t_max);

	for (unsigned int ir = 0; ir < n_repetitions; ++ir)
	{
		// generate biased extrapolation
		TVectorD delta(V.GetNrows()), rdm(V.GetNrows());
		for (int i = 0; i < V.GetNrows(); i++)
			rdm(i) = gRandom->Gaus();
		delta = m_gen * rdm;

		TVectorD a_mod(n_ff_parameters), b_mod(n_ff_parameters);
		for (int i = 0; i < n_ff_parameters; ++i)
		{
			a_mod(i) = par[i] + delta(i);
			b_mod(i) = par[n_ff_parameters + i] + delta(n_ff_parameters + i);
		}

		TVectorD p_mod = EvaluateParameters(a_mod, b_mod, ds.sqrt_s, extrapolationModel);

		// set extrapolated parameters to fit function
		TF1 *ff_mod = new TF1(*ds.ff);

		for (int i = 0; i < n_ff_parameters; ++i)
			ff_mod->SetParameter(i, p_mod(i));

		// acceptable configuration ??
		const double dsdt_mod_at_t_min = ff_mod->Eval(ds.t_min);
		const double dsdt_mod_at_t_dip = ff_mod->Eval(ds.t_dip);
		const double dsdt_mod_at_t_bmp = ff_mod->Eval(ds.t_bmp);
		const double dsdt_mod_at_t_max = ff_mod->Eval(ds.t_max);

		bool accept = (
				fabs(dsdt_mod_at_t_min - dsdt_at_t_min) / dsdt_at_t_min < 2.0
				&& fabs(dsdt_mod_at_t_dip - dsdt_at_t_dip) / dsdt_at_t_dip < 2.0
				&& fabs(dsdt_mod_at_t_bmp - dsdt_at_t_bmp) / dsdt_at_t_bmp < 2.0
				&& fabs(dsdt_mod_at_t_max - dsdt_at_t_max) / dsdt_at_t_max < 2.0
			);

		if (!accept)
			continue;

		h_dsdt_mod_at_t_min->Fill(dsdt_mod_at_t_min);
		h_dsdt_mod_at_t_dip->Fill(dsdt_mod_at_t_dip);
		h_dsdt_mod_at_t_bmp->Fill(dsdt_mod_at_t_bmp);
		h_dsdt_mod_at_t_max->Fill(dsdt_mod_at_t_max);

		// sample biased extrapolation
		for (int i = 0; i < g_ref->GetN(); ++i)
		{
			const double t = g_ref->GetX()[i];
			const double f_dsdt_mod = ff_mod->Eval(t);
			const double f_B_mod = (log(ff_mod->Eval(t + ep)) - log(f_dsdt_mod)) / ep;

			stat[i].Fill(f_dsdt_mod);
			stat_B[i].Fill(f_B_mod);
		}

		// clean up
		delete ff_mod;
	}

	// write debug plots
	h_dsdt_mod_at_t_min->Write("h_dsdt_mod_at_t_min");
	h_dsdt_mod_at_t_dip->Write("h_dsdt_mod_at_t_dip");
	h_dsdt_mod_at_t_bmp->Write("h_dsdt_mod_at_t_bmp");
	h_dsdt_mod_at_t_max->Write("h_dsdt_mod_at_t_max");

	// build graphs
	for (int i = 0; i < g_ref->GetN(); ++i)
	{
		const double t = g_ref->GetX()[i];

		const double f_dsdt_unc = stat[i].GetStdDev(0);
		g_dsdt_ext_unc->SetPoint(i, t, f_dsdt_unc);

		const double f_B_unc = stat_B[i].GetStdDev(0);
		g_B_ext_unc->SetPoint(i, t, f_B_unc);
	}
}

//----------------------------------------------------------------------------------------------------

void SaveTValueUncertaintyBand(const ROOT::Fit::FitResult &result, const string &extrapolationModel, const Dataset &ds,
		const double *par, const double ep, const TGraph *g_ref, TGraph *g_dsdt_ext_unc, TGraph *g_B_ext_unc)
{
	// prepare parameter values
	TVectorD a(n_ff_parameters), b(n_ff_parameters);
	for (int i = 0; i < n_ff_parameters; ++i)
	{
		a(i) = par[i];
		b(i) = par[n_ff_parameters + i];
	}

	TVectorD p = EvaluateParameters(a, b, ds.sqrt_s, extrapolationModel);

	// prepare pertubation generator matrix
	TMatrixDSymEigen eig_decomp(ds.m_t_value_cov);
	TVectorD eig_values(eig_decomp.GetEigenValues());
	TMatrixDSym S(ds.m_t_value_cov.GetNrows());
	for (int i = 0; i < ds.m_t_value_cov.GetNrows(); i++)
		S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
	TMatrixD m_gen = eig_decomp.GetEigenVectors() * S;

	// uncertainty band
	vector<Stat> stat(n_dsdt_samples, Stat(1)), stat_B(n_dsdt_samples, Stat(1));

	TH1D *h_t_min = new TH1D("", ";t_{min}", 1000, 0.3, 1.3);
	TH1D *h_t_dip = new TH1D("", ";t_{dip}", 1000, 0.3, 1.3);
	TH1D *h_t_bmp = new TH1D("", ";t_{bmp}", 1000, 0.3, 1.3);
	TH1D *h_t_max = new TH1D("", ";t_{max}", 1000, 0.3, 1.3);

	for (unsigned int ir = 0; ir < n_repetitions; ++ir)
	{
		// generate biased dataset
		Dataset ds_mod(ds);

		// generate bias in t-values
		TVectorD delta(ds_mod.m_t_value_cov.GetNrows()), rdm(ds_mod.m_t_value_cov.GetNrows());
		for (int i = 0; i < rdm.GetNrows(); i++)
			rdm(i) = gRandom->Gaus();
		delta = m_gen * rdm;

		ds_mod.t_min += delta(0);
		ds_mod.t_dip += delta(1);
		ds_mod.t_bmp += delta(2);
		ds_mod.t_max += delta(3);

		h_t_min->Fill(ds_mod.t_min);
		h_t_dip->Fill(ds_mod.t_dip);
		h_t_bmp->Fill(ds_mod.t_bmp);
		h_t_max->Fill(ds_mod.t_max);

		// (re)-build fit function
		BuildFitFunction(fitModel, ds_mod);

		// set extrapolated parameters to fit function
		for (int i = 0; i < n_ff_parameters; ++i)
			ds_mod.ff->SetParameter(i, p(i));

		// sample biased extrapolation
		for (int i = 0; i < g_ref->GetN(); ++i)
		{
			const double t = g_ref->GetX()[i];
			const double f_dsdt_mod = ds_mod.ff->Eval(t);
			const double f_B_mod = (log(ds_mod.ff->Eval(t + ep)) - log(f_dsdt_mod)) / ep;

			stat[i].Fill(f_dsdt_mod);
			stat_B[i].Fill(f_B_mod);
		}

		// clean up
		delete ds_mod.ff;
	}

	h_t_min->Write("h_t_min");
	h_t_dip->Write("h_t_dip");
	h_t_bmp->Write("h_t_bmp");
	h_t_max->Write("h_t_max");

	// build graphs
	for (int i = 0; i < g_ref->GetN(); ++i)
	{
		const double t = g_ref->GetX()[i];

		const double f_dsdt_unc = stat[i].GetStdDev(0);
		g_dsdt_ext_unc->SetPoint(i, t, f_dsdt_unc);

		const double f_B_unc = stat_B[i].GetStdDev(0);
		g_B_ext_unc->SetPoint(i, t, f_B_unc);
	}
}

//----------------------------------------------------------------------------------------------------

void SaveExtrapolation(const ROOT::Fit::FitResult &result, const string &extrapolationModel, const Dataset &ds, bool calculateTUnc = false)
{
	// calculate extrapolated parameters
	const double *par = result.GetParams();

	TVectorD a(n_ff_parameters), b(n_ff_parameters);
	for (int i = 0; i < n_ff_parameters; ++i)
	{
		a(i) = par[i];
		b(i) = par[n_ff_parameters + i];
	}

	TVectorD p = EvaluateParameters(a, b, ds.sqrt_s, extrapolationModel);

	// set extrapolated parameters to fit function
	for (int i = 0; i < n_ff_parameters; ++i)
		ds.ff->SetParameter(i, p(i));

	// evaluate extrapolated graph
	map<string, TGraph *> m_g_ext;
	m_g_ext["dsdt"] = new TGraph(); m_g_ext["dsdt"]->SetName("g_dsdt_ext");
	m_g_ext["B"] = new TGraph(); m_g_ext["B"]->SetName("g_B_ext");

	TGraph *g_dsdt_ext = new TGraph(); g_dsdt_ext->SetName("g_dsdt_ext");
	TGraph *g_B_ext = new TGraph(); g_B_ext->SetName("g_B_ext");

	const double ep = 1E-4;

	for (unsigned int i = 0; i < n_dsdt_samples; ++i)
	{
		const double t_min_eff = ds.t_min * 0.8;
		const double t = t_min_eff + (ds.t_max - t_min_eff) / (n_dsdt_samples - 1) * i;
		const double f_dsdt = ds.ff->Eval(t);
		const double f_B = (log(ds.ff->Eval(t + ep)) - log(f_dsdt)) / ep;

		int idx = m_g_ext["dsdt"]->GetN();
		m_g_ext["dsdt"]->SetPoint(idx, t, f_dsdt);
		m_g_ext["B"]->SetPoint(idx, t, f_B);
	}

	for (const auto &p : m_g_ext)
		p.second->Write();

	// evaluate uncertainty bands
	map<string, TGraph *> m_g_ext_unc_p;
	m_g_ext_unc_p["dsdt"] = new TGraph();
	m_g_ext_unc_p["B"] = new TGraph();
	SaveParameterUncertaintyBand(result, extrapolationModel, ds, par, ep, m_g_ext["dsdt"], m_g_ext_unc_p["dsdt"], m_g_ext_unc_p["B"]);

	map<string, TGraph *> m_g_ext_unc_t;
	m_g_ext_unc_t["dsdt"] = new TGraph();
	m_g_ext_unc_t["B"] = new TGraph();
	if (calculateTUnc)
		SaveTValueUncertaintyBand(result, extrapolationModel, ds, par, ep, m_g_ext["dsdt"], m_g_ext_unc_t["dsdt"], m_g_ext_unc_t["B"]);

	// make plots
	TDirectory *d_top = gDirectory;

	for (const string &uncType : { "p", "t", "c" })
	{
		gDirectory = d_top->mkdir(("unc " + uncType).c_str());

		for (const string &quantity : { "dsdt", "B" } )
		{
			const TGraph *g_ref = m_g_ext[quantity];

			TGraph *g_unc = new TGraph(); g_unc->SetName(("g_"+quantity+"_ext_unc").c_str());
			TGraph *g_pl_unc = new TGraph(); g_pl_unc->SetName(("g_"+quantity+"_ext_pl_unc").c_str());
			TGraph *g_mi_unc = new TGraph(); g_mi_unc->SetName(("g_"+quantity+"_ext_mi_unc").c_str());

			for (int i = 0; i < g_ref->GetN(); ++i)
			{
				const double t = g_ref->GetX()[i];
				const double v_ref = g_ref->GetY()[i];

				const double unc_p = m_g_ext_unc_p[quantity]->GetY()[i];
				const double unc_t = (calculateTUnc) ? m_g_ext_unc_t[quantity]->GetY()[i] : 0.;

				double unc = 0.;
				if (uncType == "p") unc = unc_p;
				if (uncType == "t") unc = unc_t;
				if (uncType == "c") unc = sqrt(unc_p*unc_p + unc_t*unc_t);

				g_unc->SetPoint(i, t, unc);
				g_pl_unc->SetPoint(i, t, v_ref + unc);
				g_mi_unc->SetPoint(i, t, v_ref - unc);
			}

			g_unc->Write();
			g_pl_unc->Write();
			g_mi_unc->Write();
		}
	}

	gDirectory = d_top;
}

//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: s_extrapolation [option] [option] ...\n");
	printf("OPTIONS:\n");
	printf("    -model <string>\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string tRange = "";

	// parse command line
	for (int argi = 1; (argi < argc) && (cl_error == 0); ++argi)
	{
		if (strcmp(argv[argi], "-h") == 0 || strcmp(argv[argi], "--help") == 0)
		{
			cl_error = 1;
			continue;
		}

		if (TestStringParameter(argc, argv, argi, "-range", tRange)) continue;
		if (TestStringParameter(argc, argv, argi, "-model", fitModel)) continue;

		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// input validation
	if (tRange == "" || fitModel == "")
	{
		printf("ERROR: some mandatory input not provided.\n");
		PrintUsage();
	}

	// prepare datasets
	InitDatasets();

	for (auto &ds : datasets)
	{
		BuildTRanges(tRange, ds);
		BuildFitFunction(fitModel, ds);
	}

	BuildTRanges(tRange, ds_ext);
	BuildFitFunction(fitModel, ds_ext);

	n_ff_parameters = datasets.front().ff->GetNpar();

	// get input
	TFile *f_in = TFile::Open("do_fits.root");

	TGraph *g_settings = (TGraph *) f_in->Get("g_settings");

	for (const auto ds : datasets)
	{
		if (!ds.use_for_extrapolation)
			continue;

		TVectorD *par = (TVectorD *) f_in->Get((ds.name + "/par").c_str());
		TMatrixDSym *par_V = (TMatrixDSym *) f_in->Get((ds.name + "/par_V").c_str());

		constraints.emplace_back(ds.sqrt_s, *par, *par_V);
	}

	// prepare output
	TFile *f_out = TFile::Open("s_extrapolation.root", "recreate");

	g_settings->Write("g_settings");

	// do extrapolation
	for (string extrapolationModel : { "sqrt_s", "log_sqrt_s" })
	{
		TDirectory *d_ext_mod = f_out->mkdir(extrapolationModel.c_str());

		for (bool useCorrelationMatrix : { false, true })
		{
			printf("--------------------------------------------------\n");
			printf("ext model = %s, use correlation matrix = %i\n", extrapolationModel.c_str(), useCorrelationMatrix);
			printf("--------------------------------------------------\n");


			TDirectory *d_ucm = d_ext_mod->mkdir( (useCorrelationMatrix) ? "corr" : "uncorr");
			gDirectory = d_ucm;

			// initialise fitter
			ROOT::Fit::Fitter fitter;

			S2_FCN s2Fcn;
			s2Fcn.extrapolationModel = extrapolationModel;
			s2Fcn.useCorrelationMatrix = useCorrelationMatrix;

			double pStart[2*n_ff_parameters];
  			fitter.SetFCN(2*n_ff_parameters, s2Fcn, pStart, 0, true);

			for (int i = 0; i < n_ff_parameters; i++)
			{
				char buf[100];

				sprintf(buf, "a%i", i);
				fitter.Config().ParSettings(i).Set(buf, 1., 1.);

				sprintf(buf, "b%i", i);
				fitter.Config().ParSettings(n_ff_parameters + i).Set(buf, 0., 1.);
			}

			// run fit
			fitter.FitFCN();

			// process fit
			const ROOT::Fit::FitResult &result = fitter.Result();

			SaveFitPlots(result, extrapolationModel);

			SaveExtrapolation(result, extrapolationModel, ds_ext, true);

			for (const auto &ds : datasets)
			{
				gDirectory = d_ucm->mkdir(ds.name.c_str());

				SaveExtrapolation(result, extrapolationModel, ds);
			}
		}
	}

	// clean up
	delete f_out;

	return 0;
}
