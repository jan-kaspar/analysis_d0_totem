#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TMatrixDSymEigen.h"
#include "Fit/Fitter.h"
#include "TMinuitMinimizer.h"

#include <vector>
#include <map>
#include <string>

#include "command_line_tools.h"
#include "datasets.h"
#include "stat.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

struct BinData
{
	double t, t_unc, dsdt, dsdt_unc_stat;
	double dsdt_unc_syst_t_dep;
	double dsdt_unc_syst_full;
};

//----------------------------------------------------------------------------------------------------

struct InputData
{
	vector<BinData> binData;

	TMatrixDSym m_cov;
	TMatrixDSym m_cov_inv;

	double eta_unc;
};

//----------------------------------------------------------------------------------------------------

bool useStatUnc = true;
bool useSystUnc = true;
bool useNormUnc = false;

void LoadInput(const Dataset &ds, InputData &id, bool useCheckFile = false)
{
	TFile *f_in = TFile::Open((useCheckFile) ? ds.f_in_check.c_str() : ds.f_in.c_str());

	// get input
	TGraph *g_dsdt = (TGraph *) f_in->Get("g_dsdt");
	TMatrixDSym *m_dsdt_cov_syst_t_dep = (TMatrixDSym *) f_in->Get("m_dsdt_cov_syst_t_dep");
	TVectorD *v_dsdt_syst_t_indep = (TVectorD *) f_in->Get("v_dsdt_syst_t_indep");
	TVectorD *v_rel_syst_t_indep = (TVectorD *) f_in->Get("v_rel_syst_t_indep");

	id.eta_unc = (*v_rel_syst_t_indep)(0);

	// get bin data
	double t_min_eff = ds.t_min;

	int idx_min = 100000, idx_max = -1000000;
	for (int i = 0; i < g_dsdt->GetN(); ++i)
	{
		BinData bd;
		bd.t = g_dsdt->GetX()[i];
		bd.t_unc = g_dsdt->GetErrorX(i);
		bd.dsdt = g_dsdt->GetY()[i];
		bd.dsdt_unc_stat = g_dsdt->GetErrorY(i);

		if (bd.t < t_min_eff || bd.t > ds.t_max)
			continue;

		idx_min = min(i, idx_min);
		idx_max = max(i, idx_max);

		id.binData.push_back(std::move(bd));
	}

	// prepare convariance matrices
	int dim = idx_max - idx_min + 1;

	id.m_cov.ResizeTo(dim, dim);
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			const double con_stat = (i == j) ? pow(id.binData[i].dsdt_unc_stat, 2.) : 0.;
			const double con_syst_t_dep = (*m_dsdt_cov_syst_t_dep)(idx_min + i, idx_min + j);
			const double con_syst_t_indep = (*v_dsdt_syst_t_indep)(idx_min + i) * (*v_dsdt_syst_t_indep)(idx_min + j);

			id.m_cov(i, j) = 0.;
			if (useStatUnc) id.m_cov(i, j) += con_stat;
			if (useSystUnc) id.m_cov(i, j) += con_syst_t_dep;

			// normalisation uncertainty handled with its nuissance parameter (eta)
			//if (useNormUnc) id.m_cov(i, j) += con_syst_t_indep;

			if (i == j)
			{
				id.binData[i].dsdt_unc_syst_t_dep = sqrt(con_syst_t_dep);
				id.binData[i].dsdt_unc_syst_full = sqrt(con_syst_t_dep + con_syst_t_indep);

				if (!useCheckFile)
				{
					const auto &bd = id.binData[i];
					const double dsdt_full_unc = sqrt(bd.dsdt_unc_stat*bd.dsdt_unc_stat + bd.dsdt_unc_syst_full*bd.dsdt_unc_syst_full);
					printf("%2i | %.4f, %.4f | %.4f, %.4f\n", i, bd.t, bd.t_unc, bd.dsdt, dsdt_full_unc);
				}
			}
		}
	}

	id.m_cov_inv.ResizeTo(id.m_cov);
	id.m_cov_inv = id.m_cov;
	id.m_cov_inv.Invert();

	delete f_in;
}

//----------------------------------------------------------------------------------------------------

class S2_FCN
{
	public:
		S2_FCN() {}

  		double operator() (const double *par) const;

		const Dataset *ds;
		const InputData *id;
};

//----------------------------------------------------------------------------------------------------

double S2_FCN::operator() (const double *par) const
{
	// transfer parameters to function
	const double eta = par[0];

	for (int i = 0; i < ds->ff->GetNpar(); ++i)
		ds->ff->SetParameter(i, par[i+1]);

	// evaluate difference for each bin
	vector<double> diff(id->binData.size());
	for (unsigned int i = 0; i < diff.size(); ++i)
		diff[i] = id->binData[i].dsdt - eta * ds->ff->Eval(id->binData[i].t);

	// calculate chi^2 from bins
	double S2 = 0.;
	for (unsigned int i = 0; i < diff.size(); ++i)
	{
		for (unsigned int j = 0; j < diff.size(); ++j)
		{
			S2 += diff[i] * id->m_cov_inv(i, j) * diff[j];
		}
	}

	// add normalisation constraint
	const double eta_unc_eff = (useNormUnc) ? id->eta_unc : 1E-6;
	const double ed = (eta - 1.) / eta_unc_eff;
	S2 += ed*ed;

	// add constraints
	for (const auto &p : ds->constraints)
	{
		double rd = (par[p.first + 1] - p.second.mean) / p.second.sigma;
		S2 += rd*rd;
	}

	return S2;
}

//----------------------------------------------------------------------------------------------------

void SaveInputPlots(const Dataset &ds, const InputData &id)
{
	TGraphErrors *g_dsdt = new TGraphErrors();
	TGraph *g_dsdt_unc_syst_t_dep = new TGraph();
	TGraph *g_dsdt_unc_syst_full = new TGraph();

	for (unsigned int i = 0; i < id.binData.size(); ++i)
	{
		const auto &d = id.binData[i];

		int idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, d.t, d.dsdt);
		g_dsdt->SetPointError(idx, d.t_unc, d.dsdt_unc_stat);

		g_dsdt_unc_syst_t_dep->SetPoint(idx, d.t, d.dsdt_unc_syst_t_dep);
		g_dsdt_unc_syst_full->SetPoint(idx, d.t, d.dsdt_unc_syst_full);
	}

	g_dsdt->Write("g_dsdt");
	g_dsdt_unc_syst_t_dep->Write("g_dsdt_unc_syst_t_dep");
	g_dsdt_unc_syst_full->Write("g_dsdt_unc_syst_full");
}

//----------------------------------------------------------------------------------------------------

void SaveFitResults(const ROOT::Fit::FitResult &result)
{
	// number of fit-function parameters
	int dim = result.NPar() - 1;

	TVectorD p(dim);
	TMatrixDSym V(dim);
	TMatrixDSym C(dim);

	TH2D* h2_C = new TH2D("h2_C", "", dim, -0.5, dim-0.5, dim, -0.5, dim-0.5);

	TGraphErrors *g_par = new TGraphErrors();

	for (int i = 0; i < dim; ++i)
	{
		p(i) = result.Parameter(i+1);

		for (int j = 0; j < dim; ++j)
		{
			V(i, j) = result.CovMatrix(i+1, j+1);
			C(i, j) = result.Correlation(i+1, j+1);
			h2_C->SetBinContent(i+1, j+1, C(i, j));
		}

		g_par->SetPoint(i, i, p(i));
		g_par->SetPointError(i, i, sqrt(V(i, i)));
	}

	p.Write("par");
	V.Write("par_V");
	C.Write("par_C");

	g_par->Write("g_par");
	h2_C->Write();
}

//----------------------------------------------------------------------------------------------------

void AnalyzeFit(const TF1 *ff, const Dataset &ds, double &t_dip, double &dsdt_dip, double &t_bmp, double &dsdt_bmp)
{
	dsdt_dip = 1E100;
	for (double t = (ds.t_min + ds.t_dip)/2.; t < (ds.t_dip + ds.t_bmp)/2.; t += 0.001)
	{
		const double dsdt = ff->Eval(t);
		if (dsdt < dsdt_dip)
		{
			dsdt_dip = dsdt;
			t_dip = t;
		}
	}

	dsdt_bmp = -1E100;
	for (double t = (ds.t_dip + ds.t_bmp)/2.; t < (ds.t_bmp + ds.t_max)/2.; t += 0.001)
	{
		const double dsdt = ff->Eval(t);
		if (dsdt > dsdt_bmp)
		{
			dsdt_bmp = dsdt;
			t_bmp = t;
		}
	}
}

//----------------------------------------------------------------------------------------------------

const unsigned int n_points = 200;

void BuildUncertaintyBand(const Dataset &ds, const ROOT::Fit::FitResult &result,
		double &t_dip, double &t_dip_unc, double &dsdt_dip, double &dsdt_dip_unc,
		double &t_bmp, double &t_bmp_unc, double &dsdt_bmp, double &dsdt_bmp_unc,
		double &R, double &R_unc)
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

	// get central values
	const double dsdt_at_t_min = ds.ff->Eval(ds.t_min);
	const double dsdt_at_t_max = ds.ff->Eval(ds.t_max);

	AnalyzeFit(ds.ff, ds, t_dip, dsdt_dip, t_bmp, dsdt_bmp);
	R = dsdt_bmp / dsdt_dip;

	// uncertainty band
	const unsigned int n_repetitions = 1000;
	vector<Stat> stat(n_points, Stat(1));
	
	TH1D *h_dsdt_mod_pf = new TH1D("h_dsdt_mod_pf", "", 100, 0., 0.);
	TH1D *h_dsdt_mod_pl = new TH1D("h_dsdt_mod_pl", "", 100, 0., 0.);

	Stat s_t_dip(1), s_dsdt_dip(1), s_t_bmp(1), s_dsdt_bmp(1), s_R(1);

	for (unsigned int ir = 0; ir < n_repetitions; ++ir)
	{
		// generate biased fit
		TVectorD delta(V.GetNrows()), rdm(V.GetNrows());
		for (int i = 0; i < V.GetNrows(); i++)
			rdm(i) = gRandom->Gaus();
		delta = m_gen * rdm;

		TF1 *ff_mod = new TF1(*ds.ff);

		const double eta = 1.; // want to check the effect of uncertainty of the fit function only

		for (int i = 0; i < ff_mod->GetNpar(); i++)
			ff_mod->SetParameter(i, ff_mod->GetParameter(i) + delta(i+1));
		
		// acceptable configuration ?
		const double dsdt_mod_at_t_min = eta * ff_mod->Eval(ds.t_min);
		const double dsdt_mod_at_t_max = eta * ff_mod->Eval(ds.t_max);
		const bool accept = (
			fabs(dsdt_mod_at_t_min - dsdt_at_t_min) / dsdt_at_t_min < 2.0
			&& fabs(dsdt_mod_at_t_max - dsdt_at_t_max) / dsdt_at_t_max < 2.0
		);

		if (!accept)
			continue;

		// sample biased fit
		for (unsigned int i = 0; i < n_points; ++i)
		{
			const double t = ds.t_min + (ds.t_max - ds.t_min) / (n_points - 1) * i;
			const double f_dsdt = eta * ff_mod->Eval(t);

			stat[i].Fill(f_dsdt);

			if (i == 0) h_dsdt_mod_pf->Fill(f_dsdt);
			if (i == n_points - 1) h_dsdt_mod_pl->Fill(f_dsdt);
		}

		// analyze biased fit
		double t_dip_mod, dsdt_dip_mod, t_bmp_mod, dsdt_bmp_mod;
		AnalyzeFit(ff_mod, ds, t_dip_mod, dsdt_dip_mod, t_bmp_mod, dsdt_bmp_mod);

		s_t_dip.Fill(t_dip_mod);
		s_dsdt_dip.Fill(dsdt_dip_mod);
		s_t_bmp.Fill(t_bmp_mod);
		s_dsdt_bmp.Fill(dsdt_bmp_mod);
		s_R.Fill(dsdt_bmp_mod / dsdt_dip_mod);

		// clean up
		delete ff_mod;
	}

	h_dsdt_mod_pf->Write();
	h_dsdt_mod_pl->Write();

	// store results
	t_dip_unc = s_t_dip.GetStdDev(0);
	dsdt_dip_unc = s_dsdt_dip.GetStdDev(0);
	t_bmp_unc = s_t_bmp.GetStdDev(0);
	dsdt_bmp_unc = s_dsdt_bmp.GetStdDev(0);
	R_unc = s_R.GetStdDev(0);

	// build graphs
	TGraph *g_fit = new TGraph();
	TGraph *g_fit_unc = new TGraph();
	TGraph *g_fit_pl_unc = new TGraph();
	TGraph *g_fit_mi_unc = new TGraph();

	for (unsigned int i = 0; i < n_points; ++i)
	{
		const double t = ds.t_min + (ds.t_max - ds.t_min) / (n_points - 1) * i;

		const double f_dsdt = ds.ff->Eval(t);
		const double f_unc = stat[i].GetStdDev(0);

		int idx = g_fit->GetN();
		g_fit->SetPoint(idx, t, f_dsdt);
		g_fit_unc->SetPoint(idx, t, f_unc);
		g_fit_pl_unc->SetPoint(idx, t, f_dsdt + f_unc);
		g_fit_mi_unc->SetPoint(idx, t, f_dsdt - f_unc);
	}

	g_fit->Write("g_fit");
	g_fit_unc->Write("g_fit_unc");
	g_fit_pl_unc->Write("g_fit_pl_unc");
	g_fit_mi_unc->Write("g_fit_mi_unc");
}

//----------------------------------------------------------------------------------------------------

void BuildSensitivtyPlot(const Dataset &ds)
{
	TDirectory *d_top = gDirectory;

	gDirectory = d_top->mkdir("sensitivity");

	for (int pi = 0; pi < ds.ff->GetNpar(); ++pi)
	{
		for (const string &modif : { "mi", "pl" })
		{
			TF1 *ff_mod = new TF1(*ds.ff);

			double rel_corr = 1.;
			if (modif == "mi") rel_corr = 0.9;
			if (modif == "pl") rel_corr = 1.1;

			ff_mod->SetParameter(pi, ff_mod->GetParameter(pi) * rel_corr);

			TGraph *g = new TGraph();
			char buf[100];
			sprintf(buf, "%i_%s", pi, modif.c_str());
			g->SetName(buf);

			for (unsigned int i = 0; i < n_points; ++i)
			{
				const double t = ds.t_min + (ds.t_max - ds.t_min) / (n_points - 1) * i;

				int idx = g->GetN();
				g->SetPoint(idx, t, ff_mod->Eval(t));
			}

			g->Write();

			delete ff_mod;
		}
	}

	gDirectory = d_top;
}

//----------------------------------------------------------------------------------------------------

void BuildComponentPlots(const Dataset &ds, const string fitModel)
{
	TDirectory *d_top = gDirectory;

	gDirectory = d_top->mkdir("components");

	// define components
	vector<TF1*> components;

	if (fitModel.find("e012+e012") == 0 || fitModel.find("e012+e023") == 0)
	{
		TF1 *ff_comp1 = new TF1(*ds.ff);
		ff_comp1->SetName("g_comp1");
		ff_comp1->SetParameter(3, -1E10);
		ff_comp1->SetParameter(4, 0.);
		ff_comp1->SetParameter(5, 0.);
		components.push_back(ff_comp1);

		TF1 *ff_comp2 = new TF1(*ds.ff);
		ff_comp2->SetName("g_comp2");
		ff_comp2->SetParameter(0, -1E10);
		ff_comp2->SetParameter(1, 0.);
		ff_comp2->SetParameter(2, 0.);
		components.push_back(ff_comp2);
	}

	if (fitModel.find("e012+e02") == 0)
	{
		TF1 *ff_comp1 = new TF1(*ds.ff);
		ff_comp1->SetName("g_comp1");
		ff_comp1->SetParameter(3, -1E10);
		ff_comp1->SetParameter(4, 0.);
		components.push_back(ff_comp1);

		TF1 *ff_comp2 = new TF1(*ds.ff);
		ff_comp2->SetName("g_comp2");
		ff_comp2->SetParameter(0, -1E10);
		ff_comp2->SetParameter(1, 0.);
		ff_comp2->SetParameter(2, 0.);
		components.push_back(ff_comp2);
	}

	if (fitModel.find("e01+e023") == 0)
	{
		TF1 *ff_comp1 = new TF1(*ds.ff);
		ff_comp1->SetName("g_comp1");
		ff_comp1->SetParameter(2, -1E10);
		ff_comp1->SetParameter(3, 0.);
		ff_comp1->SetParameter(4, 0.);
		components.push_back(ff_comp1);

		TF1 *ff_comp2 = new TF1(*ds.ff);
		ff_comp2->SetName("g_comp2");
		ff_comp2->SetParameter(0, -1E10);
		ff_comp2->SetParameter(1, 0.);
		components.push_back(ff_comp2);
	}

	if (fitModel.find("e01-int-e01") == 0)
	{
		TF1 *ff_comp1 = new TF1(*ds.ff);
		ff_comp1->SetName("g_comp1");
		ff_comp1->SetParameter(2, -1E10);
		ff_comp1->SetParameter(3, 0.);
		components.push_back(ff_comp1);

		TF1 *ff_comp2 = new TF1(*ds.ff);
		ff_comp2->SetName("g_comp2");
		ff_comp2->SetParameter(0, -1E10);
		ff_comp2->SetParameter(1, 0.);
		components.push_back(ff_comp2);
	}

	// sample components
	for (const auto f : components)
	{
		TGraph *g = new TGraph();
		g->SetName(f->GetName());

		for (unsigned int i = 0; i < n_points; ++i)
		{
			const double t = ds.t_min + (ds.t_max - ds.t_min) / (n_points - 1) * i;

			int idx = g->GetN();
			g->SetPoint(idx, t, f->Eval(t));
		}

		g->Write();
	}

	// clean up
	gDirectory = d_top;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void PrintUsage()
{
	printf("USAGE: do_fits [option] [option] ...\n");
	printf("OPTIONS:\n");
	printf("    -range <string>\n");
	printf("    -model <string>\n");
	printf("    -unc <string>\n");
}

//----------------------------------------------------------------------------------------------------

int main(int argc, const char **argv)
{
	// defaults
	string tRange = "";
	string fitModel = "";
	string uncChoice = "";

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
		if (TestStringParameter(argc, argv, argi, "-unc", uncChoice)) continue;
		
		printf("ERROR: unknown option '%s'.\n", argv[argi]);
		cl_error = 1;
	}

	if (cl_error)
	{
		PrintUsage();
		return 1;
	}

	// input validation
	if (tRange == "" || fitModel == "" || uncChoice == "")
	{
		printf("ERROR: some mandatory input not provided.\n");
		PrintUsage();
	}

	useStatUnc = (uncChoice.find("st") != string::npos);
	useSystUnc = (uncChoice.find("sy") != string::npos);
	useNormUnc = (uncChoice.find("no") != string::npos);

	printf("useStatUnc = %i\n", useStatUnc);
	printf("useSystUnc = %i\n", useSystUnc);
	printf("useNormUnc = %i\n", useNormUnc);

	// prepare datasets
	InitDatasets();

	for (auto &ds : datasets)
	{
		BuildTRanges(tRange, ds);
		BuildFitFunction(fitModel, ds);
	}

	// prepare output
	TFile *f_out = TFile::Open("do_fits.root", "recreate");

	TGraph *g_settings = new TGraph();
	g_settings->SetPoint(0, datasets.front().ff->GetNpar(), 0.);
	g_settings->Write("g_settings");

	// start processing
	for (const auto ds : datasets)
	{
		printf("\n\n----------- %s -----------\n", ds.name.c_str());

		TDirectory *d_ds = f_out->mkdir(ds.name.c_str());

		// get input
		InputData id;
		LoadInput(ds, id);

		// initialise fitter
		ROOT::Fit::Fitter fitter;

		S2_FCN s2Fcn;
		s2Fcn.ds = &ds;
		s2Fcn.id = &id;

		int n_ff_parameters = ds.ff->GetNpar();

		double pStart[n_ff_parameters+1];
		fitter.SetFCN(n_ff_parameters+1, s2Fcn, pStart, 0, true);

		fitter.Config().ParSettings(0).Set("eta", 1., (useNormUnc) ? id.eta_unc : 1E-6);

		for (int i = 0; i < n_ff_parameters; ++i)
		{
			double lim_min, lim_max;
			ds.ff->GetParLimits(i, lim_min, lim_max);

			if (lim_max > lim_min)
				fitter.Config().ParSettings(i+1).Set(ds.ff->GetParName(i), ds.ff->GetParameter(i), fabs(ds.ff->GetParameter(i)) * 0.01, lim_min, lim_max);
			else
				fitter.Config().ParSettings(i+1).Set(ds.ff->GetParName(i), ds.ff->GetParameter(i), fabs(ds.ff->GetParameter(i)) * 0.01);
		}

		// do fit
		fitter.FitFCN();
		fitter.FitFCN();

		// print results
		const ROOT::Fit::FitResult &result = fitter.Result();
		for (unsigned int i = 0; i < result.NPar(); ++i)
		{
			printf("idx %u [%3s]: %+.3E +- %.3E\n", i, result.ParName(i).c_str(), result.Parameter(i), sqrt(result.CovMatrix(i, i)));
		}

		// result check
		if (fabs(result.Parameter(0) - 1.) > 5E-4)
			printf("ERROR: eta is significantly different from 1: eta-1 = %.1E\n", result.Parameter(0)-1.);

		// update fit function with final fit parameters
		for (int i = 0; i < ds.ff->GetNpar(); ++i)
			ds.ff->SetParameter(i, result.Parameter(i+1));

		// save fit results
		gDirectory = d_ds;

		SaveInputPlots(ds, id);

		ds.ff->Write("f_fit");

		SaveFitResults(result);

		double t_dip, t_dip_unc;
		double dsdt_dip, dsdt_dip_unc;
		double t_bmp, t_bmp_unc;
		double dsdt_bmp, dsdt_bmp_unc;
		double R, R_unc;
		BuildUncertaintyBand(ds, result, t_dip, t_dip_unc, dsdt_dip, dsdt_dip_unc, t_bmp, t_bmp_unc, dsdt_bmp, dsdt_bmp_unc, R, R_unc);

		BuildSensitivtyPlot(ds);

		BuildComponentPlots(ds, fitModel);

		int data_points = id.binData.size() + 1;
		int fit_parameters = n_ff_parameters + 1;
		int ndf = data_points - fit_parameters;

		TGraph *g_data = new TGraph();
		g_data->SetPoint(0, ds.t_min, ds.t_max);
		g_data->SetPoint(1, data_points, ndf);
		g_data->SetPoint(2, result.Chi2(), result.Chi2() / ndf);
		g_data->SetPoint(3, TMath::Prob(result.Chi2(), ndf), 0.);
		g_data->SetPoint(4, t_dip, t_dip_unc);
		g_data->SetPoint(5, dsdt_dip, dsdt_dip_unc);
		g_data->SetPoint(6, t_bmp, t_bmp_unc);
		g_data->SetPoint(7, dsdt_bmp, dsdt_bmp_unc);
		g_data->SetPoint(8, R, R_unc);
		g_data->Write("g_data");

		// test with another "check" dataset
		InputData id_check;
		LoadInput(ds, id_check, true);
		s2Fcn.id = &id_check;

		double fpar[ds.ff->GetNpar() + 1];
		for (int i = 0; i < ds.ff->GetNpar() + 1; ++i)
			fpar[i] = result.Parameter(i);

		double chi2_check = s2Fcn(fpar);

		gDirectory = d_ds->mkdir("check");

		SaveInputPlots(ds, id_check);

		g_data = new TGraph();
		g_data->SetPoint(0, ds.t_min, ds.t_max);
		g_data->SetPoint(1, data_points, ndf);
		g_data->SetPoint(2, chi2_check, chi2_check / ndf);
		g_data->SetPoint(3, TMath::Prob(chi2_check, ndf), 0.);
		g_data->SetPoint(4, t_dip, t_dip_unc);
		g_data->SetPoint(5, dsdt_dip, dsdt_dip_unc);
		g_data->SetPoint(6, t_bmp, t_bmp_unc);
		g_data->SetPoint(7, dsdt_bmp, dsdt_bmp_unc);
		g_data->SetPoint(8, R, R_unc);
		g_data->Write("g_data");
	}

	// clean up
	delete f_out;

	return 0;
}
