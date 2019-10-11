// TODO: clean
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TFitResult.h"
#include "TMatrixDSymEigen.h"
#include "Fit/Fitter.h"
#include "TMinuitMinimizer.h"
#include "TCanvas.h"

#include <vector>
#include <map>
#include <string>

#include "command_line_tools.h"
#include "datasets.h"
#include "stat.h"

using namespace std;

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

	TVectorD a(n_parameters), b(n_parameters);
	for (int i = 0; i < n_parameters; ++i)
	{
		a(i) = par[i];
		b(i) = par[n_parameters + i];
	}

	for (const auto &cnt : constraints)
	{
		TVectorD p = EvaluateParameters(a, b, cnt.sqrt_s, extrapolationModel);

		TVectorD diff = p - cnt.par;

		for (int i = 0; i < n_parameters; i++)
		{
			for (int j = 0; j < n_parameters; j++)
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
	const double *par = result.GetParams();

	TVectorD a(n_parameters), b(n_parameters);
	for (int i = 0; i < n_parameters; ++i)
	{
		a(i) = par[i];
		b(i) = par[n_parameters + i];
	}

	for (int i = 0; i < n_parameters; ++i)
	{
		char buf[100];

		TGraphErrors *g_cnt = new TGraphErrors();
		g_cnt->SetName("g_cnt");

		for (const auto &cnt : constraints)
		{
			int idx = g_cnt->GetN();
			g_cnt->SetPoint(idx, cnt.sqrt_s, cnt.par(i));
			g_cnt->SetPointError(idx, 0., sqrt(cnt.par_V(i, i)));
		}

		TGraph *g_fit = new TGraph(); g_fit->SetName("g_fit");
		TGraph *g_fit_pl_unc = new TGraph(); g_fit_pl_unc->SetName("g_fit_pl_unc"); g_fit_pl_unc->SetLineStyle(2);
		TGraph *g_fit_mi_unc = new TGraph(); g_fit_mi_unc->SetName("g_fit_mi_unc"); g_fit_pl_unc->SetLineStyle(2);

		for (double W = 1.; W < 14.; W += 0.5)
		{
			TVectorD p = EvaluateParameters(a, b, W, extrapolationModel);

			// TODO
			double unc = 0;

			int idx = g_fit->GetN();
			g_fit->SetPoint(idx, W, p(i));
			g_fit_pl_unc->SetPoint(idx, W, p(i) + unc);
			g_fit_mi_unc->SetPoint(idx, W, p(i) - unc);
		}

		sprintf(buf, "c_par%i", i);
		TCanvas *c = new TCanvas();
		c->SetName(buf);
		g_cnt->Draw("ap");
		g_fit->Draw("l");
		g_fit_pl_unc->Draw("l");
		g_fit_mi_unc->Draw("l");

		c->Write();
	}
}

//----------------------------------------------------------------------------------------------------

const unsigned int n_points = 200;

void SaveExtrapolation(const ROOT::Fit::FitResult &result, const string &extrapolationModel, const Dataset &ds)
{
	// calculate extrapolated parameters
	const double *par = result.GetParams();

	TVectorD a(n_parameters), b(n_parameters);
	for (int i = 0; i < n_parameters; ++i)
	{
		a(i) = par[i];
		b(i) = par[n_parameters + i];
	}

	TVectorD p = EvaluateParameters(a, b, ds_ext.sqrt_s, extrapolationModel);

	// set extrapolated parameters to fit function
	for (int i = 0; i < n_parameters; ++i)
		ds.ff->SetParameter(i, p(i));

	// evaluate extrapolated graph
	TGraph *g_dsdt_ext = new TGraph(); g_dsdt_ext->SetName("g_dsdt_ext");

	for (unsigned int i = 0; i < n_points; ++i)
	{
		const double t = ds.t_min + (ds.t_max - ds.t_min) / (n_points - 1) * i;
		const double f_dsdt = ds.ff->Eval(t);

		int idx = g_dsdt_ext->GetN();
		g_dsdt_ext->SetPoint(idx, t, f_dsdt);
	}

	g_dsdt_ext->Write();

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
	const unsigned int n_repetitions = 100;
	vector<Stat> stat(n_points, Stat(1));

	for (unsigned int ir = 0; ir < n_repetitions; ++ir)
	{
		// generate biased extrapolation
		TVectorD delta(V.GetNrows()), rdm(V.GetNrows());
		for (int i = 0; i < V.GetNrows(); i++)
			rdm(i) = gRandom->Gaus();
		delta = m_gen * rdm;

		TVectorD a_mod(n_parameters), b_mod(n_parameters);
		for (int i = 0; i < n_parameters; ++i)
		{
			a_mod(i) = par[i] + delta(i);
			b_mod(i) = par[n_parameters + i] + delta(n_parameters + i);
		}

		TVectorD p_mod = EvaluateParameters(a_mod, b_mod, ds_ext.sqrt_s, extrapolationModel);

		// set extrapolated parameters to fit function
		TF1 *ff_mod = new TF1(*ds.ff);

		for (int i = 0; i < n_parameters; ++i)
			ff_mod->SetParameter(i, p_mod(i));

		// sample biased extrapolation
		for (int i = 0; i < g_dsdt_ext->GetN(); ++i)
		{
			const double t = g_dsdt_ext->GetX()[i];
			const double f_dsdt_mod = ff_mod->Eval(t);

			stat[i].Fill(f_dsdt_mod);
		}

		// clean up
		delete ff_mod;
	}

	// build graphs
	TGraph *g_dsdt_ext_unc = new TGraph();
	TGraph *g_dsdt_ext_pl_unc = new TGraph();
	TGraph *g_dsdt_ext_mi_unc = new TGraph();

	for (int i = 0; i < g_dsdt_ext->GetN(); ++i)
	{
		const double t = g_dsdt_ext->GetX()[i];
		const double f_dsdt = g_dsdt_ext->GetY()[i];

		const double f_unc = stat[i].GetStdDev(0);

		int idx = g_dsdt_ext_unc->GetN();
		g_dsdt_ext_unc->SetPoint(idx, t, f_unc);
		g_dsdt_ext_pl_unc->SetPoint(idx, t, f_dsdt + f_unc);
		g_dsdt_ext_mi_unc->SetPoint(idx, t, f_dsdt - f_unc);
	}

	g_dsdt_ext_unc->Write("g_dsdt_ext_unc");
	g_dsdt_ext_pl_unc->Write("g_dsdt_ext_pl_unc");
	g_dsdt_ext_mi_unc->Write("g_dsdt_ext_mi_unc");
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
	string fitModel = "";

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

	BuildTRanges(tRange, ds_ext);
	BuildFitFunction(fitModel, ds_ext);

	// get input
	TFile *f_in = TFile::Open("do_fits.root");

	TGraph *g_settings = (TGraph *) f_in->Get("g_settings");

	for (const auto ds : datasets)
	{
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
			TDirectory *d_ucm = d_ext_mod->mkdir( (useCorrelationMatrix) ? "corr" : "uncorr");
			gDirectory = d_ucm;

			// initialise fitter
			ROOT::Fit::Fitter fitter;

			S2_FCN s2Fcn;
			s2Fcn.extrapolationModel = extrapolationModel;
			s2Fcn.useCorrelationMatrix = useCorrelationMatrix;
			
			double pStart[n_parameters];
  			fitter.SetFCN(2*n_parameters, s2Fcn, pStart, 0, true);

			for (int i = 0; i < n_parameters; i++)
			{
				char buf[100];

				sprintf(buf, "a%i", i);
				fitter.Config().ParSettings(i).Set(buf, 1., 1.);

				sprintf(buf, "b%i", i);
				fitter.Config().ParSettings(n_parameters + i).Set(buf, 0., 1.);
			}

			// run fit
			fitter.FitFCN();

			// process fit
			const ROOT::Fit::FitResult &result = fitter.Result();

			SaveFitPlots(result, extrapolationModel);

			SaveExtrapolation(result, extrapolationModel, ds_ext);
		}
	}

	// clean up
	delete f_out;

	return 0;
}
