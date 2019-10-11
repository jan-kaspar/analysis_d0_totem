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
};

//----------------------------------------------------------------------------------------------------

bool useStatUnc = true;
bool useSystUnc = true;
bool useNormUnc = false;

void LoadInput(const Dataset &ds, InputData &id)
{
	TFile *f_in = TFile::Open(ds.f_in.c_str());

	// get input
	TGraph *g_dsdt = (TGraph *) f_in->Get("g_dsdt");
	TMatrixDSym *m_dsdt_cov_syst_t_dep = (TMatrixDSym *) f_in->Get("m_dsdt_cov_syst_t_dep");
	TVectorD *v_dsdt_syst_t_indep = (TVectorD *) f_in->Get("v_dsdt_syst_t_indep");

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
			if (useNormUnc) id.m_cov(i, j) += con_syst_t_indep;

			if (i == j)
			{
				id.binData[i].dsdt_unc_syst_t_dep = sqrt(con_syst_t_dep);
				id.binData[i].dsdt_unc_syst_full = sqrt(con_syst_t_dep + con_syst_t_indep);
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
	ds->ff->SetParameters(par);

	vector<double> diff(id->binData.size());
	for (unsigned int i = 0; i < diff.size(); ++i)
		diff[i] = id->binData[i].dsdt - ds->ff->Eval(id->binData[i].t);

	double S2 = 0.;
	for (unsigned int i = 0; i < diff.size(); ++i)
	{
		for (unsigned int j = 0; j < diff.size(); ++j)
		{
			S2 += diff[i] * id->m_cov_inv(i, j) * diff[j];
		}
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
	int dim = result.NPar();

	TVectorD p(dim);
	TMatrixDSym V(dim);
	TMatrixDSym C(dim);

	TH2D* h2_C = new TH2D("h2_C", "", dim, -0.5, dim-0.5, dim, -0.5, dim-0.5);

	for (int i = 0; i < dim; ++i)
	{
		p(i) = result.Parameter(i);

		for (int j = 0; j < dim; ++j)
		{
			V(i, j) = result.CovMatrix(i, j);
			C(i, j) = result.Correlation(i, j);
			h2_C->SetBinContent(i+1, j+1, C(i, j));
		}
	}

	p.Write("par");
	V.Write("par_V");
	C.Write("par_C");
	h2_C->Write();
}

//----------------------------------------------------------------------------------------------------

const unsigned int n_points = 200;

void BuildUncertaintyBand(const Dataset &ds, const ROOT::Fit::FitResult &result)
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
	const unsigned int n_repetitions = 1000;
	vector<Stat> stat(n_points, Stat(1));
	
	TH1D *h_dsdt_mod_pf = new TH1D("h_dsdt_mod_pf", "", 100, 0., 0.);
	TH1D *h_dsdt_mod_pl = new TH1D("h_dsdt_mod_pl", "", 100, 0., 0.);

	const double dsdt_at_t_min = ds.ff->Eval(ds.t_min);
	const double dsdt_at_t_max = ds.ff->Eval(ds.t_max);

	for (unsigned int ir = 0; ir < n_repetitions; ++ir)
	{
		// generate biased fit
		TVectorD delta(V.GetNrows()), rdm(V.GetNrows());
		for (int i = 0; i < V.GetNrows(); i++)
			rdm(i) = gRandom->Gaus();
		delta = m_gen * rdm;

		TF1 *ff_mod = new TF1(*ds.ff);

		for (int i = 0; i < V.GetNrows(); i++)
			ff_mod->SetParameter(i, ff_mod->GetParameter(i) + delta(i));
		
		// acceptable configuration ?
		const double dsdt_mod_at_t_min = ff_mod->Eval(ds.t_min);
		const double dsdt_mod_at_t_max = ff_mod->Eval(ds.t_max);
		const bool accept = (
			fabs(dsdt_mod_at_t_min - dsdt_at_t_min) / dsdt_at_t_min < 2.0
			&& fabs(dsdt_mod_at_t_max - dsdt_at_t_max) / dsdt_at_t_max < 5.0
		);

		if (!accept)
			continue;

		// sample biased fit
		for (unsigned int i = 0; i < n_points; ++i)
		{
			const double t = ds.t_min + (ds.t_max - ds.t_min) / (n_points - 1) * i;
			const double f_dsdt = ff_mod->Eval(t);

			stat[i].Fill(f_dsdt);

			if (i == 0) h_dsdt_mod_pf->Fill(f_dsdt);
			if (i == n_points - 1) h_dsdt_mod_pl->Fill(f_dsdt);
		}

		// clean up
		delete ff_mod;
	}

	h_dsdt_mod_pf->Write();
	h_dsdt_mod_pl->Write();

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
	g_settings->SetPoint(0, n_parameters, 0.);
	g_settings->Write("g_settings");

	// start processing
	for (const auto ds : datasets)
	{
		printf("\n\n----------- %s -----------\n", ds.name.c_str());

		TDirectory *d_ds = f_out->mkdir(ds.name.c_str());

		// get input
		InputData id;
		LoadInput(ds, id);

		// do fit
		ROOT::Fit::Fitter fitter;

		S2_FCN s2Fcn;
		s2Fcn.ds = &ds;
		s2Fcn.id = &id;

		int n_parameters = ds.ff->GetNpar();

		double pStart[n_parameters];
		fitter.SetFCN(n_parameters, s2Fcn, pStart, 0, true);

		for (int i = 0; i < n_parameters; ++i)
		{
			fitter.Config().ParSettings(i).Set(ds.ff->GetParName(i), ds.ff->GetParameter(i), fabs(ds.ff->GetParameter(i)) * 0.05);
		}

		fitter.FitFCN();
		fitter.FitFCN();

		// save fit results
		gDirectory = d_ds;

		SaveInputPlots(ds, id);

		const ROOT::Fit::FitResult &result = fitter.Result();
		ds.ff->SetParameters(result.GetParams());

		ds.ff->Write("f_fit");

		SaveFitResults(result);

		int ndf = id.binData.size() - n_parameters;
		TGraph *g_data = new TGraph();
		g_data->SetPoint(0, ds.t_min, ds.t_max);
		g_data->SetPoint(1, id.binData.size(), ndf);
		g_data->SetPoint(2, result.Chi2(), result.Chi2() / ndf);
		g_data->Write("g_data");

		BuildUncertaintyBand(ds, result);

		BuildSensitivtyPlot(ds);
	}

	// clean up
	delete f_out;

	return 0;
}
