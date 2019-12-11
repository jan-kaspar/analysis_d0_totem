#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"
#include "TMatrixDSym.h"
#include "TRandom3.h"
#include "TMatrixDSymEigen.h"

#include "../datasets.h"
#include "../stat.h"

#include <string>
#include <vector>

using namespace std;

//----------------------------------------------------------------------------------------------------

enum { i_t_min, i_dsdt_min, i_t_dip, i_dsdt_dip, i_t_bmp, i_dsdt_bmp, i_t_max, i_dsdt_max };

//----------------------------------------------------------------------------------------------------

void Analyze(const TF1 *f_fit, const Dataset &ds, const string &method, vector<double> &par)
{
	const double dt = 0.001;

	// find dip
	double t_dip = 0., dsdt_dip = 1E100;
	for (double t = ds.t_dip - 0.05; t < ds.t_dip + 0.05; t += dt)
	{
		const double dsdt = f_fit->Eval(t);

		if (dsdt < dsdt_dip)
		{
			t_dip = t;
			dsdt_dip = dsdt;
		}
	}

	// find bump
	double t_bmp = 0., dsdt_bmp = -1E100;
	for (double t = ds.t_bmp - 0.05; t < ds.t_bmp + 0.05; t += dt)
	{
		const double dsdt = f_fit->Eval(t);

		if (dsdt > dsdt_bmp)
		{
			t_bmp = t;
			dsdt_bmp = dsdt;
		}
	}

	// find t_min
	double dsdt_min = dsdt_bmp * (dsdt_bmp / dsdt_dip);
	if (method.find("low_t") != string::npos) dsdt_min = dsdt_dip * 17.;

	double t_min = 0., diff_min = 1E100;
	for (double t = 0.25; t < ds.t_dip; t += dt)
	{
		const double diff = fabs(dsdt_min - f_fit->Eval(t));

		if (diff < diff_min)
		{
			t_min = t;
			diff_min = diff;
		}
	}

	// find t_max
	double dsdt_max = 0.9 * dsdt_bmp;
	if (method.find("high_t") != string::npos) dsdt_max = dsdt_dip;

	double t_max = 0.;
	diff_min = 1E100;
	for (double t = t_bmp; t < 1.0; t += dt)
	{
		const double diff = fabs(dsdt_max - f_fit->Eval(t));

		if (diff < diff_min)
		{
			t_max = t;
			diff_min = diff;
		}
	}

	// export results
	par[i_t_min] = t_min;
	par[i_dsdt_min] = dsdt_min;

	par[i_t_dip] = t_dip;
	par[i_dsdt_dip] = dsdt_dip;

	par[i_t_bmp] = t_bmp;
	par[i_dsdt_bmp] = dsdt_bmp;

	par[i_t_max] = t_max;
	par[i_dsdt_max] = dsdt_max;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	// config
	struct Method
	{
		string name;
	};

	vector<Method> methods = {
		{ "minimal" },
		{ "low_t,high_t" },
	};

	InitDatasets();
	for (auto &ds : datasets)
		BuildTRanges("bootstrap", ds);

	// prepare input
	TFile *f_in = TFile::Open("../fits/bootstrap/bootstrap/st+sy+no/do_fits.root");

	// prepare output
	TFile *f_out = TFile::Open("t_investigation_new.root", "recreate");

	// process each energy
	for (const auto &method : methods)
	{
		printf("* %s\n", method.name.c_str());

		TDirectory *d_method = f_out->mkdir(method.name.c_str());

		TGraphErrors *g_t_min_vs_sqrt_s = new TGraphErrors();
		TGraphErrors *g_t_dip_vs_sqrt_s = new TGraphErrors();
		TGraphErrors *g_t_bmp_vs_sqrt_s = new TGraphErrors();
		TGraphErrors *g_t_max_vs_sqrt_s = new TGraphErrors();

		// use the same seed for each method
		gRandom->SetSeed(123);

		for (const auto &ds : datasets)
		{
			printf("  - %s\n", ds.name.c_str());

			//printf("    dataset: %.3f, %.3f\n", ds.t_dip, ds.t_bmp);

			// get input
			TGraphErrors *g_dsdt = (TGraphErrors *) f_in->Get((ds.name + "/g_dsdt").c_str());
			TF1 *ff = (TF1 *) f_in->Get((ds.name + "/f_fit").c_str());
			TMatrixDSym *par_cov = (TMatrixDSym *) f_in->Get((ds.name + "/par_V").c_str());
			TGraph *g_fit = (TGraph *) f_in->Get((ds.name + "/g_fit").c_str());

			// prepare data
			vector<double> pars(8);
			Stat st(8);

			// get central values
			Analyze(ff, ds, method.name, pars);

			// prepare perturbation generator matrix
			int dim = par_cov->GetNrows();
			TMatrixDSymEigen eig_decomp(*par_cov);
			TVectorD eig_values(eig_decomp.GetEigenValues());
			TMatrixDSym S(dim);
			for (int i = 0; i < dim; i++)
				S(i, i) = (eig_values(i) >= 0.) ? sqrt(eig_values(i)) : 0.;
			TMatrixD m_gen = eig_decomp.GetEigenVectors() * S;

			// get uncertainties
			unsigned int n_repetitions = 100;
			for (unsigned int n = 0; n < n_repetitions; ++n)
			{
				// build biased ff
				TVectorD delta(dim), rdm(dim);
				for (int i = 0; i < dim; i++)
					rdm(i) = gRandom->Gaus();
				delta = m_gen * rdm;

				TF1 *ff_mod = new TF1(*ff);

				for (int i = 0; i < dim; i++)
					ff_mod->SetParameter(i, ff_mod->GetParameter(i) + delta(i));
			
				// analyze biased ff
				vector<double> p(8);
				Analyze(ff_mod, ds, method.name, p);
				st.Fill(p);

				// clean up
				delete ff_mod;
			}

			// extract parameters
			double t_min = pars[i_t_min], t_min_unc = st.GetStdDev(i_t_min);
			double dsdt_min = pars[i_dsdt_min], dsdt_min_unc = st.GetStdDev(i_dsdt_min);

			double t_dip = pars[i_t_dip], t_dip_unc = st.GetStdDev(i_t_dip);
			double dsdt_dip = pars[i_dsdt_dip], dsdt_dip_unc = st.GetStdDev(i_dsdt_dip);

			double t_bmp = pars[i_t_bmp], t_bmp_unc = st.GetStdDev(i_t_bmp);
			double dsdt_bmp = pars[i_dsdt_bmp], dsdt_bmp_unc = st.GetStdDev(i_dsdt_bmp);

			double t_max = pars[i_t_max], t_max_unc = st.GetStdDev(i_t_max);
			double dsdt_max = pars[i_dsdt_max], dsdt_max_unc = st.GetStdDev(i_dsdt_max);

			// to somewhat account for uncertainty due to fitting function choice
			t_min_unc *= 2.;
			t_dip_unc *= 2.;
			t_bmp_unc *= 2.;
			t_max_unc *= 2.;

			/*
			// print covariance matrix for t-value uncertainty
			st.PrintCorrelation();
			const vector<unsigned int> idx_sel = { i_t_min, i_t_dip, i_t_bmp, i_t_max };
			for (int i : idx_sel)
			{
				for (int j : idx_sel)
					printf("%.3E, ", st.GetCovariance(i, j));
				printf("\n");
			}
			*/

			// manual correction
			if (ds.name == "2.76TeV")
			{
				t_bmp = 0.79;
				t_bmp_unc = 0.15;

				if (method.name == "minimal") { t_max = 0.85; t_max_unc = 0.1; }
				if (method.name == "low_t,high_t") { t_max = 0.95; t_max_unc = 0.1; }
			}

			if (ds.name == "13TeV")
			{
				t_min_unc = max(t_min_unc, 0.003);
				t_dip_unc = max(t_dip_unc, 0.003);
				t_bmp_unc = max(t_bmp_unc, 0.003);
				t_max_unc = max(t_bmp_unc, 0.005);
			}

			// print results
			printf("    t_min = %.3f +- %.3f\n", t_min, t_min_unc);
			printf("    t_dip = %.3f +- %.3f\n", t_dip, t_dip_unc);
			printf("    t_bmp = %.3f +- %.3f\n", t_bmp, t_bmp_unc);
			printf("    t_max = %.3f +- %.3f\n", t_max, t_max_unc);

			// save results
			gDirectory = d_method->mkdir(ds.name.c_str());

			g_dsdt->Write("g_dsdt");

			TGraph *g_data = new TGraph();
			g_data->SetPoint(0, t_min, t_min_unc);
			g_data->SetPoint(1, dsdt_min, dsdt_min_unc);

			g_data->SetPoint(2, t_dip, t_dip_unc);
			g_data->SetPoint(3, dsdt_dip, dsdt_dip_unc);

			g_data->SetPoint(4, t_bmp, t_bmp_unc);
			g_data->SetPoint(5, dsdt_bmp, dsdt_bmp_unc);

			g_data->SetPoint(6, t_max, t_max_unc);
			g_data->SetPoint(7, dsdt_max, dsdt_max_unc);

			g_data->Write("g_data");

			g_fit->Write("g_fit");

			// update plots
			int idx = g_t_min_vs_sqrt_s->GetN();

			g_t_min_vs_sqrt_s->SetPoint(idx, ds.sqrt_s, t_min);
			g_t_min_vs_sqrt_s->SetPointError(idx, 0., t_min_unc);

			g_t_dip_vs_sqrt_s->SetPoint(idx, ds.sqrt_s, t_dip);
			g_t_dip_vs_sqrt_s->SetPointError(idx, 0., t_dip_unc);

			g_t_bmp_vs_sqrt_s->SetPoint(idx, ds.sqrt_s, t_bmp);
			g_t_bmp_vs_sqrt_s->SetPointError(idx, 0., t_bmp_unc);

			g_t_max_vs_sqrt_s->SetPoint(idx, ds.sqrt_s, t_max);
			g_t_max_vs_sqrt_s->SetPointError(idx, 0., t_max_unc);
		}
		
		// fit s-dependence
		TF1 *ff_lin = new TF1("ff_lin", "[0] + [1]*(x - [2])");
		ff_lin->FixParameter(2, ds_ext.sqrt_s);

		TF1 *ff_log = new TF1("ff_log", "[0] + [1]*(log(x) - log([2]))");
		ff_log->FixParameter(2, ds_ext.sqrt_s);

		g_t_min_vs_sqrt_s->Fit(ff_lin, "Q");
		g_t_dip_vs_sqrt_s->Fit(ff_lin, "Q");
		g_t_bmp_vs_sqrt_s->Fit(ff_lin, "Q");
		g_t_max_vs_sqrt_s->Fit(ff_lin, "Q");

		double t_min_ext, t_min_ext_unc;
		double t_dip_ext, t_dip_ext_unc;
		double t_bmp_ext, t_bmp_ext_unc;
		double t_max_ext, t_max_ext_unc;

		g_t_min_vs_sqrt_s->Fit(ff_log, "Q+"); t_min_ext = ff_log->GetParameter(0); t_min_ext_unc = ff_log->GetParError(0);
		g_t_dip_vs_sqrt_s->Fit(ff_log, "Q+"); t_dip_ext = ff_log->GetParameter(0); t_dip_ext_unc = ff_log->GetParError(0);
		g_t_bmp_vs_sqrt_s->Fit(ff_log, "Q+"); t_bmp_ext = ff_log->GetParameter(0); t_bmp_ext_unc = ff_log->GetParError(0);
		g_t_max_vs_sqrt_s->Fit(ff_log, "Q+"); t_max_ext = ff_log->GetParameter(0); t_max_ext_unc = ff_log->GetParError(0);

		printf("  - extrapolation (%.3f)\n", ds_ext.sqrt_s);
		printf("    t_min = %.3f +- %.3f\n", t_min_ext, t_min_ext_unc);
		printf("    t_dip = %.3f +- %.3f\n", t_dip_ext, t_dip_ext_unc);
		printf("    t_bmp = %.3f +- %.3f\n", t_bmp_ext, t_bmp_ext_unc);
		printf("    t_max = %.3f +- %.3f\n", t_max_ext, t_max_ext_unc);

		printf("    %.3E, 0., 0., 0.,\n", t_min_ext_unc * t_min_ext_unc);
		printf("    0., %.3E, 0., 0.,\n", t_dip_ext_unc * t_dip_ext_unc);
		printf("    0., 0., %.3E, 0.,\n", t_bmp_ext_unc * t_bmp_ext_unc);
		printf("    0., 0., 0., %.3E\n", t_max_ext_unc * t_max_ext_unc);

		// save plots
		gDirectory = d_method;

		g_t_min_vs_sqrt_s->Write("g_t_min_vs_sqrt_s");
		g_t_dip_vs_sqrt_s->Write("g_t_dip_vs_sqrt_s");
		g_t_bmp_vs_sqrt_s->Write("g_t_bmp_vs_sqrt_s");
		g_t_max_vs_sqrt_s->Write("g_t_max_vs_sqrt_s");
	}

	// clean up
	delete f_in;
	delete f_out;

	return 0;
}
