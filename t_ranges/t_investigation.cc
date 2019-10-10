#include "TFile.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1D.h"

#include <string>
#include <vector>

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// config
	struct Dataset
	{
		double sqrt_s;
		string name;
		double t1, t2, t3;
	};

	vector<Dataset> datasets = {
		{ 2.76, "2.76TeV", 0.55, 0.65, 0.8 },
		{ 7., "7TeV", 0.45, 0.6, 0.8 },
		{ 8., "8TeV", 0.45, 0.6, 0.8 },
		{ 13., "13TeV", 0.4, 0.55, 0.75 }
	};

	Dataset ds_ext{ 1.96, "1.96TeV", 0., 0., 0. };

	struct Method
	{
		string name;
	};

	vector<Method> methods = {
		{ "first" }
	};

	// prepare output
	TFile *f_out = TFile::Open("t_investigation.root", "recreate");

	for (const auto &method : methods)
	{
		TDirectory *d_method = f_out->mkdir(method.name.c_str());

		TGraphErrors *g_t_min_vs_sqrt_s = new TGraphErrors();
		TGraphErrors *g_t_dip_vs_sqrt_s = new TGraphErrors();
		TGraphErrors *g_t_bmp_vs_sqrt_s = new TGraphErrors();
		TGraphErrors *g_t_max_vs_sqrt_s = new TGraphErrors();

		// process
		for (const auto &ds : datasets)
		{
			printf("* %s\n", ds.name.c_str());
	
			// get input
			TFile *f_in = TFile::Open(("../data/TOTEM_" + ds.name + "/data.root").c_str());

			TF1 *f_fit = (TF1 *) f_in->Get("f_dsdt");

			const double dt = 0.001;

			// find dip
			double t_dip = 0., dsdt_dip = 1E100;
			for (double t = ds.t1; t < ds.t2; t += dt)
			{
				const double dsdt = f_fit->Eval(t);

				if (dsdt < dsdt_dip)
				{
					t_dip = t;
					dsdt_dip = dsdt;
				}
			}

			double t_dip_unc = 0.;
			if (ds.name == "2.76TeV") t_dip_unc = 0.01;
			if (ds.name == "7TeV") t_dip_unc = 0.01;
			if (ds.name == "8TeV") t_dip_unc = 0.01;
			if (ds.name == "13TeV") t_dip_unc = 0.005;

			// find bump
			double t_bmp = 0., dsdt_bmp = -1E100;
			for (double t = ds.t2; t < ds.t3; t += dt)
			{
				const double dsdt = f_fit->Eval(t);

				if (dsdt > dsdt_bmp)
				{
					t_bmp = t;
					dsdt_bmp = dsdt;
				}
			}

			double t_bmp_unc = 0.;
			if (ds.name == "2.76TeV") t_bmp_unc = 0.05;
			if (ds.name == "7TeV") t_bmp_unc = 0.015;
			if (ds.name == "8TeV") t_bmp_unc = 0.02;
			if (ds.name == "13TeV") t_bmp_unc = 0.01;

			// find t_min
			const double dsdt_min = dsdt_bmp * (dsdt_bmp / dsdt_dip);
			double t_min = 0., diff_min = 1E100;
			for (double t = 0.3; t < ds.t2; t += dt)
			{
				const double diff = fabs(dsdt_min - f_fit->Eval(t));

				if (diff < diff_min)
				{
					t_min = t;
					diff_min = diff;
				}
			}

			// find t_max
			const double dsdt_max = 0.89 * dsdt_bmp;
			double t_max = 0.;
			diff_min = 1E100;
			for (double t = t_bmp; t < ds.t3 + 0.2; t += dt)
			{
				const double diff = fabs(dsdt_max - f_fit->Eval(t));

				if (diff < diff_min)
				{
					t_max = t;
					diff_min = diff;
				}
			}

			// manual correction
			if (ds.name == "2.76TeV")
			{
				t_bmp = 0.75; t_bmp_unc = 0.05;
				t_max = 0.8;
			}

			// print results
			printf("    min: t = %.3f\n", t_min);
			printf("    dip: t = %.3f +- %.3f, dsdt = %.3E\n", t_dip, t_dip_unc, dsdt_dip);
			printf("    bmp: t = %.3f +- %.3f, dsdt = %.3E\n", t_bmp, t_bmp_unc, dsdt_bmp);
			printf("    max: t = %.3f\n", t_max);

			// save results
			gDirectory = d_method->mkdir(ds.name.c_str());

			TGraph *g_data = new TGraph();
			g_data->SetPoint(0, t_min, 0.);
			g_data->SetPoint(1, dsdt_min, 0.);

			g_data->SetPoint(2, t_dip, t_dip_unc);
			g_data->SetPoint(3, dsdt_dip, 0.);

			g_data->SetPoint(4, t_bmp, t_bmp_unc);
			g_data->SetPoint(5, dsdt_bmp, 0.);

			g_data->SetPoint(6, t_max, 0.);
			g_data->SetPoint(7, dsdt_max, 0.);

			g_data->Write("g_data");

			// fill graphs
			int idx = g_t_min_vs_sqrt_s->GetN();
			g_t_min_vs_sqrt_s->SetPoint(idx, ds.sqrt_s, t_min);

			g_t_dip_vs_sqrt_s->SetPoint(idx, ds.sqrt_s, t_dip);
			g_t_dip_vs_sqrt_s->SetPointError(idx, 0., t_dip_unc);

			g_t_bmp_vs_sqrt_s->SetPoint(idx, ds.sqrt_s, t_bmp);
			g_t_bmp_vs_sqrt_s->SetPointError(idx, 0., t_bmp_unc);

			g_t_max_vs_sqrt_s->SetPoint(idx, ds.sqrt_s, t_max);

			// clean up
			delete f_in;
		}

		// fit
		TF1 *ff = new TF1("ff", "[0] + [1]*(x-[2])");
		ff->FixParameter(2, ds_ext.sqrt_s);

		g_t_min_vs_sqrt_s->Fit(ff);
		g_t_dip_vs_sqrt_s->Fit(ff);
		g_t_bmp_vs_sqrt_s->Fit(ff);
		g_t_max_vs_sqrt_s->Fit(ff);

		// save results
		gDirectory = d_method;

		g_t_min_vs_sqrt_s->Write("g_t_min_vs_sqrt_s");
		g_t_dip_vs_sqrt_s->Write("g_t_dip_vs_sqrt_s");
		g_t_bmp_vs_sqrt_s->Write("g_t_bmp_vs_sqrt_s");
		g_t_max_vs_sqrt_s->Write("g_t_max_vs_sqrt_s");
	}

	// clean up
	delete f_out;

	return 0;
}
