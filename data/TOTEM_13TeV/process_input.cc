#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include <cstdio>

#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

void LoadData(FILE *f_in, double addRelSystUnc, bool attract, const string &rebin, TGraphErrors *g_dsdt, TGraph *g_dsdt_syst_unc_t_dep)
{
	// load points
	struct Point
	{
		double t_min, t_max, dsdt, dsdt_stat_unc, dsdt_syst_unc_t_dep;
	};

	vector<Point> points;

	while (!feof(f_in))
	{
		char line[500];
		char *ret = fgets(line, 500, f_in);

		if (ret == NULL)
			break;

		if (line[0] == '#')
			continue;

		float t_min, t_max, t_repr, dsdt, dsdt_stat_unc, dsdt_syst_unc_t_dep;
		int r = sscanf(line, "%f %f %f %f %f %f", &t_min, &t_max, &t_repr, &dsdt, &dsdt_stat_unc, &dsdt_syst_unc_t_dep);

		if (r == 6)
		{
			// to account for the visible oscillations at low |t|
			const double ac = dsdt * addRelSystUnc;
			dsdt_syst_unc_t_dep = sqrt(dsdt_syst_unc_t_dep*dsdt_syst_unc_t_dep + ac*ac);

			points.push_back({t_min, t_max, dsdt, dsdt_stat_unc, dsdt_syst_unc_t_dep});
		}
	}

	// process input
	for (unsigned int pi = 0; pi < points.size();)
	{
		// determine how many points to merge
		unsigned int n_merge = 1;

		if (rebin == "rebinA")
		{
			const auto &p = points[pi];

			if (p.t_min < 0.42)
				n_merge = 2;
			if (p.t_min < 0.32)
				n_merge = 3;
		}

		if (rebin == "rebinB")
		{
			const auto &p = points[pi];

			if (p.t_min < 0.43)
				n_merge = 2;
			if (p.t_min < 0.41)
				n_merge = 3;
			if (p.t_min < 0.40)
				n_merge = 4;
		}

		// calculate merged quantities
		double t_min = +1E100, t_max = -1E100;

		double s_v_r = 0., s_stu2_r2 = 0., s_syu_r = 0.;

		for (unsigned int mi = 0; mi < n_merge; ++mi)
		{
			const auto &p = points[pi + mi];

			t_min = min(t_min, p.t_min);
			t_max = max(t_max, p.t_max);

			const double r = p.t_max - p.t_min;

			s_v_r += p.dsdt * r;
			s_stu2_r2 += pow(p.dsdt_stat_unc * r, 2);
			s_syu_r += p.dsdt_syst_unc_t_dep * r;
		}

		const double R = t_max - t_min;
		const double t_unc = R / 2.;
		const double t_mean = (t_max + t_min) / 2.;
		const double dsdt = s_v_r / R;
		const double dsdt_stat_unc = sqrt(s_stu2_r2) / R;
		const double dsdt_syst_unc_t_dep = s_syu_r / R;

		double f = 1.;

		if (attract && fabs(t_mean - 0.468) < 0.001)
			f = 0.1;

		// store merged point
		int idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, t_mean, dsdt);
		g_dsdt->SetPointError(idx, t_unc, f * dsdt_stat_unc);

		g_dsdt_syst_unc_t_dep->SetPoint(idx, 0., f * dsdt_syst_unc_t_dep);

		//printf("%.3f, %.3f: %.3f, %.3f, %.3f\n", t_min, t_max, dsdt, dsdt_stat_unc, dsdt_syst_unc_t_dep);

		// advance
		pi += n_merge;
	}
}

//----------------------------------------------------------------------------------------------------

void ProcessOne(const string &inputFile, const string &outputFile, double addRelSystUnc, bool attract, const string &rebin = "")
{
	printf("--------------------------------------------------------------------\n");
	printf("%s\n\n", outputFile.c_str());

	// prepare output
	TFile *f_out = TFile::Open(outputFile.c_str(), "recreate");

	// get input
	FILE *f_in = fopen(inputFile.c_str(), "r");

	TGraphErrors *g_dsdt = new TGraphErrors();
	TGraphErrors *g_dsdt_syst_unc_t_dep = new TGraphErrors();

	LoadData(f_in, addRelSystUnc, attract, rebin, g_dsdt, g_dsdt_syst_unc_t_dep);

	// make fit
	TF1 *ff = new TF1("ff", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)");
	ff->SetParameters(
		14.6, -100.0, 262.6, -299.4,
		-20.56, 69.92, -90.0, 36.79
	);
	ff->SetRange(0.25, 0.9);

	g_dsdt->Fit(ff, "", "", 0.25, 0.9);

	g_dsdt->Write("g_dsdt");
	ff->Write("f_dsdt");

	// make systematic-uncertainty objects
	const int dim = g_dsdt_syst_unc_t_dep->GetN();

	TMatrixDSym m_dsdt_cov_syst_t_dep(dim);

	for (int i = 0; i < dim; ++i)
	{
		const double unc = g_dsdt_syst_unc_t_dep->GetY()[i];

		m_dsdt_cov_syst_t_dep(i, i) = unc * unc; // for the moment uncorrelated approximation
	}

	m_dsdt_cov_syst_t_dep.Write("m_dsdt_cov_syst_t_dep");

	const double rel_syst_t_indep = 0.055;

	TVectorD v_rel_syst_t_indep(1);
	v_rel_syst_t_indep(0) = rel_syst_t_indep;
	v_rel_syst_t_indep.Write("v_rel_syst_t_indep");

	TVectorD v_dsdt_syst_t_indep(dim);

	TGraph *g_dsdt_syst_t_dep = new TGraph();
	TGraph *g_dsdt_syst_t_indep = new TGraph();
	TGraph *g_dsdt_syst_full = new TGraph();

	for (int i = 0; i < dim; ++i)
	{
		double t = g_dsdt->GetX()[i];

		v_dsdt_syst_t_indep(i) = (t < 1.0) ? rel_syst_t_indep * ff->Eval(t) : 0.;

		g_dsdt_syst_t_dep->SetPoint(i, t, sqrt(m_dsdt_cov_syst_t_dep(i, i)));
		g_dsdt_syst_t_indep->SetPoint(i, t, v_dsdt_syst_t_indep(i));
		g_dsdt_syst_full->SetPoint(i, t, sqrt(m_dsdt_cov_syst_t_dep(i, i) + pow(v_dsdt_syst_t_indep(i), 2)));
	}

	v_dsdt_syst_t_indep.Write("v_dsdt_syst_t_indep");

	g_dsdt_syst_t_dep->Write("g_dsdt_syst_t_dep");
	g_dsdt_syst_t_indep->Write("g_dsdt_syst_t_indep");
	g_dsdt_syst_full->Write("g_dsdt_syst_full");

	// clean up
	fclose(f_in);

	delete f_out;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	ProcessOne("dsigma_dt_13TeV_90m_10sigma_complete_01062018.txt", "data.root", 0.000, false);

	//ProcessOne("dsigma_dt_13TeV_90m_10sigma_complete_01062018.txt", "data_addUnc0.5.root", 0.005, false);
	//ProcessOne("dsigma_dt_13TeV_90m_10sigma_complete_01062018.txt", "data_addUnc0.8.root", 0.008, false);
	//ProcessOne("dsigma_dt_13TeV_90m_10sigma_complete_01062018.txt", "data_addUnc1.0.root", 0.010, false);
	//ProcessOne("dsigma_dt_13TeV_90m_10sigma_complete_01062018.txt", "data_addUnc1.5.root", 0.015, false);

	//ProcessOne("dsigma_dt_13TeV_90m_10sigma_complete_01062018.txt", "data_addUnc0.8_att.root", 0.008, true);

	ProcessOne("dsigma_dt_13TeV_90m_10sigma_complete_01062018.txt", "data_addUnc0.8_rebinB.root", 0.008, false, "rebinB");

	//ProcessOne("dsigma_dt_13TeV_90m_10sigma_complete_01062018.txt", "data_addUnc0.8_rebinA_att.root", 0.008, true, "rebinA");
	ProcessOne("dsigma_dt_13TeV_90m_10sigma_complete_01062018.txt", "data_addUnc0.8_rebinB_att.root", 0.008, true, "rebinB");

	return 0;
}
