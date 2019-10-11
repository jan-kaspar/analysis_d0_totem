#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include <cstdio>

#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// get input
	FILE *f_in = fopen("dsigma_dt_13TeV_90m_10sigma_complete_01062018.txt", "r");

	// prepare output
	TFile *f_out = TFile::Open("data.root", "recreate");

	TGraphErrors *g_dsdt = new TGraphErrors();

	TGraphErrors *g_dsdt_syst_unc_t_dep = new TGraphErrors();

	// process input
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
			const double t_unc = (t_max - t_min) / 2.;

			int idx = g_dsdt->GetN();
			g_dsdt->SetPoint(idx, t_repr, dsdt);
			g_dsdt->SetPointError(idx, t_unc, dsdt_stat_unc);

			g_dsdt_syst_unc_t_dep->SetPoint(idx, 0., dsdt_syst_unc_t_dep);
		}
	}

	// make fit
	TF1 *ff = new TF1("ff", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)");
	ff->SetParameters(
		4.57, -5.24, -31.3, 0.,
		-25.2, 89.4, -117., 49.2
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

	return 0;
}
