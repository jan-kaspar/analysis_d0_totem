#include "TFile.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

#include <cstdio>

#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

void ProcessOne(const string &inputFile, const string &outputFile, bool highT, bool attract)
{
	// get input
	FILE *f_in = fopen(inputFile.c_str(), "r");

	// prepare output
	TFile *f_out = TFile::Open(outputFile.c_str(), "recreate");

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
			if (attract && fabs(t_repr - 0.61) < 0.01)
			{
				dsdt_stat_unc /= 5.;
				dsdt_syst_unc_t_dep /= 5.;
			}

			const double t_unc = (t_max - t_min) / 2.;

			int idx = g_dsdt->GetN();
			g_dsdt->SetPoint(idx, t_repr, dsdt);
			g_dsdt->SetPointError(idx, t_unc, dsdt_stat_unc);

			g_dsdt_syst_unc_t_dep->SetPoint(idx, 0., dsdt_syst_unc_t_dep);

			//printf("%i, %f, %f\n", idx, t_repr, dsdt);
		}
	}

	// additional steering points
	/*
	{
		int idx;

	   	idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, 0.61, 0.85E-2);
		g_dsdt->SetPointError(idx, 0.01, 0.05E-2);
		g_dsdt_syst_unc_t_dep->SetPoint(idx, 0., 0.05E-2);

	   	idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, 0.70, 1.30E-2);
		g_dsdt->SetPointError(idx, 0.01, 0.05E-2);
		g_dsdt_syst_unc_t_dep->SetPoint(idx, 0., 0.05E-2);

	   	idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, 0.78, 1.55E-2);
		g_dsdt->SetPointError(idx, 0.01, 0.05E-2);
		g_dsdt_syst_unc_t_dep->SetPoint(idx, 0., 0.05E-2);

	   	idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, 0.90, 1.15E-2);
		g_dsdt->SetPointError(idx, 0.01, 0.05E-2);
		g_dsdt_syst_unc_t_dep->SetPoint(idx, 0., 0.05E-2);

		idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, 0.98, 0.85E-2);
		g_dsdt->SetPointError(idx, 0.01, 0.05E-2);
		g_dsdt_syst_unc_t_dep->SetPoint(idx, 0., 0.05E-2);
	}
	*/

	// make fit
	TF1 *ff = NULL;

	if (highT)
	{
		ff = new TF1("ff", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)");
		ff->SetParameters(
			5.14, -11.89, -8.326, 0.,
			-79.00, 264.0, -306.5, 116.7
		);

		ff->FixParameter(3, 0.);

		ff->FixParameter(4, -79.00);
		ff->FixParameter(5, 264.0);
		ff->FixParameter(6, -306.5);
		ff->FixParameter(7, 116.7);
		ff->SetRange(0.35, 1.0);

		g_dsdt->Fit(ff, "", "", 0.35, 1.0);
	} else {
		ff = new TF1("ff", "exp([0] + [1]*x + [2]*x*x)");
		ff->SetParameters(
			6.75, -19.3, 0.
		);
		//ff->FixParameter(2, 0.);
		ff->SetRange(0.06, 0.47);

		g_dsdt->Fit(ff, "", "", 0.06, 0.47);
	}

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

	const double rel_syst_t_indep = 0.060;

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
	ProcessOne("from_preprint_4sigma.txt", "data_4sigma.root", false, false);

	ProcessOne("from_preprint_13sigma.txt", "data.root", true, false);

	ProcessOne("from_preprint_13sigma.txt", "data_att.root", true, true);

	return 0;
}
