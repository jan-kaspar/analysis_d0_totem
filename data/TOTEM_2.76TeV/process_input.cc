#include "TFile.h"
#include "TGraphErrors.h"
#include "TMatrixDSym.h"

#include <cstdio>

#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// get input
	FILE *f_in = fopen("final_TOTEM_2p76TeV_13sigma_elastic_cross_section.txt", "r");

	// prepare output
	TFile *f_out = TFile::Open("data.root", "recreate");

	TGraphErrors *g_dsdt = new TGraphErrors();

	TGraphErrors *g_dsdt_syst_unc = new TGraphErrors();

	// process input
	while (!feof(f_in))
	{
		float t_min, t_max, t_cen, dsdt, dsdt_unc;
		int r = fscanf(f_in, "%f %f %f %f %f", &t_min, &t_max, &t_cen, &dsdt, &dsdt_unc);

		if (r == 5)
		{
			double t_unc = (t_max - t_min) / 2.;

			int idx = g_dsdt->GetN();
			g_dsdt->SetPoint(idx, t_cen, dsdt);
			g_dsdt->SetPointError(idx, t_unc, dsdt_unc);

			g_dsdt_syst_unc->SetPoint(idx, 0., 0.);
		}
	}

	g_dsdt->Write("g_dsdt");

	// make systematic-uncertainty matrices
	const int dim = g_dsdt_syst_unc->GetN();

	TMatrixDSym m_dsdt_cov_syst_full(dim), m_dsdt_cov_syst_no_norm(dim);

	// TODO: improve, for the moment just an uncorrelated simplification
	for (int i = 0; i < dim; ++i)
	{
		const double unc = g_dsdt_syst_unc->GetY()[i];

		m_dsdt_cov_syst_no_norm(i, i) = 0.;
		m_dsdt_cov_syst_full(i, i) = unc * unc;
	}

	m_dsdt_cov_syst_full.Write("m_dsdt_cov_syst_full");
	m_dsdt_cov_syst_no_norm.Write("m_dsdt_cov_syst_no_norm");

	// clean up
	fclose(f_in);

	delete f_out;

	return 0;
}
