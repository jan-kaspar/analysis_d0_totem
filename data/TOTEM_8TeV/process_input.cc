#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMatrixDSym.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// get input
	TFile *f_in = TFile::Open("/eos/user/j/jkaspar/work/analyses/elastic/4000GeV/beta90/high_t/DS-merged/merged.root");
	TH1D *h_in = (TH1D *) f_in->Get("ob-1-30-0.10/DS4-sc/combined/h_dsdt");

	// prepare output
	TFile *f_out = TFile::Open("data.root", "recreate");

	TGraphErrors *g_dsdt = new TGraphErrors();

	TGraphErrors *g_dsdt_syst_unc = new TGraphErrors();

	for (int bi = 1; bi <= h_in->GetNbinsX(); bi++)
	{
		const double t = h_in->GetBinCenter(bi);
		const double t_unc = h_in->GetBinWidth(bi) / 2.;
		const double dsdt = h_in->GetBinContent(bi);
		const double dsdt_unc = h_in->GetBinError(bi);

		int idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, t, dsdt);
		g_dsdt->SetPointError(idx, t_unc, dsdt_unc);

		g_dsdt_syst_unc->SetPoint(idx, 0., 0.);
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
	delete f_out;

	delete f_in;

	return 0;
}
