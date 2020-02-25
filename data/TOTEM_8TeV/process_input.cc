#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

void ProcessOne(const string &input, const string &output)
{
	printf("\n----------------------------------------------------------\n");
	printf("%s\n\n", output.c_str());

	// get input
	TFile *f_in = TFile::Open("/afs/cern.ch/work/j/jkaspar/work/analyses/elastic/4000GeV/beta90/high_t/DS-merged/merged.root");
	TH1D *h_in = (TH1D *) f_in->Get(input.c_str());

	// prepare output
	TFile *f_out = TFile::Open(output.c_str(), "recreate");

	TGraphErrors *g_dsdt = new TGraphErrors();

	// process input
	for (int bi = 1; bi <= h_in->GetNbinsX(); bi++)
	{
		const double t = h_in->GetBinCenter(bi);
		const double t_unc = h_in->GetBinWidth(bi) / 2.;
		const double dsdt = h_in->GetBinContent(bi);
		const double dsdt_stat_unc = h_in->GetBinError(bi);

		if (t < 0.03 || t > 1.8)
			continue;

		int idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, t, dsdt);
		g_dsdt->SetPointError(idx, t_unc, dsdt_stat_unc);

	}

	// make fit
	TF1 *ff = new TF1("ff", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)");
	ff->SetParameters(
		//3.39, 0.475, -34.5, 0.
		//-22.1, 65.4, -73.6, 25.8

		3.85, -2.58, -29.8, 0.,
		-25.6, 80.5, -93.8, 34.6
	);
	ff->FixParameter(3, 0.);
	ff->SetRange(0.25, 1.1);

	//ff->FixParameter(7, 25.8);
	g_dsdt->Fit(ff, "", "", 0.25, 1.1);

	g_dsdt->Write("g_dsdt");
	ff->Write("f_dsdt");

	// make systematic-uncertainty objects
	const int dim = g_dsdt->GetN();

	TMatrixDSym m_dsdt_cov_syst_t_dep(dim);

	for (int i = 0; i < dim; ++i)
	{
		const double t = g_dsdt->GetX()[i];
		const double dsdt_ref = (t < 1.01) ? ff->Eval(t) : g_dsdt->Eval(t);
		const double unc = 0.03 * dsdt_ref; // just a crude approximation

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

		v_dsdt_syst_t_indep(i) = (t < 1.01) ? rel_syst_t_indep * ff->Eval(t) : 0.;

		g_dsdt_syst_t_dep->SetPoint(i, t, sqrt(m_dsdt_cov_syst_t_dep(i, i)));
		g_dsdt_syst_t_indep->SetPoint(i, t, v_dsdt_syst_t_indep(i));
		g_dsdt_syst_full->SetPoint(i, t, sqrt(m_dsdt_cov_syst_t_dep(i, i) + pow(v_dsdt_syst_t_indep(i), 2)));
	}

	v_dsdt_syst_t_indep.Write("v_dsdt_syst_t_indep");

	g_dsdt_syst_t_dep->Write("g_dsdt_syst_t_dep");
	g_dsdt_syst_t_indep->Write("g_dsdt_syst_t_indep");
	g_dsdt_syst_full->Write("g_dsdt_syst_full");

	// clean up
	delete f_in;

	delete f_out;
}

//----------------------------------------------------------------------------------------------------

int main()
{
	ProcessOne("ob-1-30-0.10/DS4-sc/combined/h_dsdt", "data.root");

	ProcessOne("ob-1-30-0.10/DS4-sc/combined/h_dsdt", "data_b1.root");
	ProcessOne("ob-2-20-0.20/DS4-sc/combined/h_dsdt", "data_b2.root");
	ProcessOne("ob-3-10-0.30/DS4-sc/combined/h_dsdt", "data_b3.root");
	ProcessOne("ob-1-30-0.10-S/DS4-sc/combined/h_dsdt", "data_b1s.root");
	ProcessOne("ob-2-20-0.20-S/DS4-sc/combined/h_dsdt", "data_b2s.root");
	ProcessOne("bt1/DS4-sc/combined/h_dsdt", "data_bt1.root");

	ProcessOne("ob-1-30-0.10/DS4-sc/45b_56t/h_dsdt", "data_b1_45b_56t.root");
	ProcessOne("ob-1-30-0.10/DS4-sc/45t_56b/h_dsdt", "data_b1_45t_56b.root");

	return 0;
}
