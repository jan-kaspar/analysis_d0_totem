#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMatrixDSym.h"
#include "TVectorD.h"


using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// get input data
	TFile *f_in = new TFile("publication1_graph.root");
	TGraphErrors *g = (TGraphErrors *) f_in->Get("g1");
	
	// systematic error data
	double syst_unc_t[] = { 0.4, 0.5, 1.5 };
	//double syst_unc_val[] = { +29.0, +31.6, +26.3 };
	double syst_unc_val[] = { +13.0, +13.0, +13.0 }; // temporarily use more reasonable numbers
	TGraph *g_rel_syst_unc_t_dep = new TGraph(3, syst_unc_t, syst_unc_val);

	// prepare output
	TFile *f_out = new TFile("data.root", "recreate");

	TGraphErrors *g_dsdt = new TGraphErrors();

	// process input
	for (int i = 0; i < g->GetN(); i++)
	{
		double x, x_e, y, y_stat_unc;
		g->GetPoint(i, x, y);
		x_e = g->GetErrorX(i);
		y_stat_unc = g->GetErrorY(i);

		if (x < 0.377)
			continue;

		int idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, x, y);
		g_dsdt->SetPointError(idx, x_e, y_stat_unc);
	}

	// make fit
	TF1 *ff = new TF1("ff", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)");
	ff->SetParameters(
		10.5, -37.7, 15.8, 0.,
		-33.9, 113., -139., 54.7
	);
	ff->FixParameter(3, 0.);
	ff->SetRange(0.25, 0.95);

	g_dsdt->Fit(ff, "", "", 0.25, 0.95);

	g_dsdt->Write("g_dsdt");
	ff->Write("f_dsdt");

	// make systematic-uncertainty objects
	const int dim = g_dsdt->GetN();

	TMatrixDSym m_dsdt_cov_syst_t_dep(dim);

	for (int i = 0; i < dim; ++i)
	{
		const double t = g_dsdt->GetX()[i];
		const double dsdt_ref = (t < 0.95) ? ff->Eval(t) : g_dsdt->Eval(t);
		const double unc = g_rel_syst_unc_t_dep->Eval(t) / 100. * dsdt_ref;

		m_dsdt_cov_syst_t_dep(i, i) = unc * unc; // for the moment uncorrelated approximation
	}

	m_dsdt_cov_syst_t_dep.Write("m_dsdt_cov_syst_t_dep");

	//const double rel_syst_t_indep = 0.11;
	const double rel_syst_t_indep = 0.27; // temporarily use more reasonable numbers

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
	delete f_in;

	delete f_out;

	return 0;
}
