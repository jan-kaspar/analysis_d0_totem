#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TMatrixDSym.h"

#include <map>

using namespace std;

//----------------------------------------------------------------------------------------------------

int main()
{
	// get input data
	TFile *f_in = new TFile("publication1_graph.root");
	TGraphErrors *g = (TGraphErrors *) f_in->Get("g1");
	
	//TFile *f_in_h = new TFile("publication1.root");
	//TH1D *h = (TH1D *) f_in_h->Get("h_dsdt");

	// systematic error data
	double syst_err_t[] = { 0.4, 0.5, 1.5 };
	double syst_err_up[] = { +25., +28., +27. };
	double syst_err_down[] = { -37., -39., -30. };
	TGraph *g_rel_syst_err_up = new TGraph(3, syst_err_t, syst_err_up);
	TGraph *g_rel_syst_err_down = new TGraph(3, syst_err_t, syst_err_down);

	// prepare output
	TFile *f_out = new TFile("data.root", "recreate");

	TGraphErrors *g_dsdt = new TGraphErrors();

	TGraphErrors *g_dsdt_syst_unc = new TGraphErrors();

	// process input
	for (int i = 0; i < g->GetN(); i++)
	{
		double x, x_e, y, y_stat_e;
		g->GetPoint(i, x, y);
		x_e = g->GetErrorX(i);
		y_stat_e = g->GetErrorY(i);

		if (x < 0.377)
			continue;

		const double y_syst_up = g_rel_syst_err_up->Eval(x) / 100. * y;
		const double y_syst_down = -g_rel_syst_err_down->Eval(x) / 100. * y;
		const double y_syst = (y_syst_up + y_syst_down) / 2.;

		//const int h_idx = h->FindBin(x);
		//const double le = h->GetBinLowEdge(h_idx);
		//const double he = le + h->GetBinWidth(h_idx);

		/*
		printf("%f, %f, %f, %f | %.3f, %.3f, %.3f, %.3f\n",
			le, he, x, x_e,
			y*1E3, y_stat_e*1E3, y_syst_up*1E3, y_syst_down*1E3);
		*/

		int idx = g_dsdt->GetN();
		g_dsdt->SetPoint(idx, x, y);
		g_dsdt->SetPointError(idx, x_e, y_stat_e);

		g_dsdt_syst_unc->SetPoint(idx, 0., y_syst);
	}
	
	g_dsdt->Write("g_dsdt");

	// make systematic-uncertainty matrices
	const int dim = g_dsdt_syst_unc->GetN();

	TMatrixDSym m_dsdt_cov_syst_full(dim), m_dsdt_cov_syst_no_norm(dim);

	// TODO: improve, for the moment just an fixed-correlation approximation
	for (int i = 0; i < dim; ++i)
	{
		for (int j = 0; j < dim; ++j)
		{
			const double unc_i = g_dsdt_syst_unc->GetY()[i];
			const double unc_j = g_dsdt_syst_unc->GetY()[j];

			m_dsdt_cov_syst_no_norm(i, i) = 0.;
			m_dsdt_cov_syst_full(i, j) = (i == j) ? unc_i * unc_i : 0.01 * unc_i * unc_j;
		}
	}

	m_dsdt_cov_syst_full.Write("m_dsdt_cov_syst_full");
	m_dsdt_cov_syst_no_norm.Write("m_dsdt_cov_syst_no_norm");

	// clean up
	delete f_in;
	//delete f_in_h;
	delete f_out;
	return 0;
}
