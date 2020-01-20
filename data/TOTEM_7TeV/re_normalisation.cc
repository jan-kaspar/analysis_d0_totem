#include "TFile.h"
#include "TH1D.h"
#include "TGraphErrors.h"
#include "TF1.h"

using namespace std;

int main()
{
	// get original input
	TFile *f_in_3p5 = TFile::Open("publication_3.5m.root");
	TGraphErrors *g_3p5 = (TGraphErrors *) f_in_3p5->Get("g1");

	TFile *f_in_90 = TFile::Open("publication_90m.root");
	TH1D *h_90 = (TH1D *) f_in_90->Get("h_avg");

	TFile *f_out = TFile::Open("re_normalisation.root", "recreate");

	TF1 *ff = new TF1("ff", "exp([0] + [1]*x)");

	printf("%p\n", ff);

	ff->SetParameters(20., -22.);
	g_3p5->Fit(ff, "", "", 0.36, 0.40);
	g_3p5->Write("g_3p5");

	const double a_3p5 = ff->GetParameter(0);
	const double b_3p5 = ff->GetParameter(1);

	ff->SetParameters(20., -22.);
	h_90->Fit(ff, "", "", 0.34, 0.43);
	h_90->Write("h_90");

	const double a_90 = ff->GetParameter(0);
	const double b_90 = ff->GetParameter(1);

	const double t1 = 0.36, t2 = 0.40;

	const double I_3p5 = exp(a_3p5) / b_3p5 * (exp(b_3p5 * t2) - exp(b_3p5 * t1));
	const double I_90 = exp(a_90) / b_90 * (exp(b_90 * t2) - exp(b_90 * t1));
	printf("I_3p5 = %.3E, I_90 = %.3E, I_90 / I_3p5 = %.3f\n", I_3p5, I_90, I_90 / I_3p5);
	
	delete f_out;

	return 0;
}
