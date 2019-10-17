import root;
import pad_layout;

string topDir = "../../";

string dir = "fits/minimal/e01+e023/";
//string dir = "fits/low_t,high_t/e012+e023/";

int n_parameters = 5;

string uncSpecs[] = {
	"st",
	"st+sy",
	"st+sy+no",
};

string extModel = "sqrt_s";
string extFit = "corr";

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

for (int pari = 0; pari < n_parameters; ++pari)
{
	NewRow();

	NewPadLabel(format("parameter %i", pari));
	
	NewPad("$\sqrt s$", "");

	for (int usi : uncSpecs.keys)
	{
		string f = topDir + dir + uncSpecs[usi] + "/s_extrapolation.root";
		string base = extModel + "/" + extFit + "/" + format("c_par%i", pari);

		RootObject g_cnt = RootGetObject(f, base + "|g_cnt");
		RootObject g_fit = RootGetObject(f, base + "|g_fit");
		RootObject g_fit_pl_unc = RootGetObject(f, base + "|g_fit_pl_unc");
		RootObject g_fit_mi_unc = RootGetObject(f, base + "|g_fit_mi_unc");

		pen p = StdPen(usi+1);

		real fsh = 0., fsh_a = 0.2;
		if (uncSpecs.length > 1)
			fsh = fsh_a * (2*usi / (uncSpecs.length-1) - 1);

		draw(shift(fsh, 0), g_cnt, "p", p, mCi+2pt+p);
		//draw(g_fit, "l", blue);
		//draw(g_fit_pl_unc, "l", blue+dashed);
		//draw(g_fit_mi_unc, "l", blue+dashed);

		//yaxis(XEquals(sqrt_s_ext, false), heavygreen);
	}

	if (pari == 0)
	{
		NewPad(false);
		for (int usi : uncSpecs.keys)
			AddToLegend(uncSpecs[usi], StdPen(usi+1));
		AttachLegend();
	}
}

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=1mm, vSkip=1mm);
