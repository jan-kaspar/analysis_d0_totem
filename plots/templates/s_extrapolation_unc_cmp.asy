import root;
import pad_layout;

include common_code;

string f = topDir + "s_extrapolation.root";

string base = "sqrt_s/corr";

TGraph_errorBar = None;

string quantities[];
quantities.push("dsdt");
quantities.push("B");

//----------------------------------------------------------------------------------------------------

for (int qi : quantities.keys)
{
	string q = quantities[qi];

	if (q == "dsdt")
	{
		NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$");
		scale(Linear, Log);
	}

	if (q == "B")
	{
		NewPad("$|t|\ung{GeV^2}$", "$B(t)\ung{GeV^{-2}}$");
	}

	AddToLegend("<TOTEM extrapolation:");

	draw(RootGetObject(f, base + "/g_"+q+"_ext"), "l", black, "central value");

	draw(RootGetObject(f, base + "/unc p/g_"+q+"_ext_pl_unc"), "l", red, "uncertainty -- parameters");
	draw(RootGetObject(f, base + "/unc p/g_"+q+"_ext_mi_unc"), "l", red);

	draw(RootGetObject(f, base + "/unc t/g_"+q+"_ext_pl_unc"), "l", blue, "uncertainty -- $t$-values");
	draw(RootGetObject(f, base + "/unc t/g_"+q+"_ext_mi_unc"), "l", blue);

	draw(RootGetObject(f, base + "/unc c/g_"+q+"_ext_pl_unc"), "l", heavygreen, "uncertainty -- combined");
	draw(RootGetObject(f, base + "/unc c/g_"+q+"_ext_mi_unc"), "l", heavygreen);

	if (q == "dsdt")
		limits((0.4, 4e-3), (0.9, 4e-1), Crop);
	if (q == "B")
		limits((0.4, -22.), (0.9, +8), Crop);
}

AttachLegend(NW, NE);
