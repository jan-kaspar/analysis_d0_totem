import root;
import pad_layout;

include "common_code.asy";

string topDir = "../";

string results[];
results.push(".");

string extModel = "sqrt_s";
string extFit = "corr";
TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

for (int ri : results.keys)
{
	NewRow();

	NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$");
	scale(Linear, Log);

	AddToLegend("<predictions of $s$-dependent model");
	
	string f = topDir + results[ri] + "/s_extrapolation.root";

	for (int dsi : datasets.keys)
	{
		string base = extModel + "/" + extFit + "/" + datasets[dsi];
		pen p = StdPen(dsi + 1);
		draw(RootGetObject(f, base + "/g_dsdt_ext"), "l", p, d_labels[dsi]);
	}

	string base = extModel + "/" + extFit;
	draw(RootGetObject(f, base + "/g_dsdt_ext"), "l", black, "D0");
	draw(RootGetObject(f, base + "/unc c/g_dsdt_ext_pl_unc"), "l", black+dashed);
	draw(RootGetObject(f, base + "/unc c/g_dsdt_ext_mi_unc"), "l", black+dashed);
	
	limits((0.3, 4e-3), (1.0, 5e-1), Crop);

	AttachLegend(BuildLegend(NW), NE);
}

GShipout(hSkip=3mm, vSkip=1mm);
