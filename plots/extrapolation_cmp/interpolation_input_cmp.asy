import root;
import pad_layout;

include "../templates/common_code.asy";

string topDir = "../../";

string results[];
results.push("minimal/e01+e023:f=0.05/st+sy");
results.push("low_t,high_t/e012+e023:f=0.2/st+sy");

string extModel = "sqrt_s";
string extFit = "corr";

TGraph_errorBar = None;


//----------------------------------------------------------------------------------------------------

NewPad(false);

for (int dsi : datasets.keys)
	NewPadLabel(d_labels[dsi]);

//----------------------------------------------------------------------------------------------------

for (int ri : results.keys)
{
	NewRow();

	string l = "\vbox{\SetFontSizesXX";
	for (string e : split(results[ri], "/"))
		l += "\hbox{" + replace(e, "_", "\_") + "}";
	l += "}";
	NewPadLabel(l);

	for (int dsi : datasets.keys)
	{
		NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$");
		scale(Linear, Log);

		string f = topDir + "fits/" + results[ri] + "/do_fits.root";
		draw(RootGetObject(f, datasets[dsi] + "/g_fit"), "l", blue, "$\d\si/\d t$ fit");
		draw(RootGetObject(f, datasets[dsi] + "/g_fit_pl_unc"), "l", blue+dashed);
		draw(RootGetObject(f, datasets[dsi] + "/g_fit_mi_unc"), "l", blue+dashed);

		string f = topDir + "fits/" + results[ri] + "/s_extrapolation.root";
		string base = extModel + "/" + extFit + "/" + datasets[dsi];
		draw(RootGetObject(f, base + "/g_dsdt_ext"), "l", red, "parameter $s$-evolution");

		limits((0.3, 4e-3), (1.0, 5e-1), Crop);
	}

	AttachLegend(BuildLegend(NW), NE);
}

GShipout(hSkip=3mm, vSkip=1mm);
