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

for (int ri : results.keys)
{
	NewRow();

	NewPadLabel(replace(results[ri], "_", "\_"));

	for (int mode = 0; mode < 2; ++mode)
	{
		NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$");
		scale(Linear, Log);

		// draw input
		string f = topDir + "fits/" + results[ri] + "/do_fits.root";
		for (int dsi : datasets.keys)
		{
			pen p = StdPen(dsi + 1);

			string obj_pth = datasets[dsi] + "/" + ((mode == 0) ? "g_dsdt" : "g_fit");

			if (mode == 0)
				draw(RootGetObject(f, obj_pth), "p", p, mCi+2pt+p);

			if (mode == 1)
				draw(RootGetObject(f, obj_pth), "l", p, d_labels[dsi]);
		}

		// draw extrapolation
		string f = topDir + "fits/" + results[ri] + "/s_extrapolation.root";
		string base = extModel + "/" + extFit;
		draw(RootGetObject(f, base + "/g_dsdt_ext"), "l", black, "D0 energy");

		limits((0.3, 4e-3), (0.9, 4e-1), Crop);

		if (mode == 1)
			AttachLegend(NW, NE);
	}
}

GShipout(hSkip=3mm, vSkip=1mm);
