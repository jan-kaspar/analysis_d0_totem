import root;
import pad_layout;

include common_code;

string f = topDir + "s_extrapolation.root";

//----------------------------------------------------------------------------------------------------

RootObject g_settings = RootGetObject(f, "g_settings");

real ax[] = {0.};
real ay[] = {0.};
g_settings.vExec("GetPoint", 0, ax, ay); int n_parameters = (int)ax[0];

//----------------------------------------------------------------------------------------------------

NewPad(false);

for (int emi : extModels.keys)
	for (int efi : extFits.keys)
		NewPadLabel(replace(extFits[efi] + ", " + extModels[emi], "_", "\_"));

//----------------------------------------------------------------------------------------------------

for (int pari = 0; pari < n_parameters; ++pari)
{
	NewRow();

	NewPadLabel(format("parameter %i", pari));

	for (int emi : extModels.keys)
	{
		for (int efi : extFits.keys)
		{
			NewPad("$\sqrt s$", "");
			
			if (find(extModels[emi], "log") >= 0)
				scale(Log, Linear);

			if (extModels[emi] == "sqrt_s")
				currentpad.xTicks = LeftTicks(2., 1.);

			string base = extModels[emi] + "/" + extFits[efi] + "/" + format("c_par%i", pari);

			RootObject g_cnt = RootGetObject(f, base + "|g_cnt");
			RootObject g_cnt_orig = RootGetObject(f, base + "|g_cnt_orig");
			RootObject g_fit = RootGetObject(f, base + "|g_fit");
			RootObject g_fit_pl_unc = RootGetObject(f, base + "|g_fit_pl_unc");
			RootObject g_fit_mi_unc = RootGetObject(f, base + "|g_fit_mi_unc");

			TGraph_errorBar = Bars;
			draw(g_cnt_orig, "p", black, mSq+false+4pt+black);

			TGraph_errorBar = None;
			draw(g_cnt, "p", blue, mCi+2pt+blue);

			draw(g_fit, "l", red);
			draw(g_fit_pl_unc, "l", red+dashed);
			draw(g_fit_mi_unc, "l", red+dashed);

			/*
			if (pari == 5)
				limits((1., -50.), (14., 0.), Crop);
			else
				xlimits(1., 14., Crop);
			*/

			/*
			if (pari == 4)
				currentpad.yTicks = RightTicks(10., 2.);
			*/

			//AttachLegend(format("%.3f", g_fit.rExec("Eval", 2.76)));

			yaxis(XEquals(sqrt_s_ext, false), heavygreen);
		}
	}
}

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=1mm, vSkip=1mm);
