import root;
import pad_layout;

include "../dip/common_code.asy";

string topDir = "../../";

TGraph_errorBar = None;

xSizeDef = 7cm;

//----------------------------------------------------------------------------------------------------

void DrawBand(RootObject fit, RootObject unc, pen p)
{
	guide g_u, g_l;

	for (int i = 0; i < unc.iExec("GetN"); ++i)
	{
		real ax[] = {0.};
		real ay[] = {0.};
		unc.vExec("GetPoint", i, ax, ay);

		real t = ax[0];

		if (t < 0.3 || t > 0.9)
			continue;

		real dsdt_ref = fit.rExec("Eval", t);
		real dsdt_unc = ay[0];

		//write(t, dsdt_ref);

		g_u = g_u -- Scale((t, dsdt_ref + dsdt_unc));
		g_l = g_l -- Scale((t, dsdt_ref - dsdt_unc));
	}

	filldraw(g_u -- reverse(g_l) -- cycle, p, nullpen);
}

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
	NewPadLabel(d_labels[dsi]);

//----------------------------------------------------------------------------------------------------

NewRow();

for (int dsi : datasets.keys)
{
	write("* " + datasets[dsi]);

	string f = topDir + "data/TOTEM_" + datasets[dsi] + "/data.root";

	RootObject g_dsdt = RootGetObject(f, "g_dsdt");
	RootObject f_dsdt = RootGetObject(f, "f_dsdt");
	RootObject g_dsdt_syst_t_dep = RootGetObject(f, "g_dsdt_syst_t_dep");
	RootObject g_dsdt_syst_full = RootGetObject(f, "g_dsdt_syst_full");

	NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
	scale(Linear, Log);

	DrawBand(f_dsdt, g_dsdt_syst_full, yellow);
	DrawBand(f_dsdt, g_dsdt_syst_t_dep, heavygreen);

	draw(f_dsdt, "l", black);

	draw(g_dsdt, "p", red, mCi+0.5pt+red);

	limits((0.3, 4e-3), (0.9, 1e0), Crop);

	if (dsi == 0)
	{
		AddToLegend("data with stat.~unc.", mPl+red+4pt);
		AddToLegend("$t$-dependent syst.~unc.", mSq+heavygreen+4pt);
		AddToLegend("full syst.~unc.", mSq+yellow+4pt);
		AttachLegend();
	}
}

GShipout(hSkip=1mm);
