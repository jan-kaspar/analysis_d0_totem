import root;
import pad_layout;

include common_code;

TGraph_errorBar = None;

xSizeDef = 7cm;

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
	NewPadLabel(d_labels[dsi]);

//----------------------------------------------------------------------------------------------------

NewRow();

for (int dsi : datasets.keys)
{
	RootObject g_data = RootGetObject(f, datasets[dsi] + "/g_data");

	RootObject g_dsdt = RootGetObject(f, datasets[dsi] + "/g_dsdt");

	RootObject g_fit = RootGetObject(f, datasets[dsi] + "/g_fit");

	RootObject g_fit_c1 = RootGetObject(f, datasets[dsi] + "/components/g_comp1");
	RootObject g_fit_c2 = RootGetObject(f, datasets[dsi] + "/components/g_comp2");

	real ax[] = {0.};
	real ay[] = {0.};
	g_data.vExec("GetPoint", 0, ax, ay); real t_min = ax[0], t_max = ay[0];
	g_data.vExec("GetPoint", 1, ax, ay); real n_fit_points = ax[0], ndf = ay[0];
	g_data.vExec("GetPoint", 2, ax, ay); real chi2 = ax[0], chi2_ndf = ay[0];

	NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
	scale(Linear, Log);

	TGraph_x_min = t_min;
	TGraph_x_max = t_max;

	draw(g_fit, "l", red+1pt);

	draw(g_dsdt, "p", black, mCi+0.5pt);

	draw(g_fit_c1, "l", blue);
	draw(g_fit_c2, "l", heavygreen);

	limits((0.3, 4e-3), (0.9, 4e-1), Crop);

	//AttachLegend();
}

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=1mm, vSkip=1mm);
