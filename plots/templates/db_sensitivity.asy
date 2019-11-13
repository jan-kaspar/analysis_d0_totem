import root;
import pad_layout;

include common_code;

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

RootObject g_settings = RootGetObject(f, "g_settings");

real ax[] = {0.};
real ay[] = {0.};
g_settings.vExec("GetPoint", 0, ax, ay); int n_parameters = (int)ax[0];

//----------------------------------------------------------------------------------------------------

NewPad(false);

for (int dsi : datasets.keys)
	NewPadLabel(d_labels[dsi]);

//----------------------------------------------------------------------------------------------------
NewRow();

for (int pari = 0; pari < n_parameters; ++pari)
{
	write("* parameter ", pari);

	NewRow();

	NewPadLabel(format("parameter %i", pari));

	for (int dsi : datasets.keys)
	{
		RootObject g_data = RootGetObject(f, datasets[dsi] + "/g_data");
		RootObject g_fit = RootGetObject(f, datasets[dsi] + "/g_fit");

		RootObject g_mi = RootGetObject(f, datasets[dsi] + format("/sensitivity/%i_mi", pari));
		RootObject g_pl = RootGetObject(f, datasets[dsi] + format("/sensitivity/%i_pl", pari));

		real ax[] = {0.};
		real ay[] = {0.};
		g_data.vExec("GetPoint", 0, ax, ay); real t_min = ax[0], t_max = ay[0];
		g_data.vExec("GetPoint", 1, ax, ay); real n_fit_points = ax[0], ndf = ay[0];
		g_data.vExec("GetPoint", 2, ax, ay); real chi2 = ax[0], chi2_ndf = ay[0];

		NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
		scale(Linear, Log);

		draw(g_mi, "l", blue);
		draw(g_pl, "l", heavygreen);

		draw(g_fit, "l", red);

		//AttachLegend();
	}
}

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=1mm, vSkip=1mm);
