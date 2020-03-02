import root;
import pad_layout;

include common_code;

TGraph_errorBar = None;

xSizeDef = 7cm;

f = "../..//fits/low_t,high_t/e012+e023:f=0.2/st+sy+no/do_fits.root";

//----------------------------------------------------------------------------------------------------

void DrawBinSize(RootObject g, pen p, mark m)
{
	guide gd;

	for (int i = 0; i < g.iExec("GetN"); ++i)
	{
		real ax[] = {0};
		real ay[] = {0};

		g.vExec("GetPoint", i, ax, ay); real x = ax[0], y = ay[0];
		real bs = 2. * g.rExec("GetErrorX", i);

		if (ax[0] < 0.3 || ax[0] > 0.6)
			continue;

		draw((ax[0], bs), m);

		gd = gd -- (ax[0], bs);
	}

	draw(gd, p);
}

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
	NewPadLabel(datasets[dsi]);
	//NewPadLabel(d_labels[dsi]);

//----------------------------------------------------------------------------------------------------

NewRow();

for (int dsi : datasets.keys)
{
	RootObject g_dsdt = RootGetObject(f, datasets[dsi] + "/g_dsdt");

	NewPad("$|t|\ung{GeV^2}$", "bin size $\ung{GeV^2}$");
	//scale(Linear, Log);

	//TGraph_x_min = t_min;
	//TGraph_x_max = t_max;

	DrawBinSize(g_dsdt, blue, mCi+blue+1pt);

	//limits((0.3, 4e-3), (0.95, 5e-1), Crop);

	//AttachLegend(BuildLegend(vSkip=-1mm));
}


GShipout(hSkip=1mm, vSkip=1mm);
