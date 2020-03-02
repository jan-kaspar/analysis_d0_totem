import root;
import pad_layout;

include common_code;

TGraph_errorBar = None;

xSizeDef = 7cm;

f = "../..//fits/low_t,high_t/e012+e023:f=0.2/st+sy+no/do_fits.root";

//----------------------------------------------------------------------------------------------------

void DrawUncBand(RootObject g_cen, RootObject g_unc, pen p = yellow)
{
	guide g_up, g_dw;

	for (int i = 0; i < g_cen.iExec("GetN"); ++i)
	{
		real ax[] = {0};
		real ay[] = {0};

		g_cen.vExec("GetPoint", i, ax, ay); real x = ax[0], y = ay[0];
		real x_w = g_cen.rExec("GetErrorX", i);
		real x_l = x - x_w, x_h = x + x_w;

		g_unc.vExec("GetPoint", i, ax, ay); real y_unc = ay[0];

		g_up = g_up -- Scale((x_l, y + y_unc)) -- Scale((x_h, y + y_unc));
		g_dw = g_dw -- Scale((x_l, y - y_unc)) -- Scale((x_h, y - y_unc));
	}

	filldraw(g_up -- reverse(g_dw) -- cycle, p, nullpen);
}

//----------------------------------------------------------------------------------------------------

void DrawRelUncBand(RootObject g_cen, RootObject g_unc, RootObject g_fit, pen p = yellow)
{
	guide g_up, g_dw;

	for (int i = 0; i < g_cen.iExec("GetN"); ++i)
	{
		real ax[] = {0};
		real ay[] = {0};

		g_cen.vExec("GetPoint", i, ax, ay); real x = ax[0], y = ay[0];
		real x_w = g_cen.rExec("GetErrorX", i);
		real x_l = x - x_w, x_h = x + x_w;

		g_unc.vExec("GetPoint", i, ax, ay); real y_unc = ay[0];

		real f = g_fit.rExec("Eval", x);

		real r_up = (y + y_unc - f) / f;
		real r_dw = (y - y_unc - f) / f;

		g_up = g_up -- Scale((x_l, r_up)) -- Scale((x_h, r_up));
		g_dw = g_dw -- Scale((x_l, r_dw)) -- Scale((x_h, r_dw));
	}

	filldraw(g_up -- reverse(g_dw) -- cycle, p, nullpen);
}

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
	NewPadLabel(datasets[dsi]);
	//NewPadLabel(d_labels[dsi]);

//----------------------------------------------------------------------------------------------------

NewRow();

for (int dsi : datasets.keys)
{
	RootObject g_data = RootGetObject(f, datasets[dsi] + "/g_data");

	RootObject g_dsdt = RootGetObject(f, datasets[dsi] + "/g_dsdt");
	RootObject g_dsdt_syst_t_dep = RootGetObject(f, datasets[dsi] + "/g_dsdt_unc_syst_t_dep");
	RootObject g_dsdt_syst_full = RootGetObject(f, datasets[dsi] + "/g_dsdt_unc_syst_full");

	RootObject g_fit = RootGetObject(f, datasets[dsi] + "/g_fit");
	RootObject g_fit_pl_unc = RootGetObject(f, datasets[dsi] + "/g_fit_pl_unc");
	RootObject g_fit_mi_unc = RootGetObject(f, datasets[dsi] + "/g_fit_mi_unc");

	real ax[] = {0.};
	real ay[] = {0.};
	g_data.vExec("GetPoint", 0, ax, ay); real t_min = ax[0], t_max = ay[0];
	g_data.vExec("GetPoint", 1, ax, ay); real n_fit_points = ax[0], ndf = ay[0];
	g_data.vExec("GetPoint", 2, ax, ay); real chi2 = ax[0], chi2_ndf = ay[0];
	g_data.vExec("GetPoint", 3, ax, ay); real p_value = ax[0];
	g_data.vExec("GetPoint", 4, ax, ay); real t_dip = ax[0], t_dip_unc = ay[0];
	g_data.vExec("GetPoint", 5, ax, ay); real dsdt_dip = ax[0], dsdt_dip_unc = ay[0];
	g_data.vExec("GetPoint", 6, ax, ay); real t_bmp = ax[0], t_bmp_unc = ay[0];
	g_data.vExec("GetPoint", 7, ax, ay); real dsdt_bmp = ax[0], dsdt_bmp_unc = ay[0];
	g_data.vExec("GetPoint", 8, ax, ay); real R = ax[0], R_unc = ay[0];

	NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
	scale(Linear, Log);

	TGraph_x_min = t_min;
	TGraph_x_max = t_max;

	DrawUncBand(g_dsdt, g_dsdt_syst_full, yellow);
	DrawUncBand(g_dsdt, g_dsdt_syst_t_dep, heavygreen);

	draw(g_fit_pl_unc, "l", red+dashed);
	draw(g_fit_mi_unc, "l", red+dashed);
	draw(g_fit, "l", red+1pt);

	draw(g_dsdt, "p", blue, mCi+0.5pt);

	limits((0.3, 4e-3), (0.95, 5e-1), Crop);

	AddToLegend("<$\ch^2/\hbox{ndf} = " + format("%#.1f / (", chi2) + format("%.0f", n_fit_points) + format(" - %.0f)", n_fit_points - ndf) + format(" = %#.1f$", chi2_ndf));
	AddToLegend(format("<$\hbox{p-value} = %#.1f\un{\%}$", p_value * 100.));
	AddToLegend(format("<$R = %#.2f", R) + format("\pm %#.2f$", R_unc));
	AddToLegend(format("<$t_{\rm dip} = (%#.3f", t_dip) + format("\pm %#.3f)\un{GeV^2}$", t_dip_unc));

	AttachLegend(BuildLegend(vSkip=-1mm));
}

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dsi : datasets.keys)
{
	RootObject g_data = RootGetObject(f, datasets[dsi] + "/g_data");

	RootObject g_dsdt = RootGetObject(f, datasets[dsi] + "/g_dsdt");
	RootObject g_dsdt_syst_t_dep = RootGetObject(f, datasets[dsi] + "/g_dsdt_unc_syst_t_dep");
	RootObject g_dsdt_syst_full = RootGetObject(f, datasets[dsi] + "/g_dsdt_unc_syst_full");

	RootObject g_fit = RootGetObject(f, datasets[dsi] + "/g_fit");
	RootObject g_fit_pl_unc = RootGetObject(f, datasets[dsi] + "/g_fit_pl_unc");
	RootObject g_fit_mi_unc = RootGetObject(f, datasets[dsi] + "/g_fit_mi_unc");

	real ax[] = {0.};
	real ay[] = {0.};
	g_data.vExec("GetPoint", 0, ax, ay); real t_min = ax[0], t_max = ay[0];

	NewPad("$|t|\ung{GeV^2}$", "$(\d\si/\d t - \hbox{fit}) / \hbox{fit}$");

	DrawRelUncBand(g_dsdt, g_dsdt_syst_full, g_fit, yellow);
	DrawRelUncBand(g_dsdt, g_dsdt_syst_t_dep, g_fit, heavygreen);

	for (int i = 0; i < g_dsdt.iExec("GetN"); ++i)
	{
		g_dsdt.vExec("GetPoint", i, ax, ay); real t = ax[0], dsdt = ay[0];
		real t_unc = g_dsdt.rExec("GetErrorX", i);
		real dsdt_unc = g_dsdt.rExec("GetErrorY", i);

		if (t < t_min || t > t_max)
			continue;

		real fit = g_fit.rExec("Eval", t);
		real r = (dsdt - fit) / fit;
		real r_unc = dsdt_unc / fit;

		draw((t-t_unc, r)--(t+t_unc, r), blue);
		draw((t, r-r_unc)--(t, r+r_unc), blue);
	}

	guide gr_fit, gr_fit_pl_unc, gr_fit_mi_unc;

	for (int i = 0; i < g_fit.iExec("GetN"); ++i)
	{
		g_fit.vExec("GetPoint", i, ax, ay); real t = ax[0], fit = ay[0];
		g_fit_pl_unc.vExec("GetPoint", i, ax, ay); real fit_pl_unc = ay[0];
		g_fit_mi_unc.vExec("GetPoint", i, ax, ay); real fit_mi_unc = ay[0];

		real r_fit = 0.;
		real r_fit_pl_unc = (fit_pl_unc - fit) / fit;
		real r_fit_mi_unc = (fit_mi_unc - fit) / fit;

		gr_fit = gr_fit -- (t, r_fit);
		gr_fit_pl_unc = gr_fit_pl_unc -- (t, r_fit_pl_unc);
		gr_fit_mi_unc = gr_fit_mi_unc -- (t, r_fit_mi_unc);
	}
	
	draw(gr_fit, red+1pt);
	draw(gr_fit_pl_unc, red+dashed);
	draw(gr_fit_mi_unc, red+dashed);

	xlimits(0.3, 0.95, Crop);
}

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=1mm, vSkip=1mm);
