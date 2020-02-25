import root;
import pad_layout;

include common_code;

string f = topDir + "s_extrapolation.root";

string extModel = "sqrt_s";
string extFit = "corr";

//----------------------------------------------------------------------------------------------------

RootObject g_settings = RootGetObject(f, "g_settings");

real ax[] = {0.};
real ay[] = {0.};
g_settings.vExec("GetPoint", 0, ax, ay); int n_parameters = (int)ax[0];

//----------------------------------------------------------------------------------------------------

NewPad(false);

NewPadLabel(replace(extFit + ", " + extModel, "_", "\_"));

//----------------------------------------------------------------------------------------------------

void DrawValues(RootObject g, pen p)
{
	for (int i = 0; i < g.iExec("GetN"); ++i)
	{
		real ax[] = {0.};
		real ay[] = {0.};
		g.vExec("GetPoint", i, ax, ay);

		label(format("$%.2f$", ay[0]), (ax[0], ay[0]), E, p);
	}
}

//----------------------------------------------------------------------------------------------------

for (int pari = 0; pari < n_parameters; ++pari)
{
	if ((pari % 3) == 0)
		NewRow();

	NewPad("$\sqrt s$", "");
	
	if (find(extModel, "log") >= 0)
		scale(Log, Linear);

	if (extModel == "sqrt_s")
		currentpad.xTicks = LeftTicks(2., 1.);

	string base = extModel + "/" + extFit + "/" + format("c_par%i", pari);

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

	DrawValues(g_cnt, blue);

	write(
		format("parameter %i: ", pari)
		+ format("%.3f, ", g_fit.rExec("Eval", 7.))
		+ format("%.3f, ", g_fit.rExec("Eval", 8.))
	);

	if (pari == 0) currentpad.yTicks = RightTicks(0.2, 0.1);
	if (pari == 1) currentpad.yTicks = RightTicks(2., 1.);
	if (pari == 2) currentpad.yTicks = RightTicks(20., 10.);

	if (pari == 3) currentpad.yTicks = RightTicks(0.2, 0.1);
	if (pari == 4) currentpad.yTicks = RightTicks(2., 1.);
	if (pari == 5) currentpad.yTicks = RightTicks(10., 5.);

	xlimits(0, 16, Crop);

	yaxis(XEquals(sqrt_s_ext, false), heavygreen);

	AttachLegend(BuildLegend(format("parameter %i", pari), S), N);
}

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=5mm, vSkip=1mm);
