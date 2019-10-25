import root;
import pad_layout;

include "../templates/common_code.asy";

string top_dir = "../../";

string methods[];
methods.push("minimal");
methods.push("low_t,high_t");

string quantities[], q_labels[];
quantities.push("t_min"); q_labels.push("$t_{\rm min}\ung{GeV^2}$");
quantities.push("t_dip"); q_labels.push("$t_{\rm dip}\ung{GeV^2}$");
quantities.push("t_bmp"); q_labels.push("$t_{\rm bmp}\ung{GeV^2}$");
quantities.push("t_max"); q_labels.push("$t_{\rm max}\ung{GeV^2}$");

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

for (int mi : methods.keys)
{
	NewRow();

	NewPadLabel(replace(methods[mi], "_", "\_"));

	for (int qi : quantities.keys)
	{
		NewPad("$\sqrt s\ung{GeV}$", q_labels[qi]);
		scale(Log, Linear);

		string f = top_dir + "t_ranges/t_investigation_new.root";
		RootObject graph = RootGetObject(f, methods[mi] + "/g_" + quantities[qi] + "_vs_sqrt_s");
		RootObject fit_lin = RootGetObject(f, methods[mi] + "/g_" + quantities[qi] + "_vs_sqrt_s|ff_lin");
		RootObject fit_log = RootGetObject(f, methods[mi] + "/g_" + quantities[qi] + "_vs_sqrt_s|ff_log");

		draw(graph, "p", black, mCi+2pt+black);
		draw(fit_log, "l", red+1pt);
		draw(fit_lin, "l", red+dashed);

		RootObject fit = fit_log;
		real ext_sqrt_s = fit.rExec("GetParameter", 2);
		real ext_q = fit.rExec("GetParameter", 0);
		real ext_q_unc = fit.rExec("GetParError", 0);

		string l = format("$%.3f", ext_q) + format("\pm %.3f$", ext_q_unc);

		draw(Scale((ext_sqrt_s, ext_q)), heavygreen, l, mCi+3pt+heavygreen);
		draw(Scale((ext_sqrt_s, ext_q-ext_q_unc))--Scale((ext_sqrt_s, ext_q+ext_q_unc)), heavygreen);

		xlimits(1, 14, Crop);

		AttachLegend();
	}
}

GShipout(hSkip=1mm, vSkip=1mm);
