import root;
import pad_layout;

include "../dip/common_code.asy";

string top_dir = "../../";

string methods[];
methods.push("first");

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

	NewPadLabel(methods[mi]);

	for (int qi : quantities.keys)
	{
		NewPad("$\sqrt s\ung{GeV}$", q_labels[qi]);

		string f = top_dir + "t_ranges/t_investigation.root";
		RootObject graph = RootGetObject(f, methods[mi] + "/g_" + quantities[qi] + "_vs_sqrt_s");
		RootObject fit = RootGetObject(f, methods[mi] + "/g_" + quantities[qi] + "_vs_sqrt_s|ff");

		draw(graph, "p", black, mCi+2pt+black);
		draw(fit, "l", red+1pt);

		real ext_sqrt_s = fit.rExec("GetParameter", 2);
		real ext_q = fit.rExec("GetParameter", 0);
		real ext_q_unc = fit.rExec("GetParError", 0);

		string l = format("$%.3f", ext_q) + format("\pm %.3f$", ext_q_unc);

		draw((ext_sqrt_s, ext_q), heavygreen, l, mCi+3pt+heavygreen);
		draw((ext_sqrt_s, ext_q-ext_q_unc)--(ext_sqrt_s, ext_q+ext_q_unc), heavygreen);

		AttachLegend();
	}
}
