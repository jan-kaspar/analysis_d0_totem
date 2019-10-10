import root;
import pad_layout;

include "../dip/common_code.asy";

string top_dir = "../../";

string methods[];
methods.push("first");

TGraph_errorBar = None;

xSizeDef = 7cm;

//----------------------------------------------------------------------------------------------------

NewPad(false);

for (int dsi : datasets.keys)
	NewPadLabel(d_labels[dsi]);

//----------------------------------------------------------------------------------------------------
for (int mi : methods.keys)
{
	NewRow();

	NewPadLabel(methods[mi]);

	for (int dsi : datasets.keys)
	{

		string f = top_dir + "data/TOTEM_" + datasets[dsi] + "/data.root";
		RootObject g_dsdt = RootGetObject(f, "g_dsdt");
		RootObject f_fit = RootGetObject(f, "f_dsdt");

		string f = top_dir + "t_ranges/t_investigation.root";
		RootObject g_data = RootGetObject(f, methods[mi] + "/" + datasets[dsi] + "/g_data");

		real ax[] = {0.};
		real ay[] = {0.};
		g_data.vExec("GetPoint", 0, ax, ay); real t_min = ax[0], t_min_unc = ay[0];
		g_data.vExec("GetPoint", 1, ax, ay); real dsdt_min = ax[0], dsdt_min_unc = ay[0];

		g_data.vExec("GetPoint", 2, ax, ay); real t_dip = ax[0], t_dip_unc = ay[0];
		g_data.vExec("GetPoint", 3, ax, ay); real dsdt_dip = ax[0], dsdt_dip_unc = ay[0];

		g_data.vExec("GetPoint", 4, ax, ay); real t_bmp = ax[0], t_bmp_unc = ay[0];
		g_data.vExec("GetPoint", 5, ax, ay); real dsdt_bmp = ax[0], dsdt_bmp_unc = ay[0];

		g_data.vExec("GetPoint", 6, ax, ay); real t_max = ax[0], t_max_unc = ay[0];
		g_data.vExec("GetPoint", 7, ax, ay); real dsdt_max = ax[0], dsdt_max_unc = ay[0];

		NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV^2}$");
		scale(Linear, Log);

		real dt_dw = -0.01, dt_up = 0.01;
		real f_y_up = 1.2, f_y_dw = 0.8;

		draw(Scale((t_min, dsdt_min * f_y_dw))--Scale((t_min, dsdt_min * f_y_up)), cyan+1.3pt);
		draw(Scale((t_min + dt_dw, dsdt_min))--Scale((t_dip + dt_up, dsdt_min)), cyan+1.3pt);

		filldraw(Scale((t_dip - t_dip_unc, dsdt_dip * f_y_dw))--Scale((t_dip + t_dip_unc, dsdt_dip * f_y_dw))
			--Scale((t_dip + t_dip_unc, dsdt_dip * f_y_up))--Scale((t_dip - t_dip_unc, dsdt_dip * f_y_up))--cycle, red+opacity(0.3), nullpen);
		draw(Scale((t_dip, dsdt_dip * f_y_dw))--Scale((t_dip, dsdt_dip * f_y_up)), red+1.3pt);
		draw(Scale((t_dip + dt_dw, dsdt_dip))--Scale((t_dip + dt_up, dsdt_dip)), red+1.3pt);

		filldraw(Scale((t_bmp - t_bmp_unc, dsdt_bmp * f_y_dw))--Scale((t_bmp + t_bmp_unc, dsdt_bmp * f_y_dw))
			--Scale((t_bmp + t_bmp_unc, dsdt_bmp * f_y_up))--Scale((t_bmp - t_bmp_unc, dsdt_bmp * f_y_up))--cycle, heavygreen+opacity(0.3), nullpen);
		draw(Scale((t_bmp, dsdt_bmp * f_y_dw))--Scale((t_bmp, dsdt_bmp * f_y_up)), heavygreen+1.3pt);
		draw(Scale((t_dip + dt_dw, dsdt_bmp))--Scale((t_max + dt_up, dsdt_bmp)), heavygreen+1.3pt);

		draw(Scale((t_max, dsdt_max * f_y_dw))--Scale((t_max, dsdt_max * f_y_up)), magenta+1.3pt);
		draw(Scale((t_bmp + dt_dw, dsdt_max))--Scale((t_max + dt_up, dsdt_max)), magenta+1.3pt);

		draw(f_fit, "l", black+1pt);
		draw(g_dsdt, "p", blue, mCi+0.5pt);

		limits((t_min-0.05, dsdt_dip*0.7), (t_max+0.05, dsdt_min*1.5), Crop);
	}
}
