import root;
import pad_layout;

include "../templates/common_code.asy";

string topDir = "../../";

string results[];
//results.push("low_t,high_t/e012+e023:f=0.2/st+sy+no");
results.push("lts-100,ht/e012+e023:f=0.4/st+sy+no");

string extModel = "sqrt_s";
string extFit = "corr";

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

void DrawPoint(real t, real unc_t, real dsdt, real unc_dsdt)
{
	pen p = black;

	draw(Scale((t - unc_t, dsdt))--Scale((t + unc_t, dsdt)), p);
	draw(Scale((t, dsdt - unc_dsdt))--Scale((t, dsdt + unc_dsdt)), p);
}

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$");
scale(Linear, Log);

// draw extrapolation
AddToLegend("<extraplation to D0 energy:");

AddToLegend("Simone's results", mPl+5pt+black);

DrawPoint(0.487460758520051, 0.009103328134492, 0.067951413556073, 0.016344625232073);
DrawPoint(0.517979343255533, 0.009390427618439, 0.045478605998104, 0.006638457964865);
DrawPoint(0.5940506617077, 0.009801516571181, 0.013631988759905, 0.002228824230835);
DrawPoint(0.651560332791308, 0.010339399250097, 0.006242439599361, 0.001752638959779);
DrawPoint(0.718039806233958, 0.016232148258169, 0.009497358286628, 0.002044128939274);
DrawPoint(0.859522341120077, 0.020610904129264, 0.012356544618783, 0.001804992995596);
DrawPoint(0.91503290799268, 0.025366543203122, 0.00905649562718, 0.003716851593585);
DrawPoint(0.97000255036621, 0.030409684526518, 0.006547145796336, 0.001400686048484);

AddToLegend("Jan's results:");

for (int ri : results.keys)
{
	string f = topDir + "fits/" + results[ri] + "/s_extrapolation.root";
	string base = extModel + "/" + extFit;

	pen p = StdPen(ri + 1);
	string lbl = replace(results[ri], "_", "\_");

	draw(RootGetObject(f, base + "/g_dsdt_ext"), "l", p, lbl);
	draw(RootGetObject(f, base + "/unc c/g_dsdt_ext_pl_unc"), "l", p+dashed);
	draw(RootGetObject(f, base + "/unc c/g_dsdt_ext_mi_unc"), "l", p+dashed);
}

AttachLegend(BuildLegend(NW), NE);
