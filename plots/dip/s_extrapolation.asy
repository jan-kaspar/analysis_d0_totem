import root;
import pad_layout;

include common_code;

string f = topDir + "s_extrapolation.root";

string base = "sqrt_s/corr";

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$");
scale(Linear, Log);

draw(RootGetObject("../../../../../data/D0_1.9TeV/d0_dsdt_1.96TeV.root", "g_D0_dsdt_1.96TeV"), "l,p", black+dashed, mCi+1pt, "D0 measurement");

AddToLegend("<TOTEM extrapolation:");

draw(RootGetObject(f, base + "/g_dsdt_ext"), "l", red, "central value");

draw(RootGetObject(f, base + "/g_dsdt_ext_pl_unc"), "l", red+dashed, "uncertainty band");
draw(RootGetObject(f, base + "/g_dsdt_ext_mi_unc"), "l", red+dashed);

limits((0.4, 4e-3), (0.9, 1e-1), Crop);

AttachLegend(NW, NE);
