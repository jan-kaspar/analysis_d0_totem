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

draw(RootGetObject("../../../../../data/TOTEM_2.76TeV/data.root", "g_dsdt"), "l,p", blue+dashed, mCi+1pt+blue, "TOTEM at 2.76");

AddToLegend("<TOTEM extrapolation:");

draw(RootGetObject(f, base + "/g_dsdt_ext"), "l", red, "central value");

draw(RootGetObject(f, base + "/unc c/g_dsdt_ext_pl_unc"), "l", red+dashed, "uncertainty band");
draw(RootGetObject(f, base + "/unc c/g_dsdt_ext_mi_unc"), "l", red+dashed);

limits((0.4, 4e-3), (0.9, 4e-1), Crop);

AttachLegend(NW, NE);

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$B(t)\ung{GeV^{-2}}$");
//scale(Linear, Log);

AddToLegend("<TOTEM extrapolation:");

draw(RootGetObject(f, base + "/g_B_ext"), "l", red, "central value");

draw(RootGetObject(f, base + "/unc c/g_B_ext_pl_unc"), "l", red+dashed);
draw(RootGetObject(f, base + "/unc c/g_B_ext_mi_unc"), "l", red+dashed);

limits((0.4, -22.), (0.9, +8), Crop);
