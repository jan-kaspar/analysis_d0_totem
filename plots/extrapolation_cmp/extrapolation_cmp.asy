import root;
import pad_layout;

string topDir = "../../";

string results[];
pen r_pens[];

results.push("minimal/e01+e023:f=0.5/st+sy"); r_pens.push(blue+dashed);
results.push("minimal/e01+e023:f=0.05/st+sy"); r_pens.push(blue);

results.push("low_t,high_t/e012+e023:f=0.5/st+sy"); r_pens.push(red+dashed);
results.push("low_t,high_t/e012+e023:f=0.2/st+sy"); r_pens.push(red);

string extModel = "sqrt_s";
string extFit = "corr";

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$");
scale(Linear, Log);

draw(RootGetObject("../../data/D0_1.9TeV/d0_dsdt_1.96TeV.root", "g_D0_dsdt_1.96TeV"), "l,p", black+dashed, mCi+1pt, "D0 measurement");

AddToLegend("<TOTEM extrapolations:");
//string l = replace(extFit + ", " + extModel, "_", "\_");

for (int ri : results.keys)
{
	string f = topDir + "fits/" + results[ri] + "/s_extrapolation.root";

	string base = extModel + "/" + extFit;

	//pen p = StdPen(ri + 1);
	pen p = r_pens[ri];

	string l = replace(results[ri], "_","\_");
	draw(RootGetObject(f, base + "/g_dsdt_ext"), "l", p, l);
}

limits((0.4, 4e-3), (0.9, 4e-1), Crop);

AttachLegend(NW, NE);
