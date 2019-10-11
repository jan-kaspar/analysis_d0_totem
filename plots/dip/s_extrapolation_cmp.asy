import root;
import pad_layout;

include common_code;

string f = topDir + "s_extrapolation.root";

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$");
scale(Linear, Log);

for (int emi : extModels.keys)
{
	for (int efi : extFits.keys)
	{
		string base = extModels[emi] + "/" + extFits[efi];
		string l = replace(extFits[efi] + ", " + extModels[emi], "_", "\_");

		pen p = StdPen(emi + 1);
		if (extFits[efi] == "uncorr")
			p += dashed;

		draw(RootGetObject(f, base + "/g_dsdt_ext"), "l", p, l);
	}
}

draw(RootGetObject("../../../../../data/D0_1.9TeV/d0_dsdt_1.96TeV.root", "g_D0_dsdt_1.96TeV"), "l,p", black+dashed, mCi+1pt, "D0 measurement");

limits((0.4, 4e-3), (0.9, 1e-1), Crop);

AttachLegend(NW, NE);
