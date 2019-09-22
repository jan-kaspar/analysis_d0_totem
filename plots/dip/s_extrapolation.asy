import root;
import pad_layout;

include common_code;

string f = topDir + "s_extrapolation.root";

NewPad("$|t|\ung{GeV^2}$", "$\d\si/\d t\ung{mb/GeV}$");
scale(Linear, Log);

int idx = 0;

for (int emi : extModels.keys)
{
	for (int efi : extFits.keys)
	{
		string base = extModels[emi] + "/" + extFits[efi];
		string l = replace(extFits[efi] + ", " + extModels[emi], "_", "\_");
		draw(RootGetObject(f, base + "/g_dsdt_ext"), "l", StdPen(++idx), l);
	}
}

limits((0.3, 1e-3), (1.2, 1e1), Crop);

AttachLegend(NW, NE);
