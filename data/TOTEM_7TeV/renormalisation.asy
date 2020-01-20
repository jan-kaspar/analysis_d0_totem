import root;
import pad_layout;

string f = "re_normalisation.root";

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

void DrawOne(real fac, string l)
{
	NewPad();
	scale(Linear, Log);

	transform ct = shift(0, log10(fac));
	draw(ct, RootGetObject(f, "g_3p5"), "p", red, "$3.5\un{m}$");
	draw(ct, RootGetObject(f, "g_3p5|ff"), red+dashed);

	draw(RootGetObject(f, "h_90"), "eb", blue, "$90\un{m}$");
	draw(RootGetObject(f, "h_90|ff"), blue+dashed);

	limits((0.3, 1e-2), (0.5, 1e0), Crop);

	AttachLegend(l);
}

//----------------------------------------------------------------------------------------------------

DrawOne(1.0, "before");

DrawOne(0.848, "after");
