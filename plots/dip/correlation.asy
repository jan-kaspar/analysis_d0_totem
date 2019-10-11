import root;
import pad_layout;

include common_code;

TH2_palette = Gradient(blue, white, red);

TH2_text_format = "%+0.2f";

//----------------------------------------------------------------------------------------------------

for (int dsi : datasets.keys)
	NewPadLabel(d_labels[dsi]);

//----------------------------------------------------------------------------------------------------
NewRow();

for (int dsi : datasets.keys)
{
	NewPad();

	RootObject h2_C = RootGetObject(f, datasets[dsi] + "/h2_C");

	TH2_z_min = -1;
	TH2_z_max = +1;

	draw(h2_C, "p,t,bar");
}

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=1mm, vSkip=1mm);
