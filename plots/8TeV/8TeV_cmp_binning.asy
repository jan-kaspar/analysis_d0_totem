import root;
import pad_layout;

string topDir = "../../";

//string dir = "fits/minimal/e01+e023/";
string dir = "fits/low_t,high_t/e012+e023:f=0.2/st+sy+no/";

int n_parameters = 6;

string datasets[], d_labels[];
datasets.push("8TeV-b1-45b"); d_labels.push("8TeV, 45b-56t, binning 1");
datasets.push("8TeV-b1-45t"); d_labels.push("8TeV, 45t-56b, binning 1");

datasets.push("8TeV-b1"); d_labels.push("8TeV, both diagonals, binning 1");

datasets.push("8TeV-b2"); d_labels.push("8TeV, both diagonals, binning 2");
datasets.push("8TeV-b3"); d_labels.push("8TeV, both diagonals, binning 3");

TGraph_errorBar = None;

//----------------------------------------------------------------------------------------------------

for (int pari = 0; pari < n_parameters; ++pari)
{
	NewRow();

	NewPadLabel(format("parameter %i", pari));
	
	NewPad("dataset", "");

	for (int dsi : datasets.keys)
	{
		string f = topDir + dir + "/do_fits.root";
		RootObject g_par = RootGetObject(f, datasets[dsi] + "/g_par");

		real ax[] = {0.};
		real ay[] = {0.};
		g_par.vExec("GetPoint", pari, ax, ay);
		real v = ay[0];
		real u = g_par.rExec("GetErrorY", pari);

		//write(pari, v);

		pen p = StdPen(dsi);

		draw((dsi, v), mCi+3pt+p);
		draw((dsi, v-u)--(dsi, v+u), p);
	}

	xlimits(-1, datasets.length);

	if (pari == 0)
	{
		NewPad(false);
		for (int dsi : datasets.keys)
			AddToLegend(d_labels[dsi], StdPen(dsi));
		AttachLegend();
	}
}

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=1mm, vSkip=1mm);
