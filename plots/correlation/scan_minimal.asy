import root;
import pad_layout;

string topDir = "../../";

string energies[];
pen e_pens[];
energies.push("2.76TeV"); e_pens.push(black);
energies.push("7TeV"); e_pens.push(red);
energies.push("8TeV"); e_pens.push(blue);
energies.push("13TeV"); e_pens.push(heavygreen);

string dirs[];
real d_fs[];
dirs.push("fits/minimal/e01+e023:f=0.0/st+sy"); d_fs.push(0.0);
dirs.push("fits/minimal/e01+e023:f=0.05/st+sy"); d_fs.push(0.05);
dirs.push("fits/minimal/e01+e023:f=0.1/st+sy"); d_fs.push(0.1);
dirs.push("fits/minimal/e01+e023:f=0.15/st+sy"); d_fs.push(0.15);
dirs.push("fits/minimal/e01+e023:f=0.2/st+sy"); d_fs.push(0.2);
dirs.push("fits/minimal/e01+e023:f=0.3/st+sy"); d_fs.push(0.3);
dirs.push("fits/minimal/e01+e023:f=0.4/st+sy"); d_fs.push(0.4);
dirs.push("fits/minimal/e01+e023:f=0.5/st+sy"); d_fs.push(0.5);
dirs.push("fits/minimal/e01+e023:f=0.6/st+sy"); d_fs.push(0.6);
dirs.push("fits/minimal/e01+e023:f=0.7/st+sy"); d_fs.push(0.7);
dirs.push("fits/minimal/e01+e023:f=0.8/st+sy"); d_fs.push(0.8);
dirs.push("fits/minimal/e01+e023:f=0.9/st+sy"); d_fs.push(0.9);
dirs.push("fits/minimal/e01+e023:f=1.0/st+sy"); d_fs.push(1.0);

string corrs[]; int c_1[], c_2[];
corrs.push("(0, 1)"); c_1.push(0); c_2.push(1);

corrs.push("(0, 2)"); c_1.push(0); c_2.push(2);
corrs.push("(0, 3)"); c_1.push(0); c_2.push(3);
corrs.push("(0, 4)"); c_1.push(0); c_2.push(4);

corrs.push("(2, 3)"); c_1.push(2); c_2.push(3);
corrs.push("(2, 4)"); c_1.push(2); c_2.push(4);
corrs.push("(3, 4)"); c_1.push(3); c_2.push(4);

//----------------------------------------------------------------------------------------------------

for (int cri : corrs.keys)
{
	if (cri == 1 || cri == 4)
		NewRow();

	NewPad("$f$", "correlation factor");

	for (int eni : energies.keys)
	{
		guide g;

		for (int dri : dirs.keys)
		{
			RootObject hist = RootGetObject(topDir + dirs[dri] + "/do_fits.root", energies[eni] + "/h2_C");
			real c = hist.rExec("GetBinContent", c_1[cri]+1, c_2[cri]+1);

			g = g--(d_fs[dri], c);

			draw((d_fs[dri], c), mCi+1pt + e_pens[eni]);
		}

		draw(g, e_pens[eni]);
	}

	limits((0, -1.), (1., +1), Crop);

	AttachLegend(BuildLegend(corrs[cri], S), N);
}

//----------------------------------------------------------------------------------------------------

NewPad(false);

for (int eni : energies.keys)
	AddToLegend(energies[eni], e_pens[eni]);

AttachLegend();

//----------------------------------------------------------------------------------------------------

GShipout(hSkip=5mm, vSkip=1mm);
