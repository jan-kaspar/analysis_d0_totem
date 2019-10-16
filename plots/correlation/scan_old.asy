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
dirs.push("fits/e012_0.0+e012/first/st+sy"); d_fs.push(0.0);
dirs.push("fits/e012_0.1+e012/first/st+sy"); d_fs.push(0.1);
dirs.push("fits/e012_0.2+e012/first/st+sy"); d_fs.push(0.2);
dirs.push("fits/e012_0.3+e012/first/st+sy"); d_fs.push(0.3);
dirs.push("fits/e012_0.4+e012/first/st+sy"); d_fs.push(0.4);
dirs.push("fits/e012_0.5+e012/first/st+sy"); d_fs.push(0.5);
dirs.push("fits/e012_0.6+e012/first/st+sy"); d_fs.push(0.6);
dirs.push("fits/e012_0.7+e012/first/st+sy"); d_fs.push(0.7);
dirs.push("fits/e012_0.8+e012/first/st+sy"); d_fs.push(0.8);
dirs.push("fits/e012_0.9+e012/first/st+sy"); d_fs.push(0.9);
dirs.push("fits/e012_1.0+e012/first/st+sy"); d_fs.push(1.0);

string corrs[]; int c_1[], c_2[];
corrs.push("(0, 1)"); c_1.push(0); c_2.push(1);
corrs.push("(0, 2)"); c_1.push(0); c_2.push(2);
corrs.push("(1, 2)"); c_1.push(1); c_2.push(2);

corrs.push("(0, 3)"); c_1.push(0); c_2.push(3);
corrs.push("(0, 5)"); c_1.push(0); c_2.push(5);
corrs.push("(1, 3)"); c_1.push(1); c_2.push(3);

corrs.push("(3, 4)"); c_1.push(2); c_2.push(4);
corrs.push("(3, 5)"); c_1.push(3); c_2.push(5);
corrs.push("(4, 5)"); c_1.push(4); c_2.push(5);

//----------------------------------------------------------------------------------------------------

for (int cri : corrs.keys)
{
	if (cri % 3 == 0)
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
