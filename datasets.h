#include "TF1.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct Dataset
{
	double sqrt_s;
	string name;
	string f_in;
	double t_min, t_max;
	double t_dip, t_bmp;
	TF1 *ff;
};

//----------------------------------------------------------------------------------------------------

vector<Dataset> datasets;

Dataset ds_ext;

void InitDatasets()
{
	Dataset d2;
	d2.sqrt_s = 2.76;
	d2.name = "2.76TeV";
	d2.f_in = "../../data/TOTEM_2.76TeV/data.root";
	d2.t_min = 0.48;
	d2.t_dip = 0.615;
	d2.t_bmp = 0.7;
	d2.t_max = 0.77;
	datasets.push_back(d2);

	Dataset d7;
	d7.sqrt_s = 7.;
	d7.name = "7TeV";
	d7.f_in = "../../data/TOTEM_7TeV/data.root";
	d7.t_min = 0.395;
	d7.t_dip = 0.529;
	d7.t_bmp = 0.701;
	d7.t_max = 0.8;
	datasets.push_back(d7);

	Dataset d8;
	d8.sqrt_s = 8.;
	d8.name = "8TeV";
	d8.f_in = "../../data/TOTEM_8TeV/data.root";
	d8.t_min = 0.4;
	d8.t_dip = 0.526;
	d8.t_bmp = 0.720;
	d8.t_max = 0.8;
	datasets.push_back(d8);

	Dataset d13;
	d13.sqrt_s = 13.;
	d13.name = "13TeV";
	d13.f_in = "../../data/TOTEM_13TeV/data.root";
	d13.t_min = 0.39;	// TODO: should be 0.363
	d13.t_dip = 0.467;
	d13.t_bmp = 0.636;
	d13.t_max = 0.70;	// TODO: should be 0.76
	datasets.push_back(d13);

	//--------------------

	ds_ext.sqrt_s = 1.96;
	ds_ext.name = "1.96TeV";
	ds_ext.f_in = "";
	ds_ext.t_min = 0.4;
	ds_ext.t_dip = 0.64;
	ds_ext.t_bmp = 0.75;
	ds_ext.t_max = 1.1;
}

//----------------------------------------------------------------------------------------------------

int n_parameters = 0;

void BuildFitFunction(const string fitModel, Dataset &ds)
{
	if (fitModel == "e012+e012")
	{
		const double t1 = (ds.t_dip + ds.t_min) / 2.;
		//const double t2 = (ds.t_max + ds.t_dip) / 2.;
		const double t2 = ds.t_bmp;

		n_parameters = 6;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f) + [5]*(x-%.3f)*(x-%.3f))", t1, t1, t1, t2, t2, t2);
		ds.ff = new TF1("ff", buf);

		if (ds.name == "2.76TeV") ds.ff->SetParameters(-3.78, -19.3, -0.3, -4.19, -1.68, -264.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, 3.29, -27.3);
		if (ds.name == "8TeV") ds.ff->SetParameters(-3.97, -34.4, -46.5, -3.53, 0.174, -15.31);
		if (ds.name == "13TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, 3.29, -27.3);

		return;
	}

	if (fitModel == "e012+e02")
	{
		const double t1 = (ds.t_dip + ds.t_min) / 2.;
		//const double t2 = (ds.t_max + ds.t_dip) / 2.;
		const double t2 = ds.t_bmp;

		n_parameters = 5;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f)*(x-%.3f))", t1, t1, t1, t2, t2);
		ds.ff = new TF1("ff", buf);

		if (ds.name == "2.76TeV") ds.ff->SetParameters(-3.78, -19.3, -0.3, -4.19, -264.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, -27.3);
		if (ds.name == "8TeV") ds.ff->SetParameters(-3.97, -34.4, -46.5, -3.53, -15.31);
		if (ds.name == "13TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, -27.3);

		return;
	}

	if (fitModel == "e012+e012_t_dip")
	{
		n_parameters = 6;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f) + [5]*(x-%.3f)*(x-%.3f))",
			ds.t_dip, ds.t_dip, ds.t_dip, ds.t_dip, ds.t_dip, ds.t_dip);
		ds.ff = new TF1("ff", buf);

		ds.ff->SetParameters(-5.58, -53.4, -135, -3.74, 8.26, -24.1);

		return;
	}

	if (fitModel == "e012+e012_t_0")
	{
		n_parameters = 6;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*x + [2]*x*x) + exp([3] + [4]*x + [5]*x*x)");
		ds.ff = new TF1("ff", buf);

		ds.ff->SetParameters(-3.68, 40.5, -94.5, -12.9, 30.9, -24.1);

		return;
	}

	if (fitModel == "e01-int-e01")
	{
		n_parameters = 5;

		char buf[500];
		sprintf(buf, "exp(2*[0] + 2*[1]*(x-%.3f)) + exp(2*[2] + 2*[3]*(x-%.3f)) - 2 *cos([4]) * exp( [0] + [1]*(x-%.3f) + [2] + [3]*(x-%.3f) )",
			ds.t_dip, ds.t_dip, ds.t_dip, ds.t_dip);
		ds.ff = new TF1("ff", buf);

		ds.ff->SetParameters(-0.816, -8.24, -0.713, -2.535, 0.341);

		if (ds.name == "7TeV") ds.ff->SetParameters(0.48, -4.96, 0.51, -3.97, -0.072);

		return;
	}

	printf("ERROR in BuildFitFunction\n");
}
