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
	string topDir = "../../../../data/";

	Dataset d2;
	d2.sqrt_s = 2.76;
	d2.name = "2.76TeV";
	d2.f_in = topDir + "TOTEM_2.76TeV/data.root";
	datasets.push_back(d2);

	Dataset d7;
	d7.sqrt_s = 7.;
	d7.name = "7TeV";
	d7.f_in = topDir + "TOTEM_7TeV/data.root";
	datasets.push_back(d7);

	Dataset d8;
	d8.sqrt_s = 8.;
	d8.name = "8TeV";
	d8.f_in = topDir + "TOTEM_8TeV/data.root";
	datasets.push_back(d8);

	Dataset d13;
	d13.sqrt_s = 13.;
	d13.name = "13TeV";
	d13.f_in = topDir + "TOTEM_13TeV/data.root";
	datasets.push_back(d13);

	//--------------------

	ds_ext.sqrt_s = 1.96;
	ds_ext.name = "1.96TeV";
	ds_ext.f_in = "";
}

//----------------------------------------------------------------------------------------------------

void BuildTRanges(const string tRangeModel, Dataset &ds)
{
	if (tRangeModel == "minimal")
	{
		if (ds.name == "2.76TeV")	{ ds.t_min = 0.530; ds.t_dip = 0.616; ds.t_bmp = 0.750; ds.t_max = 0.800; }
		if (ds.name == "7TeV")		{ ds.t_min = 0.442; ds.t_dip = 0.530; ds.t_bmp = 0.694; ds.t_max = 0.780; }
		if (ds.name == "8TeV")		{ ds.t_min = 0.433; ds.t_dip = 0.518; ds.t_bmp = 0.687; ds.t_max = 0.800; }
		if (ds.name == "13TeV")		{ ds.t_min = 0.400; ds.t_dip = 0.470; ds.t_bmp = 0.638; ds.t_max = 0.707; }

		if (ds.name == "1.96TeV")	{ ds.t_min = 0.524; ds.t_dip = 0.611; ds.t_bmp = 0.745; ds.t_max = 0.808; }
	}

	if (tRangeModel == "low_t")
	{
		if (ds.name == "2.76TeV")	{ ds.t_min = 0.450; ds.t_dip = 0.616; ds.t_bmp = 0.750; ds.t_max = 0.800; }
		if (ds.name == "7TeV")		{ ds.t_min = 0.368; ds.t_dip = 0.530; ds.t_bmp = 0.694; ds.t_max = 0.780; }
		if (ds.name == "8TeV")		{ ds.t_min = 0.375; ds.t_dip = 0.518; ds.t_bmp = 0.687; ds.t_max = 0.800; }
		if (ds.name == "13TeV")		{ ds.t_min = 0.338; ds.t_dip = 0.470; ds.t_bmp = 0.638; ds.t_max = 0.707; }

		if (ds.name == "1.96TeV")	{ ds.t_min = 0.443; ds.t_dip = 0.611; ds.t_bmp = 0.745; ds.t_max = 0.808; }
	}

	if (tRangeModel == "high_t")
	{
		if (ds.name == "2.76TeV")	{ ds.t_min = 0.530; ds.t_dip = 0.616; ds.t_bmp = 0.750; ds.t_max = 0.990; }
		if (ds.name == "7TeV")		{ ds.t_min = 0.442; ds.t_dip = 0.530; ds.t_bmp = 0.694; ds.t_max = 0.878; }
		if (ds.name == "8TeV")		{ ds.t_min = 0.433; ds.t_dip = 0.518; ds.t_bmp = 0.687; ds.t_max = 0.900; }
		if (ds.name == "13TeV")		{ ds.t_min = 0.400; ds.t_dip = 0.470; ds.t_bmp = 0.638; ds.t_max = 0.858; }

		if (ds.name == "1.96TeV")	{ ds.t_min = 0.524; ds.t_dip = 0.611; ds.t_bmp = 0.745; ds.t_max = 0.838; }
	}

	if (tRangeModel == "low_t,high_t")
	{
		if (ds.name == "2.76TeV")	{ ds.t_min = 0.450; ds.t_dip = 0.616; ds.t_bmp = 0.750; ds.t_max = 0.990; }
		if (ds.name == "7TeV")		{ ds.t_min = 0.368; ds.t_dip = 0.530; ds.t_bmp = 0.694; ds.t_max = 0.878; }
		if (ds.name == "8TeV")		{ ds.t_min = 0.375; ds.t_dip = 0.518; ds.t_bmp = 0.687; ds.t_max = 0.900; }
		if (ds.name == "13TeV")		{ ds.t_min = 0.338; ds.t_dip = 0.470; ds.t_bmp = 0.638; ds.t_max = 0.858; }

		if (ds.name == "1.96TeV")	{ ds.t_min = 0.443; ds.t_dip = 0.611; ds.t_bmp = 0.745; ds.t_max = 0.838; }
	}
}

//----------------------------------------------------------------------------------------------------

int n_parameters = 0;

void BuildFitFunction(const string fitModel, Dataset &ds)
{
	printf("\n>> BuildFitFunction(%s)\n", fitModel.c_str());

	if (fitModel == "e012+e012")
	{
		const double t1 = (ds.t_dip + ds.t_min) / 2.;
		const double t2 = ds.t_bmp;

		n_parameters = 6;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f) + [5]*(x-%.3f)*(x-%.3f))", t1, t1, t1, t2, t2, t2);
		ds.ff = new TF1("ff", buf);

		if (ds.name == "2.76TeV") ds.ff->SetParameters(-4.35, -27.3, -147., -4.58, -16.8, -160.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-4.45, -30.4, -65.5, -3.62, 0.34, -21.7);
		if (ds.name == "8TeV") ds.ff->SetParameters(-4.27, -41.8, -297., -3.50, 0.59, -22.4);
		if (ds.name == "13TeV") ds.ff->SetParameters(-4.00, -44.7, -140., -3.03, 0.1, -24.2);

		return;
	}

	// e012_x.y+e012
	if (fitModel.find("e012_") == 0 && fitModel.find("e012", 8) == 9)
	{
		double f = atof(fitModel.substr(5, 3).c_str());

		printf("f = %.2f\n", f);

		const double t1 = ds.t_min * (1. - f) + ds.t_dip * f;
		const double t1_ref = (ds.t_min + ds.t_dip) / 2.;
		const double t2 = ds.t_bmp;

		printf("t1 = %.3f\n", t1);

		n_parameters = 6;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f) + [5]*(x-%.3f)*(x-%.3f))", t1, t1, t1, t2, t2, t2);
		ds.ff = new TF1("ff", buf);

		// these are the parameters at t1_ref
		if (ds.name == "2.76TeV") ds.ff->SetParameters(-4.35, -27.3, -147., -4.58, -16.8, -160.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-4.45, -30.4, -65.5, -3.62, 0.34, -21.7);
		if (ds.name == "8TeV") ds.ff->SetParameters(-4.27, -41.8, -297., -3.50, 0.59, -22.4);
		if (ds.name == "13TeV") ds.ff->SetParameters(-4.00, -44.7, -140., -3.03, 0.1, -24.2);

		// do evolution from t1_ref to t1
		const double p0 = ds.ff->GetParameter(0);
		const double p1 = ds.ff->GetParameter(1);
		const double p2 = ds.ff->GetParameter(2);

		const double pp2 = p2;
		const double pp1 = p1 - 2.*p2 * (t1_ref - t1);
		const double pp0 = p0 - p1*t1_ref + pp1*t1 + p2 * (t1_ref*t1_ref - t1*t1);

		ds.ff->SetParameter(0, pp0);
		ds.ff->SetParameter(1, pp1);
		ds.ff->SetParameter(2, pp2);

		for (double t : { 0.40, 0.45, 0.50 })
		{
			printf("t=%.2f --> f=%.4f\n", t, ds.ff->Eval(t));
		}

		return;
	}

	if (fitModel == "e012+e0123")
	{
		const double t1 = (ds.t_dip + ds.t_min) / 2.;
		//const double t2 = (ds.t_max + ds.t_dip) / 2.;
		const double t2 = ds.t_bmp;

		n_parameters = 7;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f) + [5]*(x-%.3f)*(x-%.3f) + [6]*(x-%.3f)*(x-%.3f)*(x-%.3f))", t1, t1, t1, t2, t2, t2, t2, t2, t2);
		ds.ff = new TF1("ff", buf);

		if (ds.name == "2.76TeV") ds.ff->SetParameters(-3.78, -19.3, -0.3, -4.19, -1.68, -264., 0.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, 3.29, -24.3, 40.);
		if (ds.name == "8TeV") ds.ff->SetParameters(-3.97, -34.4, -46.5, -3.53, 0.174, -12.31, 40.);
		if (ds.name == "13TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, 3.29, -24.3, 40.);

		return;
	}

	if (fitModel == "e012+e023")
	{
		const double t1 = (ds.t_dip + ds.t_min) / 2.;
		//const double t2 = (ds.t_max + ds.t_dip) / 2.;
		const double t2 = ds.t_bmp;

		n_parameters = 6;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f)*(x-%.3f) + [5]*(x-%.3f)*(x-%.3f)*(x-%.3f))", t1, t1, t1, t2, t2, t2, t2, t2);
		ds.ff = new TF1("ff", buf);

		if (ds.name == "2.76TeV") ds.ff->SetParameters(-3.78, -19.3, -0.3, -4.19, -264., 0.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, -24.3, 40.);
		if (ds.name == "8TeV") ds.ff->SetParameters(-3.97, -34.4, -46.5, -3.53, -12.31, 40.);
		if (ds.name == "13TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, -24.3, 40.);

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
