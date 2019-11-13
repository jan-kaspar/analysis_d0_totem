#include "TF1.h"
#include "TMatrixDSym.h"

#include <vector>
#include <string>

using namespace std;

//----------------------------------------------------------------------------------------------------

struct ParameterConstraint
{
	double mean, sigma;
};

struct Dataset
{
	double sqrt_s;
	string name;
	string f_in;
	double t_min, t_max, t_dip, t_bmp;
	TMatrixDSym m_t_value_cov;
	TF1 *ff = NULL;
	map<unsigned int, ParameterConstraint> constraints;

	void AddConstraint(unsigned int pi, double mean, double sigma)
	{
		constraints.insert({pi, {mean, sigma}});
	}
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
	printf("\n>> BuildTRanges(%s, %s)\n", tRangeModel.c_str(), ds.name.c_str());

	if (tRangeModel == "bootstrap")
	{
		if (ds.name == "2.76TeV")	{ ds.t_min = 0.350; ds.t_dip = 0.616; ds.t_bmp = 0.700; ds.t_max = 0.800; }
		if (ds.name == "7TeV")		{ ds.t_min = 0.250; ds.t_dip = 0.530; ds.t_bmp = 0.694; ds.t_max = 0.950; }
		if (ds.name == "8TeV")		{ ds.t_min = 0.250; ds.t_dip = 0.518; ds.t_bmp = 0.687; ds.t_max = 1.100; }
		if (ds.name == "13TeV")		{ ds.t_min = 0.320; ds.t_dip = 0.470; ds.t_bmp = 0.638; ds.t_max = 0.870; }
	}

	if (tRangeModel == "minimal")
	{
		if (ds.name == "2.76TeV")	{ ds.t_min = 0.530; ds.t_dip = 0.616; ds.t_bmp = 0.790; ds.t_max = 0.850; }
		if (ds.name == "7TeV")		{ ds.t_min = 0.442; ds.t_dip = 0.529; ds.t_bmp = 0.693; ds.t_max = 0.780; }
		if (ds.name == "8TeV")		{ ds.t_min = 0.433; ds.t_dip = 0.521; ds.t_bmp = 0.701; ds.t_max = 0.800; }
		if (ds.name == "13TeV")		{ ds.t_min = 0.400; ds.t_dip = 0.468; ds.t_bmp = 0.638; ds.t_max = 0.719; }

		if (ds.name == "1.96TeV")	{ ds.t_min = 0.548; ds.t_dip = 0.653; ds.t_bmp = 0.819; ds.t_max = 0.875;
			double cov_data[16] = {
				2.904E-04, 0., 0., 0.,
				0., 8.527E-05, 0., 0.,
				0., 0., 5.252E-04, 0.,
				0., 0., 0., 9.390E-04
			};
			ds.m_t_value_cov.ResizeTo(4, 4);
			ds.m_t_value_cov.SetMatrixArray(cov_data);
		}
	}

	if (tRangeModel == "low_t,high_t")
	{
		if (ds.name == "2.76TeV")	{ ds.t_min = 0.448; ds.t_dip = 0.616; ds.t_bmp = 0.790; ds.t_max = 0.950; }
		if (ds.name == "7TeV")		{ ds.t_min = 0.368; ds.t_dip = 0.529; ds.t_bmp = 0.693; ds.t_max = 0.878; }
		if (ds.name == "8TeV")		{ ds.t_min = 0.374; ds.t_dip = 0.521; ds.t_bmp = 0.701; ds.t_max = 0.909; }
		if (ds.name == "13TeV")		{ ds.t_min = 0.338; ds.t_dip = 0.468; ds.t_bmp = 0.638; ds.t_max = 0.863; }

		if (ds.name == "1.96TeV")	{ ds.t_min = 0.467; ds.t_dip = 0.653; ds.t_bmp = 0.819; ds.t_max = 0.971;
			double cov_data[16] = {
				6.656E-05, 0., 0., 0.,
				0., 8.527E-05, 0., 0.,
				0., 0., 5.252E-04, 0.,
				0., 0., 0., 4.593E-03
			};
			ds.m_t_value_cov.ResizeTo(4, 4);
			ds.m_t_value_cov.SetMatrixArray(cov_data);
		}
	}
}

//----------------------------------------------------------------------------------------------------

void BuildFitFunction(const string fitModel, Dataset &ds)
{
	printf("\n>> BuildFitFunction(%s, %s)\n", fitModel.c_str(), ds.name.c_str());

	if (fitModel == "bootstrap")
	{
		if (ds.name == "2.76TeV")
		{
			ds.ff = new TF1("ff", "exp([0] + [1]*x) + exp([2] + [3]*x + [4]*x*x)");
			ds.ff->SetParameters(6.75, -19.3, -116., 319., -228.);
		}

		if (ds.name == "7TeV")
		{
			ds.ff = new TF1("ff", "exp([0] + [1]*x + [2]*x*x) + exp([3] + [4]*x + [5]*x*x + [6]*x*x*x)");
			ds.ff->SetParameters(10.5, -37.7, 15.8, -33.9, 113., -139., 54.7);
		}

		if (ds.name == "8TeV")
		{
			ds.ff = new TF1("ff", "exp([0] + [1]*x + [2]*x*x) + exp([3] + [4]*x + [5]*x*x + [6]*x*x*x)");
			ds.ff->SetParameters(10.5, -37.7, 15.8, -33.9, 113., -139., 54.7);
		}

		if (ds.name == "13TeV")
		{
			ds.ff = new TF1("ff", "exp([0] + [1]*x + [2]*x*x + [3]*x*x*x) + exp([4] + [5]*x + [6]*x*x + [7]*x*x*x)");
			//ds.ff->SetParameters(4.57, -5.24, -31.3, 1., -25.2, 89.4, -117., 49.2);
			ds.ff->SetParameters(10.9, -65.9, 160., -198., -21.4, 73.5, -95.1, +39.2);
		}

		return;
	}

	if (fitModel.find("e012+e012") == 0)
	{
		double f = 0.5; // default

		const auto pos = fitModel.find(":f=");
		if (pos != string::npos)
		{
			f = atof(fitModel.substr(pos+3).c_str());
			printf("f = %.2f\n", f);
		}

		const double t1 = ds.t_min * (1. - f) + ds.t_dip * f;
		const double t1_def = (ds.t_min + ds.t_dip) / 2.;
		const double t2 = ds.t_bmp;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f) + [5]*(x-%.3f)*(x-%.3f))", t1, t1, t1, t2, t2, t2);
		ds.ff = new TF1("ff", buf);

		// these are the parameters at t1_def
		if (ds.name == "2.76TeV") ds.ff->SetParameters(-4.35, -27.3, -147., -4.58, -16.8, -160.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-4.45, -30.4, -65.5, -3.62, 0.34, -21.7);
		if (ds.name == "8TeV") ds.ff->SetParameters(-4.27, -41.8, -297., -3.50, 0.59, -22.4);
		if (ds.name == "13TeV") ds.ff->SetParameters(-4.00, -44.7, -140., -3.03, 0.1, -24.2);

		// do evolution from t1_def to t1
		const double p0 = ds.ff->GetParameter(0);
		const double p1 = ds.ff->GetParameter(1);
		const double p2 = ds.ff->GetParameter(2);

		const double pp2 = p2;
		const double pp1 = p1 - 2.*p2 * (t1_def - t1);
		const double pp0 = p0 - p1*t1_def + pp1*t1 + p2 * (t1_def*t1_def - t1*t1);

		ds.ff->SetParameter(0, pp0);
		ds.ff->SetParameter(1, pp1);
		ds.ff->SetParameter(2, pp2);

		//for (double t : { 0.40, 0.45, 0.50 })
		//	printf("t=%.2f --> f=%.4f\n", t, ds.ff->Eval(t));

		return;
	}

	if (fitModel == "e012+e0123")
	{
		const double t1 = (ds.t_dip + ds.t_min) / 2.;
		const double t2 = ds.t_bmp;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f) + [5]*(x-%.3f)*(x-%.3f) + [6]*(x-%.3f)*(x-%.3f)*(x-%.3f))", t1, t1, t1, t2, t2, t2, t2, t2, t2);
		ds.ff = new TF1("ff", buf);

		if (ds.name == "2.76TeV") ds.ff->SetParameters(-3.78, -19.3, -0.3, -4.19, -1.68, -264., 1.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, 3.29, -24.3, 40.);
		if (ds.name == "8TeV") ds.ff->SetParameters(-3.97, -34.4, -46.5, -3.53, 0.174, -12.31, 40.);
		if (ds.name == "13TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, 3.29, -24.3, 40.);

		return;
	}

	if (fitModel.find("e012+e023") == 0)
	{
		double f = 0.5; // default

		const auto pos = fitModel.find(":f=");
		if (pos != string::npos)
		{
			f = atof(fitModel.substr(pos+3).c_str());
			printf("f = %.2f\n", f);
		}

		const double t1 = ds.t_min * (1. - f) + ds.t_dip * f;
		const double t1_def = (ds.t_min + ds.t_dip) / 2.;
		const double t2 = ds.t_bmp;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f)*(x-%.3f) + [5]*(x-%.3f)*(x-%.3f)*(x-%.3f))", t1, t1, t1, t2, t2, t2, t2, t2);
		ds.ff = new TF1("ff", buf);

		// these are the parameters at f=0.5
		if (ds.name == "2.76TeV") { ds.ff->SetParameters(-3.47, -18.8, -14.9, -4., -26., 40.); ds.AddConstraint(4, -26., 3.); ds.AddConstraint(5, 40., 15.); }
		if (ds.name == "7TeV") ds.ff->SetParameters(-3.31, -25.3, -6.27, -3.62, -22.6, 39.0);
		if (ds.name == "8TeV") ds.ff->SetParameters(-3.24, -30.6, -56.7, -3.50, -22.8, 42.0);
		if (ds.name == "13TeV") ds.ff->SetParameters(-2.70, -33.2, -67.9, -3.03, -20.6, 43.1);

		// do evolution from t1_def to t1
		const double p0 = ds.ff->GetParameter(0);
		const double p1 = ds.ff->GetParameter(1);
		const double p2 = ds.ff->GetParameter(2);

		const double pp2 = p2;
		const double pp1 = p1 - 2.*p2 * (t1_def - t1);
		const double pp0 = p0 - p1*t1_def + pp1*t1 + p2 * (t1_def*t1_def - t1*t1);

		ds.ff->SetParameter(0, pp0);
		ds.ff->SetParameter(1, pp1);
		ds.ff->SetParameter(2, pp2);

		//for (double t : { 0.40, 0.45, 0.50 })
		//	printf("t=%.2f --> f=%.4f\n", t, ds.ff->Eval(t));

		return;
	}

	if (fitModel == "e012+e02")
	{
		const double t1 = (ds.t_dip + ds.t_min) / 2.;
		//const double t2 = (ds.t_max + ds.t_dip) / 2.;
		const double t2 = ds.t_bmp;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f) + [2]*(x-%.3f)*(x-%.3f)) + exp([3] + [4]*(x-%.3f)*(x-%.3f))", t1, t1, t1, t2, t2);
		ds.ff = new TF1("ff", buf);

		if (ds.name == "2.76TeV") ds.ff->SetParameters(-3.78, -19.3, -0.3, -4.19, -264.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, -27.3);
		if (ds.name == "8TeV") ds.ff->SetParameters(-3.97, -34.4, -46.5, -3.53, -15.31);
		if (ds.name == "13TeV") ds.ff->SetParameters(-3.6, -36.3, -82.4, -3.13, -27.3);

		return;
	}

	if (fitModel == "e01+e012")
	{
		const double t1 = (ds.t_dip + ds.t_min) / 2.;
		const double t2 = ds.t_bmp;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f)) + exp([2] + [3]*(x-%.3f) + [4]*(x-%.3f)*(x-%.3f))", t1, t2, t2, t2);
		ds.ff = new TF1("ff", buf);

		if (ds.name == "2.76TeV") ds.ff->SetParameters(-4.35, -27.3, -4.58, -16.8, -160.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-4.45, -30.4, -3.62, 0.34, -21.7);
		if (ds.name == "8TeV") ds.ff->SetParameters(-4.27, -41.8, -3.50, 0.59, -22.4);
		if (ds.name == "13TeV") ds.ff->SetParameters(-4.00, -44.7, -3.03, 0.1, -24.2);

		return;
	}

	if (fitModel == "e01+e02")
	{
		const double t1 = (ds.t_dip + ds.t_min) / 2.;
		const double t2 = ds.t_bmp;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f)) + exp([2] + [4]*(x-%.3f)*(x-%.3f))", t1, t2, t2);
		ds.ff = new TF1("ff", buf);

		if (ds.name == "2.76TeV") ds.ff->SetParameters(-4.35, -27.3, -4.58, -160.);
		if (ds.name == "7TeV") ds.ff->SetParameters(-4.45, -30.4, -3.62, -21.7);
		if (ds.name == "8TeV") ds.ff->SetParameters(-4.27, -41.8, -3.50, -22.4);
		if (ds.name == "13TeV") ds.ff->SetParameters(-4.00, -44.7, -3.03, -24.2);

		return;
	}

	if (fitModel.find("e01+e023") == 0)
	{
		double f = 0.5; // default

		const auto pos = fitModel.find(":f=");
		if (pos != string::npos)
		{
			f = atof(fitModel.substr(pos+3).c_str());
			printf("f = %.2f\n", f);
		}

		const double t1 = ds.t_min * (1. - f) + ds.t_dip * f;
		const double t1_def = (ds.t_dip + ds.t_min) / 2.;
		const double t2 = ds.t_bmp;

		char buf[500];
		sprintf(buf, "exp([0] + [1]*(x-%.3f)) + exp([2] + [3]*(x-%.3f)*(x-%.3f) + [4]*(x-%.3f)*(x-%.3f)*(x-%.3f))", t1, t2, t2, t2, t2, t2);
		ds.ff = new TF1("ff", buf);

		// these are the parameters for f=0.5
		if (ds.name == "2.76TeV") { ds.ff->SetParameters(-3.5, -17., -4.0, -20., 40.); ds.AddConstraint(2, -4.0, 0.3); ds.AddConstraint(3, -20., 10.); ds.AddConstraint(4, 40., 15.); }
		if (ds.name == "7TeV") ds.ff->SetParameters(-4.08, -21.6, -3.64, -19.4, 99.4);
		if (ds.name == "8TeV") ds.ff->SetParameters(-3.98, -26.1, -3.51, -23.5, 87.5);
		if (ds.name == "13TeV") ds.ff->SetParameters(-3.76, -35.1, -3.03, -20.2, 55.0);

		// do evolution from t1_def to t1
		const double p0 = ds.ff->GetParameter(0);
		const double p1 = ds.ff->GetParameter(1);

		const double pp1 = p1;
		const double pp0 = p0 - p1*t1_def + pp1*t1;

		ds.ff->SetParameter(0, pp0);
		ds.ff->SetParameter(1, pp1);

		//for (double t : { 0.40, 0.45, 0.50 })
		//	printf("t=%.2f --> f=%.4f\n", t, ds.ff->Eval(t));

		return;
	}

	if (fitModel == "e01-int-e01")
	{
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
