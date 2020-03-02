string topDir = "../";

string f = topDir + "do_fits.root";

string datasets[], d_labels[];
real d_sqrt_s[];
//datasets.push("13TeV"); d_labels.push(""); d_sqrt_s.push(13);
//datasets.push("13TeV-addUnc0.5"); d_labels.push(""); d_sqrt_s.push(13);
//datasets.push("13TeV-addUnc0.8"); d_labels.push(""); d_sqrt_s.push(13);
datasets.push("13TeV-addUnc0.8-att"); d_labels.push(""); d_sqrt_s.push(13);
datasets.push("13TeV-addUnc0.8-rebinA-att"); d_labels.push(""); d_sqrt_s.push(13);
datasets.push("13TeV-addUnc0.8-rebinB-att"); d_labels.push(""); d_sqrt_s.push(13);
//datasets.push("13TeV-addUnc1.0"); d_labels.push(""); d_sqrt_s.push(13);
//datasets.push("13TeV-addUnc1.5"); d_labels.push(""); d_sqrt_s.push(13);
//datasets.push("13TeV-addUnc-att"); d_labels.push("13TeV-addUnc-att"); d_sqrt_s.push(13);

real sqrt_s_ext = 1.96;

string extModels[];
extModels.push("sqrt_s");
extModels.push("log_sqrt_s");

string extFits[];
extFits.push("corr");
//extFits.push("uncorr");
