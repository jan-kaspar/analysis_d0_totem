string topDir = "../";

string f = topDir + "do_fits.root";

string datasets[], d_labels[];
real d_sqrt_s[];
datasets.push("2.76TeV"); d_labels.push("$2.76\un{TeV}$"); d_sqrt_s.push(2.76);
datasets.push("7TeV"); d_labels.push("$7\un{TeV}$"); d_sqrt_s.push(7);
datasets.push("8TeV"); d_labels.push("$8\un{TeV}$"); d_sqrt_s.push(8);
datasets.push("13TeV"); d_labels.push("$13\un{TeV}$"); d_sqrt_s.push(13);

real sqrt_s_ext = 1.96;

string extModels[];
extModels.push("sqrt_s");
//extModels.push("log_sqrt_s");

string extFits[];
extFits.push("corr");
//extFits.push("uncorr");
