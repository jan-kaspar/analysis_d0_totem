string topDir = "../";

string f = topDir + "do_fits.root";

string datasets[], d_labels[];
real d_sqrt_s[];
//datasets.push("8TeV"); d_labels.push("$8\un{TeV}$"); d_sqrt_s.push(8);

//datasets.push("8TeV-b1-45b"); d_labels.push("$8\un{TeV}$"); d_sqrt_s.push(8);
//datasets.push("8TeV-b1-45t"); d_labels.push("$8\un{TeV}$"); d_sqrt_s.push(8);

//datasets.push("8TeV-b1"); d_labels.push("$8\un{TeV}$"); d_sqrt_s.push(8);

//datasets.push("8TeV-b2"); d_labels.push("$8\un{TeV}$"); d_sqrt_s.push(8);
//datasets.push("8TeV-b3"); d_labels.push("$8\un{TeV}$"); d_sqrt_s.push(8);

datasets.push("8TeV-bt1"); d_labels.push("$8\un{TeV}$"); d_sqrt_s.push(8);

real sqrt_s_ext = 1.96;

string extModels[];
extModels.push("sqrt_s");
extModels.push("log_sqrt_s");

string extFits[];
extFits.push("corr");
//extFits.push("uncorr");
