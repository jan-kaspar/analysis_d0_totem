import root;
import pad_layout;

xSizeDef = 10cm;

void DrawPoint(real W, real R, real R_unc, pen p)
{
	draw((W, R), mCi+2pt+p);
	draw((W, R-R_unc)--(W, R+R_unc), p+1pt);
}

NewPad("$\sqrt s\ung{TeV}$", "$\displaystyle R = {\d\si/\d t|_{\rm bump} \over \d\si/\d t|_{\rm dip}}$");

DrawPoint(1.96, 1.190, 0.226, red);

DrawPoint(7, 1.6, 0.3, blue);
DrawPoint(8, 2.0, 0.3, magenta);
DrawPoint(13, 1.775, 0.05325, blue);

real a = (1.775 - 1.6) / (13 - 7);
real b = 1.775 - a*13.;

real W_min = 1, W_max = 15;

draw((W_min, a*W_min + b)--(W_max, a*W_max + b), blue+dashed);

limits((W_min, 0.8), (W_max, 2.4), Crop);

yaxis(XEquals(2.76, false), heavygreen);
