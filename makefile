all: do_fits s_extrapolation

do_fits: do_fits.cc datasets.h command_line_tools.h
	g++ -g -O3 -Wall `root-config --libs` -lMinuit2 `root-config --cflags`\
		do_fits.cc -o do_fits

s_extrapolation: s_extrapolation.cc datasets.h command_line_tools.h
	g++ -g -O3 -Wall `root-config --libs` -lMinuit2 `root-config --cflags`\
		s_extrapolation.cc -o s_extrapolation
