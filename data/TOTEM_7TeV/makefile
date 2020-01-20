all: process_input re_normalisation

process_input: process_input.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` --std=c++17 -Wall -Werror \
		process_input.cc -o process_input

re_normalisation: re_normalisation.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` --std=c++17 -Wall -Werror \
		re_normalisation.cc -o re_normalisation
