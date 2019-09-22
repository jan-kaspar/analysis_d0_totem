all: process_input

process_input: process_input.cc
	g++ `root-config --libs` -lMinuit `root-config --cflags` --std=c++17 -Wall -Werror \
		process_input.cc -o process_input
