all: t_investigation t_investigation_new

t_investigation: t_investigation.cc
	g++ -g -O3 -Wall `root-config --libs` `root-config --cflags`\
		t_investigation.cc -o t_investigation

t_investigation_new: t_investigation_new.cc ../datasets.h ../stat.h
	g++ -g -O3 -Wall `root-config --libs` `root-config --cflags`\
		t_investigation_new.cc -o t_investigation_new
