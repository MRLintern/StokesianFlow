CXXFLAGS = $(shell pkg-config --cflags eigen3) -O2

all : main

clean :
	rm -f *~ *.dat main
