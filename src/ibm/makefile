all : xreinforcedRT xreadhisto

xreinforcedRT : reinforcedRT_ExpEnhPerf_stepsize.cpp
	g++ -Wall -O3 -o xreinforcedRT reinforcedRT_ExpEnhPerf_stepsize.cpp -fopenmp -lgsl -lgslcblas 

xreadhisto : read_histograms.cpp
	g++ -Wall -ggdb -O3 -o xreadhisto read_histograms.cpp -lm -lrt -lgsl -lgslcblas

clean :
	rm -rf xreinforcedRT
	rm -rf xreadhisto

cleanout:
	rm -f data_work_alloc*.txt
	rm -f f_dist*.txt
	rm -f header*.txt
	rm -f lastgen*.txt
	rm -f allele_distrib*
	rm -f branch.txt
