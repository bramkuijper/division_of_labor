all : xfixed_response xreinforcedRT xreadhisto

xfixed_response : fixed_response_threshold.cpp
	g++ -Wall -O3 -o xfixed_response fixed_response_threshold.cpp -lgsl -lgslcblas

xreinforcedRT : reinforcedRT_ExpEnhPerf_stepsize.cpp
	g++ -Wall -O3 -o xreinforcedRT reinforcedRT_ExpEnhPerf_stepsize.cpp -fopenmp -lgsl -lgslcblas 

xreadhisto : read_histograms.cpp
	g++ -Wall -ggdb -O3 -o xreadhisto read_histograms.cpp -lm -lrt -lgsl -lgslcblas

clean :
	rm -rf xfixed_response
	rm -rf xreinforcedRT
	rm -rf xreadhisto

cleanout:
	rm -f data_work_alloc*.txt
	rm -f f_dist*.txt
	rm -f header*.txt
	rm -f lastgen*.txt
	rm -f allele_distrib*
	rm -f branch.txt
	rm -f threshold.txt
	rm -f thresholds*.txt
	rm -f ant_beh*.txt
	rm -f data_1gen*.txt
