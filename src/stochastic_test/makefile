all : xstochastic_sine

xstochastic_sine : stochastic_sine_test.cpp
	g++ -Wall -O3 -o xstochastic_sine stochastic_sine_test.cpp -lgsl -lgslcblas

clean :
	rm -rf xstochastic_sine
