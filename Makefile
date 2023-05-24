# This makefile is to run local .cpp files to test the library
# The local .cpp files are in the main directory.
# Simply run "make" in the terminal


.PHONY: all

all:

	mkdir -p bin 
	mkdir -p build

	g++ -Iinclude -I./ -MMD -MP -Wall -g -O3 -std=c++17 -c ./src/exp_math.cpp -o ./build/exp_math.o
	g++ -Iinclude -I./ -MMD -MP -Wall -g -O3 -std=c++17 -c ./src/exp_robots.cpp -o ./build/exp_robots.o
	g++ -Iinclude -I./ -MMD -MP -Wall -g -O3 -std=c++17 -c ./src/exp_utils.cpp -o ./build/exp_utils.o

	g++ -Iinclude -I./ -MMD -MP -Wall -g -O3 -std=c++17 -c ./main/test.cpp 	 -o ./build/test.o
	g++ -Iinclude -I./ -MMD -MP -Wall -g -O3 -std=c++17 -c ./main/timing.cpp -o ./build/timing.o

	g++ -Iinclude -I./ -MMD -MP -Wall -g -O3 -std=c++17 -g ./build/test.o   ./build/exp_math.o ./build/exp_robots.o ./build/exp_utils.o -o ./bin/test
	g++ -Iinclude -I./ -MMD -MP -Wall -g -O3 -std=c++17 -g ./build/timing.o ./build/exp_math.o ./build/exp_robots.o ./build/exp_utils.o -o ./bin/timing
	