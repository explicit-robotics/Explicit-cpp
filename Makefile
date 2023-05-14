
.PHONY: all

all:

	g++ -Iinclude -I./ -MMD -MP -Wall -g -std=c++17 -c ./src/exp_math.cpp -o ./build/exp_math.o
	g++ -Iinclude -I./ -MMD -MP -Wall -g -std=c++17 -c ./src/exp_robots.cpp -o ./build/exp_robots.o
	g++ -Iinclude -I./ -MMD -MP -Wall -g -std=c++17 -c ./src/exp_utils.cpp -o ./build/exp_utils.o

	g++ -Iinclude -I./ -MMD -MP -Wall -g -std=c++17 -c ./main/test.cpp 	 -o ./build/test.o
	g++ -Iinclude -I./ -MMD -MP -Wall -g -std=c++17 -c ./main/timing.cpp -o ./build/timing.o

	g++ -Iinclude -I./ -MMD -MP -Wall -g -std=c++17 -g ./build/test.o   ./build/exp_math.o ./build/exp_robots.o ./build/exp_utils.o -o ./bin/test
	g++ -Iinclude -I./ -MMD -MP -Wall -g -std=c++17 -g ./build/timing.o ./build/exp_math.o ./build/exp_robots.o ./build/exp_utils.o -o ./bin/timing
	