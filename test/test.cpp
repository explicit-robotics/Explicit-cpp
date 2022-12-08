#include <iostream>

#include <Eigen/Dense>

#include "exp_utils.h"
#include "exp_math.h"

int main( )
{
	Eigen::Matrix2d tmp { { 0, 1 }, {1 , 0} };

	std::cout << tmp << std::endl;

	return 0;
}