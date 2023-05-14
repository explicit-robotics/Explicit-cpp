#include <iostream>
#include <chrono>
#include <Eigen/Dense>

#include "exp_utils.h"
#include "exp_math.h"
#include "exp_robots.h"
#include "exp_constants.h"

int main( )
{
	Eigen::Matrix2d tmp { { 0, 1 }, {1 , 0} };

	std::cout << tmp  << std::endl;

	iiwa14 myLBR( 1, "iiwa1" );

	Eigen::VectorXd q_arr_iiwa ( myLBR.nq );
	q_arr_iiwa << 0.2, 0.1, 0.1, 0.1, 0.2, 0.3, 0.2;

	auto start = std::chrono::steady_clock::now( );
	Eigen::MatrixXd H1 = myLBR.getMassMatrix( q_arr_iiwa );
	auto end = std::chrono::steady_clock::now( );

    std::cout << "Elapsed time for M1: "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
        << " us" << std::endl;

	start = std::chrono::steady_clock::now( );
	Eigen::MatrixXd H2 = myLBR.getMassMatrix2( q_arr_iiwa );
	end = std::chrono::steady_clock::now( );

	std::cout << "Elapsed time for M2: "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
        << " us" << std::endl;



	

	return 0;
}