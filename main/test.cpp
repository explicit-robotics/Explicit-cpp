#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <cassert>  

#include "exp_utils.h"
#include "exp_math.h"
#include "exp_robots.h"
#include "exp_constants.h"

int main( )
{

	
	iiwa14 myLBR( 1, "iiwa1" );
	myLBR.init( );

	Eigen::VectorXd q_arr_iiwa ( myLBR.nq );
	q_arr_iiwa << 0.2, 0.1, 0.1, 0.1, 0.2, 0.3, 0.2;

	auto start = std::chrono::steady_clock::now( );
	auto end   = std::chrono::steady_clock::now( );

	start = std::chrono::steady_clock::now( );
	Eigen::MatrixXd H3 = myLBR.getForwardKinematics( q_arr_iiwa );
	Eigen::MatrixXd H4 = myLBR.getHybridJacobian( q_arr_iiwa );
	Eigen::MatrixXd H1 = myLBR.getMassMatrix( q_arr_iiwa );
	end = std::chrono::steady_clock::now( );

	std::cout << "Elapsed time for FK, HJ, M1 "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
        << " us" << std::endl;		

	start = std::chrono::steady_clock::now( );
	H3 = myLBR.getForwardKinematics( q_arr_iiwa );
	H4 = myLBR.getHybridJacobian( q_arr_iiwa );
	Eigen::MatrixXd H2 = myLBR.getMassMatrix2( q_arr_iiwa );
	end = std::chrono::steady_clock::now( );
	
	std::cout << "Elapsed time for FK, HJ, M2 "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
        << " us" << std::endl;		

	std::cout << H1 << std::endl;
	std::cout << "========================" << std::endl;
	std::cout << H2 << std::endl;
	std::cout << "========================" << std::endl;
	std::cout << H3 << std::endl;	
	std::cout << "========================" << std::endl;
	std::cout << H4 << std::endl;		

	

	return 0;
}