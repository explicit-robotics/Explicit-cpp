#include <iostream>
#include <chrono>
#include <Eigen/Dense>

#include "exp_utils.h"
#include "exp_math.h"
#include "exp_robots.h"
#include "exp_constants.h"

int main( )
{

	/* ************************************** */
	/* *************** iiwa14 *************** */
	/* ************************************** */
	iiwa14 myLBR( 1, "iiwa1" );

	std::cout << std::endl << "Robot Name: ";
	std::cout << myLBR.Name << std::endl;

	std::cout << std::endl << "DOF: ";
	std::cout << myLBR.nq << std::endl;
	
	Eigen::VectorXd q_arr_iiwa ( myLBR.nq );
	q_arr_iiwa << 0.2, 0.1, 0.1, 0.1, 0.2, 0.3, 0.2;

	auto start = std::chrono::steady_clock::now( );
	Eigen::MatrixXd H = myLBR.getForwardKinematics( q_arr_iiwa );
	auto end = std::chrono::steady_clock::now( );

    std::cout << "Elapsed time for Forward Kinematics: "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
        << " us" << std::endl;

	start = std::chrono::steady_clock::now( );
	Eigen::MatrixXd J_hy = myLBR.getHybridJacobian( q_arr_iiwa );	
	end = std::chrono::steady_clock::now( );

    std::cout << "Elapsed time for Hybrid Jacobian: "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
        << " us" << std::endl;

	start = std::chrono::steady_clock::now( );
	Eigen::MatrixXd M_iiwa = myLBR.getMassMatrix( q_arr_iiwa );	
	end = std::chrono::steady_clock::now( );

	std::cout << "Elapsed time for Mass Matrix: "
        << std::chrono::duration_cast<std::chrono::microseconds>(end - start).count()
        << " us" << std::endl;

	return 0;
}