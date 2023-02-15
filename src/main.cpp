#include <iostream>

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

	Eigen::MatrixXd H = myLBR.getForwardKinematics( q_arr_iiwa );

	std::cout << std::endl << "H: " << std::endl;
	std::cout << H << std::endl;

	std::cout << std::endl << "Axis Origin " << std::endl;
	std::cout << myLBR.AxisOrigins << std::endl;	

	Eigen::MatrixXd J_sp = myLBR.getSpatialJacobian( q_arr_iiwa );	

	std::cout << std::endl << "J_s: " << std::endl;
	std::cout << J_sp << std::endl;

	Eigen::MatrixXd J_hy = myLBR.getHybridJacobian( q_arr_iiwa );	

	std::cout << std::endl << "J_hy: " << std::endl;
	std::cout << J_hy << std::endl;	

	Eigen::MatrixXd J_b = myLBR.getBodyJacobian( q_arr_iiwa, 7, TYPE_COM );	

	std::cout << std::endl << "J_b: " << std::endl;
	std::cout << J_b << std::endl;		

	Eigen::MatrixXd M_iiwa = myLBR.getMassMatrix( q_arr_iiwa );	

	std::cout << std::endl << "M_iiwa: " << std::endl;
	std::cout << M_iiwa << std::endl;	

	Eigen::VectorXd tau ( myLBR.nq );
	tau << 20, 20, 20, 20, 20, 20, 20;
	Eigen::VectorXd qDot ( myLBR.nq );
	qDot << 10, 10, 10, 10, 10, 10, 10;
	Eigen::MatrixXd Minv  = M_iiwa;
	
	Eigen::VectorXd tauSat = myLBR.addIIWALimits( &myLBR, q_arr_iiwa, qDot, Minv, tau, 0.005 );
	std::cout << std::endl << "tauSat: " << std::endl;
	std::cout << tauSat << std::endl;

	return 0;
}