#include <iostream>

#include <Eigen/Dense>

#include "exp_utils.h"
#include "exp_math.h"
#include "exp_robots.h"
#include "exp_constants.h"

int main( )
{
	/* ************************************** */
	/* ************** SnakeBot ************** */
	/* ************************************** */

	// The number of degrees of freedom of the robot
	int nq = 8;

	SnakeBot robot( 1, "first_robot" , nq, 1.0 , 1.0  );

	// std::cout << robot.JointTwists << std::endl;

	Eigen::VectorXd q_arr( nq ); 
	q_arr << 1, 1, 1, 1, 1, 1, 1, 1;

	Eigen::MatrixXd tmp1 = robot.getForwardKinematics( q_arr );
	Eigen::MatrixXd tmp2 = robot.getForwardKinematics( q_arr, 1, TYPE_COM );

	// std::cout << tmp1 << std::endl;
	// std::cout << tmp2 << std::endl;

	Eigen::MatrixXd JS = robot.getSpatialJacobian( q_arr );	
	
	// std::cout << JS << std::endl;

	Eigen::MatrixXd JH = robot.getHybridJacobian( q_arr );	
	
	// std::cout << JH << std::endl;	
	// std::cout << robot.M_Mat << std::endl;	

	Eigen::MatrixXd JB = robot.getBodyJacobian( q_arr, 3, TYPE_COM );	
	
	// std::cout << JB << std::endl;	

	Eigen::MatrixXd M = robot.getMassMatrix( q_arr );	

	// std::cout << M << std::endl;

	/* ************************************** */
	/* *************** iiwa14 *************** */
	/* ************************************** */
	iiwa14 myLBR( 1, "iiwa1" );

	std::cout << std::endl << "Robot Name: ";
	std::cout << myLBR.Name << std::endl;

	std::cout << std::endl << "DOF: ";
	std::cout << myLBR.nq << std::endl;

	Eigen::VectorXd q_arr_iiwa ( myLBR.nq );
	q_arr_iiwa << 2, 2, 2, 2, 2, 2, 2;

	Eigen::MatrixXd H = myLBR.getForwardKinematics( q_arr_iiwa );

	std::cout << std::endl << "H: ";
	std::cout << H << std::endl;

	Eigen::MatrixXd J_sp = myLBR.getSpatialJacobian( q_arr_iiwa );	

	std::cout << std::endl << "J_s: ";
	std::cout << J_sp << std::endl;

	Eigen::MatrixXd J_hy = myLBR.getHybridJacobian( q_arr_iiwa );	

	std::cout << std::endl << "J_hy: ";
	std::cout << J_hy << std::endl;	

	Eigen::MatrixXd J_b = myLBR.getBodyJacobian( q_arr_iiwa, 7, TYPE_COM );	

	std::cout << std::endl << "J_b: ";
	std::cout << J_b << std::endl;		

	Eigen::MatrixXd M_iiwa = myLBR.getMassMatrix( q_arr_iiwa );	

	std::cout << std::endl << "M_iiwa: ";
	std::cout << M_iiwa << std::endl;	

	return 0;
}