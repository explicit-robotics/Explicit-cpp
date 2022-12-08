#include <iostream>

#include <Eigen/Dense>

#include "exp_utils.h"
#include "exp_math.h"
#include "exp_robots.h"
#include "exp_constants.h"

int main( )
{

	// The number of degrees of freedom of the robot
	int nq = 8;

	SnakeBot robot( 1, "first_robot" , nq, 1.0 , 1.0  );

	std::cout << robot.JointTwists << std::endl;

	Eigen::VectorXd q_arr( nq ); 
	q_arr << 1, 1, 1, 1, 1, 1, 1, 1;

	Eigen::MatrixXd tmp1 = robot.getForwardKinematics( q_arr );
	Eigen::MatrixXd tmp2 = robot.getForwardKinematics( q_arr, 1, TYPE_COM );

	std::cout << tmp1 << std::endl;
	std::cout << tmp2 << std::endl;

	Eigen::MatrixXd JS = robot.getSpatialJacobian( q_arr );	
	
	std::cout << JS << std::endl;

	Eigen::MatrixXd JH = robot.getHybridJacobian( q_arr );	
	
	std::cout << JH << std::endl;	
	std::cout << robot.M_Mat << std::endl;	

	Eigen::MatrixXd JB = robot.getBodyJacobian( q_arr, 3, TYPE_COM );	
	
	std::cout << JB << std::endl;	

	Eigen::MatrixXd M = robot.getMassMatrix( q_arr );	

	std::cout << M << std::endl;

	return 0;
}