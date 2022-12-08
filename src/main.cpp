#include <iostream>

#include <Eigen/Dense>

#include "exp_utils.h"
#include "exp_math.h"
#include "exp_robots.h"

int main( )
{

	// The number of degrees of freedom of the robot
	int nq = 3;

	Eigen::VectorXd JointTypes( nq ); 
	Eigen::MatrixXd AxisOrigins( 3, nq ); 
	Eigen::MatrixXd AxisDirections( 3, nq ); 

	JointTypes << 1, 1, 1;
	
	AxisOrigins << 0, 1, 2, 
				   0, 0, 0,
				   0, 0, 0;	

	AxisDirections << 0, 0, 0, 
				      0, 0, 0,
				      1, 1, 1;

	RobotPrimitive robot( 1, "first_robot" , JointTypes, AxisOrigins, AxisDirections );

	robot.setJointTwists( );

	Eigen::VectorXd q_arr( 3 ); 
	q_arr << 1, 1, 1;

	Eigen::MatrixXd tmp = robot.getForwardKinematics( q_arr );

	std::cout << tmp << std::endl;

	return 0;
}