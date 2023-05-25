#include <iostream>
#include <chrono>
#include <Eigen/Dense>
#include <cassert>  
#include <fstream>


#include "exp_utils.h"
#include "exp_math.h"
#include "exp_robots.h"
#include "exp_constants.h"
#include "exp_trajs.h"


int main( )
{
	std::ofstream file1("./data/test1.txt");
	std::ofstream file2("./data/test2.txt");
	std::ofstream file3("./data/test3.txt");

	Eigen::VectorXd q0i = Eigen::VectorXd::Zero( 7 );
	q0i << 0.2, 0.1, 0.1, 0.1, 0.2, 0.3, 0.2;

	Eigen::VectorXd q0f = Eigen::VectorXd::Zero( 7 );
	q0f << 0.5, 0.4, 0.2, 0.8, 1.0, 0.7, 0.4;

	MinimumJerkTrajectory traj1 = MinimumJerkTrajectory( 7, q0i, q0f, 6, 2 );
	double dt = 0.01;

	for( int i = 0; i<= 1000; i++)
	{	
		file1 << traj1.getPosition( dt * i ) << std::endl;
		file2 << traj1.getVelocity( dt * i ) << std::endl;
		file3 << traj1.getAcceleration( dt * i ) << std::endl;
	}		
	
	return 0;
}