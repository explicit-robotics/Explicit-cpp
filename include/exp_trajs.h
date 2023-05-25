/*
 * EXPlicit - A robotics toolbox based on the product of exponentials formula.
 *
 * Copyright (c) 2022 MIT
 * Authors
 * 			Johannes Lachner  	<jlachner@mit.edu>	
 * 			Moses C. Nah 		<mosesnah@mit.edu>
 * 
 * Classes for Trajectories
 */

#ifndef EXP_TRAJECTORIES
#define EXP_TRAJECTORIES

#include <Eigen/Dense>


class MinimumJerkTrajectory
{
	public:
		int n;					// The number of degrees-of-freedom
		Eigen::VectorXd q0i;	// Final   Posture 
		Eigen::VectorXd q0f;	// Initial Posture
		double D;				// Duration of Movement
		double ti;				// Initial Time of the Movement

		MinimumJerkTrajectory( const int n, const Eigen::VectorXd &q0i, const Eigen::VectorXd &q0f, const double D, const double ti );
		~MinimumJerkTrajectory(  ){ };

		Eigen::VectorXd getPosition( const double t );		
		Eigen::VectorXd getVelocity( const double t );		
		Eigen::VectorXd getAcceleration( const double t );		

};

#endif // EXP_TRAJECTORIES