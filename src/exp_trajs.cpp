#include <cmath>
#include <cassert>
#include <iostream>
#include <Eigen/Dense>

#include "exp_math.h"
#include "exp_utils.h"
#include "exp_trajs.h"
#include "exp_robots.h"
#include "exp_constants.h"

MinimumJerkTrajectory::MinimumJerkTrajectory( const int n, const Eigen::VectorXd &q0i, const Eigen::VectorXd &q0f, const double D, const double ti )
{
	// Some Assertions
	assert( n >= 1 && D > 0 && ti > 0 );

	// Set the Member Variables
	this->n   = n;
	this->ti  = ti;
	this->D   = D;
	this->q0i = q0i;
	this->q0f = q0f;

}

Eigen::VectorXd MinimumJerkTrajectory::getPosition( const double t )
{	
	assert( t >= 0 ); 
	if( t <= this->ti )
	{
		return this->q0i;
	}
	else if( t >= this->ti && t <= ( this->ti + this->D ) )
	{
		return this->q0i + ( this->q0f - this->q0i ) * ( 10.0 * pow( ( t - this->ti )/this->D, 3 ) -  15.0 * pow( ( t - this->ti )/this->D, 4 ) + 6.0 * pow( ( t - this->ti )/this->D, 5 ) );
	}
	else
	{
		return this->q0f;
	}
}

Eigen::VectorXd MinimumJerkTrajectory::getVelocity( const double t )
{
	assert( t >= 0 ); 
	if( this->ti <= t && t <= this->ti + this->D )
	{
		return 1.0/this->D * ( this->q0f - this->q0i ) * ( 30 * pow( ( t - this->ti )/this->D, 2 ) - 60 * pow( ( t - this->ti )/this->D, 3 ) + 30 * pow( ( t - this->ti )/this->D, 4 ) );
	}
	else
	{
		return Eigen::VectorXd::Zero( this->n );
	}
}

Eigen::VectorXd MinimumJerkTrajectory::getAcceleration( const double t )
{
	assert( t >= 0 ); 
	if( this->ti <= t && t <= this->ti + this->D )
	{
		return 1.0/( this->D * this->D ) * ( this->q0f - this->q0i ) * ( 60 * pow( ( t - this->ti )/this->D, 1 ) - 180 * pow( ( t - this->ti )/this->D, 2 ) + 120 * pow( ( t - this->ti )/this->D, 3 ) );
	}
	else
	{
		return Eigen::VectorXd::Zero( this->n );
	}

}



