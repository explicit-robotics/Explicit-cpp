#include <cmath>
#include <iostream>
#include <Eigen/Dense>

#include "exp_math.h"
#include "exp_utils.h"
#include "exp_robots.h"
#include "exp_constants.h"

Eigen::Matrix3d vec2SkewSym( const Eigen::Vector3d &v )
{
	Eigen::Matrix3d M { 
		{ 	    0, -v( 2 ),  v( 1 ) },
		{  v( 2 ),       0, -v( 0 ) },
		{ -v( 1 ),  v( 0 ),       0 }
	  };

	return M;
}

Eigen::Vector3d skewSym2Vec( const Eigen::Matrix3d &M )
{
	// Check whether the given matrix is a skew-symmetric matrix
	assert( isSkewSymmetric( M ) );

	Eigen::Vector3d v( -M( 1, 2 ), M( 0, 2 ), -M( 0, 1 ) );
		
	return v;
}

Eigen::Matrix3d getExpSO3( const Eigen::Vector3d &w, const double theta )
{
	// Assertion that 
	// (1) w should be a non-zero vector 
	// (2) w should be a unit vector
	assert( !w.isZero( 0 ) && w.norm( ) == 1. );

	Eigen::Matrix3d w_hat = vec2SkewSym( w );

	// The Rodriguez's Formula
	return Eigen::Matrix3d::Identity( 3, 3 ) + sin( theta ) * w_hat + ( 1 - cos( theta ) ) * w_hat * w_hat;
}

Eigen::Matrix4d getExpSE3( const Eigen::Vector3d &w, const Eigen::Vector3d &v, const double theta )
{
	Eigen::Matrix4d H = Eigen::Matrix4d::Identity( 4, 4 );
	H.block< 3, 3 >( 0, 0 ) = getExpSO3( w, theta );

	// Check whether the w is a zero vector ( i.e., prismatic), if not, then revolute 
	// The input value of isZero method is the tolerance to identify the zero element
	if( w.isZero( 0 ) ) 
	{
		H.block< 3, 1 >( 0, 3 ) = v * theta;
	}
	else
	{
		H.block< 3, 1 >( 0, 3 ) = ( Eigen::Matrix3d::Identity( 3, 3 ) - getExpSO3( w, theta ) ) * ( vec2SkewSym( w ) * v );
	}
	
	return H;

}

Eigen::MatrixXd getAdjoint( const Eigen::Matrix4d &H )
{
	Eigen::MatrixXd Adj = Eigen::MatrixXd::Zero( 6, 6 );

	Adj.block< 3, 3 >( 0, 0 ) = H.block< 3, 3 >( 0, 0 );
	Adj.block< 3, 3 >( 3, 3 ) = H.block< 3, 3 >( 0, 0 );
	Adj.block< 3, 3 >( 0, 3 ) = vec2SkewSym( H.block< 3, 1 >( 0, 3 )  ) * H.block< 3, 3 >( 0, 0 );

	return Adj;

}

Eigen::MatrixXd getInvAdjoint( const Eigen::Matrix4d &H )
{
	Eigen::MatrixXd Adj = Eigen::MatrixXd::Zero( 6, 6 );

	Adj.block< 3, 3 >( 0, 0 ) =  H.block< 3, 3 >( 0, 0 ).transpose( );
	Adj.block< 3, 3 >( 3, 3 ) =  H.block< 3, 3 >( 0, 0 ).transpose( );
	Adj.block< 3, 3 >( 0, 3 ) = -H.block< 3, 3 >( 0, 0 ).transpose( ) * vec2SkewSym( H.block< 3, 1 >( 0, 3 )  );

	return Adj;

}

Eigen::MatrixXd getExpProd( const Eigen::MatrixXd &H )
{	
	// Given a T : ( 4 x ( 4 x nq ) ) 2D array, doing the exponential product
	// Check whether the given H Matrix's column is divided by 4.
	assert( H.rows( ) == 4 && H.cols( ) % 4 == 0 );

	// Get the number of H matrix for iteration
	int n = H.cols( ) / 4;

	Eigen::MatrixXd H_calc( 4, 4 );
	H_calc = Eigen::MatrixXd::Identity( 4, 4 );

	for( int i = 0; i < n; i ++ )
	{
		H_calc *= H.block< 4, 4 >( 0, 4*i );
	}

	return H_calc;

}
