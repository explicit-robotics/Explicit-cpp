#include <iostream>
#include <Eigen/Dense>

#include "exp_math.h"
#include "exp_utils.h"
#include "exp_robots.h"
#include "exp_constants.h"

bool isSymmetric( const Eigen::MatrixXd &M )
{
	int nr = M.rows( );
	int nc = M.cols( );

	// Rows and Columns must be the same size
	assert( nr == nc );

	for ( int i = 0; i < nr; i++ )
	{
		for( int j = 0; j < nr; j ++ )
		{
			if ( M( i, j ) != M( j, i ) )
			{
				return false;
			}
		}
	}
		
	return true;
}

bool isSkewSymmetric( const Eigen::MatrixXd &M )
{
	int nr = M.rows( );
	int nc = M.cols( );

	// Rows and Columns must be the same size
	assert( nr == nc );

	for ( int i = 0; i < nr; i++ )
	{
		for( int j = 0; j < nr; j ++ )
		{
			if ( M( i, j ) != -M( j, i ) )
			{
				return false;
			}
		}
	}
	return true;
}

