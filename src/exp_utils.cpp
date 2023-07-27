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

double getMaxValue( const Eigen::VectorXd &myVector )
{
    double auxMax = myVector[ 0 ];

    for ( int i = 0; i < myVector.size(); i++ )
    {
        if ( myVector[ i ] > auxMax )
            auxMax = myVector[ i ];
    }
    return auxMax;
}

double getMinValue( const Eigen::VectorXd &myVector )
{
    double auxMin = myVector[ 0 ];
	 
    for ( int i = 0; i < myVector.size(); i++ )
    {
        if ( myVector[ i ] < auxMin )
            auxMin = myVector[ i ];
    }
    return auxMin;
}

