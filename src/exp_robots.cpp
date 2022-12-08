#include <cmath>
#include <iostream>
#include <Eigen/Dense>

#include "exp_math.h"
#include "exp_utils.h"
#include "exp_robots.h"

RobotPrimitive::RobotPrimitive( )
{
	
}

RobotPrimitive::RobotPrimitive( int ID, const char* Name, const Eigen::VectorXd &JointTypes, const Eigen::MatrixXd &AxisOrigins, const Eigen::MatrixXd &AxisDirections )
{
	this->ID   = ID;
	this->Name = Name;

	// Length of JointTypes, AxisOrigins and AxisDirections must be the same.
	// The length of Joint Types   : ( 1 x nq )
	// The shape of Axis Origins   : ( 3 x nq )
	// The shape of Axis Directions: ( 3 x nq )
	assert( JointTypes.size( ) == AxisOrigins.cols( ) && JointTypes.size( ) == AxisDirections.cols( ) );

	// If passed assertion, save the nq size 
	this->nq = JointTypes.size( );

	// Saving the passed parameters to member variables
	this->JointTypes  	 = JointTypes;
	this->AxisOrigins 	 = AxisOrigins;
	this->AxisDirections = AxisDirections;
}

void RobotPrimitive::setJointTwists( )
{
	// The Joint Twist Matrix of the Robot: ( 6 x nq )
	// Initializing the Joint Twist Matrix that will be later saved as a member variable
	Eigen::MatrixXd JointTwists( 6, this->nq );
	JointTwists = Eigen::MatrixXd::Zero( 6, this->nq );

	for( int i = 0; i < this->nq; i++ )
	{
		// If Rotational Joint
		
		if( this->JointTypes( i ) == 1 )
		{	
			JointTwists.block< 3, 1 >( 0, i ) = -vec2SkewSym( this->AxisDirections.col( i ) ) * this->AxisOrigins.col( i );
			JointTwists.block< 3, 1 >( 3, i ) = this->AxisDirections.col( i );
		}
		else
		{
			JointTwists.block< 3, 1 >( i, 0 ) = this->AxisDirections.col( i );
		}
		
	}

	this->JointTwists = JointTwists;
}

Eigen::Matrix4d RobotPrimitive::getForwardKinematics( const Eigen::VectorXd &q_arr )
{
	// Currently, we are interested only at the end-effector of the robot
	// [2022.12.07] [TODO] [Moses C. Nah] 
	// Defining more "getForwardKinematics" to get the forward kineamtics of any points.

	// [2022.12.07] [TODO] [Moses C. Nah] 
	// Currently, the H_init of the robot is not initially defined, hence defined manually 
	// in the function. We should be aware of this and fix soon.

	// Assertion that q_arr length must be same with nq
	assert( this->nq == q_arr.size( ) );

	// The end-effector H_init matrix
	Eigen::MatrixXd H_init( 4, 4 );
	H_init = Eigen::MatrixXd::Identity( 4, 4 );

	// Produce the ( 4 x ( 4 x nq ) ) 2D Matrix for the forward kinematics
	Eigen::MatrixXd H_arr( 4, 4 * this->nq );

	for( int i = 0; i < this->nq; i ++ )
	{
		H_arr.block< 4, 4 >( 0, 4 * i ) = getExpSE3( this->JointTwists.block< 3, 1 >( 0, i ), this->JointTwists.block< 3, 1 >( 3, i ), q_arr( i ) );
	}

	// The Return the Forward Kinematics of the robot
	return getExpProd( H_arr ) * H_init;

} 

Eigen::MatrixXd RobotPrimitive::getSpatialJacobian( const Eigen::VectorXd &q_arr );		
{
	// Assertion that q_arr length must be same with nq
	assert( this->nq == q_arr.size( ) );

	// Initialize the Spatial Jacobian Matrix, which is a ( 6 x nq ) matrix
	Eigen::MatrixXd JS( 6, this->nq );

	// The First column of the Spatial Jacobian is the first joint twist array
	// Hence, saving the first joint twist values to the Spatial Jacobian
	JS.block< 6, 1 >( 0, 0 ) = this->JointTwists.col( 0 )

	Eigen::MatrixXd H_tmp( 4, 4 );
	H_tmp = Eigen::MatrixXd::Identity( 4, 4 );

	for( int i = 0; i < ( this->nq - 1 ); i ++ )
	{
		H_tmp *= getExpSE3( this->JointTwists.block< 3, 1 >( 0, i ), this->JointTwists.block< 3, 1 >( 3, i ), q_arr( i ) );

		JS.block< 6, 1 >( 0, i+1 ) = getAdjoint( H_tmp ) * this->JointTwists( i+1 );
	}

	return JS;
}

Eigen::MatrixXd RobotPrimitive::getHybridJacobian( const Eigen::VectorXd &q_arr )
{

}	

Eigen::MatrixXd RobotPrimitive::getBodyJacobian( const Eigen::VectorXd &q_arr )
{
	
}	

