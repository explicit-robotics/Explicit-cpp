#include <cmath>
#include <iostream>
#include <Eigen/Dense>

#include "exp_math.h"
#include "exp_utils.h"
#include "exp_robots.h"
#include "exp_constants.h"

/*************************************************************/
/********************* ROBOT PRIMITIVES **********************/
/*************************************************************/

RobotPrimitive::RobotPrimitive( const int ID, const char* Name )
{
	this->ID   = ID;
	this->Name = Name;
}

void RobotPrimitive::init( )
{
	
	// Initialize the Jacobian Matrices for the Robot
	this->JS = Eigen::MatrixXd::Zero( 6,  this->nq );
	this->JH = Eigen::MatrixXd::Zero( 6,  this->nq );
	this->JB = Eigen::MatrixXd::Zero( 6,  this->nq );

	// Initialize the Joint Twists and A Matrices
	this->JointTwists = Eigen::MatrixXd::Zero( 	   		  6,  this->nq );
	this->A_Mat1 	  = Eigen::MatrixXd::Zero( 		      6,  this->nq );
	this->A_Mat2 	  = Eigen::MatrixXd::Zero( 6 * this->nq,  this->nq );

	// Initialize the Mass Matrices
	this->M_Mat1 	  = Eigen::MatrixXd::Zero( 6 * this->nq,  6 * this->nq );
	this->M_Mat2 	  = Eigen::MatrixXd::Zero( 6 * this->nq,  6 * this->nq );

	// Initialize the L Matrix used for the Mass Matrix Computation
	this->L_Mat = Eigen::MatrixXd::Identity( 6 * this->nq, 6 * this->nq );	

	// Once Initialization Complete, set the joint twists and generalized mass matrix
	this->setJointTwists( );
	this->setGeneralizedMassMatrix( );

}

void RobotPrimitive::setJointTwists( )
{
	// The Joint Twist Matrix of the Robot: ( 6 x nq )
	// Initializing the Joint Twist Matrix, and A matrices
	Eigen::VectorXd Ai( 6, 1 );	

	// Iteration over the degrees of freedom.
	for( int i = 0; i < this->nq; i++ )
	{
		// Before calculation, normalize the axis direction in case if it wasn't
		this->AxisDirections.col( i ).normalize( );

		// Create the joint twist of the i-th joint with respect to the {i} frame
		if( this->JointTypes( i ) == REVOLUTE_JOINT )
		{	
			
			this->JointTwists.block< 3, 1 >( 0, i ) = -vec2SkewSym( this->AxisDirections.col( i ) ) * this->AxisOrigins.col( i );
			this->JointTwists.block< 3, 1 >( 3, i ) = this->AxisDirections.col( i );
		}
		else if( this->JointTypes( i ) == PRISMATIC_JOINT )
		{
			this->JointTwists.block< 3, 1 >( i, 0 ) = this->AxisDirections.col( i );
		}
		else 
		{	// Something wrong, raise assertion
			assert( false );
		}

		Ai = getInvAdjoint( this->H_COM_init.block< 4,4 >( 0,4*i ) ) * this->JointTwists.block< 6, 1 >( 0, i );
	
		this->A_Mat1.block<6,1>(0, i)	= Ai; 		// Stacking Ai horizontally
		this->A_Mat2.block<6,1>(6*i, i) = Ai;		// Stacking Ai horizontally with an offset

	}

}

void RobotPrimitive::setGeneralizedMassMatrix( )
{
	for( int i = 0; i < this->nq; i++)
	{
		// Generalized inertia matrix
		this->M_Mat1.block< 3, 3 >( 0, 6 * i 	 ) = this->Masses( i ) * Eigen::Matrix3d::Identity( 3, 3 );
		this->M_Mat1.block< 3, 3 >( 3, 6 * i + 3 ) = this->Inertias.block< 3, 3 >( 0, 3 * i );

		// Generalized inertia matrix, in diagonal form
		this->M_Mat2.block< 6, 6 >( 6*i, 6*i ) = this->M_Mat1.block< 6, 6 >( 0, 6*i );
	}
}

Eigen::Matrix4d RobotPrimitive::getForwardKinematics( const Eigen::VectorXd &q_arr )
{
	// Assertion that q_arr length must be same with nq
	assert( this->nq == q_arr.size( ) );

	// Use the H_arrs for the calculation
	for( int i = 0; i < this->nq; i ++ )
	{
		this->H_arrs.block< 4, 4 >( 0, 4 * i ) = getExpSE3( this->JointTwists.block< 3, 1 >( 3, i ), this->JointTwists.block< 3, 1 >( 0, i ), q_arr( i ) );
	}

	// The Return the Forward Kinematics of the end-effector of robot
	return getExpProd( this->H_arrs ) * this->H_init.block< 4, 4 >( 0, 4*this->nq );

} 

Eigen::Matrix4d RobotPrimitive::getForwardKinematics( const Eigen::VectorXd &q_arr, const int bodyID, const int type )
{
	// Assertion that q_arr length must be same with nq
	assert( this->nq == q_arr.size( ) );

	// Assertion that bodyID must be smaller or equal to nq
	assert( bodyID <= this->nq && bodyID >= 1  );

	// Assertion that type is either JOINT or COM (Center of Mass)
	assert( type == TYPE_JOINT || type == TYPE_COM );

	// Assign the H_arr matrix
	for( int i = 0; i < bodyID; i ++ )
	{
		this->H_arrs.block< 4, 4 >( 0, 4 * i ) = getExpSE3( this->JointTwists.block< 3, 1 >( 3, i ), this->JointTwists.block< 3, 1 >( 0, i ), q_arr( i ) );
	}

	// If either JOINT or COM
	if( type == TYPE_JOINT )
	{
		return getExpProd( this->H_arrs.block( 0, 0, 4, 4 * bodyID ) ) * this->H_init.block< 4, 4 >( 0, 4 * ( bodyID - 1 ) );
	}
	else if ( type == TYPE_COM )
	{	
		return getExpProd( this->H_arrs.block( 0, 0, 4, 4 * bodyID ) ) * this->H_COM_init.block< 4, 4 >( 0, 4 * ( bodyID - 1 ) );
	}
	else 
	{
		assert( false );
	}

}

Eigen::MatrixXd RobotPrimitive::getSpatialJacobian( const Eigen::VectorXd &q_arr )		
{
	// Assertion that q_arr length must be same with nq
	assert( this->nq == q_arr.size( ) );

	// The First column of the Spatial Jacobian is the first joint twist array
	// Hence, saving the first joint twist values to the Spatial Jacobian
	this->JS = Eigen::MatrixXd::Zero( 6,  this->nq );
	this->JS.col( 0 ) = this->JointTwists.col( 0 );

	Eigen::MatrixXd H_tmp = Eigen::MatrixXd::Identity( 4, 4 );

	for( int i = 0; i < ( this->nq - 1 ); i ++ )
	{
		H_tmp *= getExpSE3( this->JointTwists.block< 3, 1 >( 3, i ), this->JointTwists.block< 3, 1 >( 0, i ), q_arr( i ) );
		this->JS.col( i + 1 ) = getAdjoint( H_tmp ) * this->JointTwists.col( i + 1 );
	}

	return this->JS;
}

Eigen::MatrixXd RobotPrimitive::getSpatialJacobian( const Eigen::VectorXd &q_arr, const int bodyID )		
{
	// Assertion that q_arr length must be same with nq
	assert( this->nq == q_arr.size( ) );

	// Assertion that bodyID must be smaller or equal to nq
	assert( bodyID <= this->nq && bodyID >= 1  );

	// The First column of the Spatial Jacobian is the first joint twist array
	// Initialization
	this->JS = Eigen::MatrixXd::Zero( 6,  this->nq );
	this->JS.col( 0 ) = this->JointTwists.col( 0 );

	Eigen::MatrixXd H_tmp = Eigen::MatrixXd::Identity( 4, 4 );

	for( int i = 0; i < bodyID-1; i ++ )
	{
		H_tmp *= getExpSE3( this->JointTwists.block< 3, 1 >( 3, i ), this->JointTwists.block< 3, 1 >( 0, i ), q_arr( i ) );
		this->JS.col( i + 1 ) = getAdjoint( H_tmp ) * this->JointTwists.col( i + 1 );
	}

	return this->JS;
}


Eigen::MatrixXd RobotPrimitive::getHybridJacobian( const Eigen::VectorXd &q_arr  )		
{
	// Assertion that q_arr length must be same with nq
	assert( this->nq == q_arr.size( ) );

	// Getting the end-effector's Hybrid Jacobian.
	// Note that "Hybrid Jacobian" is by definition, the Jacobian of the end-effector coordinate frame
	// with respect to the base coordinate frame.

	// First, get the Spatial Jacobian
	// Get the position of the end-effector
	Eigen::Matrix4d H_EE = this->getForwardKinematics( q_arr );

	Eigen::MatrixXd A = Eigen::MatrixXd::Identity( 6, 6 );

	A.block< 3, 3 >( 0, 3 ) = -vec2SkewSym( H_EE.block< 3, 1 >( 0, 3 ) );
	
	return A * this->getSpatialJacobian( q_arr );

}

Eigen::MatrixXd RobotPrimitive::getBodyJacobian( const Eigen::VectorXd &q_arr, const int bodyID, const int type )
{
	return getInvAdjoint( this->getForwardKinematics( q_arr, bodyID, type )  ) * this->getSpatialJacobian( q_arr, bodyID );
}

Eigen::MatrixXd RobotPrimitive::getMassMatrix( const Eigen::VectorXd &q_arr )
{

	Eigen::MatrixXd Ai( 6, 1 );
	Eigen::MatrixXd Hi( 4, 4 );
	Eigen::MatrixXd Hj( 4, 4 );

	// Initialize
	this->L_Mat = Eigen::MatrixXd::Identity( 6 * this->nq, 6 * this->nq );	

	for( int i = 1; i < this->nq; i++ )
	{
		Ai = this->A_Mat1.block< 6, 1 >( 0, i );               
		Hi = this->H_COM_init.block< 4, 4 >( 0, 4 * i );
		Hj = this->H_COM_init.block< 4, 4 >( 0, 4 * (i-1) );

		this->L_Mat.block( 6*i, 0, 6, 6*i ) = getAdjoint( getExpSE3( -Ai.block<3,1>(3,0), -Ai.block<3,1>(0,0), q_arr( i ) ) * Hi.inverse( ) * Hj  ) *  L_Mat.block( 6*( i - 1 ), 0, 6, 6*i  );
	}

	return this->A_Mat2.transpose( ) * this->L_Mat.transpose( ) * this->M_Mat2 *  this->L_Mat * this->A_Mat2;
}


/*************************************************************/
/********************** KUKA LBR iiwa 14 *********************/
/*************************************************************/
iiwa14::iiwa14( const int ID, const char* name ) : RobotPrimitive( ID , name )
{

	// ======================================================== //
	// =================== BASIC PROPERTIES  ================== //
	// ======================================================== //

	// There are 7 joints for iiwa14
	this->nq = 7;

	// ======================================================== //
	// ================= JOINT LIMITS  ======================== //
	// ======================================================== //
	
	// Joint limits, including some safe distance to physical limit
	this->q_max = Eigen::VectorXd( this->nq );
	this->q_min = Eigen::VectorXd( this->nq );
	this->q_max <<  163,  113,  163,  115,  160,  110,  165;
	this->q_min << -163, -113, -163, -115, -160, -110, -165;

	// Changing degrees to radian
	this->q_max *= M_PI/180;
	this->q_min *= M_PI/180;

	// Velocity and acceleration limits
	this->dq_max  =  150 * Eigen::VectorXd::Ones( this->nq );
	this->dq_min  = -150 * Eigen::VectorXd::Ones( this->nq );
	this->ddq_max =  300 * Eigen::VectorXd::Ones( this->nq );
	this->ddq_min = -300 * Eigen::VectorXd::Ones( this->nq );

	// ======================================================== //
	// ================= GEOMETRIC PROPERTIES  ================ //
	// ======================================================== //

	// The joint types are all 1 (i.e., revolute joints)
	this->JointTypes = REVOLUTE_JOINT * Eigen::VectorXd::Ones( this->nq );

	// Initialization of the Axis Origins Matrix
	this->AxisOrigins = Eigen::MatrixXd::Zero( 3, this->nq );

	this->AxisOrigins.col( 0 ) = Eigen::Vector3d( 0.0,      0.0, 152.5e-3 );  						  
	this->AxisOrigins.col( 1 ) = this->AxisOrigins.col( 0 ) + Eigen::Vector3d( 0.0, -13.0e-3, 207.5e-3 );
	this->AxisOrigins.col( 2 ) = this->AxisOrigins.col( 1 ) + Eigen::Vector3d( 0.0, +13.0e-3, 232.5e-3 );
	this->AxisOrigins.col( 3 ) = this->AxisOrigins.col( 2 ) + Eigen::Vector3d( 0.0, +11.0e-3, 187.5e-3 );
	this->AxisOrigins.col( 4 ) = this->AxisOrigins.col( 3 ) + Eigen::Vector3d( 0.0, -11.0e-3, 212.5e-3 );
	this->AxisOrigins.col( 5 ) = this->AxisOrigins.col( 4 ) + Eigen::Vector3d( 0.0, -62.0e-3, 187.5e-3 );
	this->AxisOrigins.col( 6 ) = this->AxisOrigins.col( 5 ) + Eigen::Vector3d( 0.0, +62.0e-3,  79.6e-3 );

	// Axis Directions of the robot, Defining it for joints 1~7
	this->AxisDirections = Eigen::MatrixXd::Zero( 3, this->nq );
	this->AxisDirections.col( 0 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );
	this->AxisDirections.col( 1 ) = Eigen::Vector3d( 0.0,  1.0,  0.0 );
	this->AxisDirections.col( 2 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );
	this->AxisDirections.col( 3 ) = Eigen::Vector3d( 0.0, -1.0,  0.0 );
	this->AxisDirections.col( 4 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );
	this->AxisDirections.col( 5 ) = Eigen::Vector3d( 0.0,  1.0,  0.0 );
	this->AxisDirections.col( 6 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );

	// ======================================================== //
	// ================== INERTIAL PROPERTIES  ================ //
	// ======================================================== //
	// Define the masses of the links
	this->Masses = Eigen::VectorXd( this->nq );
	this->Masses << 6.404, 7.89, 2.54, 4.82, 1.76, 2.50, 0.42;
	
	// Eigen::Vector3d com[ nq ];
	this->COM = Eigen::MatrixXd( 3, this->nq );
	this->COM.col( 0 ) = this->AxisOrigins.col( 0 ) + Eigen::Vector3d( 0.0,  -14.0e-3,  102.0e-3);  
	this->COM.col( 1 ) = this->AxisOrigins.col( 1 ) + Eigen::Vector3d( 0.0,   16.0e-3,   64.0e-3);  
	this->COM.col( 2 ) = this->AxisOrigins.col( 2 ) + Eigen::Vector3d( 0.0,   19.0e-3,   98.0e-3);  
	this->COM.col( 3 ) = this->AxisOrigins.col( 3 ) + Eigen::Vector3d( 0.0,  -20.0e-3,   86.0e-3);
	this->COM.col( 4 ) = this->AxisOrigins.col( 4 ) + Eigen::Vector3d( 0.0,  -13.0e-3,   66.0e-3);  
	this->COM.col( 5 ) = this->AxisOrigins.col( 5 ) + Eigen::Vector3d( 0.0,   60.0e-3,   16.0e-3);  
	this->COM.col( 6 ) = this->AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0,    0.0e-3,   11.0e-3); 

	// Define the Principal axes/moments of inertia of the link
	this->Inertias = Eigen::MatrixXd( 3, 3 * this->nq );	
	this->Inertias.block< 3, 3 >( 0, 3 * 0 ) << 0.069, 0.000, 0.00,
										  		0.000, 0.071, 0.00,
										  		0.000, 0.000, 0.02;

	this->Inertias.block< 3, 3 >( 0, 3 * 1 ) << 0.08, 0.00, 0.00,
			  			            	 		0.00, 0.08, 0.00,
								 		  		0.00, 0.00, 0.01;

	this->Inertias.block< 3, 3 >( 0, 3 * 2 ) << 0.02, 0.00, 0.00,
				            	 		  		0.00, 0.02, 0.00,
								 		  		0.00, 0.00, 0.06;
	
	this->Inertias.block< 3, 3 >( 0, 3 * 3 ) << 0.04, 0.00, 0.00, 
				            	 		  		0.00, 0.03, 0.00,
								 		  		0.00, 0.00, 0.01;
	
	this->Inertias.block< 3, 3 >( 0, 3 * 4 ) << 0.01, 0.00, 0.00, 
				            	 		  		0.00, 0.01, 0.00,
								 		  		0.00, 0.00, 0.01;

	this->Inertias.block< 3, 3 >( 0, 3 * 5 ) << 0.007, 0.000, 0.000, 
				        		 		  		0.000, 0.006, 0.000,
								 		  		0.000, 0.000, 0.005;

	this->Inertias.block< 3, 3 >( 0, 3 * 6 ) << 0.0003, 0.0000, 0.0000, 
				            	 		  		0.0000, 0.0003, 0.0000,
								 		  		0.0000, 0.0000, 0.0005;

	// Concatenate matrices for transformations and Generalized inertia matrix
	this->H_init 	 = Eigen::MatrixXd( 4, 4 * ( this->nq + 1 ) );
	this->H_COM_init = Eigen::MatrixXd( 4, 4 * ( this->nq     ) );
	this->H_arrs	 = Eigen::MatrixXd( 4, 4 * ( this->nq     ) );

	for( int i = 0; i < this->nq; i++)
	{
		// The Rotation matrix is an identity matrix
		this->H_init.block< 4, 4 >( 0, 4 * i ) = Eigen::Matrix4d::Identity( 4, 4 );

		// The translational part of H_init
		this->H_init.block< 3, 1 >( 0, 4 * i + 3 ) = this->AxisOrigins.col( i );

		this->H_COM_init.block< 4, 4 >( 0, 4 * i ) = Eigen::Matrix4d::Identity( 4, 4 );
		this->H_COM_init.block< 3, 1 >( 0, 4 * i + 3 ) = this->COM.col( i );

		// Initializing the H_arrs matrix with identity matrices
		this->H_arrs.block< 4, 4 >( 0, 4 * i ) = Eigen::Matrix4d::Identity( 4, 4 );
	}

	// Setup the final H_init Matrix for Media Flange Touch
	Eigen::Vector3d FlangePos = AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0, 0.0, 0.071 );
	// For all other Media Flanges:                    
	// Eigen::Vector3d FlangePos = this->AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0, 0.0, 0.0314 );                    
	this->H_init.block< 4, 4 >( 0, 4 * this->nq	 	) = Eigen::Matrix4d::Identity( 4, 4 );
	this->H_init.block< 3, 1 >( 0, 4 * this->nq + 3 ) = FlangePos;

}
