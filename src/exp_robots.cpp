#include <cmath>
#include <iostream>
#include <Eigen/Dense>
#include <vector>

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


	// Initialize the JS, JH, JB to avoid the re-construction of the matrix 
	Eigen::MatrixXd JS( 6, this->nq );
	Eigen::MatrixXd JH( 6, this->nq );
	Eigen::MatrixXd JB( 6, this->nq );

	this->JS = JS;
	this->JH = JH;
	this->JB = JB;

	// Initialize the L Matrix used for the 
	Eigen::MatrixXd L_Mat( 6 * this->nq, 6 * this->nq );
	L_Mat = Eigen::MatrixXd::Identity( 6 * this->nq, 6 * this->nq );	
	this->L_Mat = L_Mat;

	Eigen::MatrixXd JointTwists( 6, this->nq );
	JointTwists = Eigen::MatrixXd::Zero( 6, this->nq );
	this->JointTwists = JointTwists;

	Eigen::MatrixXd A_Mat1( 		   6,  this->nq );
	Eigen::MatrixXd A_Mat2( 6 * this->nq,  this->nq );

	A_Mat1 = Eigen::MatrixXd::Zero( 		   6,  this->nq );
	A_Mat2 = Eigen::MatrixXd::Zero( 6 * this->nq,  this->nq );

	this->A_Mat1 = A_Mat1;
	this->A_Mat2 = A_Mat2;

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
		this->A_Mat2.block<6,1>(6*i, i) = Ai;			// Stacking Ai horizontally with an offset

	}

}

void RobotPrimitive::setGeneralizedMassMatrix( )
{

	Eigen::MatrixXd M_Mat1( 6, 6 * this->nq );
	Eigen::MatrixXd M_Mat2( 6 * this->nq, 6 * this->nq );

	for( int i = 0; i < this->nq; i++)
	{
		// Generalized inertia matrix
		M_Mat1.block< 3, 3 >( 0, 6 * i 	) 	 = this->Masses ( i ) * Eigen::Matrix3d::Identity( 3, 3 );
		M_Mat1.block< 3, 3 >( 3, 6 * i + 3 ) = this->Inertias.block< 3, 3 >( 0, 3 * i );

		// Generalized inertia matrix, in diagonal form
		M_Mat2.block< 6, 6 >( 6*i, 6*i ) = M_Mat1.block< 6, 6 >( 0, 6*i );
	}

	this->M_Mat1 = M_Mat1;
	this->M_Mat2 = M_Mat2;
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
	this->JS.col( 0 ) = this->JointTwists.col( 0 );

	Eigen::MatrixXd H_tmp( 4, 4 );
	H_tmp = Eigen::MatrixXd::Identity( 4, 4 );

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
	// Hence, saving the first joint twist values to the Spatial Jacobian
	this->JS.col( 0 ) = this->JointTwists.col( 0 );

	Eigen::MatrixXd H_tmp( 4, 4 );
	H_tmp = Eigen::MatrixXd::Identity( 4, 4 );

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
	// Note that "Hybrid Jacobian" is by definition, the Jacobian of the end-effector

	// First, get the Spatial Jacobian
	// Get the position of the end-effector
	Eigen::Matrix4d H_EE = this->getForwardKinematics( q_arr );

	Eigen::MatrixXd A;
	A = Eigen::MatrixXd::Identity( 6, 6 );

	A.block< 3, 3 >( 0, 3 ) = -vec2SkewSym( H_EE.block< 3, 1 >( 0, 3 ) );
	
	return A * this->getSpatialJacobian( q_arr );

}

Eigen::MatrixXd RobotPrimitive::getBodyJacobian( const Eigen::VectorXd &q_arr, const int bodyID, const int type )
{
	return getInvAdjoint( this->getForwardKinematics( q_arr, bodyID, type )  ) * this->getSpatialJacobian( q_arr, bodyID );
}

Eigen::MatrixXd RobotPrimitive::getMassMatrix( const Eigen::VectorXd &q_arr )
{
	Eigen::MatrixXd M( this->nq, this->nq );
	M = Eigen::MatrixXd::Zero( this->nq, this->nq );

	for( int i = 0; i < this->nq; i++ )
	{
		M += this->getBodyJacobian( q_arr, i + 1, TYPE_COM ).transpose( ) * this->M_Mat1.block< 6, 6 >( 0, 6 * i ) * this->getBodyJacobian( q_arr, i+ 1 , TYPE_COM );
	} 

	return M;
}

Eigen::MatrixXd RobotPrimitive::getMassMatrix2( const Eigen::VectorXd &q_arr )
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
	Eigen::VectorXd q_max = Eigen::VectorXd::Zero( this->nq );
	Eigen::VectorXd q_min = Eigen::VectorXd::Zero( this->nq );

	// q_max/q_min in degrees
	// Currently, q_min is -q_max, but the values can change, thus manually defining the values.
	q_max <<  163,  113,  163,  115,  160,  110,  165;
	q_min << -163, -113, -163, -115, -160, -110, -165;

	// Changing degrees to radian
	q_max *= M_PI/180;
	q_min *= M_PI/180;

	// Velocity and acceleration limits
	Eigen::VectorXd dq_max  =  150 * Eigen::VectorXd::Ones( this->nq );
	Eigen::VectorXd dq_min  = -150 * Eigen::VectorXd::Ones( this->nq );
	Eigen::VectorXd ddq_max =  300 * Eigen::VectorXd::Ones( this->nq );
	Eigen::VectorXd ddq_min = -300 * Eigen::VectorXd::Ones( this->nq );

	// Saving all defined limits
	this->q_max   =   q_max;
	this->q_min   =   q_min;
	this->dq_max  =  dq_max;
	this->dq_min  =  dq_min;		
	this->ddq_max = ddq_max;
	this->ddq_min = ddq_min;	

	// ======================================================== //
	// ================= GEOMETRIC PROPERTIES  ================ //
	// ======================================================== //

	// The joint types are all 1 (i.e., revolute joints)
	Eigen::VectorXd JointTypes = REVOLUTE_JOINT * Eigen::VectorXd::Ones( this->nq );

	// Define the joint axes positions and directions 
	// These are relative displacements, hence should add the values cumulatively 
	// [2023.02.15] [Moses C. Nah] It might be good to have cumulative sum function 
	//							   To simplify the code.
	Eigen::MatrixXd AxisOrigins( 3, this->nq );
	AxisOrigins.col( 0 ) = 						  Eigen::Vector3d( 0.0,      0.0, 152.5e-3 );  
	AxisOrigins.col( 1 ) = AxisOrigins.col( 0 ) + Eigen::Vector3d( 0.0, -13.0e-3, 207.5e-3 );
	AxisOrigins.col( 2 ) = AxisOrigins.col( 1 ) + Eigen::Vector3d( 0.0, +13.0e-3, 232.5e-3 );
	AxisOrigins.col( 3 ) = AxisOrigins.col( 2 ) + Eigen::Vector3d( 0.0, +11.0e-3, 187.5e-3 );
	AxisOrigins.col( 4 ) = AxisOrigins.col( 3 ) + Eigen::Vector3d( 0.0, -11.0e-3, 212.5e-3 );
	AxisOrigins.col( 5 ) = AxisOrigins.col( 4 ) + Eigen::Vector3d( 0.0, -62.0e-3, 187.5e-3 );
	AxisOrigins.col( 6 ) = AxisOrigins.col( 5 ) + Eigen::Vector3d( 0.0, +62.0e-3,  79.6e-3 );

	// Axis Directions of the robot, Defining it for joints 1~7
	Eigen::MatrixXd AxisDirections( 3, this->nq );
	AxisDirections.col( 0 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );
	AxisDirections.col( 1 ) = Eigen::Vector3d( 0.0,  1.0,  0.0 );
	AxisDirections.col( 2 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );
	AxisDirections.col( 3 ) = Eigen::Vector3d( 0.0, -1.0,  0.0 );
	AxisDirections.col( 4 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );
	AxisDirections.col( 5 ) = Eigen::Vector3d( 0.0,  1.0,  0.0 );
	AxisDirections.col( 6 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );

	// ======================================================== //
	// ================== INERTIAL PROPERTIES  ================ //
	// ======================================================== //
	// Define the masses of the links
	Eigen::VectorXd Masses( this->nq );	
	Masses << 6.404, 7.89, 2.54, 4.82, 1.76, 2.50, 0.42;
	
	// Eigen::Vector3d com[ nq ];
	Eigen::MatrixXd COM( 3, this->nq );
	COM.col( 0 ) = AxisOrigins.col( 0 ) + Eigen::Vector3d( 0.0,  -14.0e-3,  102.0e-3);  
	COM.col( 1 ) = AxisOrigins.col( 1 ) + Eigen::Vector3d( 0.0,   16.0e-3,   64.0e-3);  
	COM.col( 2 ) = AxisOrigins.col( 2 ) + Eigen::Vector3d( 0.0,   19.0e-3,   98.0e-3);  
	COM.col( 3 ) = AxisOrigins.col( 3 ) + Eigen::Vector3d( 0.0,  -20.0e-3,   86.0e-3);
	COM.col( 4 ) = AxisOrigins.col( 4 ) + Eigen::Vector3d( 0.0,  -13.0e-3,   66.0e-3);  
	COM.col( 5 ) = AxisOrigins.col( 5 ) + Eigen::Vector3d( 0.0,   60.0e-3,   16.0e-3);  
	COM.col( 6 ) = AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0,    0.0e-3,   11.0e-3); 

	// Define the Principal axes/moments of inertia of the link
	Eigen::MatrixXd Inertias( 3, 3 * this->nq );	
	Inertias.block< 3, 3 >( 0, 3 * 0 ) << 0.069, 0.000, 0.00,
										  0.000, 0.071, 0.00,
										  0.000, 0.000, 0.02;

	Inertias.block< 3, 3 >( 0, 3 * 1 ) << 0.08, 0.00, 0.00,
				            	 		  0.00, 0.08, 0.00,
								 		  0.00, 0.00, 0.01;

	Inertias.block< 3, 3 >( 0, 3 * 2 ) << 0.02, 0.00, 0.00,
				            	 		  0.00, 0.02, 0.00,
								 		  0.00, 0.00, 0.06;
	
	Inertias.block< 3, 3 >( 0, 3 * 3 ) << 0.04, 0.00, 0.00, 
				            	 		  0.00, 0.03, 0.00,
								 		  0.00, 0.00, 0.01;
	
	Inertias.block< 3, 3 >( 0, 3 * 4 ) << 0.01, 0.00, 0.00, 
				            	 		  0.00, 0.01, 0.00,
								 		  0.00, 0.00, 0.01;

	Inertias.block< 3, 3 >( 0, 3 * 5 ) << 0.007, 0.000, 0.000, 
				        		 		  0.000, 0.006, 0.000,
								 		  0.000, 0.000, 0.005;

	Inertias.block< 3, 3 >( 0, 3 * 6 ) << 0.0003, 0.0000, 0.0000, 
				            	 		  0.0000, 0.0003, 0.0000,
								 		  0.0000, 0.0000, 0.0005;

	// Concatenate matrices for transformations and Generalized inertia matrix
	Eigen::MatrixXd H_init( 4, 4 * ( this->nq + 1 ) );
	Eigen::MatrixXd H_COM_init( 4, 4 * ( this->nq     ) );
	Eigen::MatrixXd H_arrs( 4, 4 * ( this->nq     ) );

	for( int i = 0; i < this->nq; i++)
	{
		// The Rotation matrix is an identity matrix
		H_init.block< 4, 4 >( 0, 4 * i ) = Eigen::Matrix4d::Identity( 4, 4 );

		// The translational part of H_init
		H_init.block< 3, 1 >( 0, 4 * i + 3 ) = AxisOrigins.col( i );

		H_COM_init.block< 4, 4 >( 0, 4 * i ) = Eigen::Matrix4d::Identity( 4, 4 );
		H_COM_init.block< 3, 1 >( 0, 4 * i + 3 ) = COM.col( i );

		H_arrs.block< 4, 4 >( 0, 4 * i ) = Eigen::Matrix4d::Identity( 4, 4 );
	}

	// Setup the final H_init Matrix for Media Flange Touch
	// Eigen::Vector3d FlangePos = AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0, 0.0, 0.071 );                    
	Eigen::Vector3d FlangePos = AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0, 0.0, 0.0314 );                    
	H_init.block< 4, 4 >( 0, 4 * this->nq 	  ) = Eigen::Matrix4d::Identity( 4, 4 );
	H_init.block< 3, 1 >( 0, 4 * this->nq + 3 ) = FlangePos;

	// Saving all of the calculated matrices
	this->H_init 		 = H_init;
	this->H_arrs 	 	 = H_arrs;
	this->H_COM_init 	 = H_COM_init;
	this->JointTypes  	 = JointTypes;
	this->AxisOrigins 	 = AxisOrigins;
	this->AxisDirections = AxisDirections;
	this->Masses		 = Masses;
	this->Inertias		 = Inertias;

}
