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

void RobotPrimitive::setJointTwists( )
{
	// [2022.12.07] [Moses C. Nah] [Backup]
	// The length of Joint Types   : ( 1 x nq )
	// The shape of Axis Origins   : ( 3 x nq )
	// The shape of Axis Directions: ( 3 x nq )
	// assert( JointTypes.size( ) == AxisOrigins.cols( ) && JointTypes.size( ) == AxisDirections.cols( ) );

	// The Joint Twist Matrix of the Robot: ( 6 x nq )
	// Initializing the Joint Twist Matrix that will be later saved as a member variable
	Eigen::MatrixXd JointTwists( 6, this->nq );
	JointTwists = Eigen::MatrixXd::Zero( 6, this->nq );

	for( int i = 0; i < this->nq; i++ )
	{
		// If Rotational Joint
		if( this->JointTypes( i ) == REVOLUTE_JOINT )
		{	
			JointTwists.block< 3, 1 >( 0, i ) = -vec2SkewSym( this->AxisDirections.col( i ) ) * this->AxisOrigins.col( i );
			JointTwists.block< 3, 1 >( 3, i ) = this->AxisDirections.col( i );
		}
		else if( this->JointTypes( i ) == PRISMATIC_JOINT )
		{
			JointTwists.block< 3, 1 >( i, 0 ) = this->AxisDirections.col( i );
		}
		else 
		{	// Something wrong, raise assertion
			assert( false );
		}
		
	}

	this->JointTwists = JointTwists;

	// Set the A Matrices
	Eigen::MatrixXd A_mat1( 		   6,  this->nq );
	Eigen::MatrixXd A_mat2( 6 * this->nq,  this->nq );

	Eigen::VectorXd Ai( 6, 1 );
	for( int i = 0; i < this->nq; i++ )
	{
		Ai = getInvAdjoint( this->H_COM_init.block< 4,4 >( 0,4*i ) ) * this->JointTwists.block< 6, 1 >( 0, i );
		A_mat1.block<6,1>(0, i)	  = Ai; 		// Stacking Ai horizontally
		A_mat2.block<6,1>(6*i, i) = Ai;			// Stacking Ai horizontally with an offset
	}
	this->A_mat1 = A_mat1;
	this->A_mat2 = A_mat2;
}

void RobotPrimitive::setGeneralizedMassMatrix( )
{
	
	this->M_mat1 = M_mat1;
	this->M_mat2 = M_mat2;
}

Eigen::Matrix4d RobotPrimitive::getForwardKinematics( const Eigen::VectorXd &q_arr )
{
	// Assertion that q_arr length must be same with nq
	assert( this->nq == q_arr.size( ) );

	// Produce the ( 4 x ( 4 x nq ) ) 2D Matrix for the forward kinematics
	Eigen::MatrixXd H_arr( 4, 4 * this->nq );

	for( int i = 0; i < this->nq; i ++ )
	{
		H_arr.block< 4, 4 >( 0, 4 * i ) = getExpSE3( this->JointTwists.block< 3, 1 >( 3, i ), this->JointTwists.block< 3, 1 >( 0, i ), q_arr( i ) );
	}

	// The Return the Forward Kinematics of the end-effector of robot
	return getExpProd( H_arr ) * this->H_init.block< 4, 4 >( 0, 4*this->nq );

} 

Eigen::Matrix4d RobotPrimitive::getForwardKinematics( const Eigen::VectorXd &q_arr, const int bodyID, const int type )
{
	// Assertion that q_arr length must be same with nq
	assert( this->nq == q_arr.size( ) );

	// Assertion that bodyID must be smaller or equal to nq
	assert( bodyID <= this->nq && bodyID >= 1  );

	// Assertion that type is either JOINT or COM (Center of Mass)
	assert( type == TYPE_JOINT || type == TYPE_COM );

	// Produce the ( 4 x ( 4 x bodyID ) ) 2D Matrix for the forward kinematics
	Eigen::MatrixXd H_arr( 4, 4 * bodyID );

	for( int i = 0; i < bodyID; i ++ )
	{
		H_arr.block< 4, 4 >( 0, 4 * i ) = getExpSE3( this->JointTwists.block< 3, 1 >( 3, i ), this->JointTwists.block< 3, 1 >( 0, i ), q_arr( i ) );
	}

	// If either JOINT or COM
	if( type == TYPE_JOINT )
	{
		return getExpProd( H_arr ) * this->H_init.block< 4, 4 >( 0, 4 * ( bodyID - 1 ) );
	}
	else if ( type == TYPE_COM )
	{	
		return getExpProd( H_arr ) * this->H_COM_init.block< 4, 4 >( 0, 4 * ( bodyID - 1 ) );
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

	// Initialize the Spatial Jacobian Matrix, which is a ( 6 x nq ) matrix
	Eigen::MatrixXd J_Spatial( 6, this->nq );

	// The First column of the Spatial Jacobian is the first joint twist array
	// Hence, saving the first joint twist values to the Spatial Jacobian
	J_Spatial.col( 0 ) = this->JointTwists.col( 0 );

	Eigen::MatrixXd H_tmp( 4, 4 );
	H_tmp = Eigen::MatrixXd::Identity( 4, 4 );

	for( int i = 0; i < ( this->nq - 1 ); i ++ )
	{
		H_tmp *= getExpSE3( this->JointTwists.block< 3, 1 >( 3, i ), this->JointTwists.block< 3, 1 >( 0, i ), q_arr( i ) );
		J_Spatial.col( i + 1 ) = getAdjoint( H_tmp ) * this->JointTwists.col( i + 1 );
	}

	return J_Spatial;
}

Eigen::MatrixXd RobotPrimitive::getSpatialJacobian( const Eigen::VectorXd &q_arr, const int bodyID )		
{
	// Assertion that q_arr length must be same with nq
	assert( this->nq == q_arr.size( ) );

	// Assertion that bodyID must be smaller or equal to nq
	assert( bodyID <= this->nq && bodyID >= 1  );

	// Initialize the Spatial Jacobian Matrix, which is a ( 6 x nq ) matrix
	Eigen::MatrixXd J_Spatial( 6, this->nq );
	J_Spatial = Eigen::MatrixXd::Zero( 6, this->nq );

	// The First column of the Spatial Jacobian is the first joint twist array
	// Hence, saving the first joint twist values to the Spatial Jacobian
	J_Spatial.col( 0 ) = this->JointTwists.col( 0 );

	Eigen::MatrixXd H_tmp( 4, 4 );
	H_tmp = Eigen::MatrixXd::Identity( 4, 4 );

	for( int i = 0; i < bodyID-1; i ++ )
	{
		H_tmp *= getExpSE3( this->JointTwists.block< 3, 1 >( 3, i ), this->JointTwists.block< 3, 1 >( 0, i ), q_arr( i ) );
		J_Spatial.col( i + 1 ) = getAdjoint( H_tmp ) * this->JointTwists.col( i + 1 );
	}

	return J_Spatial;
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

	A.block< 3, 3 >( 0, 3 ) = - vec2SkewSym( H_EE.block< 3, 1 >( 0, 3 ) );
	
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
		M += this->getBodyJacobian( q_arr, i + 1, TYPE_COM ).transpose( ) * this->M_Mat.block< 6, 6 >( 0, 6 * i ) * this->getBodyJacobian( q_arr, i+ 1 , TYPE_COM );
	} 

	return M;
}

Eigen::MatrixXd RobotPrimitive::getMassMatrix2( const Eigen::VectorXd &q_arr )
{
	Eigen::MatrixXd M( this->nq, this->nq );
	Eigen::MatrixXd L_Mat( 6 * this->nq, 6 * this->nq );

	Eigen::MatrixXd Ai( 6, 1 );
	Eigen::MatrixXd Hi( 4, 4 );
	Eigen::MatrixXd Hj( 4, 4 );

	L_Mat = Eigen::MatrixXd::Identity( 6 * this->nq, 6 * this->nq );
	

	for( int i = 1; i < this->nq; i++ )
	{
		Ai = this->A_mat.block< 6, 1 >( 6*i, i );               
		Hi = this->H_COM_init.block< 4, 4 >( 0, 4 * i );
		Hj = this->H_COM_init.block< 4, 4 >( 0, 4 * (i-1) );

		L_Mat.block( 6*i, 0, 6, 6*i ) = getAdjoint( getExpSE3( -Ai.block<3,1>(3,0), -Ai.block<3,1>(0,0), q_arr( i ) ) * Hi.inverse( ) * Hj  ) *  L_Mat.block( 6*( i - 1 ), 0, 6, 6*i  );
	}

	M = this->A_mat.transpose( ) * L_Mat.transpose( ) * this->M_Mat2 *  L_Mat * this->A_mat;

	return M;
}


/*************************************************************/
/************************* SNAKE BOT *************************/
/*************************************************************/
SnakeBot::SnakeBot( const int ID, const char* name, const int nq, const double m, const double l ): RobotPrimitive( ID, name )
{
	// Assertion of nq, m and l, which should be positive values
	assert( nq >= 1 && m > 0 && l > 0 );
	this->nq = nq;

	// ======================================================== //
	// ================== JOINT TWIST SET-UP ================== //
	// ======================================================== //
	// Construct the H_init and H_COM_init matrices
	// Sizes are ( 4 x { 4 * (nq + 1) } ) and ( 4 x { 4 * nq } ), respectively.
	// The final one includes the end-effector's SE(3) Matrix
	// Due to the end-effector, matrix is 4 * nq "+1"
	Eigen::MatrixXd     H_init( 4, 4 * ( this->nq + 1 ) );
	Eigen::MatrixXd H_COM_init( 4, 4 * ( this->nq     ) );

	// Define the Axis Origin, Axis Direction and Joint Types
	// The Size of each Matrix:
	// Axis Origin:    3 x nq
	// Axis Direction: 3 x nq
	// Joint Types:    1 x nq
	Eigen::MatrixXd AxisOrigins( 3, this->nq );
	Eigen::MatrixXd AxisDirections( 3, this->nq );
	Eigen::VectorXd JointTypes( this->nq );

	// The Joint Types are all 1 (i.e., revolute joint)
	AxisOrigins    = Eigen::MatrixXd::Zero( 3, this->nq );
	AxisDirections = Eigen::MatrixXd::Zero( 3, this->nq );
	JointTypes 	   = REVOLUTE_JOINT * Eigen::VectorXd::Ones( this-> nq );
	
	// The inertial parameters of the robot.
	Eigen::MatrixXd M_Mat(  6, 6 * this->nq );
	Eigen::MatrixXd Inertias( 3, 3 * this->nq );
	Eigen::VectorXd Masses( this->nq );	

	Eigen::MatrixXd M_Mat2( 6 * this->nq, 6 * this->nq );

	for( int i = 0; i < this->nq; i++ )
	{
		// The Rotation matrix is an identity matrix
		H_init.block< 4, 4 >( 0, 4*i ) = Eigen::Matrix4d::Identity( 4, 4 );

		// The x-position of H_init, others are zeros
		H_init( 0, 4*i+3 ) = l * i;

		H_COM_init.block< 4, 4 >( 0, 4*i ) = Eigen::Matrix4d::Identity( 4, 4 );
		H_COM_init( 0, 4*i+3 ) = l * ( i + 0.5 );

		// The x-position of the axis, others are zeros
		AxisOrigins( 0, i ) = l * i;

		// +1 along Z axis
		AxisDirections( 2, i ) = 1;
		
		// The mass of the (i+1)-th segment 
		Masses( i ) = m;

		// The inertia along the +z axis
		Inertias( 2, 3*i + 2 ) = 1.0/12. * m * l * l;			
		M_Mat.block< 3, 3 >( 0, 6*i 	) = m * Eigen::Matrix3d::Identity( 3, 3 );
		M_Mat.block< 3, 3 >( 3, 6*i + 3 ) = Inertias.block< 3, 3 >( 0, 3*i );

	}

	// Setup the final H_init Matrix
	H_init.block< 4, 4 >( 0, 4 * this->nq ) = Eigen::Matrix4d::Identity( 4, 4 );
	H_init( 0, 4*this->nq + 3 ) = l * this->nq;

	// Saving all of the calculated matrices
	this->H_init 		 = H_init;
	this->H_COM_init 	 = H_COM_init;
	this->JointTypes  	 = JointTypes;
	this->AxisOrigins 	 = AxisOrigins;
	this->AxisDirections = AxisDirections;
	this->Masses		 = Masses;
	this->Inertias		 = Inertias;
	this->M_Mat		 	 = M_Mat;

	// Once the values are assigned, set Joint Twists
	RobotPrimitive::setJointTwists(  );

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
	Eigen::MatrixXd M_Mat1( 6, 6 * this->nq );
	Eigen::MatrixXd M_Mat2( 6 * this->nq, 6 * this->nq );

	// Setup the final H_init Matrix for Media Flange Touch
	// Eigen::Vector3d FlangePos = AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0, 0.0, 0.071 );                    
	Eigen::Vector3d FlangePos = AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0, 0.0, 0.0314 );                    
	H_init.block< 4, 4 >( 0, 4 * this->nq 	  ) = Eigen::Matrix4d::Identity( 4, 4 );
	H_init.block< 3, 1 >( 0, 4 * this->nq + 3 ) = FlangePos;

	// Saving all of the calculated matrices
	this->H_init 		 = H_init;
	this->H_COM_init 	 = H_COM_init;
	this->JointTypes  	 = JointTypes;
	this->AxisOrigins 	 = AxisOrigins;
	this->AxisDirections = AxisDirections;
	this->Masses		 = Masses;
	this->Inertias		 = Inertias;

}
