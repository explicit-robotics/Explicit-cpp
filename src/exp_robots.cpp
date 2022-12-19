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
	Eigen::MatrixXd M_Mat( 6, 6 * this->nq );
	Eigen::MatrixXd Inertias( 3, 3 * this->nq );
	Eigen::VectorXd Masses( this->nq );	

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
	const int nq = 7;
	this->nq = nq;
	
	// TODO: Joint limits

	// ======================================================== //
	// ================= GEOMETRIC PROPERTIES  ================ //
	// ======================================================== //

	// The joint types are all 1 (i.e., revolute joints)
	Eigen::VectorXd JointTypes = REVOLUTE_JOINT * Eigen::VectorXd::Ones( this-> nq );

	// Define the joint axes positions and directions 
	Eigen::MatrixXd AxisOrigins( 3, this->nq );
	AxisOrigins.col( 0 ) = Eigen::Vector3d( 0.0,   0.0,    152.5e-3 );  
	AxisOrigins.col( 1 ) = Eigen::Vector3d( 0.0, -13.0e-3, 207.5e-3 );
	AxisOrigins.col( 2 ) = Eigen::Vector3d( 0.0, +13.0e-3, 232.5e-3 );
	AxisOrigins.col( 3 ) = Eigen::Vector3d( 0.0, +11.0e-3, 187.5e-3 );
	AxisOrigins.col( 4 ) = Eigen::Vector3d( 0.0, -11.0e-3, 212.5e-3 );
	AxisOrigins.col( 5 ) = Eigen::Vector3d( 0.0, -62.0e-3, 187.5e-3 );
	AxisOrigins.col( 6 ) = Eigen::Vector3d( 0.0, +62.0e-3,  79.6e-3 );

	Eigen::MatrixXd AxisDirections( 3, this->nq );
	AxisDirections.col( 0 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );
	AxisDirections.col( 1 ) = Eigen::Vector3d( 0.0,  1.0,  0.0 );
	AxisDirections.col( 2 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );
	AxisDirections.col( 3 ) = Eigen::Vector3d( 0.0, -1.0,  0.0 );
	AxisDirections.col( 4 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );
	AxisDirections.col( 5 ) = Eigen::Vector3d( 0.0,  1.0,  0.0 );
	AxisDirections.col( 6 ) = Eigen::Vector3d( 0.0,  0.0,  1.0 );

	// Define inertial parameters of the robot
	Eigen::VectorXd Masses( this->nq );	
	Masses( 0 ) = 6.404;
	Masses( 1 ) = 7.89;
	Masses( 2 ) = 2.54;
	Masses( 3 ) = 4.82;
	Masses( 4 ) = 1.76;
	Masses( 5 ) = 2.5;
	Masses( 6 ) = 0.42;

	// Eigen::Vector3d com[ nq ];
	Eigen::MatrixXd COM( 3, this->nq );
	COM.col( 0 ) = Eigen::Vector3d( 0.0,	-14.0e-3,  102.0e-3);  
	COM.col( 1 ) = Eigen::Vector3d( 0.0,  16.0e-3,   64.0e-3);  
	COM.col( 2 ) = Eigen::Vector3d( 0.0,  19.0e-3,   98.0e-3);  
	COM.col( 3 ) = Eigen::Vector3d( 0.0, -20.0e-3,   86.0e-3);
	COM.col( 4 ) = Eigen::Vector3d( 0.0, -13.0e-3,   66.0e-3);  
	COM.col( 5 ) = Eigen::Vector3d( 0.0,  60.0e-3,   16.0e-3);  
	COM.col( 6 ) = Eigen::Vector3d( 0.0,   0.0e-3,   11.0e-3); 

	Eigen::MatrixXd Inertias( 3, 3 * this->nq );	
	Inertias.block< 3, 3 >( 0, 3 * 0 ) << 0.069, 0.0,   0.0,
										  0.0 ,  0.071, 0.0,
										  0.0,   0.0,   0.02;

	Inertias.block< 3, 3 >( 0, 3 * 1 ) << 0.08, 0.0,  0.0,
				            	 		  0.0,  0.08, 0.0,
								 		  0.0,  0.0,  0.01;

	Inertias.block< 3, 3 >( 0, 3 * 2 ) << 0.02, 0.0,  0.0,
				            	 		  0.0,  0.02, 0.0,
								 		  0.0,  0.0,  0.06;
	
	Inertias.block< 3, 3 >( 0, 3 * 3 ) << 0.04, 0.0,  0.0, 
				            	 		  0.0,  0.03, 0.0,
								 		  0.0,  0.0,  0.01;
	
	Inertias.block< 3, 3 >( 0, 3 * 4 ) << 0.01, 0.0,  0.0, 
				            	 		  0.0,  0.01, 0.0,
								 		  0.0,  0.0,  0.01;

	Inertias.block< 3, 3 >( 0, 3 * 5 ) << 0.007, 0.0,   0.0, 
				        		 		  0.0,   0.006, 0.0,
								 		  0.0,   0.0,   0.005;

	Inertias.block< 3, 3 >( 0, 3 * 6 ) << 0.0003, 0.0,    0.0, 
				            	 		  0.0,    0.0003, 0.0,
								 		  0.0,    0.0,    0.0005;

	// Concatenate matrices for transformations and Generalized inertia matrix
	Eigen::MatrixXd H_init( 4, 4 * ( this->nq + 1 ) );
	Eigen::MatrixXd H_COM_init( 4, 4 * ( this->nq     ) );
	Eigen::MatrixXd M_Mat( 6, 6 * this->nq );
	for( int i = 0; i < this->nq; i++)
	{
		// The Rotation matrix is an identity matrix
		H_init.block< 4, 4 >( 0, 4 * i ) = Eigen::Matrix4d::Identity( 4, 4 );

		// The translational part of H_init
		H_init.block< 3, 1 >( 0, 4 * i + 3 ) = AxisOrigins.col( i );

		H_COM_init.block< 4, 4 >( 0, 4 * i ) = Eigen::Matrix4d::Identity( 4, 4 );
		H_COM_init.block< 3, 1 >( 0, 4 * i + 3 ) = COM.col( i );

		// Generalized inertia matrix
		M_Mat.block< 3, 3 >( 0, 6 * i 	) = Masses ( i ) * Eigen::Matrix3d::Identity( 3, 3 );
		M_Mat.block< 3, 3 >( 3, 6 * i + 3 ) = Inertias.block< 3, 3 >( 0, 3 * i );
	}

	// Setup the final H_init Matrix for Media Flange Touch
	Eigen::Vector3d FlangePos = Eigen::Vector3d( 0.0, 0.0, 0.071 );                    
	H_init.block< 4, 4 >( 0, 4 * this->nq ) = Eigen::Matrix4d::Identity( 4, 4 );
	H_init.block< 3, 1 >( 0, 4 * this->nq + 3 ) = FlangePos;

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