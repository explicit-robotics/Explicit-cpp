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
	// this->q_max <<   70, 70, 70, 70, 70, 70, 70;
	// this->q_min <<  -70,-70,-70,-70,-70,-70,-70;

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
	// The Center of Mass Locations are based on the Data from Drake Software of Toyota Research Institute.
	// [REF] https://github.com/RobotLocomotion/drake/blob/master/manipulation/models/iiwa_description/urdf/iiwa14_spheres_dense_collision.urdf 
	this->COM = Eigen::MatrixXd( 3, this->nq );
	this->COM.col( 0 ) = this->AxisOrigins.col( 0 ) + Eigen::Vector3d( 0.0000, -0.0300,  0.1200 );  
	this->COM.col( 1 ) = this->AxisOrigins.col( 1 ) + Eigen::Vector3d( 0.0003,  0.0590,  0.0420 );  
	this->COM.col( 2 ) = this->AxisOrigins.col( 2 ) + Eigen::Vector3d( 0.0000,  0.0300,  0.1300 );  
	this->COM.col( 3 ) = this->AxisOrigins.col( 3 ) + Eigen::Vector3d( 0.0000,  0.0670,  0.0340 );
	this->COM.col( 4 ) = this->AxisOrigins.col( 4 ) + Eigen::Vector3d( 0.0001,  0.0210,  0.0760 );  
	this->COM.col( 5 ) = this->AxisOrigins.col( 5 ) + Eigen::Vector3d( 0.0000,  0.0006,  0.0004 );  
	this->COM.col( 6 ) = this->AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0000,  0.0000,  0.0200 ); 

	// Define the Principal axes/moments of inertia of the link
	// The inertia values are based on the Data from Drake Software of Toyota Research Institute.
	// [REF] https://github.com/RobotLocomotion/drake/blob/master/manipulation/models/iiwa_description/urdf/iiwa14_spheres_dense_collision.urdf 	
	this->Inertias = Eigen::MatrixXd( 3, 3 * this->nq );	
	this->Inertias.block< 3, 3 >( 0, 3 * 0 ) << 0.0330, 0.0000, 0.0000,
										  		0.0000, 0.0333, 0.0000,
										  		0.0000, 0.0000, 0.0123;

	this->Inertias.block< 3, 3 >( 0, 3 * 1 ) << 0.0305, 0.0000, 0.0000,
			  			            	 		0.0000, 0.0304, 0.0000,
								 		  		0.0000, 0.0000, 0.0110;

	this->Inertias.block< 3, 3 >( 0, 3 * 2 ) << 0.0250, 0.0000, 0.0000,
				            	 		  		0.0000, 0.0238, 0.0000,
								 		  		0.0000, 0.0000, 0.0076;
	
	this->Inertias.block< 3, 3 >( 0, 3 * 3 ) << 0.0170, 0.0000, 0.0000, 
				            	 		  		0.0000, 0.0164, 0.0000,
								 		  		0.0000, 0.0000, 0.0060;
	
	this->Inertias.block< 3, 3 >( 0, 3 * 4 ) << 0.0100, 0.0000, 0.0000, 
				            	 		  		0.0000, 0.0087, 0.0000,
								 		  		0.0000, 0.0000, 0.0045;

	this->Inertias.block< 3, 3 >( 0, 3 * 5 ) << 0.0049, 0.0000, 0.0000, 
				        		 		  		0.0000, 0.0036, 0.0000,
								 		  		0.0000, 0.0000, 0.0047;

	this->Inertias.block< 3, 3 >( 0, 3 * 6 ) << 0.0010, 0.0000, 0.0000, 
				            	 		  		0.0000, 0.0010, 0.0000,
								 		  		0.0000, 0.0000, 0.0010;

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
	Eigen::Vector3d FlangePos = AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0, 0.0, 0.0 );
	// One can also comment out to use the following line
	// where the end-effector is defined at the center of the flange.
	// Eigen::Vector3d FlangePos = AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0, 0.0, 0.071 );
	// Also for a different flange with a shorter length.
	// Eigen::Vector3d FlangePos = this->AxisOrigins.col( 6 ) + Eigen::Vector3d( 0.0, 0.0, 0.0314 );                    

	this->H_init.block< 4, 4 >( 0, 4 * this->nq	 	) = Eigen::Matrix4d::Identity( 4, 4 );
	this->H_init.block< 3, 1 >( 0, 4 * this->nq + 3 ) = FlangePos;

}



Eigen::VectorXd iiwa14::addIIWALimits( Eigen::VectorXd q, Eigen::VectorXd dq, Eigen::MatrixXd M, Eigen::VectorXd tau, double dt )
{
	/* This method is the implementation based on the paper:
	Muñoz Osorio, Juan & Fiore, Mario & Allmendinger, Felix. (2018). 
	Operational Space Formulation Under Joint Constraints. 
	doi: 10.1115/DETC2018-86058.   
	*/

	Eigen::VectorXd dt2 	 = Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd dt_var 	 = Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd rho_down = Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd rho_up 	 = Eigen::VectorXd::Zero( this->nq, 1 );

    Eigen::VectorXd dq_max_Q 		= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd dq_min_Q 		= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd dq_max_ddQ  	= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd dq_min_ddQ  	= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::Vector3d v_max_vec   	= Eigen::Vector3d::Zero( 3 );
	Eigen::Vector3d v_min_vec 		= Eigen::Vector3d::Zero( 3 );
	Eigen::VectorXd dq_max_fin 		= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd dq_min_fin	 	= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd a_max_dq 		= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd a_min_dq 		= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd a_max_Q 		= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd a_min_Q 		= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::Vector3d a_max_vec 		= Eigen::Vector3d::Zero( 3 );
	Eigen::Vector3d a_min_vec 		= Eigen::Vector3d::Zero( 3 );
	Eigen::VectorXd ddq_max_fin		= Eigen::VectorXd::Zero( this->nq, 1 );
	Eigen::VectorXd ddq_min_fin		= Eigen::VectorXd::Zero( this->nq, 1 );
    Eigen::MatrixXd Iden 			= Eigen::MatrixXd::Identity( this->nq, this->nq );
    Eigen::VectorXd tau_jl 			= Eigen::VectorXd::Zero( this->nq ,1 );
    Eigen::VectorXd ddq_got 		= Eigen::VectorXd::Zero( this->nq, 1 );
    Eigen::MatrixXd Js 				= Eigen::MatrixXd::Zero( 3 , this->nq );

    double low_dt_fac = 10;

    rho_down = q - this->q_min;
    rho_up = this->q_max - q;

    dt_var[ 0 ] = 3 * dt;
    dt_var[ 1 ] = 3 * dt;
    dt_var[ 2 ] = 2 * dt;
    dt_var[ 3 ] = 3 * dt;
    dt_var[ 4 ] = dt;
    dt_var[ 5 ] = dt;
    dt_var[ 6 ] = dt;

    for ( int i = 0; i < this->nq; i ++ )
    {
        dt2[ i ] = dt_var[ i ];
        if ( rho_down[ i ] < 10 * M_PI / 180 )
        {

            if ( rho_down[ i ] < 0 )
                rho_down[i] = 0;

            dt2[ i ] = ( low_dt_fac + ( sqrt( low_dt_fac ) * sqrt( rho_down[ i ] * 180 / M_PI ) ) ) * dt_var[ i ];

            if ( dt2[ i ] < low_dt_fac * dt_var[ i ] )
                dt2[ i ] = low_dt_fac * dt_var[ i ];
        }
        if ( rho_up[ i ] < 10 * M_PI / 180 )
        {

            if ( rho_up[ i ] < 0 )
                rho_up[ i ] = 0;

            dt2[ i ] = ( low_dt_fac + ( sqrt( low_dt_fac ) * sqrt( rho_up[ i ] * 180 / M_PI ) ) ) * dt_var[ i ];
            if ( dt2[ i ] < low_dt_fac * dt_var[ i ] )
                dt2[ i ] = low_dt_fac * dt_var[ i ];

        }


        dq_max_Q[ i ] = ( this->q_max[ i ] - q[ i ] ) / dt2[ i ];
        dq_min_Q[ i ] = ( this->q_min[ i ] - q[ i ]) / dt2[ i ];
        dq_max_ddQ[ i ] = sqrt( 2 * this->ddq_max[ i ] * ( this->q_max[ i ] - q[ i ] ) );
        dq_min_ddQ[ i ] = -sqrt( 2 * this->ddq_max[ i ] * ( q[ i ] - this->q_min[ i ] ) );

        if( this->q_max[ i ] - q[ i ] < 0 )
            dq_max_ddQ[ i ] = 1000000;

        if( q[ i ] - this->q_min[ i ] < 0 )
            dq_min_ddQ[ i ] = -1000000;

        v_max_vec = Eigen::Vector3d( this->dq_max[ i ], dq_max_Q[ i ], dq_max_ddQ[ i ] );
        dq_max_fin[ i ] = getMinValue( v_max_vec );

        v_min_vec = Eigen::Vector3d( this->dq_min[ i ], dq_min_Q[ i ], dq_min_ddQ[ i ] );
        dq_min_fin[ i ] = getMaxValue( v_min_vec );


        a_max_dq[ i ] = ( dq_max_fin[ i ] - dq[ i ]) / dt_var[ i ];
        a_min_dq[ i ] = ( dq_min_fin[ i ] - dq[ i ]) / dt_var[ i ];
        a_max_Q[ i ] = 2 * ( this->q_max[ i ] - q[ i ] - dq[ i ] * dt2[ i ] ) / pow( dt2[ i ], 2 );
        a_min_Q[ i ] = 2 * ( this->q_min[ i ] - q[ i ] - dq[ i ] * dt2[ i ] ) / pow( dt2[ i ], 2 );

        a_max_vec = Eigen::Vector3d( a_max_Q[ i ], a_max_dq[ i ], 10000000 );
        ddq_max_fin[ i ] = getMinValue( a_max_vec );
        a_min_vec = Eigen::Vector3d( a_min_Q[ i ], a_min_dq[ i ], -10000000 );
        ddq_min_fin[ i ] = getMaxValue( a_min_vec );

        if( ddq_max_fin[ i ] < ddq_min_fin[ i ] )
        {
            v_max_vec = Eigen::Vector3d( INFINITY, dq_max_Q[ i ], dq_max_ddQ[ i ] );
            dq_max_fin[ i ] = getMinValue( v_max_vec );

            v_min_vec = Eigen::Vector3d( -INFINITY, dq_min_Q[ i ], dq_min_ddQ[ i ] );
            dq_min_fin[ i ] = getMaxValue( v_min_vec );

            a_max_dq[ i ] = ( dq_max_fin[ i ] - dq[ i ] ) / dt_var[ i ];
            a_min_dq[ i ] = ( dq_min_fin[ i ] - dq[ i ] ) / dt_var[ i ];

            a_max_vec = Eigen::Vector3d( a_max_Q[ i ], a_max_dq[ i ], 10000000 );
            ddq_max_fin[ i ] = getMinValue( a_max_vec );
            a_min_vec = Eigen::Vector3d( a_min_Q[ i ], a_min_dq[ i ], -10000000 );
            ddq_min_fin[ i ] = getMaxValue( a_min_vec );
        }
    }


    Eigen::VectorXd ddq_s = Eigen::VectorXd::Zero( this->nq );
    Eigen::VectorXd tau_s = Eigen::VectorXd::Zero( this->nq );
    Eigen::MatrixXd P_sat = Iden;
    bool LimitedExceeded = true;
    bool CreateTaskSat = false;
    int NumSatJoints = 0;
    Eigen::VectorXd theMostCriticalOld = Eigen::VectorXd::Zero( this->nq );
    theMostCriticalOld.conservativeResize( 1 );
    theMostCriticalOld[ 0 ] = 100;
    bool isThere = false;
    int iO = 0;
    int cycle = 0;
    while ( LimitedExceeded == true )
    {

        LimitedExceeded = false;

        if ( CreateTaskSat == true )
        {
            Js.conservativeResize( NumSatJoints, this->nq );
            for ( int i = 0; i < NumSatJoints; i++ )
            {
                for( int k = 0; k < this->nq; k++ )
                {
                    Js( i, k ) = 0;
                }
                int m = theMostCriticalOld[ i ];
                Js( i, m ) = 1;
            }

            Eigen::MatrixXd Minv = M.inverse();
            Eigen::MatrixXd lambda_sat_inv = Js * Minv * Js.transpose();
            Eigen::MatrixXd lambda_sat_inv_aux = lambda_sat_inv * lambda_sat_inv.transpose();
            Eigen::MatrixXd lambda_sat_aux = lambda_sat_inv_aux.inverse();
            Eigen::MatrixXd lambda_sat = lambda_sat_inv.transpose() * lambda_sat_aux;

            Eigen::MatrixXd J_sat_bar = Minv * Js.transpose() * lambda_sat;
            P_sat = Iden - Js.transpose() * J_sat_bar.transpose();
            Eigen::VectorXd ddx_s = Js * ddq_s;
            tau_s = Js.transpose() * ( lambda_sat * ddx_s );

        }

        tau_jl = tau_s + P_sat * tau;
        ddq_got = M.inverse() * tau_jl;				// Assuming gravity and Coriolis are compensated

        isThere = false;
        for ( int i = 0; i < this->nq; i++ )
        {
            if ( ( ddq_max_fin[ i ] + 0.001 < ddq_got[ i ] )  || ( ddq_got[ i ] < ddq_min_fin[ i ] - 0.001 ) )
            {
                LimitedExceeded = true;
                CreateTaskSat = true;

                for ( int k = 0; k < theMostCriticalOld.size(); k++ )
                {
                    if ( i == theMostCriticalOld[ k ] )
                    {
                        isThere = true;
                    }
                }
                if ( isThere == false )
                {

                    theMostCriticalOld.conservativeResize( iO + 1 );
                    theMostCriticalOld[ iO ] = i;
                    iO += 1;
                }
            }
        }

        if ( LimitedExceeded == true )
        {
            NumSatJoints = iO;
            theMostCriticalOld.conservativeResize( iO );
            cycle += 1;
            if ( cycle > 8 )
                LimitedExceeded = false;

            for ( int i = 0; i < theMostCriticalOld.size(); i++ )
            {
                int jM = theMostCriticalOld[ i ];

                if ( ddq_got[ jM ] > ddq_max_fin[ jM ] )
                    ddq_s[ jM ] = ddq_max_fin[ jM ];

                if ( ddq_got[ jM ] < ddq_min_fin[ jM ] )
                    ddq_s[ jM ] = ddq_min_fin[ jM ];
            }
        }
    }

    return tau_jl;

}