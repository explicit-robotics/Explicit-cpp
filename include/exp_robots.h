/*
 * EXPlicit - A robotics toolbox based on the product of exponential formulae.
 *
 * Copyright (c) 2022 MIT
 * Authors
 * 			Johannes Lachner  	<jlachner@mit.edu>	
 * 			Moses C. Nah 		<mosesnah@mit.edu>
 *
 */

#ifndef EXP_ROBOT
#define EXP_ROBOT

#include "exp_constants.h"
#include <Eigen/Dense>

/*
	@brief The primitive robot class which will be inherited for all robots
*/
class RobotPrimitive
{
	public:
		int ID;
		int parentID;
		const char* Name;

		/* 
			Joint Properties
		*/
		int nq;
		Eigen::VectorXd JointTypes;

		Eigen::MatrixXd AxisOrigins;
		Eigen::MatrixXd AxisDirections;
		Eigen::MatrixXd JointTwists;
		Eigen::MatrixXd COM;

		Eigen::VectorXd q;
		Eigen::VectorXd q_init;
		Eigen::VectorXd q_max;
		Eigen::VectorXd q_min;
		Eigen::VectorXd dq_max;
		Eigen::VectorXd dq_min;		
		Eigen::VectorXd ddq_max;
		Eigen::VectorXd ddq_min;				

		/* 
			The H Matrix (SE3)
		*/
		Eigen::MatrixXd H_init;
		Eigen::MatrixXd H_COM_init;
		Eigen::MatrixXd H_ij;
		Eigen::MatrixXd H_base;

		Eigen::MatrixXd H_arrs; 	// This is a 4 x (4xnq) array which temporarily saves the H matrices for calculation

		/* 
			Inertial Properties
		*/

		Eigen::VectorXd Masses;
		Eigen::MatrixXd Inertias;
		Eigen::MatrixXd M_Mat1;	
		Eigen::MatrixXd M_Mat2;		// The Diagonalized version of M_Mat1, Eq. 8.61 of Modern Robotics. 
		Eigen::MatrixXd L_Mat;		// The Matrix used for calculating the MAss matrix

		// The Jacobian Matrices
		Eigen::MatrixXd JS;
		Eigen::MatrixXd JB;
		Eigen::MatrixXd JH;

		// The A Matrix
		Eigen::MatrixXd A_Mat1;		// (6 x nq) array, collection of joint twists expressed in {i} frame	
		Eigen::MatrixXd A_Mat2;		// The Diagonalized version of A_Mat1, Eq. 8.60 of Modern Robotics. 

		/*
			Etc. 	
		*/
		double Grav = 9.81;

	public:

		/* 
			Default Constructor
		*/
		RobotPrimitive( ) {};

		/*
			Constructor
		*/
		RobotPrimitive( const int ID, const char* Name );

		void init( );
		void setJointTwists(  );
		void setGeneralizedMassMatrix( );

		Eigen::Matrix4d getForwardKinematics( const Eigen::VectorXd &q_arr );		
		Eigen::Matrix4d getForwardKinematics( const Eigen::VectorXd &q_arr, const int bodyID, const int type );
		Eigen::Matrix4d getForwardKinematics( const Eigen::VectorXd &q_arr, const int bodyID, const Eigen::Vector3d &p_pos );		

		Eigen::MatrixXd getSpatialJacobian( const Eigen::VectorXd &q_arr );		
		Eigen::MatrixXd getSpatialJacobian( const Eigen::VectorXd &q_arr, const int bodyID );		

		Eigen::MatrixXd getHybridJacobian( const Eigen::VectorXd &q_arr );
		Eigen::MatrixXd getHybridJacobian( const Eigen::VectorXd &q_arr, const Eigen::Vector3d &p_pos );

		Eigen::MatrixXd getBodyJacobian( const Eigen::VectorXd &q_arr, const int bodyID, const int type );	

		Eigen::MatrixXd getMassMatrix( const Eigen::VectorXd &q_arr );	

};

class iiwa14 : public RobotPrimitive
{
	private:

	public:
		iiwa14( ){};
		iiwa14( const int ID, const char* name );		
		iiwa14( const int ID, const char* name, const Eigen::Vector3d &flange );				

		Eigen::VectorXd addIIWALimits( Eigen::VectorXd q, Eigen::VectorXd dq, Eigen::MatrixXd M, Eigen::VectorXd tau, double dt );	

};

#endif
