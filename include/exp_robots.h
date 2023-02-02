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

		/* 
			Inertial Properties
		*/
		Eigen::VectorXd Masses;
		Eigen::MatrixXd Inertias;
		Eigen::MatrixXd M_Mat;	

		/* 
			Animation Properties
			TODO
		*/

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

		void setJointTwists(  );

		Eigen::Matrix4d getForwardKinematics( const Eigen::VectorXd &q_arr );		

		Eigen::Matrix4d getForwardKinematics( const Eigen::VectorXd &q_arr, const int bodyID, const int type );

		// Eigen::Matrix4d getForwardKinematics( const Eigen::VectorXd &q_arr, const int bodyID, const Eigen::Vector3d position );

		Eigen::MatrixXd getSpatialJacobian( const Eigen::VectorXd &q_arr );		

		Eigen::MatrixXd getSpatialJacobian( const Eigen::VectorXd &q_arr, const int bodyID );		

		Eigen::MatrixXd getHybridJacobian( const Eigen::VectorXd &q_arr );

		Eigen::MatrixXd getBodyJacobian( const Eigen::VectorXd &q_arr, const int bodyID, const int type );	

		Eigen::MatrixXd getMassMatrix( const Eigen::VectorXd &q_arr );	

		// void switchJoint(  );		

		// void init( );			
		
		// void getCoriolisMatrix(  );		

		// void addKinematics(  );		

		

};

class SnakeBot : public RobotPrimitive
{
	private:

	public:
		SnakeBot( ) {};

		SnakeBot( const int ID, const char* name, const int nq, const double m, const double l );

		// SnakeBot( const int ID, const char* name, const int nq, const Eigen::VectorXd &m_arr, const Eigen::VectorXd &m_arr );

};

class iiwa14 : public RobotPrimitive
{
	private:

	public:
		iiwa14( ){};

		iiwa14( const int ID, const char* name );

		Eigen::VectorXd addIIWALimits( iiwa14 *myIIWA, Eigen::VectorXd q, Eigen::VectorXd qDot, Eigen::MatrixXd Minv, const Eigen::VectorXd tau, double dt );	

};




#endif