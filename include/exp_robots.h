/*
 * EXPlicit - A robotics toolbox based on produce of exponential formula.
 *
 * Copyright (c) 2022 MIT
 * Authors
 * 			Johannes Lachner  	<jlachner@mit.edu>	
 * 			Moses C. Nah 		<mosesnah@mit.edu>
 *
 */

#ifndef EXP_ROBOT
#define EXP_ROBOT

#include <Eigen/Dense>

class RobotPrimitive
{
	private:
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
			Constructors
		*/
		RobotPrimitive( );

		/*
			Constructor

		*/
		RobotPrimitive( const int ID, const char* Name, const Eigen::VectorXd &JointTypes, const Eigen::MatrixXd &AxisOrigins, const Eigen::MatrixXd &AxisDirections );

		// void init( );

		void setJointTwists(  );

		Eigen::Matrix4d getForwardKinematics( const Eigen::VectorXd &q_arr );		

		Eigen::MatrixXd getSpatialJacobian( const Eigen::VectorXd &q_arr );		

		Eigen::MatrixXd getHybridJacobian( const Eigen::VectorXd &q_arr );				

		Eigen::MatrixXd getBodyJacobian( const Eigen::VectorXd &q_arr );		

		// void getMassMatrix(  );		

		// void getCoriolisMatrix(  );		

		// void switchJoint(  );		

		// void addKinematics(  );		

		

};

#endif