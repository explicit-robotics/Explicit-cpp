/*
 * EXPlicit - A robotics toolbox based on the product of exponential formulae.
 *
 * Copyright (c) 2022 MIT
 * Authors
 * 			Johannes Lachner  	<jlachner@mit.edu>	
 * 			Moses C. Nah 		<mosesnah@mit.edu>
 *
 */

#ifndef EXP_UTILS
#define EXP_UTILS

#include <Eigen/Dense>

/**
 * @brief Check whether parameter M is a symmetric matrix
 * @param M a 3x3 matrix
 * @return True or False
 */
bool isSymmetric( const Eigen::MatrixXd &M );

/**
 * @brief Check whether parameter M is a skew-symmetric matrix
 * @param M a 3x3 matrix
 * @return True or False
 */
bool isSkewSymmetric( const Eigen::MatrixXd &M );

double getMaxValue( Eigen::VectorXd &myVector );

double getMinValue( Eigen::VectorXd &myVector );


#endif