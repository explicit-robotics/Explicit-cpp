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
 * @brief FILL IN
 * @param w FILL IN
 * @param theta FILL IN
 * @return FILL IN
 */

bool isSymmetric( const Eigen::MatrixXd &M );

bool isSkewSymmetric( const Eigen::MatrixXd &M );

double getMaxValue( Eigen::VectorXd &myVector );

double getMinValue( Eigen::VectorXd &myVector );


#endif