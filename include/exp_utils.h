/*
 * EXPlicit - A robotics toolbox based on produce of exponential formula.
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


#endif