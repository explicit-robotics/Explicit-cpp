/*
 * EXPlicit - A robotics toolbox based on product of exponential formula.
 *
 * Copyright (c) 2022 MIT
 * Authors
 * 			Johannes Lachner  	<jlachner@mit.edu>	
 * 			Moses C. Nah 		<mosesnah@mit.edu>
 *
 */

#ifndef EXP_MATH
#define EXP_MATH


#include <Eigen/Dense>

/**
 * @brief The so(3) to SO(3) exponential map, the Rodriguez's formula
 * @param w The so(3) vector or skew symmetric form
 * @param theta The rotation about w
 * @return The exp( [w] ) = I + sin(theta) [w] + ( 1 - cos( theta ) ) [w]^2
 */
Eigen::Matrix3d getExpSO3( const Eigen::Vector3d &w, const double theta );

Eigen::Matrix4d getExpSE3( const Eigen::Vector3d &w, const Eigen::Vector3d &v, const double theta );

Eigen::MatrixXd getAdjoint( const Eigen::Matrix4d &H );

Eigen::MatrixXd getInvAdjoint( const Eigen::Matrix4d &H );

Eigen::Matrix3d vec2SkewSym( const Eigen::Vector3d &v );

Eigen::Vector3d skewSym2Vec( const Eigen::Matrix3d &M );

Eigen::MatrixXd getExpProd( const Eigen::MatrixXd &H );

#endif