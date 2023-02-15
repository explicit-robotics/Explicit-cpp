/*
 * EXPlicit - A robotics toolbox based on the product of exponential formulae.
 *
 * Copyright (c) 2022 MIT
 * Authors
 * 			Johannes Lachner  	<jlachner@mit.edu>	
 * 			Moses C. Nah 		<mosesnah@mit.edu>
 * 
 * Function prototypes (function declarations) 
 * for basic math operations for the produce of exponentials formula.
 */

#ifndef EXP_MATH
#define EXP_MATH

#include <Eigen/Dense>

/**
 * @brief The so(3) to SO(3) exponential map, the Rodriguez's formula
 * @param w An element of so(3), represented as a (3x1) or (1x3) vector
 * @param theta The angle of rotation about w
 * @return 3x3 matrix: exp( [w]theta ) = I + sin(theta) [w] + ( 1 - cos( theta ) ) [w]^2
 */
Eigen::Matrix3d getExpSO3( const Eigen::Vector3d &w, const double theta );

/**
 * @brief The se(3) to SE(3) exponential map.
 * @param w The rotational    velocity
 * @param v The translational velocity
 * @param theta The angle of rotation
 * @return The exp( [ [w], v; 0, 1] ) = 
 * 				[1] if w is a     zero vector: v * theta 
 * 				[1] if w is a non-zero vector: ( I - exp( [w]theta ) ) * [w]v
 */
Eigen::Matrix4d getExpSE3( const Eigen::Vector3d &w, const Eigen::Vector3d &v, const double theta );

/**
 * @brief An adjoint map which maps between se(3) to se(3)
 * @param H A 4x4 SE(3) matrix
 * @return The exp( [ [w], v; 0, 1] ) = 
 * 				[1] if w is a     zero vector: v * theta 
 * 				[1] if w is a non-zero vector: ( I - exp( [w]theta ) ) * [w]v
 */
Eigen::MatrixXd getAdjoint( const Eigen::Matrix4d &H );

/**
 * @brief An adjoint map which maps between se(3) to se(3)
 * @param H A 4x4 SE(3) matrix, which consists of 
 * 	      H = [R, p;
 * 			   0, 1], R is an element of SO(3) and p is a R3 vector
 * @return A 6x6 matrix:
 * 			[ R, [p]R;
 * 		  	  0,    R]
 */
Eigen::MatrixXd getInvAdjoint( const Eigen::Matrix4d &H );

/**
 * @brief A transformation from a R3 to so(3) 3x3 matrix
 * @param H A 4x4 SE(3) matrix, which consists of 
 * 	      H = [R, p;
 * 			   0, 1], R is an element of SO(3) and p is a R3 vector
 * @return A 6x6 matrix:
 * 			[ R^T, -R^T[p];
 * 		  	  0,    R^T	  ]
 */
Eigen::Matrix3d vec2SkewSym( const Eigen::Vector3d &v );

/**
 * @brief A transformation from so(3) 3x3 matrix to R3 vector
 * @param M A skew-symmetric matrix, given by:
 * 		    M = [  0, -Mz,  My; 
 * 		     	  Mz,   0, -Mx;
 * 				 -My,  Mx,   0]
 * @return An R3 vector [ Mx, My, Mz ]
 */
Eigen::Vector3d skewSym2Vec( const Eigen::Matrix3d &M );

/**
 * @brief Conducting a product of exponentials formula for n arrays
 * @param H n of 4x4 matrix, expressed as a 2D matrix of dimension ( 4 x (4 x n ) )
 * 		    H = [ H1;
 * 				  H2;
 * 				  ..
 * 				  Hn], where Hi is a (4x4 matrix) 
 * @return A 4x4 matrix of H1 * H2 * ... * Hn
 */
Eigen::MatrixXd getExpProd( const Eigen::MatrixXd &H );

#endif