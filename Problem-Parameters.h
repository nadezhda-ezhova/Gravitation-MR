/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: Using BSF Skeleton for Solving Gravitational Problem
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/
#pragma once
/* ========================= Compilation Modes ======================== */
//#define NDEBUG

//=========================== Skeleton Parameters =========================
#define PP_PRECISION 4	// Decimal precision on output

//=========================== Problem Parameters =========================
#define PP_SPACE_DIMENSION 3	// Dimension of space = 3 (change is prohibited)
#define PP_FLOAT_POINT_TYPE double

#define PP_N_DIM		5		// Number of big material points PD_Y[*] per one dimension
#define PP_GRID		0		// 0 - false; 1 - true
#define PP_N  (PP_N_DIM	* PP_SPACE_DIMENSION +(PP_GRID ? PP_N_DIM ^ PP_SPACE_DIMENSION : 0))	// Number of big material points PD_Y[*]

#define PP_DIM_STEP	1E2		// Distance between coodinates of neighbour points in one dimrnsion
#define PP_dT 1				// Time step (sec.)        
#define PP_ITER_COUNT 10			
#define PP_OUTPUT_LIMIT 5	// Maximal count of points to output
#define PP_G 6.67408E-11		// Gravitational Constant
#define PP_VELOSITY 1			// Initial velosity of small material point PD_X
#define PP_SMALL_MASS 10			// Mass of small material point PD_X
#define PP_BIG_MASS 1E12			// Mass of bvig material point PD_Y[*]