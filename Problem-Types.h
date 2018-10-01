/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: Using BSF Skeleton for Solving Gravitational Problem
Module: Problem-Types.h (BSF Types)
Prefix: PT
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/
#pragma once						
#include "Problem-Parameters.h"		// Problem Parameters 

//=========================== Problem Types =========================
struct PT_smallPoint_T {
	double coordinates[PP_SPACE_DIMENSION];
	double velosity[PP_SPACE_DIMENSION];
	double acceleration[PP_SPACE_DIMENSION];
	double mass;
};
struct PT_bigPoint_T {
	double coordinates[PP_SPACE_DIMENSION];
	double mass;
};

//=========================== BSF Types =========================
struct PT_bsf_data_T {				// Data for workers
	// Pointers are not allowed
	double coordinates[PP_SPACE_DIMENSION];
	double velosity[PP_SPACE_DIMENSION];
};

struct PT_bsf_mapElem_T {			// Element of map list
	// Pointers are not allowed
	int materialPointNo;
};

struct PT_bsf_reduceElem_T {		// Element of reduce list	
	// Pointers are not allowed
	double acceleration[PP_SPACE_DIMENSION];
};