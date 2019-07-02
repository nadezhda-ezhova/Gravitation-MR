/*==============================================================================
Project: Bulk Synchronous Farm (BSF)
Theme: Using BSF Skeleton for Solving Gravitational Problem
Module: Problem-Implementation.h (Implementation of the Problem)
Prefix: PI
Authors: Nadezhda A. Ezhova, Leonid B. Sokolinsky
Creation Date: 09.04.2017
==============================================================================*/
#include <math.h>
#include "Problem-Include.h"
#include "Problem-Types.h"			// Problem Types 
#include "Problem-Data.h"			// Problem Data 
#include "Problem-Parameters.h"		// Problem Parameters 
#include "Problem-Forwards.h"		// Function Forwards
using namespace std;

#define MIN(x,y) (x<y?x:y)

void PI_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, PT_bsf_data_T* data, int* success) {
	double sub[PP_SPACE_DIMENSION], sum = 0, F[PP_SPACE_DIMENSION];
	for (int i = 0; i < PP_SPACE_DIMENSION; i++) {
		sub[i] = PD_Y[mapElem->materialPointNo].coordinates[i] - PD_X.coordinates[i];
		sum += pow(sub[i], 2);
	}
	double norm = sqrt(sum);
	double mm = PD_X.mass * PD_Y[mapElem->materialPointNo].mass;
	for (int i = 0; i < PP_SPACE_DIMENSION; i++) {
		F[i] = PP_G * mm * sub[i] / norm;
		reduceElem->acceleration[i] = F[i] / PD_X.mass;
	}

	/*debug
	cout << "acceleration:\t";
	cout << reduceElem->acceleration[0] << "\t";
	cout << reduceElem->acceleration[1] << "\t";
	cout << reduceElem->acceleration[2] << endl;/**/
};

void PI_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	for (int i = 0; i < PP_SPACE_DIMENSION; i++)
		z->acceleration[i] = x->acceleration[i] + y->acceleration[i];
};

void PI_bsf_ProcessResults(bool* exit, PT_bsf_reduceElem_T* reduceResult, int count, PT_bsf_data_T* data) {
	static int counter = 0;	// Iteration Counter
	counter++;
	if (counter < PP_ITER_COUNT)
		*exit = false;
	else
		*exit = true;

	for (int i = 0; i < PP_SPACE_DIMENSION; i++) {
		data->velosity[i] = data->velosity[i] + reduceResult->acceleration[i] * PP_dT;
		data->coordinates[i] = data->coordinates[i] + data->velosity[i] * PP_dT;
	}
};

void PI_bsf_CopyData(PT_bsf_data_T* dataIn, PT_bsf_data_T* dataOut) {
	for (int j = 0; j < PP_SPACE_DIMENSION; j++) {
		dataOut->coordinates[j] = dataIn->coordinates[j];
		dataOut->velosity[j] = dataIn->velosity[j];
	}
};

bool PI_bsf_Init(PT_bsf_data_T* data, char** message) {
	cout << setprecision(PP_PRECISION);
	if (PP_SPACE_DIMENSION != 3) {
		static char st[] = "PP_SPACE_DIMENSION must be equal to 3";
		*message = st;
		return false;
	};

	PD_X.mass = PP_SMALL_MASS;
	for (int j = 0; j < PP_SPACE_DIMENSION; j++) {
		PD_X.coordinates[j] = 0;
		PD_X.velosity[j] = PP_VELOSITY;
		PD_X.acceleration[j] = 0;
	};

	// Generating PD_Y[*]
	if (!PP_GRID) {
		for (int j = 0; j < PP_SPACE_DIMENSION; j++)
			for (int l = 0; l < PP_N_DIM; l++)
				PD_Y[l + j*PP_N_DIM].coordinates[j] = (l+1)*PP_DIM_STEP;
	} else {
		exit(0);
		int i = 0;
		for (int j1 = 0; j1 < PP_N_DIM; j1++)
			for (int j2 = 0; j2 < PP_N_DIM; j2++)
				for (int j3 = 0; j3 < PP_N_DIM; j3++) {
					if (j1 == j2 && j2 == j3)
						continue;
					PD_Y[i].coordinates[0] = j1 * PP_DIM_STEP;
					PD_Y[i].coordinates[1] = j2 * PP_DIM_STEP;
					PD_Y[i].coordinates[2] = j3 * PP_DIM_STEP;
					i++;
					if (i > PP_N)
						break;
				};
	};

	for (int i = 0; i < PP_N; i++)
		PD_Y[i].mass = PP_BIG_MASS;
	//PD_Y[0].mass = PP_BIG_MASS * 1000; // Удалить!!!

	for (int j = 0; j < PP_SPACE_DIMENSION; j++) {
		data->coordinates[j] = PD_X.coordinates[j];
		data->velosity[j] = PD_X.velosity[j];
	};
	return true;
};

void PI_bsf_AssignListSize(int* listSize) {
	*listSize = PP_N;
};

void PI_bsf_AssignMapList(PT_bsf_mapElem_T *mapList, int listSize) {
	for (int j = 0; j < listSize; j++) // First column has number 0
		mapList[j].materialPointNo = j;
};

void PI_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int count, PT_bsf_data_T data,
	int iterCount, double t, double t_L, double t_s_L, double t_S, double t_r_L, double t_W,
	double t_A_w, double t_A_m, double t_p) {
	cout << "---------------------" << endl;
	cout << "Iteration count = " << iterCount << endl;
	cout << "t_L = " << t_L << "\tt_s_L = " << t_s_L << "\tt_r_L = " << t_r_L << "\tt_W = " << t_W << "\tt_A_w = " << t_A_w << "\tt_A_m = " << t_A_m << "\tt_p = " << t_p << endl;
	cout << "Runtime = " << t << endl;
	cout << "---------- Solution ----------- " << endl;

	cout << "Coordinates:";
	for (int j = 0; j < PP_SPACE_DIMENSION; j++)
		cout << data.coordinates[j] << "\t";
	cout << endl;
	cout << "Velosity: ";
	for (int j = 0; j < PP_SPACE_DIMENSION; j++)
		cout << data.velosity[j] << "\t";
	cout << endl;

	cout << "---------------------" << endl;
	/**/
};

void PI_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int count, PT_bsf_data_T data,
	int iterCount, double elapsedTime) {
/*	cout << "coordinates:\t";
	for (int j = 0; j < PP_SPACE_DIMENSION; j++)
		cout << data.coordinates[j] << "\t";
	cout << "\t\tvelosity:\t";
	for (int j = 0; j < PP_SPACE_DIMENSION; j++)
		cout << data.velosity[j] << "\t";
	cout << endl;/**/
};

void PI_bsf_ParametersOutput(int numOfWorkers, PT_bsf_data_T data) {
	cout << "---------- Gravitation MRA -----------" << endl;
	cout << "Number of workers: " << numOfWorkers << endl;
	cout << "Number of big material points: " << PP_N << endl;
	cout << "Dimension of space n = " << PP_SPACE_DIMENSION << endl;

	cout << "---------- PD_X -----------" << endl;
	cout << "Initial coordinates:\t";
	for (int j = 0; j < PP_SPACE_DIMENSION; j++)
		cout << data.coordinates[j] << "\t";
	cout << endl;
	cout << "Initial velosity:\t";
	for (int j = 0; j < PP_SPACE_DIMENSION; j++)
		cout << data.velosity[j] << "\t";
	cout << endl;
	cout << "Mass:\t" << PD_X.mass << endl;

	cout << "---------- PD_Y -----------" << endl;
	for (int i = 0; i < MIN(PP_OUTPUT_LIMIT, PP_N); i++) {
		cout << "PD_Y[" << i << "].coordinates:\t";
		for (int j = 0; j < PP_SPACE_DIMENSION; j++)
			cout << PD_Y[i].coordinates[j] << "\t";
		cout << "\tmass:\t" << PD_Y[i].mass << endl;
	};
	cout << (PP_OUTPUT_LIMIT < PP_N ? "..." : "") << endl;
	/**/
	cout << "___________________________________________________________________________" << endl;
};
