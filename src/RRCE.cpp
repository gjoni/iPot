/*
 * RRCE.cpp
 *
 *  Created on: May 27, 2016
 *      Author: ivan
 */

#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>

#include "RRCE.h"

RRCE::RRCE(RRCE_TYPE type, double d, int k) :
		h(NULL), J(NULL), dmax(d) {

	char *DATADIR = getenv("TMDOCKDAT");
	if (DATADIR == NULL) {
		printf("Error: environment variable 'TMDOCKDAT' not set\n");
		exit(1);
	}

	const int STRMAX = 1000;
	char tname[STRMAX], rrce_name[STRMAX];

	switch (type) {

	case RRCE20RC:
		sprintf(tname, "%s/RRCE20RC/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(rrce_name, "RRCE20RC_%.1fA_k%d", dmax, k);
		break;

	case RRCE20SCC:
		sprintf(tname, "%s/RRCE20SCC/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(rrce_name, "RRCE20SCC_%.1fA_k%d", dmax, k);
		break;

	case RRCE20CB:
		sprintf(tname, "%s/RRCE20CB/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(rrce_name, "RRCE20CB_%.1fA_k%d", dmax, k);
		break;

	default:
		printf("Error: wrong RRCE type (%d)\n", type);
		exit(1);
		break;

	}

	Allocate(20);

	ReadCouplings(tname);

	name = rrce_name;

}

RRCE::RRCE(const RRCE &source) :
		h(NULL), J(NULL), dmax(source.dmax), name(source.name) {

	const int N = 20;

	Allocate(N);

	memcpy(h, source.h, N * sizeof(double));
	for (int i = 0; i < N; i++) {
		memcpy(J[i], source.J[i], N * sizeof(double));
	}

}

RRCE::RRCE(const char *fname, double dmax_) :
	h(NULL), J(NULL), dmax(dmax_) {

	Allocate(20);
	ReadCouplings(fname);
	name = fname;

}

RRCE::RRCE() :
		h(NULL), J(NULL), dmax(7.8), name("") {

	Allocate(20);

	name = "RRCE20_7.8A_k5";

	double h_[20] = { -8.07207e-01, -5.81121e-01, -5.10194e-01, -1.03498e+00, 
		2.02255e+00, -5.04953e-01, -1.27513e+00, -7.24711e-01, 5.08540e-01, 
		4.91374e-01, -3.37279e-01, -1.06075e+00, 1.25663e+00, 9.30134e-01, 
		-3.86274e-01, -6.55560e-01, -3.54092e-01, 1.87274e+00, 8.61632e-01, 
		2.88653e-01 };

	double J_[20][20] = { 
		{-1.13461e-01, 1.87919e-01, 2.42264e-01, 3.40469e-01, -1.01989e-01, 2.32959e-01, 4.11017e-01, -2.08627e-03, 1.16568e-01, -1.73213e-01, -1.63894e-01, 4.06361e-01, -9.11694e-02, -1.65980e-01, 9.37284e-02, 1.40787e-01, 3.63397e-02, -1.23698e-01, -7.02783e-02, -2.20705e-01 }, 
		{1.87919e-01, 4.19627e-01, 1.81807e-01, -3.74378e-01, 2.88583e-02, 1.95305e-01, -3.23708e-01, 1.65534e-02, 1.95470e-01, 1.46923e-01, 9.94067e-02, 7.96317e-01, 1.45739e-01, 3.57403e-02, 7.94297e-02, 1.09684e-01, 1.02048e-01, -1.62307e-01, -4.41038e-02, 7.63653e-02 }, 
		{2.42264e-01, 1.81807e-01, -2.23739e-01, 2.37409e-02, 2.22284e-02, 9.80413e-02, 2.38547e-01, -7.18402e-02, 5.81354e-02, 1.84437e-01, 3.37032e-01, 1.50394e-01, 1.65917e-01, 1.28305e-01, 6.06594e-02, -4.50238e-02, -2.69840e-02, 5.42123e-02, 5.04195e-03, 1.86089e-01 }, 
		{3.40469e-01, -3.74378e-01, 2.37409e-02, 4.89458e-01, 1.97153e-01, 2.87025e-01, 7.73059e-01, 2.93577e-02, -2.15499e-01, 3.41730e-01, 4.73982e-01, -2.42539e-01, 3.38528e-01, 3.21939e-01, 1.78864e-01, 3.88854e-02, 7.19048e-02, 1.05585e-01, 3.90860e-02, 2.70405e-01 }, 
		{-1.01989e-01, 2.88583e-02, 2.22284e-02, 1.97153e-01, -1.55601e+00, 3.76619e-02, 2.31575e-01, -9.96109e-02, -1.98640e-01, -2.74738e-01, -2.47791e-01, 1.20711e-01, -2.90761e-01, -3.97638e-01, -1.25118e-01, -7.16624e-02, -6.19429e-02, -4.01366e-01, -2.61906e-01, -2.92217e-01 }, 
		{2.32959e-01, 1.95305e-01, 9.80413e-02, 2.87025e-01, 3.76619e-02, 1.45838e-01, 4.19084e-01, 8.55031e-02, 1.44604e-01, 1.09701e-01, 1.74024e-01, 3.43527e-01, 1.30458e-01, 6.40932e-02, 9.39209e-02, 7.52075e-02, 3.53378e-02, -8.68396e-02, 2.26002e-02, 8.21436e-02 }, 
		{4.11017e-01, -3.23708e-01, 2.38547e-01, 7.73059e-01, 2.31575e-01, 4.19084e-01, 8.49596e-01, 2.22735e-01, -6.07080e-02, 2.35773e-01, 3.71978e-01, -2.03471e-01, 3.24733e-01, 2.46411e-01, 2.80698e-01, 1.74891e-01, 1.65182e-01, 9.59850e-02, 6.65437e-02, 1.77406e-01 }, 
		{-2.08627e-03, 1.65534e-02, -7.18402e-02, 2.93577e-02, -9.96109e-02, 8.55031e-02, 2.22735e-01, -1.96608e-01, -2.37346e-02, 1.23396e-01, 1.81726e-01, 1.52528e-01, -1.76834e-03, 9.45760e-03, -3.13937e-02, -6.70803e-02, -2.33967e-02, -5.80547e-02, -1.29298e-02, 4.71345e-02 }, 
		{1.16568e-01, 1.95470e-01, 5.81354e-02, -2.15499e-01, -1.98640e-01, 1.44604e-01, -6.07080e-02, -2.37346e-02, -3.64770e-01, 3.46150e-02, 6.34589e-02, 4.44976e-01, -6.02530e-02, -9.90441e-02, -1.52002e-02, -6.98510e-02, -4.53098e-02, -2.37447e-01, -1.73527e-01, -1.35154e-03 }, 
		{-1.73213e-01, 1.46923e-01, 1.84437e-01, 3.41730e-01, -2.74738e-01, 1.09701e-01, 2.35773e-01, 1.23396e-01, 3.46150e-02, -6.30830e-01, -4.71128e-01, 1.35460e-01, -3.28754e-01, -4.84823e-01, 1.12265e-01, 8.39153e-02, -9.35238e-02, -2.89457e-01, -3.80859e-01, -5.10251e-01 }, 
		{-1.63894e-01, 9.94067e-02, 3.37032e-01, 4.73982e-01, -2.47791e-01, 1.74024e-01, 3.71978e-01, 1.81726e-01, 6.34589e-02, -4.71128e-01, -4.80393e-01, 2.86524e-01, -2.59288e-01, -4.55632e-01, 9.57444e-02, 1.58471e-01, 3.76260e-02, -3.51207e-01, -3.09638e-01, -4.22827e-01 }, 
		{4.06361e-01, 7.96317e-01, 1.50394e-01, -2.42539e-01, 1.20711e-01, 3.43527e-01, -2.03471e-01, 1.52528e-01, 4.44976e-01, 1.35460e-01, 2.86524e-01, 6.44255e-01, 2.99044e-01, 1.75588e-01, 3.89149e-01, 1.84925e-01, 1.45488e-01, 1.17772e-01, -1.47717e-02, 1.26723e-01 }, 
		{-9.11694e-02, 1.45739e-01, 1.65917e-01, 3.38528e-01, -2.90761e-01, 1.30458e-01, 3.24733e-01, -1.76834e-03, -6.02530e-02, -3.28754e-01, -2.59288e-01, 2.99044e-01, -6.40928e-01, -4.53883e-01, 5.21250e-03, 7.11233e-02, 1.78978e-02, -4.00702e-01, -3.03968e-01, -2.81088e-01 }, 
		{-1.65980e-01, 3.57403e-02, 1.28305e-01, 3.21939e-01, -3.97638e-01, 6.40932e-02, 2.46411e-01, 9.45760e-03, -9.90441e-02, -4.84823e-01, -4.55632e-01, 1.75588e-01, -4.53883e-01, -6.58396e-01, -7.01911e-02, 3.45045e-03, -5.23787e-02, -5.11494e-01, -4.68451e-01, -4.11131e-01 }, 
		{9.37284e-02, 7.94297e-02, 6.06594e-02, 1.78864e-01, -1.25118e-01, 9.39209e-02, 2.80698e-01, -3.13937e-02, -1.52002e-02, 1.12265e-01, 9.57444e-02, 3.89149e-01, 5.21250e-03, -7.01911e-02, -8.24372e-02, 8.88011e-02, 4.88165e-02, -3.73890e-01, -2.46803e-01, 3.62246e-02 }, 
		{1.40787e-01, 1.09684e-01, -4.50238e-02, 3.88854e-02, -7.16624e-02, 7.52075e-02, 1.74891e-01, -6.70803e-02, -6.98510e-02, 8.39153e-02, 1.58471e-01, 1.84925e-01, 7.11233e-02, 3.45045e-03, 8.88011e-02, -6.93395e-02, -2.72636e-02, -4.46887e-02, 1.65519e-03, 3.40155e-02 }, 
		{3.63397e-02, 1.02048e-01, -2.69840e-02, 7.19048e-02, -6.19429e-02, 3.53378e-02, 1.65182e-01, -2.33967e-02, -4.53098e-02, -9.35238e-02, 3.76260e-02, 1.45488e-01, 1.78978e-02, -5.23787e-02, 4.88165e-02, -2.72636e-02, -1.35707e-01, -1.52535e-02, -5.43936e-02, -1.44581e-01 }, 
		{-1.23698e-01, -1.62307e-01, 5.42123e-02, 1.05585e-01, -4.01366e-01, -8.68396e-02, 9.59850e-02, -5.80547e-02, -2.37447e-01, -2.89457e-01, -3.51207e-01, 1.17772e-01, -4.00702e-01, -5.11494e-01, -3.73890e-01, -4.46887e-02, -1.52535e-02, -6.29388e-01, -3.75391e-01, -2.50661e-01 }, 
		{-7.02783e-02, -4.41038e-02, 5.04195e-03, 3.90860e-02, -2.61906e-01, 2.26002e-02, 6.65437e-02, -1.29298e-02, -1.73527e-01, -3.80859e-01, -3.09638e-01, -1.47717e-02, -3.03968e-01, -4.68451e-01, -2.46803e-01, 1.65519e-03, -5.43936e-02, -3.75391e-01, -3.53153e-01, -3.01465e-01 }, 
		{-2.20705e-01, 7.63653e-02, 1.86089e-01, 2.70405e-01, -2.92217e-01, 8.21436e-02, 1.77406e-01, 4.71345e-02, -1.35154e-03, -5.10251e-01, -4.22827e-01, 1.26723e-01, -2.81088e-01, -4.11131e-01, 3.62246e-02, 3.40155e-02, -1.44581e-01, -2.50661e-01, -3.01465e-01, -5.53687e-01 } 
		};

	for (int i = 0; i < 20; i++) {
		h[i] = h_[i];
		for (int j = 0; j < 20; j++) {
			J[i][j] = J_[i][j];
		}
	}

}

RRCE::~RRCE() {

	Free(20);

}

RRCE & RRCE::operator =(const RRCE & source) {

	assert(this != &source);

	const int N = 20;

	Free(N);

	dmax = source.dmax;
	name = source.name;

	Allocate(N);

	memcpy(h, source.h, N * sizeof(double));
	for (int i = 0; i < N; i++) {
		memcpy(J[i], source.J[i], N * sizeof(double));
	}

	return *this;

}

void RRCE::Allocate(int N) {

	h = (double*) malloc(N * sizeof(double));
	J = (double**) malloc(N * sizeof(double*));
	for (int i = 0; i < N; i++) {
		J[i] = (double*) malloc(N * sizeof(double));
	}

}

void RRCE::Free(int N) {

	if (J != NULL) {
		for (int i = 0; i < N; i++) {
			free(J[i]);
		}
		free(J);
	}
	free(h);

}

int RRCE::ReadCouplings(const char *name) {

	FILE *F = fopen(name, "r");
	if (F == NULL) {
		printf("Error: cannot open %s file\n", name);
		exit(1);
	}

	const int STRMAX = 10000;
	char buf[STRMAX];

	const int N = 20;

	for (int i = 0; i < N; i++) {
		fgets(buf, STRMAX, F);
		char * pch;
		pch = strtok(buf, " "); /* label */
		pch = strtok(NULL, " "); /* local field */
		h[i] = atof(pch);
		for (int j = 0; j < N; j++) {
			pch = strtok(NULL, " ");
			J[i][j] = atof(pch);
		}
	}
	fclose(F);

	return 0;

}

std::string RRCE::GetName() {

	return name;

}

void RRCE::GetCouplings(double **J_) {

	for (int i = 0; i < 20; i++) {
		for (int j = 0; j < 20; j++) {
			J_[i][j] = J[i][j];
		}
	}

}

void RRCE::GetFields(double *h_) {

	for (int i = 0; i < 20; i++) {
		h_[i] = h[i];
	}

}

double RRCE::GetEnergy(const Chain &C) {

	double E = 0.0;

	for (int i = 0; i < C.nRes; i++) {

		char ti = C.residue[i].type;
		if (ti < 0 || ti > 19) {
			continue;
		}

		Residue *R = &(C.residue[i]);

		kdres *res = kd_nearest_range3f(C.kdCent, R->centroid[0],
				R->centroid[1], R->centroid[2], dmax);
		double pos[3];

		while (!kd_res_end(res)) {

			Atom *Aj = (Atom*) kd_res_item(res, pos);
			int j = Aj->residue - C.residue;
			char tj = Aj->residue->type;
			if (tj >= 0 && tj < 20 && i < j) {
				E += J[(int) ti][(int) tj];
			}
			kd_res_next(res);
		}
		kd_res_free(res);
	}

	return E;

}

double RRCE::GetEnergy(const Complex &C) {

	double E = 0;

	for (int i = 0; i < C.Receptor.nRes; i++) {

		char ti = C.Receptor.residue[i].type;
		if (ti < 0 || ti > 19) {
			continue;
		}

		Residue *R = &(C.Receptor.residue[i]);
		double *ci = R->centroid;

		/* check whether the current atom is within the bounding box
		 * of the Receptor */
		kdtree *kd = C.Ligand.kdCent;
		if (ci[0] < kd->rect->min[0] - dmax || ci[1] < kd->rect->min[1] - dmax
				|| ci[2] < kd->rect->min[2] - dmax
				|| ci[0] > kd->rect->max[0] + dmax
				|| ci[1] > kd->rect->max[1] + dmax
				|| ci[2] > kd->rect->max[2] + dmax) {

			continue;

		}

		kdres *res = kd_nearest_range3f(kd, ci[0], ci[1], ci[2], dmax);
		double pos[3];

		while (!kd_res_end(res)) {

			Atom *Aj = (Atom*) kd_res_item(res, pos);
			char tj = Aj->residue->type;

			if (tj >= 0 && tj < 20) {
				E += J[(int) ti][(int) tj];
			}
			kd_res_next(res);
		}
		kd_res_free(res);
	}

	return E;

}
