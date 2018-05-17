/*
 * AACE.cpp
 *
 *  Created on: Mar 30, 2016
 *      Author: ivan
 */

#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdio>

#include "AACE.h"

AACE::AACE(const GROUPS_TYPE &GROUPS, const std::vector<double> &h_, 
			const std::vector<std::vector<double> > &J_, 
		double dmax_, const std::string &name_) : h(NULL), J(NULL), 
		dmax(dmax_), name(name_) {

	int dim = ReadGroups(GROUPS);
	Allocate(dim);

	assert((int)h_.size() == dim);
	assert((int)J_.size() == dim);
	assert((int)J_[0].size() == dim);

	for (int i = 0; i < dim; i++) {
		h[i] = h_[i];
		for (int j = 0; j < dim; j++) {
			J[i][j] = J_[i][j];
		}
	}

}

AACE::AACE(const GROUPS_TYPE &GROUPS, double dmax_, 
		const std::string &name_) : 
		h(NULL), J(NULL), dmax(dmax_), name(name_) {

	int dim = ReadGroups(GROUPS);
	Allocate(dim);

	FILE *F = fopen(name.c_str(), "r");
	if (F == NULL) {
		printf("Error: cannot open table file %s\n", name.c_str());
		exit(1);
	}

	const int STRMAX = 10000;
	char buf[STRMAX];

	int N = 0;
	
	while (fgets(buf, STRMAX, F)) {
		N++;
	}
	rewind(F);
	
	assert(N == dim);

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

}

AACE::AACE(AACE_TYPE type, double d, int k) :
		h(NULL), J(NULL), dmax(d) {

	char *DATADIR = getenv("TMDOCKDAT");
	if (DATADIR == NULL) {
		printf("Error: environment variable 'TMDOCKDAT' not set\n");
		exit(1);
	}

	const int STRMAX = 1000;
	char gname[STRMAX], tname[STRMAX], aace_name[STRMAX];

	switch (type) {

	case AACE18:
		sprintf(gname, "%s/groups18.dat", DATADIR);
		sprintf(tname, "%s/AACE18/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(aace_name, "AACE18_%.1fA_k%d", dmax, k);
		break;

	case AACE167:
		sprintf(gname, "%s/groups167.dat", DATADIR);
		sprintf(tname, "%s/AACE167/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(aace_name, "AACE167_%.1fA_k%d", dmax, k);
		break;

	case AACE20:
		sprintf(gname, "%s/groups20.dat", DATADIR);
		sprintf(tname, "%s/AACE20/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(aace_name, "AACE20_%.1fA_k%d", dmax, k);
		break;

	case AACE24:
		sprintf(gname, "%s/groups24.dat", DATADIR);
		sprintf(tname, "%s/AACE24/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(aace_name, "AACE24_%.1fA_k%d", dmax, k);
		break;

	case AACE33:
		sprintf(gname, "%s/groups33.dat", DATADIR);
		sprintf(tname, "%s/AACE33/table.%.1fA_k%d", DATADIR, dmax, k);
		sprintf(aace_name, "AACE33_%.1fA_k%d", dmax, k);
		break;

	default:
		printf("Error: wrong AACE groups file\n");
		exit(1);
		break;

	}

	ReadGroups(gname);

	Allocate(TYPE_IDX.size());

	ReadCouplings(tname);

	name = aace_name;

}

AACE::AACE(const AACE &source) :
		ATOM_IDX(source.ATOM_IDX), TYPE_IDX(source.TYPE_IDX), h(NULL), J(NULL), dmax(
				source.dmax), name(source.name) {

	int N = TYPE_IDX.size();

	Allocate(N);

	memcpy(h, source.h, N * sizeof(double));
	for (int i = 0; i < N; i++) {
		memcpy(J[i], source.J[i], N * sizeof(double));
	}

}

AACE::AACE() :
		h(NULL), J(NULL), dmax(0.0), name("") {

}

AACE::~AACE() {

	Free(TYPE_IDX.size());

}

AACE & AACE::operator =(const AACE & source) {

	assert(this != &source); /* an attempt to assign Residue to itself */

	Free(TYPE_IDX.size());

	ATOM_IDX = source.ATOM_IDX;
	TYPE_IDX = source.TYPE_IDX;
	dmax = source.dmax;
	name = source.name;

	int N = TYPE_IDX.size();
	Allocate(N);

	memcpy(h, source.h, N * sizeof(double));
	for (int i = 0; i < N; i++) {
		memcpy(J[i], source.J[i], N * sizeof(double));
	}

	return *this;

}

void AACE::Allocate(int N) {

	h = (double*) malloc(N * sizeof(double));
	J = (double**) malloc(N * sizeof(double*));
	for (int i = 0; i < N; i++) {
		J[i] = (double*) malloc(N * sizeof(double));
	}

}

void AACE::Free(int N) {

	if (J != NULL) {
		for (int i = 0; i < N; i++) {
			free(J[i]);
		}
		free(J);
	}
	free(h);

}

int AACE::ReadGroups(const char *name) {

	FILE *F = fopen(name, "r");
	if (F == NULL) {
		printf("Error: cannot open %s file\n", name);
		exit(1);
	}

	const int STRMAX = 1000;
	char atom[10], type[10], buf[STRMAX];
	int n = 0;

	while (fgets(buf, STRMAX, F)) {
		if (buf[0] != '#') {
			if (sscanf(buf, "%s %s", atom, type) == 2) {
				std::map<std::string, int>::const_iterator it;
				it = TYPE_IDX.find(type);
				int t;
				if (it == TYPE_IDX.end()) {
					TYPE_IDX[type] = n;
					//printf("%s %d\n", type, n);
					t = n;
					n++;
				} else {
					t = it->second;
				}
				ATOM_IDX[atom] = t;
			}
		}
	}
	fclose(F);

	return TYPE_IDX.size();

}

int AACE::ReadGroups(const GROUPS_TYPE &GROUPS) {

	int n = 0;

	for (auto g: GROUPS) {
		std::map<std::string, int>::const_iterator it;
		it = TYPE_IDX.find(g.second);
		int t;
		if (it == TYPE_IDX.end()) {
			TYPE_IDX[g.second] = n;
			//printf("%s %d\n", type, n);
			t = n++;
			//n++;
		} else {
			t = it->second;
		}
		ATOM_IDX[g.first] = t;
	}

	return TYPE_IDX.size();

}

int AACE::ReadCouplings(const char *name) {

	FILE *F = fopen(name, "r");
	if (F == NULL) {
		printf("Error: cannot open %s file\n", name);
		exit(1);
	}

	const int STRMAX = 10000;
	char buf[STRMAX];

	int N = TYPE_IDX.size();

	for (int i = 0; i < N; i++) {
		fgets(buf, STRMAX, F);
		char * pch;
		pch = strtok(buf, " "); /* label */
		pch = strtok(NULL, " "); /* local field */
		h[i] = atof(pch);
//		printf("%.5e, ", h[i]);
		printf("{");
		for (int j = 0; j < N; j++) {
			pch = strtok(NULL, " ");
			J[i][j] = atof(pch);
			printf("%.5e, ", J[i][j]);
		}
		printf("}, \n");
	}
//	printf("\n");
	fclose(F);

	return 0;

}

std::vector<int> AACE::GetTypes(const Chain &C) {

	int n = 0;

	std::map<std::string, int>::iterator it;
	std::vector<int> types(C.nAtoms);

	for (int i = 0; i < C.nRes; i++) {

		std::string rname(C.residue[i].name);

		for (int j = 0; j < C.residue[i].nAtoms; j++) {

			/* ACE atom type */
			std::string name = rname + "_"
					+ std::string(C.residue[i].atom[j].name);

			it = ATOM_IDX.find(name);
			if (it != ATOM_IDX.end()) {
				types[n] = it->second;
			} else {
				types[n] = -1;
			}
			n++;
		}
	}

	return types;

}

double AACE::GetEnergy(const Chain &C, const std::vector<int> &t) {

	double E = 0.0;

	for (int i = 0; i < C.nAtoms; i++) {

		if (t[i] < 0) {
			continue;
		}

		Atom *A = C.atom[i];

		kdres *res = kd_nearest_range3f(C.kd, A->x, A->y, A->z, dmax);
		double pos[3];

		while (!kd_res_end(res)) {

			int j = ((Atom**) kd_res_item(res, pos)) - C.atom;
			if (t[j] >= 0 && i > j
			/* && C.atom[i]->residue != C.atom[j]->residue */) {
				E += J[t[i]][t[j]];
			}
			kd_res_next(res);
		}
		kd_res_free(res);
	}

	return E;

}

double AACE::GetEnergy(const Chain &C) {

	return GetEnergy(C, GetTypes(C));

}

double AACE::GetEnergy(const Complex &C, const std::vector<int> &tRec,
		const std::vector<int> &tLig) {

	double E = 0;

	for (int i = 0; i < C.Receptor.nAtoms; i++) {

		if (tRec[i] < 0) {
			continue;
		}

		Atom *A = C.Receptor.atom[i];

		/* check whether the current atom is within the bounding box
		 * of the Receptor */
		kdtree *kd = C.Ligand.kd;
		if (A->x < kd->rect->min[0] - dmax || A->y < kd->rect->min[1] - dmax
				|| A->z < kd->rect->min[2] - dmax
				|| A->x > kd->rect->max[0] + dmax
				|| A->y > kd->rect->max[1] + dmax
				|| A->z > kd->rect->max[2] + dmax) {

			continue;

		}

		kdres *res = kd_nearest_range3f(kd, A->x, A->y, A->z, dmax);
		double pos[3];

		while (!kd_res_end(res)) {

			int j = ((Atom**) kd_res_item(res, pos)) - C.Ligand.atom;
			if (tLig[j] >= 0) {
				E += J[tRec[i]][tLig[j]];
			}
			kd_res_next(res);
		}
		kd_res_free(res);
	}

	return E;

}

double AACE::GetEnergy(const Complex &C) {

	return GetEnergy(C, GetTypes(C.Receptor), GetTypes(C.Ligand));

}

std::string AACE::GetName() {

	return name;

}
