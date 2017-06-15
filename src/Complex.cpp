/*
 * Complex.cpp
 *
 *  Created on: May 27, 2015
 *      Author: ivan
 */

#include "Complex.h"

#include <cstring>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <algorithm>

Complex::Complex() :
		intRec(NULL), intLig(NULL), intRecN(0), intLigN(0), Receptor(), Ligand() {

	/* nothing to do */

}

Complex::Complex(const Complex & C) :
		contacts(C.contacts), intRec(NULL), intLig(NULL), intRecN(C.intRecN), intLigN(
				C.intLigN), Receptor(C.Receptor), Ligand(C.Ligand) {

	AllocateComplex();

	/* copy interface residues */
	memcpy(intRec, C.intRec, Receptor.nRes * sizeof(bool));
	memcpy(intLig, C.intLig, Ligand.nRes * sizeof(bool));

}

Complex::Complex(const Chain & R, const Chain & L) :
		intRec(NULL), intLig(NULL), intRecN(0), intLigN(0), Receptor(R), Ligand(
				L) {

	AllocateComplex();

}

Complex::Complex(const char *nameRec, const char *nameLig) :
		intRec(NULL), intLig(NULL), intRecN(0), intLigN(0), Receptor(nameRec), Ligand(
				nameLig) {

	AllocateComplex();

}

Complex::~Complex() {

	FreeComplex();

}

Complex & Complex::operator =(const Complex & C) {

	assert(this != &C); /* error: an attempt to assign Complex to itself */

	FreeComplex();

	/* copy proteins */
	Receptor = C.Receptor;
	Ligand = C.Ligand;

	AllocateComplex();

	/* copy interface residues and contacts */
	contacts = C.contacts;
	memcpy(intRec, C.intRec, Receptor.nRes * sizeof(bool));
	memcpy(intLig, C.intLig, Ligand.nRes * sizeof(bool));
	intRecN = C.intRecN;
	intLigN = C.intLigN;

	return *this;

}

void Complex::SetInterface1(double d) {

	/* search for interface residues */
	kdres *res;
	double pos[3];
	for (int i = 0; i < Ligand.nRes; i++) {

		for (int j = 0; j < Ligand.residue[i].nAtoms; j++) {

			Atom *A = &(Ligand.residue[i].atom[j]);

			if (A->type == 'H') { /* exclude hydrogens */
				continue;
			}

			/* check whether the current atom is within the bounding box
			 * of the Receptor */
			kdtree *kd = Receptor.kd;
			if (A->x < kd->rect->min[0] - d || A->y < kd->rect->min[1] - d
					|| A->z < kd->rect->min[2] - d
					|| A->x > kd->rect->max[0] + d
					|| A->y > kd->rect->max[1] + d
					|| A->z > kd->rect->max[2] + d) {

				continue;

			}

			/* find neighbors of A in the Receptor */
			res = kd_nearest_range3f(kd, A->x, A->y, A->z, d);
			while (!kd_res_end(res)) {

				Atom *B = *((Atom**) kd_res_item(res, pos));
				if (B->type != 'H') { /* exclude hydrogens */
					int idx = B->residue - Receptor.residue;
					intRec[idx] = true;
					intLig[i] = true;
					contacts.push_back(std::make_pair(idx, i));
				}
				kd_res_next(res);
			}
			kd_res_free(res);

		}
	}

	ArrangeContatcs();

}

void Complex::SetInterface2(double d) {

	/* search for interface residues */
	kdres *res;
	double pos[3];
	for (int i = 0; i < Ligand.nRes; i++) {

		Atom *A = Ligand.residue[i].CA;

		if (A == NULL) {
			continue;
		}

		/* check whether the current atom is within the bounding box
		 * of the Receptor */
		kdtree *kdCA = Receptor.kdCA;
		if (A->x < kdCA->rect->min[0] - d || A->y < kdCA->rect->min[1] - d
				|| A->z < kdCA->rect->min[2] - d
				|| A->x > kdCA->rect->max[0] + d
				|| A->y > kdCA->rect->max[1] + d
				|| A->z > kdCA->rect->max[2] + d) {

			continue;

		}

		/* find neighbors of A in the Receptor */
		res = kd_nearest_range3f(kdCA, A->x, A->y, A->z, d);
		while (!kd_res_end(res)) {
			//int idx = (Residue*) kd_res_item(res, pos) - Receptor.residues;
			Atom *A = (Atom*) kd_res_item(res, pos);
			int idx = A->residue - Receptor.residue;
			intRec[idx] = true;
			intLig[i] = true;
			contacts.push_back(std::make_pair(idx, i));
			kd_res_next(res);
		}
		kd_res_free(res);

	}

	ArrangeContatcs();

}

void Complex::SetInterface3(double ksi) {

	/* maximal cut-off at a given ksi */
	double dmax = 0.0;
	for (int i = 0; i < 20; i++) {
		for (int j = i + 1; j < 20; j++) {
			double d = PCALIGN_DIST[i][j] + ksi * PCALIGN_STD[i][j];
			if (d > dmax) {
				dmax = d;
			}
		}
	}

	/* search for interface residues */
	kdres *res;
	double pos[3];
	double x, y, z, d;
	for (int i = 0; i < Ligand.nRes; i++) {

		Atom *A = Ligand.residue[i].CA;

		if (A == NULL) {
			continue;
		}

		/* check whether the current atom is within the bounding box
		 * of the Receptor */
		kdtree *kdCA = Receptor.kdCA;
		if (A->x < kdCA->rect->min[0] - dmax || A->y < kdCA->rect->min[1] - dmax
				|| A->z < kdCA->rect->min[2] - dmax
				|| A->x > kdCA->rect->max[0] + dmax
				|| A->y > kdCA->rect->max[1] + dmax
				|| A->z > kdCA->rect->max[2] + dmax) {

			continue;

		}

		/* find neighbors of A in the Receptor */
		res = kd_nearest_range3f(kdCA, A->x, A->y, A->z, dmax);
		while (!kd_res_end(res)) {
			Atom *A2 = (Atom*) kd_res_item(res, pos);
			int idxRec = A2->residue - Receptor.residue;
			int type1 = A2->residue->type;
			int type2 = Ligand.residue[i].type;

			if (type1 < 20 && type2 < 20) { /* process standard residues only */

				/* calculate residue-dependent cut-off */
				d = PCALIGN_DIST[type1][type2]
						+ ksi * PCALIGN_STD[type1][type2];

				/* find actual Ca-Ca distance */
				x = A->x - pos[0];
				y = A->y - pos[1];
				z = A->z - pos[2];

				if (sqrt(x * x + y * y + z * z) < d) {
					intRec[idxRec] = true;
					intLig[i] = true;
					contacts.push_back(std::make_pair(idxRec, i));
				}
			}

			kd_res_next(res);
		}
		kd_res_free(res);

	}

	ArrangeContatcs();

}

void Complex::SetInterface4(double d) {

	/* search for interface residues */
	double pos[3];
	kdtree *kd = Receptor.kdCent;
	double *min = kd->rect->min;
	double *max = kd->rect->max;
	for (int i = 0; i < Ligand.nRes; i++) {

		double *xyz = Ligand.residue[i].centroid;

		/* check whether the current atom is within the bounding box
		 * of the Receptor */
		if (xyz[0] < min[0] - d || xyz[1] < min[1] - d || xyz[2] < min[2] - d
				|| xyz[0] > max[0] + d || xyz[1] > max[1] + d
				|| xyz[2] > max[2] + d) {

			continue;

		}

		/* find neighbors of A in the Receptor */
		kdres *res = kd_nearest_range(kd, xyz, d);
		while (!kd_res_end(res)) {
			//int idxRec = (Residue*) kd_res_item(res, pos) - Receptor.residues;
			Atom *A = (Atom*) kd_res_item(res, pos);
			int idxRec = A->residue - Receptor.residue;
			intRec[idxRec] = true;
			intLig[i] = true;
			contacts.push_back(std::make_pair(idxRec, i));
			kd_res_next(res);
		}
		kd_res_free(res);

	}

	ArrangeContatcs();

}

void Complex::SetInterface5(double ksi) {

	/* maximal cut-off at a given ksi */
	double dmax = 0.0;
	for (int i = 0; i < 20; i++) {
		for (int j = i + 1; j < 20; j++) {
			double d = CENTROID_DIST[i][j] + ksi * CENTROID_STD[i][j];
			if (d > dmax) {
				dmax = d;
			}
		}
	}

	/* search for interface residues */
	double pos[3];
	kdtree *kd = Receptor.kdCent;
	double *min = kd->rect->min;
	double *max = kd->rect->max;
	for (int i = 0; i < Ligand.nRes; i++) {

		double *xyz = Ligand.residue[i].centroid;

		/* check whether the current atom is within the bounding box
		 * of the Receptor */
		if (xyz[0] < min[0] - dmax || xyz[1] < min[1] - dmax
				|| xyz[2] < min[2] - dmax || xyz[0] > max[0] + dmax
				|| xyz[1] > max[1] + dmax || xyz[2] > max[2] + dmax) {

			continue;

		}

		/* find neighbors of A in the Receptor */
		kdres *res = kd_nearest_range(kd, xyz, dmax);
		while (!kd_res_end(res)) {
			//int idxRec = (Residue*) kd_res_item(res, pos) - Receptor.residues;
			Atom *A = (Atom*) kd_res_item(res, pos);
			int idxRec = A->residue - Receptor.residue;
			int type1 = Receptor.residue[idxRec].type;
			int type2 = Ligand.residue[i].type;

			if (type1 < 20 && type2 < 20) { /* process standard residues only */

				/* calculate residue-dependent cut-off */
				double d = CENTROID_DIST[type1][type2]
						+ ksi * CENTROID_STD[type1][type2];

				/* find actual Ca-Ca distance */
				double x = xyz[0] - pos[0];
				double y = xyz[1] - pos[1];
				double z = xyz[2] - pos[2];

				if (sqrt(x * x + y * y + z * z) < d) {
					intRec[idxRec] = true;
					intLig[i] = true;
					contacts.push_back(std::make_pair(idxRec, i));
				}
			}

			kd_res_next(res);
		}
		kd_res_free(res);

	}

	ArrangeContatcs();

}

void Complex::CountIntRes() {

	intRecN = 0;
	for (int i = 0; i < Receptor.nRes; i++) {
		if (intRec[i]) {
			intRecN++;
		}
	}

	intLigN = 0;
	for (int i = 0; i < Ligand.nRes; i++) {
		if (intLig[i]) {
			intLigN++;
		}
	}

}

void Complex::ArrangeContatcs() {

	std::sort(contacts.begin(), contacts.end());
	std::vector<std::pair<int, int> >::iterator it;
	it = std::unique(contacts.begin(), contacts.end());
	contacts.resize(it - contacts.begin());

}

void Complex::AllocateComplex() {

	intRec = (bool*) malloc(Receptor.nRes * sizeof(bool));
	intLig = (bool*) malloc(Ligand.nRes * sizeof(bool));

}

void Complex::FreeComplex() {

	free(intRec);
	free(intLig);

}

void Complex::SetInterface(double dist, int mode) {

	contacts.clear();
	memset(intRec, 0, Receptor.nRes * sizeof(bool));
	memset(intLig, 0, Ligand.nRes * sizeof(bool));

	if (mode == 1) {
		SetInterface1(dist);
	} else if (mode == 2) {
		SetInterface2(dist);
	} else if (mode == 3) {
		SetInterface3(dist);
	} else if (mode == 4) {
		SetInterface4(dist);
	} else if (mode == 5) {
		SetInterface5(dist);
	} else {
		assert(0); /* error: wrong mode for interface extraction */
	}

	CountIntRes();

}

std::vector<std::pair<int, int> > Complex::GetContacts() {

	return contacts;

}

