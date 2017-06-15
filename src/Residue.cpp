/*
 * Residue.cpp
 *
 *  Created on: Jun 3, 2014
 *      Author: ivan
 */

#include "Residue.h"
#include "Info.h"

#include <cstring>
#include <cstdlib>
#include <cassert>
#include <new>

Residue::Residue() :
		type(23), seqNum(-1), insCode(' '), atom(NULL), nAtoms(0), chainId(' '), ntFlag(
				false), ctFlag(false), N(NULL), CA(NULL), C(NULL), O(NULL) {

	strcpy(name, "UNK");

	centroid[0] = centroid[1] = centroid[2] = 0;

}

Residue::Residue(const Residue &source) :
		type(source.type), seqNum(source.seqNum), insCode(source.insCode), atom(
		NULL), nAtoms(source.nAtoms), chainId(source.chainId), ntFlag(
				source.ntFlag), ctFlag(source.ctFlag), N(NULL), CA(NULL), O(
		NULL) {

	assert(nAtoms > 0); /* pedantic check */

	strcpy(name, source.name);

	atom = new Atom[nAtoms];
	for (int i = 0; i < nAtoms; i++) {
		atom[i] = source.atom[i];
		atom[i].residue = this;
	}

	SetBBAtoms(); /* N, CA, C, O are pointers - they need  to be reinitialized */

	memcpy(centroid, source.centroid, 3 * sizeof(double));

}

Residue::Residue(const std::vector<AtomRecord>& source) :
		seqNum(source[0].resNum), insCode(source[0].insCode), atom(NULL), nAtoms(
				source.size()), chainId(source[0].chainId), ntFlag(false), ctFlag(
				false), N(NULL), CA(NULL), C(
		NULL), O(NULL) {

	assert(nAtoms > 0); /* pedantic check */

	strcpy(name, source[0].resName); /* first atom is used to name the residue*/

	atom = new Atom[nAtoms];
	int idx = 0;
	for (int i = 0; i < nAtoms; i++) {
		atom[i] = Atom(source[i]);
		atom[i].residue = this;
		idx++;
	}

	SetBBAtoms(); /* N, CA, C , O are pointers - they need to be reinitialized */

	type = SetType(name);

	SetCentroid();

}

Residue::Residue(std::vector<AtomRecord>::iterator begin_it,
		std::vector<AtomRecord>::iterator end_it) :
		seqNum(begin_it->resNum), insCode(begin_it->insCode), atom(NULL), nAtoms(
				end_it - begin_it + 1), chainId(begin_it->chainId), ntFlag(
				false), ctFlag(false), N(NULL), CA(
		NULL), C(NULL), O(NULL) {

	assert(nAtoms > 0); /* pedantic check */

	strcpy(name, begin_it->resName); /* first atom is used to name the residue*/

	atom = new Atom[nAtoms];
	int idx = 0;
	std::vector<AtomRecord>::iterator it;
	for (it = begin_it; it <= end_it; it++) {
		atom[idx] = Atom(*it);
		atom[idx].residue = this;
		idx++;
	}

	SetBBAtoms(); /* N, CA, C, O are pointers - they need to be reinitialized */

	type = SetType(name);

	SetCentroid();

}

Residue::~Residue() {

	delete[] atom;
	atom = NULL;

}

Residue & Residue::operator =(const Residue & source) {

	assert(this != &source); /* an attempt to assign Residue to itself */

	type = source.type;
	strcpy(name, source.name);
	seqNum = source.seqNum;
	insCode = source.insCode;
	nAtoms = source.nAtoms;
	chainId = source.chainId;

	delete[] atom;

	atom = new Atom[nAtoms];
	for (int i = 0; i < nAtoms; i++) {
		atom[i] = source.atom[i];
		atom[i].residue = this;
	}

	SetBBAtoms();

	memcpy(centroid, source.centroid, 3 * sizeof(double));

	return *this;

}

int Residue::SetBBAtoms() {

	bool flN = false, flCA = false, flC = false, flO = false;

	for (int i = 0; i < nAtoms; i++) {
		if (strcmp(atom[i].name, "N") == 0) {
			N = &(atom[i]);
			flN = true;
		} else if (strcmp(atom[i].name, "CA") == 0) {
			CA = &(atom[i]);
			flCA = true;
		} else if (strcmp(atom[i].name, "C") == 0) {
			C = &(atom[i]);
			flC = true;
		} else if (strcmp(atom[i].name, "O") == 0
				|| strcmp(atom[i].name, "OXT") == 0) {
			O = &(atom[i]);
			flO = true;
		} else if (strcmp(atom[i].name, "H1") == 0) {
			ntFlag = true;
		} else if (strcmp(atom[i].name, "H2") == 0) {
			ntFlag = true;
		} else if (strcmp(atom[i].name, "H3") == 0) {
			ntFlag = true;
		} else if (strcmp(atom[i].name, "OXT") == 0) {
			ctFlag = true;
		}
	}

	int n = flN + flCA + flC + flO;

	return n;

}

char Residue::SetType(const char* name) {

	for (int i = 0; i < 24; i++) {
		if (strcmp(name, AAA3[i]) == 0) {
			return i;
		}
	}
	return 23;
}

void Residue::SetCentroid() {

	centroid[0] = centroid[1] = centroid[2] = 0.0;

	if (type == 7 /* GLY */) {
		centroid[0] = CA->x;
		centroid[1] = CA->y;
		centroid[2] = CA->z;
	} else /* other residues */{
		for (int j = 0; j < nAtoms; j++) {
			centroid[0] += atom[j].x;
			centroid[1] += atom[j].y;
			centroid[2] += atom[j].z;
		}
		centroid[0] /= nAtoms;
		centroid[1] /= nAtoms;
		centroid[2] /= nAtoms;
	}

}
