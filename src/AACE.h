/*****************************************************************************
 *
 * 2016, Ivan Anishchenko <anishchenko.ivan@gmail.com>
 *
 * Vakser Lab,
 * Center for Computational biology,
 * University of Kansas
 *
 * V20160331
 *
 *****************************************************************************/

#ifndef AACE_H_
#define AACE_H_

#include <map>
#include <string>
#include <vector>

#include "Complex.h"

class AACE {
private:

	std::map<std::string, int> ATOM_IDX; /* atom full name --> idx */
	std::map<std::string, int> TYPE_IDX; /* atom type --> idx */
	double *h;
	double **J;

	double dmax;
	//int kmin;

	std::string name;

	int ReadGroups(const char *name); /* sets ATOM_IDX, TYPE_IDX */
	int ReadCouplings(const char *name); /* sets h,J */

	std::vector<int> GetTypes(const Chain &C);

	/* types should be precalculated */
	double GetEnergy(const Chain &C, const std::vector<int> &t);
	double GetEnergy(const Complex &C, const std::vector<int> &tRec,
			const std::vector<int> &tLig);

	void Allocate(int n);
	void Free(int n);

public:

	enum AACE_TYPE {
		AACE18 = 1, /* 18 atom types */
		AACE167 = 2, /* 167 atom types - all heavy atoms in 20 standard AA */
		AACE20 = 3, /* 20 atom types - one per residue */
		AACE24 = 4, /* 24 atom types from Rosetta */
		AACE33 = 5 /* 33 atom types from CHARMM22 */
	};

	AACE(AACE_TYPE type, double d, int k);
	AACE(const AACE &source);
	AACE();

	~AACE();

	AACE& operator=(const AACE &source);

	/*
	 * per chain
	 */
	double GetEnergy(const Chain &C);

	/*
	 * TODO: between residues
	 */
	//double GetEnergy(const Residue &R1, const Residue &R2);

	/*
	 * per complex
	 */
	double GetEnergy(const Complex &C); /* interaction energy */

	std::string GetName();

};

#endif /* AACE_H_ */
