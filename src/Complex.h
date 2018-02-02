/*****************************************************************************
 *
 * Ivan Anishchenko <anishchenko.ivan@gmail.com>
 *
 * 2013-2016, Vakser Lab, Center for Computational biology, University of Kansas
 *
 * V20160215
 *
 *****************************************************************************/

#ifndef COMPLEX_H_
#define COMPLEX_H_

#include <vector>
#include <utility>

#include "Chain.h"
#include "Info.h"
#include "kdtree.h"

/* Base class for Target and Template */
class Complex {

	friend class TMdock;
	friend class Target;
	friend class Template;

protected:

	/* interface residue-residue contacts */
	std::vector<std::pair<int, int> > contacts;

	/* interface residues */
	bool *intRec; /* for the Receptor */
	bool *intLig; /* for the Ligand */

	/*
	 * functions to calculate interface residues and contacts
	 */

	/* (1) heavy atoms with const. distance cut-off */
	/* TODO: check this function */
	void SetInterface1(double d);

	/* (2) C-alphas with const. distance cut-off */
	void SetInterface2(double d);

	/* (3) C-alphas with residue-dependent distances
	 * [1]  Sh Cheng, Y Zhang, Ch L Brooks. PCalign: a method to quantify
	 *      physicochemical similarity of protein-protein interfaces.
	 *      BMC Bioinformatics. 2015 Feb 1;16(1):33 */
	void SetInterface3(double ksi);

	/* (4) centroids with const. distances */
	void SetInterface4(double d);

	/* (5) centroids with residue-dependent distances */
	void SetInterface5(double ksi);

	int intRecN, intLigN; /* number of interface residues */
	void CountIntRes(); /* sets intRecN and intLigN */

	/* sort and remove duplicate contacts */
	void ArrangeContatcs();

	/*
	 * memory management
	 */
	void AllocateComplex();
	void FreeComplex();

public:

	Chain Receptor;
	Chain Ligand;

	Complex();
	Complex(const Complex &C);
	Complex(const Chain &R, const Chain &L);
	Complex(const char *nameRec, const char *nameLig);

	~Complex();

	Complex& operator=(const Complex &C);

	/* mode = 1: 'dist' - distance between heavy atoms (default)
	 * mode = 2: 'dist' - distance between C-alphas, same to all residue-residues pairs
	 * mode = 3: residue-dependent distances taken from PCalign [1],
	 *           'dist' := 'ksi' (see eq. (1) in [1])
	 * mode = 4: similar to mode=2 but for centroids
	 * mode = 5: similar to mode=3 but for centroids */
	void SetInterface(double dist, int mode = 1);

	std::vector<std::pair<int, int> > GetContacts();

};

#endif /* COMPLEX_H_ */
