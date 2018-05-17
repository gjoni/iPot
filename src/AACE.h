/*****************************************************************************
 *
 * Ivan Anishchenko <anishchenko.ivan@gmail.com>
 *
 * 2013-2016, Vakser Lab, Center for Computational biology, University of Kansas
 *
 * V20160331
 *
 *****************************************************************************/

#ifndef AACE_H_
#define AACE_H_

#include <map>
#include <string>
#include <vector>
#include <utility>

#include "Complex.h"

typedef std::vector<std::pair<std::string, std::string> > GROUPS_TYPE;

class AACE {
private:

	std::map<std::string, int> ATOM_IDX; /* atom full name --> idx */
	std::map<std::string, int> TYPE_IDX; /* atom type --> idx */
	double *h;
	double **J;

	double dmax;
	//int kmin;

	std::string name;

	/* set ATOM_IDX, TYPE_IDX from variable */
	int ReadGroups(const GROUPS_TYPE&);

	/* set ATOM_IDX, TYPE_IDX  from file */
	int ReadGroups(const char *name);

	/* set h,J from variables */
	int ReadCouplings(double *h_, double **J_);

	/* set h,J from file */
	int ReadCouplings(const char *name);

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

	AACE(const GROUPS_TYPE&, const std::vector<double>&, 
			const std::vector<std::vector<double> >&, double, 
			const std::string&);
	AACE(const GROUPS_TYPE&, double, const std::string&);
	AACE(AACE_TYPE type, double d, int k);
	AACE(const AACE &source);
	AACE();

	~AACE();

	AACE& operator=(const AACE &source);

	/* intrachain energy */
	double GetEnergy(const Chain &C);

	/* interaction energy */
	double GetEnergy(const Complex &C);

	std::string GetName();

};

using namespace std;

const GROUPS_TYPE AACE167_GROUPS = {
	{ "ALA_N", "ALA_N" }, { "ALA_CA", "ALA_CA" }, { "ALA_CB", "ALA_CB" }, {
				"ALA_C", "ALA_C" }, { "ALA_O", "ALA_O" }, { "ARG_N", "ARG_N" },
		{ "ARG_CA", "ARG_CA" }, { "ARG_CB", "ARG_CB" }, { "ARG_CG", "ARG_CG" },
		{ "ARG_CD", "ARG_CD" }, { "ARG_NE", "ARG_NE" }, { "ARG_CZ", "ARG_CZ" },
		{ "ARG_NH1", "ARG_NH1" }, { "ARG_NH2", "ARG_NH2" },
		{ "ARG_C", "ARG_C" }, { "ARG_O", "ARG_O" }, { "ASN_N", "ASN_N" }, {
				"ASN_CA", "ASN_CA" }, { "ASN_CB", "ASN_CB" }, { "ASN_CG",
				"ASN_CG" }, { "ASN_OD1", "ASN_OD1" }, { "ASN_ND2", "ASN_ND2" },
		{ "ASN_C", "ASN_C" }, { "ASN_O", "ASN_O" }, { "ASP_N", "ASP_N" }, {
				"ASP_CA", "ASP_CA" }, { "ASP_CB", "ASP_CB" }, { "ASP_CG",
				"ASP_CG" }, { "ASP_OD1", "ASP_OD1" }, { "ASP_OD2", "ASP_OD2" },
		{ "ASP_C", "ASP_C" }, { "ASP_O", "ASP_O" }, { "CYS_N", "CYS_N" }, {
				"CYS_CA", "CYS_CA" }, { "CYS_CB", "CYS_CB" }, { "CYS_SG",
				"CYS_SG" }, { "CYS_C", "CYS_C" }, { "CYS_O", "CYS_O" }, {
				"GLN_N", "GLN_N" }, { "GLN_CA", "GLN_CA" },
		{ "GLN_CB", "GLN_CB" }, { "GLN_CG", "GLN_CG" }, { "GLN_CD", "GLN_CD" },
		{ "GLN_OE1", "GLN_OE1" }, { "GLN_NE2", "GLN_NE2" },
		{ "GLN_C", "GLN_C" }, { "GLN_O", "GLN_O" }, { "GLU_N", "GLU_N" }, {
				"GLU_CA", "GLU_CA" }, { "GLU_CB", "GLU_CB" }, { "GLU_CG",
				"GLU_CG" }, { "GLU_CD", "GLU_CD" }, { "GLU_OE1", "GLU_OE1" }, {
				"GLU_OE2", "GLU_OE2" }, { "GLU_C", "GLU_C" },
		{ "GLU_O", "GLU_O" }, { "GLY_N", "GLY_N" }, { "GLY_CA", "GLY_CA" }, {
				"GLY_C", "GLY_C" }, { "GLY_O", "GLY_O" }, { "HIS_N", "HIS_N" },
		{ "HIS_CA", "HIS_CA" }, { "HIS_CB", "HIS_CB" }, { "HIS_CG", "HIS_CG" },
		{ "HIS_ND1", "HIS_ND1" }, { "HIS_CD2", "HIS_CD2" }, { "HIS_NE2",
				"HIS_NE2" }, { "HIS_CE1", "HIS_CE1" }, { "HIS_C", "HIS_C" }, {
				"HIS_O", "HIS_O" }, { "ILE_N", "ILE_N" },
		{ "ILE_CA", "ILE_CA" }, { "ILE_CB", "ILE_CB" },
		{ "ILE_CG2", "ILE_CG2" }, { "ILE_CG1", "ILE_CG1" }, { "ILE_CD1",
				"ILE_CD1" }, { "ILE_C", "ILE_C" }, { "ILE_O", "ILE_O" }, {
				"LEU_N", "LEU_N" }, { "LEU_CA", "LEU_CA" },
		{ "LEU_CB", "LEU_CB" }, { "LEU_CG", "LEU_CG" },
		{ "LEU_CD1", "LEU_CD1" }, { "LEU_CD2", "LEU_CD2" },
		{ "LEU_C", "LEU_C" }, { "LEU_O", "LEU_O" }, { "LYS_N", "LYS_N" }, {
				"LYS_CA", "LYS_CA" }, { "LYS_CB", "LYS_CB" }, { "LYS_CG",
				"LYS_CG" }, { "LYS_CD", "LYS_CD" }, { "LYS_CE", "LYS_CE" }, {
				"LYS_NZ", "LYS_NZ" }, { "LYS_C", "LYS_C" },
		{ "LYS_O", "LYS_O" }, { "MET_N", "MET_N" }, { "MET_CA", "MET_CA" }, {
				"MET_CB", "MET_CB" }, { "MET_CG", "MET_CG" }, { "MET_SD",
				"MET_SD" }, { "MET_CE", "MET_CE" }, { "MET_C", "MET_C" }, {
				"MET_O", "MET_O" }, { "PHE_N", "PHE_N" },
		{ "PHE_CA", "PHE_CA" }, { "PHE_CB", "PHE_CB" }, { "PHE_CG", "PHE_CG" },
		{ "PHE_CD1", "PHE_CD1" }, { "PHE_CD2", "PHE_CD2" }, { "PHE_CE1",
				"PHE_CE1" }, { "PHE_CE2", "PHE_CE2" }, { "PHE_CZ", "PHE_CZ" }, {
				"PHE_C", "PHE_C" }, { "PHE_O", "PHE_O" }, { "PRO_N", "PRO_N" },
		{ "PRO_CD", "PRO_CD" }, { "PRO_CA", "PRO_CA" }, { "PRO_CB", "PRO_CB" },
		{ "PRO_CG", "PRO_CG" }, { "PRO_C", "PRO_C" }, { "PRO_O", "PRO_O" }, {
				"SER_N", "SER_N" }, { "SER_CA", "SER_CA" },
		{ "SER_CB", "SER_CB" }, { "SER_OG", "SER_OG" }, { "SER_C", "SER_C" }, {
				"SER_O", "SER_O" }, { "THR_N", "THR_N" },
		{ "THR_CA", "THR_CA" }, { "THR_CB", "THR_CB" },
		{ "THR_OG1", "THR_OG1" }, { "THR_CG2", "THR_CG2" },
		{ "THR_C", "THR_C" }, { "THR_O", "THR_O" }, { "TRP_N", "TRP_N" }, {
				"TRP_CA", "TRP_CA" }, { "TRP_CB", "TRP_CB" }, { "TRP_CG",
				"TRP_CG" }, { "TRP_CD2", "TRP_CD2" }, { "TRP_CE2", "TRP_CE2" },
		{ "TRP_CE3", "TRP_CE3" }, { "TRP_CD1", "TRP_CD1" }, { "TRP_NE1",
				"TRP_NE1" }, { "TRP_CZ2", "TRP_CZ2" }, { "TRP_CZ3", "TRP_CZ3" },
		{ "TRP_CH2", "TRP_CH2" }, { "TRP_C", "TRP_C" }, { "TRP_O", "TRP_O" }, {
				"TYR_N", "TYR_N" }, { "TYR_CA", "TYR_CA" },
		{ "TYR_CB", "TYR_CB" }, { "TYR_CG", "TYR_CG" },
		{ "TYR_CD1", "TYR_CD1" }, { "TYR_CE1", "TYR_CE1" }, { "TYR_CD2",
				"TYR_CD2" }, { "TYR_CE2", "TYR_CE2" }, { "TYR_CZ", "TYR_CZ" }, {
				"TYR_OH", "TYR_OH" }, { "TYR_C", "TYR_C" },
		{ "TYR_O", "TYR_O" }, { "VAL_N", "VAL_N" }, { "VAL_CA", "VAL_CA" }, {
				"VAL_CB", "VAL_CB" }, { "VAL_CG1", "VAL_CG1" }, { "VAL_CG2",
				"VAL_CG2" }, { "VAL_C", "VAL_C" }, { "VAL_O", "VAL_O" } };

const GROUPS_TYPE AACE20_GROUPS = { { "ALA_N", "ALA" }, {
		"ALA_CA", "ALA" }, { "ALA_CB", "ALA" }, { "ALA_C", "ALA" }, { "ALA_O",
		"ALA" }, { "ARG_N", "ARG" }, { "ARG_CA", "ARG" }, { "ARG_CB", "ARG" }, {
		"ARG_CG", "ARG" }, { "ARG_CD", "ARG" }, { "ARG_NE", "ARG" }, { "ARG_CZ",
		"ARG" }, { "ARG_NH1", "ARG" }, { "ARG_NH2", "ARG" }, { "ARG_C", "ARG" },
		{ "ARG_O", "ARG" }, { "ASN_N", "ASN" }, { "ASN_CA", "ASN" }, { "ASN_CB",
				"ASN" }, { "ASN_CG", "ASN" }, { "ASN_OD1", "ASN" }, { "ASN_ND2",
				"ASN" }, { "ASN_C", "ASN" }, { "ASN_O", "ASN" }, { "ASP_N",
				"ASP" }, { "ASP_CA", "ASP" }, { "ASP_CB", "ASP" }, { "ASP_CG",
				"ASP" }, { "ASP_OD1", "ASP" }, { "ASP_OD2", "ASP" }, { "ASP_C",
				"ASP" }, { "ASP_O", "ASP" }, { "CYS_N", "CYS" }, { "CYS_CA",
				"CYS" }, { "CYS_CB", "CYS" }, { "CYS_SG", "CYS" }, { "CYS_C",
				"CYS" }, { "CYS_O", "CYS" }, { "GLN_N", "GLN" }, { "GLN_CA",
				"GLN" }, { "GLN_CB", "GLN" }, { "GLN_CG", "GLN" }, { "GLN_CD",
				"GLN" }, { "GLN_OE1", "GLN" }, { "GLN_NE2", "GLN" }, { "GLN_C",
				"GLN" }, { "GLN_O", "GLN" }, { "GLU_N", "GLU" }, { "GLU_CA",
				"GLU" }, { "GLU_CB", "GLU" }, { "GLU_CG", "GLU" }, { "GLU_CD",
				"GLU" }, { "GLU_OE1", "GLU" }, { "GLU_OE2", "GLU" }, { "GLU_C",
				"GLU" }, { "GLU_O", "GLU" }, { "GLY_N", "GLY" }, { "GLY_CA",
				"GLY" }, { "GLY_C", "GLY" }, { "GLY_O", "GLY" }, { "HIS_N",
				"HIS" }, { "HIS_CA", "HIS" }, { "HIS_CB", "HIS" }, { "HIS_CG",
				"HIS" }, { "HIS_ND1", "HIS" }, { "HIS_CD2", "HIS" }, {
				"HIS_NE2", "HIS" }, { "HIS_CE1", "HIS" }, { "HIS_C", "HIS" }, {
				"HIS_O", "HIS" }, { "ILE_N", "ILE" }, { "ILE_CA", "ILE" }, {
				"ILE_CB", "ILE" }, { "ILE_CG2", "ILE" }, { "ILE_CG1", "ILE" }, {
				"ILE_CD1", "ILE" }, { "ILE_C", "ILE" }, { "ILE_O", "ILE" }, {
				"LEU_N", "LEU" }, { "LEU_CA", "LEU" }, { "LEU_CB", "LEU" }, {
				"LEU_CG", "LEU" }, { "LEU_CD1", "LEU" }, { "LEU_CD2", "LEU" }, {
				"LEU_C", "LEU" }, { "LEU_O", "LEU" }, { "LYS_N", "LYS" }, {
				"LYS_CA", "LYS" }, { "LYS_CB", "LYS" }, { "LYS_CG", "LYS" }, {
				"LYS_CD", "LYS" }, { "LYS_CE", "LYS" }, { "LYS_NZ", "LYS" }, {
				"LYS_C", "LYS" }, { "LYS_O", "LYS" }, { "MET_N", "MET" }, {
				"MET_CA", "MET" }, { "MET_CB", "MET" }, { "MET_CG", "MET" }, {
				"MET_SD", "MET" }, { "MET_CE", "MET" }, { "MET_C", "MET" }, {
				"MET_O", "MET" }, { "PHE_N", "PHE" }, { "PHE_CA", "PHE" }, {
				"PHE_CB", "PHE" }, { "PHE_CG", "PHE" }, { "PHE_CD1", "PHE" }, {
				"PHE_CD2", "PHE" }, { "PHE_CE1", "PHE" }, { "PHE_CE2", "PHE" },
		{ "PHE_CZ", "PHE" }, { "PHE_C", "PHE" }, { "PHE_O", "PHE" }, { "PRO_N",
				"PRO" }, { "PRO_CD", "PRO" }, { "PRO_CA", "PRO" }, { "PRO_CB",
				"PRO" }, { "PRO_CG", "PRO" }, { "PRO_C", "PRO" }, { "PRO_O",
				"PRO" }, { "SER_N", "SER" }, { "SER_CA", "SER" }, { "SER_CB",
				"SER" }, { "SER_OG", "SER" }, { "SER_C", "SER" }, { "SER_O",
				"SER" }, { "THR_N", "THR" }, { "THR_CA", "THR" }, { "THR_CB",
				"THR" }, { "THR_OG1", "THR" }, { "THR_CG2", "THR" }, { "THR_C",
				"THR" }, { "THR_O", "THR" }, { "TRP_N", "TRP" }, { "TRP_CA",
				"TRP" }, { "TRP_CB", "TRP" }, { "TRP_CG", "TRP" }, { "TRP_CD2",
				"TRP" }, { "TRP_CE2", "TRP" }, { "TRP_CE3", "TRP" }, {
				"TRP_CD1", "TRP" }, { "TRP_NE1", "TRP" }, { "TRP_CZ2", "TRP" },
		{ "TRP_CZ3", "TRP" }, { "TRP_CH2", "TRP" }, { "TRP_C", "TRP" }, {
				"TRP_O", "TRP" }, { "TYR_N", "TYR" }, { "TYR_CA", "TYR" }, {
				"TYR_CB", "TYR" }, { "TYR_CG", "TYR" }, { "TYR_CD1", "TYR" }, {
				"TYR_CE1", "TYR" }, { "TYR_CD2", "TYR" }, { "TYR_CE2", "TYR" },
		{ "TYR_CZ", "TYR" }, { "TYR_OH", "TYR" }, { "TYR_C", "TYR" }, { "TYR_O",
				"TYR" }, { "VAL_N", "VAL" }, { "VAL_CA", "VAL" }, { "VAL_CB",
				"VAL" }, { "VAL_CG1", "VAL" }, { "VAL_CG2", "VAL" }, { "VAL_C",
				"VAL" }, { "VAL_O", "VAL" } };

const GROUPS_TYPE AACE18_GROUPS = { { "ALA_N", "N" }, {
		"ALA_CA", "CA" }, { "ALA_CB", "CB" }, { "ALA_C", "C" },
		{ "ALA_O", "O" }, { "ARG_N", "N" }, { "ARG_CA", "CA" },
		{ "ARG_CB", "CB" }, { "ARG_CG", "FC" }, { "ARG_CD", "RNE" }, { "ARG_NE",
				"RNE" }, { "ARG_CZ", "RNZ" }, { "ARG_NH1", "RNZ" }, { "ARG_NH2",
				"RNZ" }, { "ARG_C", "C" }, { "ARG_O", "O" }, { "ASN_N", "N" }, {
				"ASN_CA", "CA" }, { "ASN_CB", "CB" }, { "ASN_CG", "NN" }, {
				"ASN_OD1", "NN" }, { "ASN_ND2", "NN" }, { "ASN_C", "C" }, {
				"ASN_O", "O" }, { "ASP_N", "N" }, { "ASP_CA", "CA" }, {
				"ASP_CB", "CB" }, { "ASP_CG", "DO" }, { "ASP_OD1", "DO" }, {
				"ASP_OD2", "DO" }, { "ASP_C", "C" }, { "ASP_O", "O" }, {
				"CYS_N", "N" }, { "CYS_CA", "CA" }, { "CYS_CB", "CB" }, {
				"CYS_SG", "CS" }, { "CYS_C", "C" }, { "CYS_O", "O" }, { "GLN_N",
				"N" }, { "GLN_CA", "CA" }, { "GLN_CB", "CB" },
		{ "GLN_CG", "FC" }, { "GLN_CD", "NN" }, { "GLN_OE1", "NN" }, {
				"GLN_NE2", "NN" }, { "GLN_C", "C" }, { "GLN_O", "O" }, {
				"GLU_N", "N" }, { "GLU_CA", "CA" }, { "GLU_CB", "CB" }, {
				"GLU_CG", "FC" }, { "GLU_CD", "DO" }, { "GLU_OE1", "DO" }, {
				"GLU_OE2", "DO" }, { "GLU_C", "C" }, { "GLU_O", "O" }, {
				"GLY_N", "N" }, { "GLY_CA", "GC" }, { "GLY_C", "C" }, { "GLY_O",
				"O" }, { "HIS_N", "N" }, { "HIS_CA", "CA" }, { "HIS_CB", "CB" },
		{ "HIS_CG", "HN" }, { "HIS_ND1", "HN" }, { "HIS_CD2", "HN" }, {
				"HIS_NE2", "HN" }, { "HIS_CE1", "HN" }, { "HIS_C", "C" }, {
				"HIS_O", "O" }, { "ILE_N", "N" }, { "ILE_CA", "CA" }, {
				"ILE_CB", "CB" }, { "ILE_CG2", "LC" }, { "ILE_CG1", "FC" }, {
				"ILE_CD", "LC" }, { "ILE_CD1", "LC" }, { "ILE_C", "C" }, {
				"ILE_O", "O" }, { "LEU_N", "N" }, { "LEU_CA", "CA" }, {
				"LEU_CB", "CB" }, { "LEU_CG", "FC" }, { "LEU_CD1", "LC" }, {
				"LEU_CD2", "LC" }, { "LEU_C", "C" }, { "LEU_O", "O" }, {
				"LYS_N", "N" }, { "LYS_CA", "CA" }, { "LYS_CB", "CB" }, {
				"LYS_CG", "FC" }, { "LYS_CD", "KC" }, { "LYS_CE", "KN" }, {
				"LYS_NZ", "KN" }, { "LYS_C", "C" }, { "LYS_O", "O" }, { "MET_N",
				"N" }, { "MET_CA", "CA" }, { "MET_CB", "CB" },
		{ "MET_CG", "FC" }, { "MET_SD", "FC" }, { "MET_CE", "LC" }, { "MET_C",
				"C" }, { "MET_O", "O" }, { "PHE_N", "N" }, { "PHE_CA", "CA" }, {
				"PHE_CB", "CB" }, { "PHE_CG", "FC" }, { "PHE_CD1", "FC" }, {
				"PHE_CD2", "FC" }, { "PHE_CE1", "FC" }, { "PHE_CE2", "FC" }, {
				"PHE_CZ", "FC" }, { "PHE_C", "C" }, { "PHE_O", "O" }, { "PRO_N",
				"N" }, { "PRO_CD", "CB" }, { "PRO_CA", "CA" },
		{ "PRO_CB", "CB" }, { "PRO_CG", "CB" }, { "PRO_C", "C" },
		{ "PRO_O", "O" }, { "SER_N", "N" }, { "SER_CA", "CA" },
		{ "SER_CB", "SO" }, { "SER_OG", "SO" }, { "SER_C", "C" },
		{ "SER_O", "O" }, { "THR_N", "N" }, { "THR_CA", "CA" },
		{ "THR_CB", "CB" }, { "THR_OG1", "SO" }, { "THR_CG2", "FC" }, { "THR_C",
				"C" }, { "THR_O", "O" }, { "TRP_N", "N" }, { "TRP_CA", "CA" }, {
				"TRP_CB", "CB" }, { "TRP_CG", "FC" }, { "TRP_CD2", "FC" }, {
				"TRP_CE2", "FC" }, { "TRP_CE3", "FC" }, { "TRP_CD1", "FC" }, {
				"TRP_NE1", "HN" }, { "TRP_CZ2", "FC" }, { "TRP_CZ3", "FC" }, {
				"TRP_CH2", "FC" }, { "TRP_C", "C" }, { "TRP_O", "O" }, {
				"TYR_N", "N" }, { "TYR_CA", "CA" }, { "TYR_CB", "CB" }, {
				"TYR_CG", "FC" }, { "TYR_CD1", "FC" }, { "TYR_CE1", "YC" }, {
				"TYR_CD2", "FC" }, { "TYR_CE2", "YC" }, { "TYR_CZ", "YC" }, {
				"TYR_OH", "SO" }, { "TYR_C", "C" }, { "TYR_O", "O" }, { "VAL_N",
				"N" }, { "VAL_CA", "CA" }, { "VAL_CB", "CB" },
		{ "VAL_CG1", "LC" }, { "VAL_CG2", "LC" }, { "VAL_C", "C" }, { "VAL_O",
				"O" }, { "ALA_OXT", "O" }, { "ARG_OXT", "O" },
		{ "ASN_OXT", "O" }, { "ASP_OXT", "O" }, { "CYS_OXT", "O" }, { "GLN_OXT",
				"O" }, { "GLU_OXT", "O" }, { "GLY_OXT", "O" },
		{ "HIS_OXT", "O" }, { "ILE_OXT", "O" }, { "LEU_OXT", "O" }, { "LYS_OXT",
				"O" }, { "MET_OXT", "O" }, { "PHE_OXT", "O" },
		{ "PRO_OXT", "O" }, { "SER_OXT", "O" }, { "THR_OXT", "O" }, { "TRP_OXT",
				"O" }, { "TYR_OXT", "O" }, { "VAL_OXT", "O" } };

#endif /* AACE_H_ */
