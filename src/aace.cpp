/*****************************************************************************
 *
 * Ivan Anishchenko <anishchenko.ivan@gmail.com>
 *
 * 2013-2016, Vakser Lab, Center for Computational biology, University of Kansas
 *
 *****************************************************************************
 *
 * TODO (2017May02): -L list.file - calculate energies for all structures in the list
 *
 * V20170306 - support of Rosetta- (AACE24) and CHARMM- (AACE33) based
 *             contact potentials
 * V20160331
 *
 *****************************************************************************
 *
 * Atom-atom contact energy (AACE) calculation program.
 *
 * Two modes of operation:
 *
 * (1) internal/folding (free) energy for the entire PDB structure
 *     ('-l' option is omitted)
 *
 * (2) binding (free) energy between the two proteins
 *     (both '-r' and '-l' options are used)
 *
 *****************************************************************************/

#include <cstring>

#include "AACE.h"

#define SLEN 1000

const char usage[] =
		"Usage: \n"
				"./aace -r receptor.pdb \\  # PDB file with receptor's coordinates\n"
				"       -l ligand.pdb \\    # (optional) PDB file with ligand's coordinates\n"
				"       -t AACE_TYPE \\     # AACE type: AACE18, AACE167, AACE20, AACE24, AACE33\n"
				"       -d dmax \\          # contact distance\n"
				"       -k kmin            # sequence separation\n";

bool GetArguments(char *rname, char *lname, AACE::AACE_TYPE &type, double &dmax,
		int &kmin, int argc, char *argv[]);

int main(int argc, char *argv[]) {

	/*
	 * (0) process arguments
	 */
	char rname[SLEN], lname[SLEN];
	AACE::AACE_TYPE type;
	double dmax;
	int kmin;

	if (!GetArguments(rname, lname, type, dmax, kmin, argc, argv)) {
		printf("%s", usage);
		exit(1);
	}

	/*
	 * (1) get AACE energy
	 */
	AACE AACE_(type, dmax, kmin);
	double E =
			lname[0] == '\0' ?
					AACE_.GetEnergy(Chain(rname)) :
					AACE_.GetEnergy(Complex(rname, lname));

	printf("E(%s)= %.5e\n", AACE_.GetName().c_str(), E);

	return 0;

}

bool GetArguments(char *rname, char *lname, AACE::AACE_TYPE &type, double &dmax,
		int &kmin, int argc, char *argv[]) {

	rname[0] = lname[0] = '\0';

	bool fld = false, flk = false;

	if (argc < 9 || argc > 11 || argc % 2 != 1) {
		printf("Error: invalid number of parameters\n");
		return false;
	}

	for (int i = 0; i < argc / 2; i++) {

		char *fl = argv[2 * i + 1];
		char *arg = argv[2 * i + 2];

		/* -r */
		if (strcmp(fl, "-r") == 0) {
			strcpy(rname, arg);
		} else

		/* -l */
		if (strcmp(fl, "-l") == 0) {
			strcpy(lname, arg);
		} else

		/* -t */
		if (strcmp(fl, "-t") == 0) {
			char tname[20];
			strcpy(tname, arg);
			if (strcmp(tname, "AACE18") == 0) {
				type = AACE::AACE18;
			} else if (strcmp(tname, "AACE167") == 0) {
				type = AACE::AACE167;
			} else if (strcmp(tname, "AACE20") == 0) {
				type = AACE::AACE20;
			} else if (strcmp(tname, "AACE24") == 0) {
				type = AACE::AACE24;
			} else if (strcmp(tname, "AACE33") == 0) {
				type = AACE::AACE33;
			} else {
				printf("Error: wrong AACE type '%s'\n", arg);
				return false;
			}
		} else

		/* -d */
		if (strcmp(fl, "-d") == 0) {
			dmax = atof(arg);
			fld = true;
		} else

		/* -k */
		if (strcmp(fl, "-k") == 0) {
			kmin = atoi(arg);
			flk = true;
		} else

		/* default */{
			printf("Error: wrong argument '%s %s'\n", fl, arg);
			return false;
		}

	}

	if (rname[0] == '\0') {
		printf("Error: 'receptor.pdb' not set\n");
		return false;
	}
	if (!fld) {
		printf("Error: 'dmax' not set\n");
		return false;
	}
	if (!flk) {
		printf("Error: 'kmin' not set\n");
		return false;
	}

	return true;

}
