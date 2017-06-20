/*****************************************************************************
 *
 * 2016, Ivan Anishchenko <anishchenko.ivan@gmail.com>
 *
 * Vakser Lab,
 * Center for Computational biology,
 * The University of Kansas
 *
 * V20160527
 *
 *****************************************************************************
 *
 * RRCE energy calculation program.
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

//#include "AACE.h"
#include "RRCE.h"
//#include "EigenRRCE.h"

#define SLEN 1000

const char usage[] =
		"Usage: \n"
				"./rrce -r receptor.pdb \\  # PDB file with receptor's coordinates\n"
				"       -l ligand.pdb \\    # (optional) PDB file with ligand's coordinates\n"
				"       -t AACE_TYPE \\     # RRCE type: RRCE20RC, RRCE20SCC (not impl.), RRCE20CB (not impl.) \n"
				"       -d dmax \\          # contact distance\n"
				"       -k kmin            # sequence separation\n";

bool GetArguments(char *rname, char *lname, RRCE::RRCE_TYPE &type, double &dmax,
		int &kmin, int argc, char *argv[]);

int main(int argc, char *argv[]) {

	/*
	 * (0) process arguments
	 */
	char rname[SLEN], lname[SLEN];
	RRCE::RRCE_TYPE type;
	double dmax;
	int kmin;

	if (!GetArguments(rname, lname, type, dmax, kmin, argc, argv)) {
		printf("%s", usage);
		exit(1);
	}

	/*
	 * (1) get AACE energy
	 */
	RRCE RRCE_(type, dmax, kmin);
	double E =
			lname[0] == '\0' ?
					RRCE_.GetEnergy(Chain(rname)) :
					RRCE_.GetEnergy(Complex(rname, lname));

	printf("#E(%s)= %.5e\n", RRCE_.GetName().c_str(), E);

	/*
	 * (2) spectral decomposition
	 */
	double **J = (double**) malloc(20 * sizeof(double*));
	for (int i = 0; i < 20; i++) {
		J[i] = (double*) malloc(20 * sizeof(double));
	}
	RRCE_.GetCouplings(J);
	//EigenRRCE Eigen(J);

	/*
	for (int i = 1; i <= 20; i++) {
		printf("# %d %f %f %f\n", i, Eigen.GetEigenvalue(i),
				1.0 - Eigen.GetReconstructionError(i),
				Eigen.GetReconstructionCorrel(i));
	}
	*/
	//Eigen.GetReconstructionCorrel(8);

	return 0;

}

bool GetArguments(char *rname, char *lname, RRCE::RRCE_TYPE &type, double &dmax,
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
			if (strcmp(tname, "RRCE20RC") == 0) {
				type = RRCE::RRCE20RC;
			} else if (strcmp(tname, "RRCE20SCC") == 0) {
				type = RRCE::RRCE20SCC;
			} else if (strcmp(tname, "RRCE20CB") == 0) {
				type = RRCE::RRCE20CB;
			} else {
				printf("Error: wrong RRCE type '%s'\n", arg);
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
