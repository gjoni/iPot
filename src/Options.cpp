
#include <unistd.h>

#include "Options.h"

void PrintOpts(const std::string &exec, const double DMAX) {

	printf("\nUsage:   ./%s [-option] [argument]\n\n", exec.c_str());
	printf("Options:  -r receptor.pdb \\  # PDB file with receptor's coordinates\n");
	printf("          -l ligand.pdb \\    # (optional) PDB file with ligand's coordinates\n");
	printf("          -t table.txt \\     # (optional) %s contact potential table \n", exec.c_str());
	printf("          -d dmax            # (optional) contact distance, default dmax=%.1fA\n", DMAX);
	printf("\n");

}

bool GetOpts(int argc, char *argv[], OPTS &opts) {

	char tmp;
	while ((tmp = getopt(argc, argv, "hs:r:l:t:d:")) != -1) {
		switch (tmp) {
		case 'h': /* help */
			printf("!!! HELP !!!\n");
			return false;
			break;
		case 'r': /* receptor file */
			opts.rec = std::string(optarg);
			break;
		case 'l': /* ligand file */
			opts.lig = std::string(optarg);
			break;
		case 't': /* table file */
			opts.table = std::string(optarg);
			break;
		case 'd': /* distance cutoff */
			opts.dmax = atof(optarg);
			break;
		default:
			return false;
			break;
		}
	}

	if (opts.rec == "") {
		printf("Error: receptor file not specified ('-r')\n");
		return false;
	}

	return true;

}

