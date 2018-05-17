/*****************************************************************************
 *
 * Ivan Anishchenko <anishchenko.ivan@gmail.com>
 *
 * 2013-2016, Vakser Lab, Center for Computational biology, University of Kansas
 *
 *****************************************************************************/

#include <cstring>

#include "RRCE.h"
#include "Options.h"

#define DMAX 7.8

int main(int argc, char *argv[]) {

	/*
	 * (1) process arguments
	 */
	OPTS opts = { "", "", "", DMAX };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts("rrce20", DMAX);
		return 1;
	}

	/*
	 * (2) set potential
	 */
	RRCE RRCE_;
	if (opts.table != "") {
		RRCE_ = RRCE(opts.table.c_str(), opts.dmax);
	}

	/*
	 * (2) get RRCE energy
	 */
	double E;
	if (opts.lig == "") {
		E = RRCE_.GetEnergy(Chain(opts.rec.c_str()));
	} else {
		E = RRCE_.GetEnergy(Complex(opts.rec.c_str(), opts.lig.c_str()));
	}

	printf("E(%s)= %.5e\n", RRCE_.GetName().c_str(), E);

	return 0;

}
