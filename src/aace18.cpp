/*****************************************************************************
 *
 * Ivan Anishchenko <anishchenko.ivan@gmail.com>
 *
 * 2013-2016, Vakser Lab, Center for Computational biology, University of Kansas
 *
 *****************************************************************************/

#include <cstring>

#include "AACE.h"
#include "Options.h"

#define DMAX 6.8

int main(int argc, char *argv[]) {

	/*
	 * (0) default energies (dmax=6.8A, kmin=5)
	 */
	const std::vector<double> h = {
		-1.60823e+00, -1.50102e+00, -1.47335e+00, -1.82202e+00, -1.75968e+00, -7.72545e-01, 6.04329e-01, 2.55978e-01, -4.19681e-01, -1.01144e+00, 4.42831e+00, 1.12513e+00, 7.39351e-01, 5.13158e-01, 7.03065e-01, -4.79007e-03, 1.23034e-01, 1.88040e+00
		};
	const std::vector<std::vector<double> > J = {
		{5.14494e-01, 2.32982e-01, -1.18825e-02, -4.79267e-01, -2.60971e-01, 7.92265e-03, 7.65882e-02, 3.00580e-02, 2.19054e-02, 5.83614e-02, -1.36975e-01, 1.96902e-01, 9.23782e-03, -2.08941e-02, 6.00891e-02, 1.23249e-01, -4.56176e-02, -8.81898e-03 }, 
		{2.32982e-01, -2.00039e-01, -6.56869e-02, 2.70771e-01, -2.92405e-01, -1.29688e-03, 2.66237e-02, 2.78053e-02, 4.34264e-02, 9.09723e-02, -3.18841e-03, -1.84093e-01, 1.75533e-02, 4.55108e-03, 4.44661e-03, 8.85454e-02, -5.36141e-03, -2.89657e-03 }, 
		{-1.18825e-02, -6.56869e-02, -3.58054e-03, 1.61754e-02, 1.01122e-01, -5.20250e-02, 1.47322e-02, -7.77877e-03, 2.93541e-02, 6.68389e-02, -3.44028e-02, -3.95784e-02, -9.05682e-03, -6.92520e-02, 8.67456e-02, 3.06997e-02, 2.63688e-02, -7.55891e-02 }, 
		{-4.79267e-01, 2.70771e-01, 1.61754e-02, -7.34460e-02, 4.22984e-01, -4.73875e-03, 1.43135e-02, -1.76282e-02, -7.93952e-03, 7.52781e-02, -1.01343e-01, 1.66995e-01, -1.17393e-02, -3.35176e-02, 4.17009e-02, 3.22659e-02, -2.99271e-02, -2.45088e-02 }, 
		{-2.60971e-01, -2.92405e-01, 1.01122e-01, 4.22984e-01, 1.36008e-01, -7.13044e-03, -7.48668e-03, -4.55776e-03, -7.20017e-03, 7.64353e-02, -1.09748e-01, -2.75419e-01, -8.12165e-03, -1.61199e-02, 6.78784e-02, 3.29076e-02, 1.26425e-02, -4.97552e-02 }, 
		{7.92265e-03, -1.29688e-03, -5.20250e-02, -4.73875e-03, -7.13044e-03, -1.05435e-01, -3.19825e-02, 1.00205e-03, 1.33486e-02, 2.53556e-02, -8.93222e-02, 2.00271e-02, -3.76764e-02, -1.29922e-01, 1.39732e-02, 4.12190e-02, -2.60661e-03, -9.08305e-02 }, 
		{7.65882e-02, 2.66237e-02, 1.47322e-02, 1.43135e-02, -7.48668e-03, -3.19825e-02, 6.23459e-02, 8.37308e-02, -4.59253e-03, -2.63784e-01, 5.51613e-02, -1.71107e-02, 2.27089e-03, 2.19795e-02, 2.54067e-01, 2.69225e-01, -1.72442e-02, -8.22190e-02 }, 
		{3.00580e-02, 2.78053e-02, -7.77877e-03, -1.76282e-02, -4.55776e-03, 1.00205e-03, 8.37308e-02, 1.56708e-02, -9.66613e-03, -2.53002e-01, 1.37543e-01, 7.54680e-03, -7.22425e-03, 7.30337e-02, 2.41662e-01, 1.82912e-01, -1.14706e-02, -7.38459e-02 }, 
		{2.19054e-02, 4.34264e-02, 2.93541e-02, -7.93952e-03, -7.20017e-03, 1.33486e-02, -4.59253e-03, -9.66613e-03, -8.89144e-02, -4.49119e-02, 1.13148e-01, -1.74796e-02, -2.83279e-02, 9.27447e-02, 7.43568e-03, -2.85799e-02, -5.86095e-02, -5.35649e-02 }, 
		{5.83614e-02, 9.09723e-02, 6.68389e-02, 7.52781e-02, 7.64353e-02, 2.53556e-02, -2.63784e-01, -2.53002e-01, -4.49119e-02, 1.21987e-01, 1.30581e-01, -5.94875e-03, -1.52840e-01, 1.10914e-01, -3.08656e-01, -3.40476e-01, -8.70755e-02, -9.98389e-02 }, 
		{-1.36975e-01, -3.18841e-03, -3.44028e-02, -1.01343e-01, -1.09748e-01, -8.93222e-02, 5.51613e-02, 1.37543e-01, 1.13148e-01, 1.30581e-01, -1.16021e+00, 1.89012e-01, -7.49936e-02, -8.24475e-02, 2.10707e-01, 1.61511e-01, 8.90165e-02, -5.08981e-02 }, 
		{1.96902e-01, -1.84093e-01, -3.95784e-02, 1.66995e-01, -2.75419e-01, 2.00271e-02, -1.71107e-02, 7.54680e-03, -1.74796e-02, -5.94875e-03, 1.89012e-01, -2.76854e-01, -1.09103e-02, 1.42503e-01, 1.03347e-01, -2.61694e-02, -1.89722e-02, -1.13523e-02 }, 
		{9.23782e-03, 1.75533e-02, -9.05682e-03, -1.17393e-02, -8.12165e-03, -3.76764e-02, 2.27089e-03, -7.22425e-03, -2.83279e-02, -1.52840e-01, -7.49936e-02, -1.09103e-02, -1.45362e-01, 1.86466e-02, 8.87609e-02, 6.51157e-02, -7.00287e-02, -1.05835e-01 }, 
		{-2.08941e-02, 4.55108e-03, -6.92520e-02, -3.35176e-02, -1.61199e-02, -1.29922e-01, 2.19795e-02, 7.30337e-02, 9.27447e-02, 1.10914e-01, -8.24475e-02, 1.42503e-01, 1.86466e-02, -2.57466e-01, 8.46647e-02, 1.38211e-01, 3.62556e-02, -1.05258e-01 }, 
		{6.00891e-02, 4.44661e-03, 8.67456e-02, 4.17009e-02, 6.78784e-02, 1.39732e-02, 2.54067e-01, 2.41662e-01, 7.43568e-03, -3.08656e-01, 2.10707e-01, 1.03347e-01, 8.87609e-02, 8.46647e-02, 3.68852e-01, 3.06581e-01, 3.14345e-02, -1.34007e-01 }, 
		{1.23249e-01, 8.85454e-02, 3.06997e-02, 3.22659e-02, 3.29076e-02, 4.12190e-02, 2.69225e-01, 1.82912e-01, -2.85799e-02, -3.40476e-01, 1.61511e-01, -2.61694e-02, 6.51157e-02, 1.38211e-01, 3.06581e-01, 1.85065e-01, -4.92545e-02, -1.02915e-01 }, 
		{-4.56176e-02, -5.36141e-03, 2.63688e-02, -2.99271e-02, 1.26425e-02, -2.60661e-03, -1.72442e-02, -1.14706e-02, -5.86095e-02, -8.70755e-02, 8.90165e-02, -1.89722e-02, -7.00287e-02, 3.62556e-02, 3.14345e-02, -4.92545e-02, -4.79142e-02, -2.83613e-02 }, 
		{-8.81898e-03, -2.89657e-03, -7.55891e-02, -2.45088e-02, -4.97552e-02, -9.08305e-02, -8.22190e-02, -7.38459e-02, -5.35649e-02, -9.98389e-02, -5.08981e-02, -1.13523e-02, -1.05835e-01, -1.05258e-01, -1.34007e-01, -1.02915e-01, -2.83613e-02, -1.01090e-01 }, 
		};

	/*
	 * (1) process arguments
	 */
	OPTS opts = { "", "", "", DMAX };
	if (!GetOpts(argc, argv, opts)) {
		PrintOpts("aace18", DMAX);
		return 1;
	}

	/*
	 * (2) set potential
	 */
	AACE AACE_;
	if (opts.table != "") {
		AACE_ = AACE(AACE18_GROUPS, opts.dmax, opts.table);
	} else {
		AACE_ = AACE(AACE18_GROUPS, h, J, opts.dmax, "AACE18_6.8A_k5");
	}

	/*
	 * (2) get AACE energy
	 */
	double E;
	if (opts.lig == "") {
		E = AACE_.GetEnergy(Chain(opts.rec.c_str()));
	} else {
		E = AACE_.GetEnergy(Complex(opts.rec.c_str(), opts.lig.c_str()));
	}

	printf("E(%s)= %.5e\n", AACE_.GetName().c_str(), E);

	return 0;

}
