#pragma once
#include <vector>

//this file contains the important thermodynamical functions used in the simulation as methods
class THDQuantities
{

public:
	 //0D limit function
	 static double phizeroD(double x);
	 //derivation of 0D limit funciton
	 static double phiZeroDPrime(double x);

	 /* calculates the canonical potential of a homogene density profile
	 note that wall_pos_1 and 2 are actually the positions after and before the barriers respectively (which have nonzero density) because of the fluid domain */
	 static double canonialPotHomogene(double rho_0, int wall_pos_1, int wall_pos_2, int rod_length);

	 //calculates the canonical potential of a density profile with a given external potential "ext_pot" and a chemical potential "mu"
	 static double canonicalPot(std::vector<double> rho_s, std::vector<double> ext_pot, double mu, int rod_length, int wall_pos, int lattice_length);

	 //returns the grand potential density of each lattice point (small omega_s) given a density profile "rho_s"
	 static std::vector<double> canonicalPotDensity(std::vector<double> rho_s, std::vector<double> ext_pot, double mu, int rod_length, int wall_pos, int lattice_length);
	
	 //calculates the pressure of a homogene density distribution
	 static double surfacePressure(double rho_0, int rod_length);

	 //calculates the exakt Adsorbtion for a 1D lattice with Rods with determined density distribution rho_s a starting density rho_0 
	 //here wall_pos_1 is actually the correct wall position, so the input must be V - 1
	 static double exactAdsorbtionOneDim(std::vector<double> rho_s, double rho_0, int wall_pos_1, int wall_pos_2);
	 //calculates the weighted density eta^{(1)}(x)
	 static double weightedDensityOne(int L, int S, std::vector<double>& density_profile);
	 //calculates weighted density eta^{(0)}(x)
	 static double weightedDensityZero(int L, int S, std::vector<double>& density_profile);
	 //calculates the local excess chemical potential mu^{Ex}_s
	 static double muExS(int L, int S, std::vector<double>& density_profile);

	 


};





