#include "THDQuantities.h"
#include <cmath>
#include <iostream>

double THDQuantities::phizeroD(double x)
{
	if (x >= 1) { //avoiding log(0)
		x = 0.9999999;
	}
	return x + (1 - x) * log(1 - x);
}



double THDQuantities::phiZeroDPrime(double x)
{
	if (x >= 1) { //avoiding log(0)
		x = 0.9999999;
	}
	return -1*log(1-x);
}


double THDQuantities::canonialPotHomogene(double rho_0, int wall_pos_1, int wall_pos_2, int rod_length)
{
	return -1 * surfacePressure(rho_0, rod_length) * (wall_pos_2 - wall_pos_1 + 1);
}

double THDQuantities::canonicalPot(std::vector<double> rho_s, std::vector<double> ext_pot, double mu, int rod_length,int wall_pos, int lattice_length){
	int M = lattice_length;
	int V = wall_pos;
	int L = rod_length;
	double sum = 0;
	for (int i = V; i < M - V - L + 3 ; i++)//provides values from V to M - V - L - 1 (we only take the points with nonzero density)
	{
		
		double rho = rho_s[i]; //assign current rho
		if(rho == 0){ //avoid log(0) just for safety
			sum += phizeroD(weightedDensityOne(L, i, rho_s)) - phizeroD(weightedDensityZero(L, i, rho_s));

		}
		else {
			double a = rho * (log(rho) - 1)
				+ phizeroD(weightedDensityOne(L, i, rho_s))
				- phizeroD(weightedDensityZero(L, i, rho_s))
				- (mu) * rho;
			
			sum += a;
		}

	}


	return sum;

}


std::vector<double> THDQuantities::canonicalPotDensity(std::vector<double> rho_s, std::vector<double> ext_pot, double mu, int rod_length, int wall_pos, int lattice_length) {
	int M = lattice_length;
	int V = wall_pos;
	int L = rod_length;
	
	std::vector<double> omega_s = {};

	
	for (int i = V; i < M - V - L + 3; i++) //provides values from V to M - V - L - 1 
	{
		double rho = rho_s[i]; //assign current rho
		if (rho == 0) { //avoid log(0) just for safety
			omega_s.push_back(phizeroD(weightedDensityOne(L, i, rho_s)) - phizeroD(weightedDensityZero(L, i, rho_s)));

		}
		else {
			double a = rho * (log(rho) - 1)
				+ phizeroD(weightedDensityOne(L, i, rho_s))
				- phizeroD(weightedDensityZero(L, i, rho_s))
				- (mu)*rho;

			omega_s.push_back(a);
		}

	}


	return omega_s;

}


double THDQuantities::surfacePressure(double rho_0 , int rod_length) {
	int L = rod_length;
	return log((1 - ((L - 1) * rho_0)) / (1 - L * rho_0));
}

double THDQuantities::exactAdsorbtionOneDim(std::vector<double> rho_s, double rho_0,  int wall_pos_1, int wall_pos_2) {
	double sum = 0;
	for (int i = wall_pos_1 + 1; i < wall_pos_2; i++) 
	{
		sum += rho_s[i] - rho_0;
	}
	return sum;
}



double THDQuantities::weightedDensityOne(int L, int S, std::vector<double> &density_profile)
{
	double rho_s = 0; //initialize sum
	for (int i = S - (L - 1); i < S + 1; i++)//sum over densities of L preceding lattice points
	{
		rho_s += density_profile[i];
	}
	return rho_s;
}

double THDQuantities::weightedDensityZero(int L ,int S, std::vector<double> &density_profile)
{
	double rho_s = 0; //initialize sum
	for (int i = S - (L - 1); i < S; i++)//sum over densities of L-1 preceding lattice points
	{
		rho_s += density_profile[i];
	}

	return rho_s;
}




double THDQuantities::muExS(int L, int S, std::vector<double> &density_profile)
{
	double sum1 = 0; //initialize sum
	for (int i = S; i < S + L; i++)
	{
		sum1 += phiZeroDPrime(weightedDensityOne(L, i, density_profile));
	}

	double sum2 = 0;
	for (int i = S + 1; i < S + L; i++)
	{
		sum2 += phiZeroDPrime(weightedDensityZero(L, i, density_profile));
	}

	return sum1 - sum2; //return the Differenz of both the sums
}



