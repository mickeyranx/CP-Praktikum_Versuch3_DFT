
#include "THDQuantities.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <tuple>



//-----------------------------------
//          global params
//-----------------------------------
double const STARTING_ALPHA = 0.01; //starting alpha for picard mixing
double const ALPHA_MIN = 0.0001;
double const ALPHA_MAX = 0.03;

int const V = 10; //technically the wall is at 9 because the profile has a 0-th value so in the code you may find a lot of the times V - 1
//call count of startSimulation() to manage data writing
static int callcount = 0;
static int callcount2 = 0;


//-----------------------------------
//          methods
//-----------------------------------

//return the next picard alpha for the next iteration
static double nextAlpha(double epsilon_old, double epsilon_new, double alpha) {
	if (epsilon_old == -1) { //do nothing if this is the first iteration (alpha stays)	
		return alpha;
	}
	else {
		//if new epsilon is smaller than old multiply by 1.1 otherwise divide by 5
		(epsilon_new < epsilon_old) ? (alpha = 1.1 * alpha) : (alpha = alpha / 5);
		
		if (alpha > ALPHA_MAX) {
			return 0.03;
		}
		else if (alpha < ALPHA_MIN) {
			return 0.0001;
		
		}
		else {
			return alpha;
		}
	}
}

//this function calculates the exact Adsorption Gamma_Ex and the an estimated Value from finite-differences and returns them as a tuple
static std::tuple<double, double> checkAdsorbtionCondition(std::vector<double> rho_s, double delta_rho, double rho_s_0, int M, int L)
{
	
	//calculate and print the exat Adsorbtion
	double G_ex = THDQuantities::exactAdsorbtionOneDim(rho_s, rho_s_0, V - 1, M - V - L) / 2;
	
	//estimate the Adsorption with a variation of rho_0 and integrating with Finite-Differences. p is "plus", m is "minus"
	double mu_id_p = log(rho_s_0 + delta_rho);
	double mu_ex_p = -L * log(1 - (rho_s_0 + delta_rho) * L) + (L - 1) * log(1 - (rho_s_0 + delta_rho) * (L - 1));
	double mu_id_m = log(rho_s_0 - delta_rho);
	double mu_ex_m = -L * log(1 - (rho_s_0 - delta_rho) * L) + (L - 1) * log(1 - (rho_s_0 - delta_rho) * (L - 1));
	double mu_p = mu_id_p + mu_ex_p;
	double mu_m = mu_id_m + mu_ex_m;
	double gamma_p = (mu_ex_p - (2 * L - 1) * THDQuantities::surfacePressure(rho_s_0 + delta_rho, L)) / 2;
	double gamma_m = (mu_ex_m - (2 * L - 1) * THDQuantities::surfacePressure(rho_s_0 - delta_rho, L)) / 2;
	double G_est = -1 * (gamma_p - gamma_m) / (mu_p - mu_m);
	return std::make_tuple(G_ex, G_est);
}

//returns the density profile and epsilon_new (needed to determine next alpha) of the next iteration
static std::tuple<std::vector<double>, double> iterate(
	std::vector<double>& wall,
	std::vector<double>& rho_prev,
	double alpha, double rho_s_0, double mu_ex, double epsilon_old ,int M, int L) {
	//---------------------------------
	//  calculate next density profile
	//---------------------------------
	std::vector<double> rho_new(M, 0);
	for (int i = L; i < M - L; i++)
	{
		rho_new[i] = rho_s_0 * exp(mu_ex - THDQuantities::muExS(L, i, rho_prev)) * wall[i];
	}

	//---------------------------------
	//         picard mixing
	//---------------------------------
	std::vector<double> rho_next(M, 0);
	double epsilon_new = 0;
	for (int i = 0; i < M; i++)
	
	{
		//prefix c = current
		double c_rho_new = rho_new[i];

		double c_rho = rho_prev[i];

		rho_next[i] = (1 - alpha) * c_rho + (alpha)*c_rho_new;

		epsilon_new += pow((c_rho_new - c_rho), 2) * 1 / M;


	}
	
	return std::make_tuple(rho_next, epsilon_new);

}


//this function initializes important lists, the first density profile and values and starts to iterate
static void startSimulation(int M ,int L , double eta_0, int max_iter, std::string filename1, std::string filename2) {

	//------------------------------------
	//           setup
	//------------------------------------
	//create wall profile (0 outside walls, 1 inside walls)
	std::vector<double> wall(M, 0);
	for (int i = 0; i < M; i++)
	{
		if ((i > V - 1) && (i < (M - V) - L)) { //if lattice point within walls 
			wall[i] = 1;
		}
		else {
			wall[i] = 0;
		}
	}

	//calculate starting values for starting density profile
	double rho_s_0 = eta_0 / L; 
	double mu_ex = -L * log(1 - rho_s_0 * L) + (L - 1) * log(1 - rho_s_0 * (L - 1));

	//calculate first (0th) density_profile
	std::vector<double> rho_0(M);
	for (int i = 0; i < M; i++)
	{
		rho_0[i] = rho_s_0 * wall[i];
	}

	//initialize values for first iteration
	double alpha = STARTING_ALPHA;
	double epsilon_old = -1; //set first epsilon to -1 to detect first iteration
	double epsilon_new = 0;
	std::vector<double> rho_prev = rho_0;
	std::vector<double> rho_next = rho_0;
	//------------------------------------
	//           iterations
	//------------------------------------
	for (int i = 0; i < max_iter; i++)
	{
		if (epsilon_old <= pow(10, -9) && epsilon_old > 0) break; //stop iterating if epsilon is very small

		//calculate next density profile and epsilon new
		std::tie(rho_next, epsilon_new) = iterate(wall, rho_prev, alpha, rho_s_0, mu_ex, epsilon_old, M, L); //unpack tuple
		//prepare for the next iteration
		rho_prev = rho_next;
		alpha = nextAlpha(epsilon_old, epsilon_new, alpha);
		epsilon_old = epsilon_new;
		
	}
	
	//------------------------------------
	//        write data to file
	//------------------------------------
	
	//------------------------------------
	//        1. density profile
	//------------------------------------


	if (callcount == 0) { //first function call
		std::ofstream file(filename1); //create new file
		file << std::fixed << std::setprecision(8);
		for (double density : rho_next) //put data into textfile
		{
			file << density << "\n";
		}
	}else{

		std::ifstream inFile(filename1);
		std::vector<std::string> lines(M, "");

		std::string line; 
		int row = 0;

		while (std::getline(inFile, line) && row < M) {
			lines[row] = line; 
			row++;
		}
		inFile.close();
		for (int i = 0; i < M; i++) //put data into textfile
		{
			std::stringstream ss; 
			ss << std::fixed << std::setprecision(8) << (rho_next[i]);
			if (!lines[i].empty()) {
				lines[i] += "\t";
			}
			lines[i] += ss.str();
		}

		std::ofstream outFile(filename1);
		for (const auto& l : lines)
		{
			outFile << l << std::endl;
		}
		outFile.close();
	}
	callcount++;

	//------------------------------------
	//    2. other quantities
	//------------------------------------
	if (callcount2 == 0) {
		std::ofstream file2(filename2);
		file2.close();
	}
	
	std::ofstream file2(filename2, std::ios::app); //open textfile in append mode
	file2 << "----------------------------------" << "\n";
	file2 << "     eta_0 = " << eta_0 << "          " << "\n";
	file2 << "----------------------------------" << "\n";
	file2 << std::fixed << std::setprecision(10);
	
	double mu = log(rho_s_0) + mu_ex;
	file2 << "mu:             " << mu << "\n";
	file2 << "pressure:       " << THDQuantities::surfacePressure(rho_s_0, L) << "\n";
	//canonical potentials
	double omega_rho = THDQuantities::canonicalPot(rho_prev, wall, mu, L, V, M);
	double omega_rho_0 = THDQuantities::canonialPotHomogene(rho_s_0, V, M - V - L - 1, L);
	file2 << "omega_rho:      " << omega_rho << "\n";
	file2 << "omega_rho_0:    " << omega_rho_0 << "\n";
	//surface pressure
	file2 << "gamma_ana:      " << (mu_ex - (2 * L - 1) * THDQuantities::surfacePressure(rho_s_0, L)) / 2 << "\n";
	file2 << "gamma:          " << (omega_rho - omega_rho_0) / 2 << "\n";
	//Adsorbtion
	double G_ex;
	double G_est;
	double delta_rho = 0.001; //deviation value of rho for finite-differences
	std::tie(G_ex, G_est) = checkAdsorbtionCondition(rho_prev, delta_rho, rho_s_0, M, L); //calculate both values
	file2 << "exakt Gamma:     " << G_ex << "\n";
	file2 << "estimated Gamma: " << G_est << "\n";
	callcount2++;

	//this part is only to check if the canonical potential density is not symmetric around the slit midpoint
	/*
	if (eta_0 > 0.8) {
		std::vector<double> omega_s = THDQuantities::canonicalPotDensity(rho_prev, wall, mu, L, V - 1, M);
		for (double omega : omega_s)
		{
			std::cout << omega << ", ";
		}
		
	}
	*/
}






int main(int argc, char* argv[])
{

	if (argc != 5) { //if the given no. of params are not correct provide a error message
		std::cerr << "Usage:\n";
		std::cerr << "Usage: sim.exe <rod_length> <number_of_lattice_points> <wall_position> <manual_max_steps>  \n";
		return 1;
	}
	//-----------------------------------------------------
	//  parse parameters for simulation from a batch file
	//-----------------------------------------------------
	int L = std::atoi(argv[1]); //rodlength
	int M = std::atoi(argv[2]); //number of lattice points
	//int V = 10; 
	int V = std::atoi(argv[3]); //wall position (2nd wall (M-V))
	
	//int manual_max_steps = 5000;
	int manual_max_steps = std::atoi(argv[4]); //max no. iterations
	

	//the filename is generated form the input params

	//--------------------------------
	//  simulate with given parameters
	//    from eta_0 = 0.1 to 0.9
	//--------------------------------
	clock_t start = clock();
	//initial density
	double eta_0 = 0.1;
	for (int i = 0; i < 9; i++) //simulate for bulk packing fractions between 0.1 and 0.9 with stepsize of 0.1 
	{
		//data is put into two seperate files being ..._rho containing the density profile and ..._gamma containing the measured quantities on these profiles
		startSimulation(M, L, eta_0, manual_max_steps, "results/L" + std::to_string(L) + "_rho.txt", "results/L" + std::to_string(L) + "_gamma.txt");
		eta_0 += 0.1;
	}

	
	clock_t end = clock();
	double elapsed = double(end - start) / CLOCKS_PER_SEC;
	printf("exectution time: %.3f sec \n", elapsed);



}


