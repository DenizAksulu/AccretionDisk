#include "stdafx.h"
#include <fstream>
#include <iostream>
#include <amp_math.h>
#include <algorithm>
#include <amp.h>
#include <ppl.h>
#include <chrono>

using namespace ::std;
using namespace ::Concurrency;
using namespace ::Concurrency::fast_math;

double NormalizationGaussian(double Mu, double Sigma, double R_isco, double R_outer, double M_disk);
void WriteGraphData(double* X, double* Y, int length, string filename, bool append);
void WriteGraphData(array_view<double, 1> X, array_view<double, 1> Y, int length, string filename, bool append);
void Inverse(array_view<double, 2> A, int length);

int main()
{
	// CONSTANTS*******************************************************************
	const double day = 86400e0;                                          // (day) s
	const double month = 2629743.83e0;                                   // (month) s
	const double year = 31556926e0;                                      // (year) s
	const double M_solar = 1.9891e33;                                    // (solar mass) g
	const double AU = 1.496e13;                                          // (astronomical unit) cm
	const double ly = 9.463e17;                                          // (light year) cm
	const double pc = 3.086e18;                                          // (parsec)cm
	const double c = 2.99792458e10;                                      // (speed of light) cm / s
	const double h = 6.6260755e-27;                                      // (planck's constant) erg s
	const double hbar = 1.05457266e-27;                                  // (reduced planck's constant) erg s
	const double G = 6.67259e-8;                                         // (gravitational constant) cm3 g-1 s-2
	const double e = 4.8032068e-10;                                      // (electron charge) esu
	const double m_e = 9.1093897e-28;                                    // (mass of electron) g
	const double m_p = 1.6726231e-24;                                    // (mass of proton) g
	const double m_n = 1.6749286e-24;                                    // (mass of neutron) g
	const double m_H = 1.6733e-24;                                       // (mass of hydrogen) g
	const double k = 1.380658e-16;                                       // (boltzmann's constant) erg K-1
	const double eV = 1.6021772e-12;                                     // (electron volt) erg
	const double a = 5.67051e-5;                                         // (stefan boltzmann constant) erg cm-2 K-4 s-1
	const double thompson = 6.6524e-25;                                  // (thompson cross-section for electron) cm-2
	const double PI = 3.14159265358979323846;
	//******************************************************************************************************************

	// VARIABLES
	double M_compact;								// Mass of the compact object
	double M_disk;									// Mass of the disk
	double M_dot_boundary;									// Accretion rate
	double R_isco, R_outer, X_isco, X_outer;		// Inner and outer radius
	int N_grids;									// Number of grids
	int N_sample;									// Number of samples
	double delta_X;									// Step size
	double alpha;									// alpha parameter for disk
	double mu_p = 0.6;								// proton ratio
	double T_max;
	double T = 0;
	double L_instant = 0;

	cout << "Accretion disk simulation with C++ AMP\n\n";
	accelerator acc;

	wcout << "Default accelerator is: " << acc.description << "\n";
	wcout << "Accelerator memory is: " << acc.get_dedicated_memory() << " kB\n";
	wcout << "Supports double precision operations: " << acc.get_supports_double_precision() << "\n\n";

	cout << "Please enter the mass of the compact object. (M_solar)\n";
	double n; cin >> n;
	M_compact = n * M_solar;
	cout << "Please enter the mass of the disk. (M_solar)\n";
	cin >> n;
	M_disk = n * M_solar;
	cout << "Please enter the accretion rate. (M_solar s-1)\n";
	cin >> n;
	M_dot_boundary = n * M_solar;
	cout << "Please enter the alpha parameter.\n";
	cin >> alpha;

	// Calculate some initial values
	R_isco = 6 * G * M_compact / pow(c, 2);
	R_outer = 1000 * R_isco;
	X_isco = sqrt(R_isco);
	X_outer = sqrt(R_outer);

	cout << "Mass of the compact object = " << M_compact << " g.\n";
	cout << "Mass of the disk           = " << M_disk << " g.\n";
	cout << "Innermost stable orbit     = " << R_isco << " cm.\n";
	cout << "Outer radius               = " << R_outer << " cm.\n";
	cout << "Accretion rate             = " << M_dot_boundary << " g s-1.\n";
	cout << "Please enter the number of grids for the radial coordinate.\n";
	cin >> N_grids;
	cout << "Creating initial conditions...\n";
	// Calculate steps in X space
	delta_X = (X_outer - X_isco) / N_grids;

	// Create vectors
	double* vE = new double[N_grids];			// Surface mass density vector
	double*  vE_new = new double[N_grids];		// Next surface mass density vector
	double*  vR = new double[N_grids];			// Radius
	double*  vS = new double[N_grids];
	double*  vS_new = new double[N_grids];
	double*  vX = new double[N_grids];
	double*  vM_dot = new double[N_grids];
	double*  vV = new double[N_grids];			// Viscosity
	double*  vV_new = new double[N_grids];
	double* vT_eff = new double[N_grids];			// Surface mass density vector
	double* vdelta_T = new double[N_grids];			// Surface mass density vector
	double dT;
	double vJ_total = 0;
	vM_dot[N_grids - 2] = M_dot_boundary;

	array_view<double, 1> R(N_grids, vR);
	array_view<double, 1> X(N_grids, vX);
	array_view<double, 1> E(N_grids, vE);
	array_view<double, 1> S(N_grids, vS);
	array_view<double, 1> V(N_grids, vV);
	array_view<double, 1> M_dot(N_grids, vM_dot);
	array_view<double, 1> E_new(N_grids, vE_new);
	array_view<double, 1> S_new(N_grids, vS_new);
	array_view<double, 1> V_new(N_grids, vV_new);
	array_view<double, 1> T_eff(N_grids, vT_eff);
	array_view<double, 1> deltaT(N_grids, vdelta_T);
	array_view<double, 1> J_total(1, &vJ_total);

	
	// Create initial conditions******************************************************************************************************************************
	//********************************************************************************************************************************************************
	auto start = chrono::high_resolution_clock::now();

	double VIS_C = pow(
		(27. / 32.) * (pow(alpha, 4) * pow(k, 4) * thompson)
		/ (pow(mu_p, 4) * pow(m_p, 5) * a * G * M_compact)
		, (1. / 3.));
	double Mu = (R_outer - R_isco) / 2;
	double Sigma = (R_outer - R_isco) / 10;

	double E_0 = NormalizationGaussian(Mu, Sigma, R_isco, R_outer, M_disk);

	parallel_for_each(R.extent, [=](index<1> idx) restrict(amp)
	{
		X(idx[0]) = X_isco + idx[0] * delta_X;
		M_dot(idx[0]) = 0;
	});
	X.synchronize();
	parallel_for_each(R.extent, [=](index<1> idx) restrict(amp) 
	{  
		R(idx) = X(idx)*X(idx);
		E(idx) = E_0 * exp(-1 * (pow((R(idx) - Mu), 2.)) / (2. * pow(Sigma, 2.)));
		S(idx) = X(idx)*E(idx);
		V(idx) = VIS_C * pow(X(idx), (4. / 3.)) * pow(S(idx), (2. / 3.));
		J_total(0) += 4 * PI * sqrt(G * M_compact) * X(idx) * S(idx) * delta_X;
	});
	Concurrency::extent<1> ext(N_grids - 2);
	parallel_for_each(ext, [=](index<1> idx) restrict(amp)
	{
		M_dot(idx[0]) = 3. * PI * (V(idx[0] + 1) * S(idx[0] + 1) - V(idx[0]) * S(idx[0])) / delta_X;
		if (M_dot(idx[0]) >= 0)
			T_eff(idx[0]) = pow(3 * G * M_compact * M_dot(idx[0]) / (8 * PI * a * pow(X(idx[0]), 6)) * (1 - sqrt(X_isco / pow(X(idx[0]), 2))), 0.25);
		else
			T_eff(idx[0]) = 0;
	});
	WriteGraphData(R, E, N_grids, "EvsR.txt", false);
	WriteGraphData(R, S, N_grids, "SvsR.txt", false);
	WriteGraphData(R, V, N_grids, "VvsR.txt", false);
	WriteGraphData(R, M_dot, N_grids-1, "MdotvsR.txt", false);
	WriteGraphData(R, T_eff, N_grids-1, "TeffvsR.txt", false);

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	cout << "Initial conditions created in " <<	elapsed.count() << " ms.\n";
	cout << "Total angular momentum is " << J_total(0) << " g cm2 s-1.\n";
	//********************************************************************
	//********************************************************************

	cout << "Please enter the evolution duration. (months)\n";
	cin >> n; T_max = n*month;

	cout << "How many graphs are needed?\n";
	cin >> N_sample;
	double* vT_sample = new double[N_sample];
	double deltaT_sample = log(T_max) / (double)N_sample;
	for (int i = 1; i <= N_sample; i++)
	{
		vT_sample[i - 1] = exp(i*deltaT_sample);
	}

	int j = 0;
	start = chrono::high_resolution_clock::now();

	FILE* file;
	file = fopen("lightcurve.txt", "w");
	while (T < T_max)
	{
		// Determine outer boundary condition*************************************************************
		//************************************************************************************************
		S[N_grids - 1] = pow(((M_dot[N_grids - 2] * delta_X / (3. * PI) + S[N_grids - 2] * V[N_grids - 2]) 
			/ (VIS_C * pow(X[N_grids - 1], (4. / 3.)))), 
			(3. / 5.));
		//************************************************************************************************
		//************************************************************************************************

		// Determine time step****************************************************************************
		//************************************************************************************************

		parallel_for_each(deltaT.extent, [=](index<1> idx) restrict(amp)
		{
			deltaT(idx) = (1. / 6.) * (pow(delta_X * X(idx), 2) / V(idx)) / 0.75;
		});
		deltaT.synchronize();
		dT = *min_element(vdelta_T, vdelta_T + N_grids);
		
		//************************************************************************************************
		//************************************************************************************************

		// Iterate in radial coordinate ******************************************************************
		//************************************************************************************************
		Concurrency::extent<1> ext2(N_grids - 2);
		parallel_for_each(ext2, [=](index<1> idx) restrict(amp)
		{
			
			S_new(idx[0] + 1) = S(idx[0] + 1) + 0.75 * dT / pow((X(idx[0] + 1) * delta_X), 2) *
				(V(idx[0] + 2) * S(idx[0] + 2) + V(idx[0]) * S(idx[0]) - 2. * V(idx[0] + 1) * S(idx[0] + 1));
			V_new(idx[0] + 1) = VIS_C * pow(X(idx[0] + 1), (4./3.)) * pow(S_new(idx[0] + 1), (2. / 3.));
		});
		//*************************************************************************************************
		//*************************************************************************************************

		// Obtain new values*******************************************************************************
		//*************************************************************************************************
		parallel_for_each(ext2, [=](index<1> idx) restrict(amp)
		{
			S(idx[0] + 1) = S_new(idx[0] + 1);
			V(idx[0] + 1) = V_new(idx[0] + 1);
		});
		//*************************************************************************************************
		//*************************************************************************************************

		// Obtain new M_dot values ************************************************************************
		//*************************************************************************************************
		Concurrency::extent<1> ext3(N_grids - 2);
		parallel_for_each(ext3, [=](index<1> idx) restrict(amp)
		{
			M_dot(idx[0]) = 3. * PI * (V(idx[0] + 1) * S(idx[0] + 1) - V(idx[0]) * S(idx[0])) / delta_X;
			if (M_dot(idx[0]) >= 0)
				T_eff(idx[0]) = pow(3 * G * M_compact * M_dot(idx[0]) / (8 * PI*a*pow(X(idx[0]), 6)) * (1 - sqrt(X_isco / pow(X(idx[0]), 2))), 0.25);
			else
				T_eff(idx[0]) = 0;
		});
		//*************************************************************************************************
		//*************************************************************************************************


		T += dT; // Increase time

		L_instant = (M_dot(0) * G * M_compact) / (2 * R_isco); // Luminosity in ergs/s
		fprintf(file, "%lf\t%lf\n", T, L_instant);

		// Take samples ***********************************************************************************
		//*************************************************************************************************
		if (T >= vT_sample[j])
		{
			j++;
			parallel_for_each(R.extent, [=](index<1> idx) restrict(amp)
			{
				E(idx) = S(idx) / X(idx);
			});
			WriteGraphData(R, E, N_grids, "EvsR.txt", true);
			WriteGraphData(R, V, N_grids, "VvsR.txt", true);
			WriteGraphData(R, M_dot, N_grids- 1, "MdotvsR.txt", true);
			WriteGraphData(R, T_eff, N_grids - 1, "TeffvsR.txt", true);
			end = std::chrono::high_resolution_clock::now();
			elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			cout << "Current time step is " << dT << " s.\n";
			cout << T/T_max * 100 << fixed << " percent completed! " << elapsed.count() << " ms have elapsed.\n";
		}
		//*************************************************************************************************
		//*************************************************************************************************
	}
	fclose(file);
	cin.get();
	return 0;
}

double NormalizationGaussian(double Mu, double Sigma, double R_isco, double R_outer, double M_disk)
{
	const double PI = 3.14159265358979323846;
	double IntegralResult = 0.5 * Sigma * (2 * (exp(-0.5 * pow((R_isco - Mu), 2) / pow(Sigma, 2)) - exp(-0.5 *pow((R_outer - Mu), 2) / pow(Sigma, 2))) * Sigma
		- Mu * sqrt(2 * PI) * erf((R_isco - Mu) / (sqrt(2) * Sigma)) + Mu * sqrt(2 * PI) * erf((R_outer - Mu) / (sqrt(2) * Sigma)));
	double E_0 = M_disk / (2 * PI * IntegralResult);
	return E_0;
}

void WriteGraphData(array_view<double, 1> X, array_view<double, 1> Y, int length, string filename, bool append)
{
	FILE* file;
	if(!append)
		file = fopen(filename.c_str(), "w");
	else
		file = fopen(filename.c_str(), "a");
	parallel_for(0, length, [=, &file](int i)
	{
		fprintf(file, "%lf\t%lf\n", X[i], Y[i]);
	});
	fprintf(file, "\n");
	fclose(file);
}

void WriteGraphData(double* X, double* Y, int length, string filename, bool append)
{
	FILE* file;
	if (!append)
		file = fopen(filename.c_str(), "w");
	else
		file = fopen(filename.c_str(), "a");
	parallel_for(0, length, [=, &file](int i)
	{
		fprintf(file, "%lf\t%lf\n", X[i], Y[i]);
	});
	fprintf(file, "\n");
	fclose(file);
}
void Inverse(array_view<double, 2> A, int length)
{
	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < length; j++)
		{
			cout << A(i, j) << "\t";
		}
		cout << "\n";
	}

	for (int i = 1; i < length; i++)
	{
		for(int j = 0; j < i; j++)
		{
			parallel_for_each(A.extent, [=](index<2> idx) restrict(amp)
			{
				if (A(j, j) != 0)
					A(i, idx[1]) -= (A(i, j) / A(j, j)) * A(j, idx[1]);

				if (A(i, idx[1]) < 1e-9 && A(i, idx[1]) > -1e-9)
					A(i, idx[1]) = 0;
			});
		}
	}
	/*parallel_for_each(A.extent, [=](index<2> idx) restrict(amp)
	{
		A(0, idx[1]) = A(0, idx[1]) / A(0, 0);
		
		for (int i = 0; i < idx[0]; i++)
		{
			if (A(i, i) != 0)
				A(idx[0], idx[1]) -= A(i, idx[1])*A(idx[0], i) / A(i, i);
		}
	});*/

	for (int i = 0; i < length; i++)
	{
		for (int j = 0; j < length; j++)
		{
			cout << A(i, j) << "\t";
		}
		cout << "\n";
	}
}
