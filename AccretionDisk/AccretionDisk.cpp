#include "stdafx.h"
#include <fstream>
#include <iostream>
//#include <amp_math.h>
#include <algorithm>
//#include <amp.h>
#include <ppl.h>
#include <chrono>

using namespace std;
using namespace ::Concurrency;
//using namespace ::Concurrency::fast_math;

double OpticalThickness(double SurfaceDensity);
double NormalizationGaussian(double Mu, double Sigma, double R_isco, double R_outer, double M_disk);
void WriteGraphData(double* X, double* Y, int length, string filename, bool append);
double eVtoHz(double eV);
void ExtractSpectrum(double* T, double* X, double delta_X, int N_grids, double minEnergyEV, double maxEnergyEV, double resolutionEV, bool append);
void IrradiationTemperature(double* T_irr, int N_grids, double nu, double epsilon, double L, double* R, double* H);
double* IrradiationTemperaturePower_4(int N_grids, double nu, double epsilon, double L, double* R, double* H);

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

	cout << "Accretion disk simulation with parallel CPU computing.\n\n";
	/*accelerator acc;

	wcout << "Default accelerator is: " << acc.description << "\n";
	wcout << "Accelerator memory is: " << acc.get_dedicated_memory() << " kB\n";
	wcout << "Supports double precision operations: " << acc.get_supports_double_precision() << "\n\n";*/

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
	double* vE_new = new double[N_grids];		// Next surface mass density vector
	double* vR = new double[N_grids];			// Radius
	double* vS = new double[N_grids];
	double* vS_new = new double[N_grids];
	double* vX = new double[N_grids];
	double* vM_dot = new double[N_grids];
	double* vV = new double[N_grids];			// Viscosity
	double* vV_new = new double[N_grids];
	double* vT_eff = new double[N_grids];			// Surface mass density vector
	double* vT_c = new double[N_grids];			// Surface mass density vector
	double* vT_irr = new double[N_grids];			// Surface mass density vector
	double* vT_sur = new double[N_grids];			// Surface mass density vector
	double* vdelta_T = new double[N_grids];			// Surface mass density vector
	double* vH = new double[N_grids];
	double dT;
	double vJ_total = 0;
	vM_dot[N_grids - 1] = M_dot_boundary;
	vM_dot[N_grids - 2] = M_dot_boundary;


	
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

	parallel_for(0, N_grids,[=](int i)
	{
		vX[i] = X_isco + i * delta_X;
		vM_dot[i] = 0;
	});

	parallel_for(0, N_grids, [=, &vJ_total](int i)
	{  
		vR[i] = vX[i]*vX[i];
		vE[i] = E_0 * exp(-1 * (pow((vR[i] - Mu), 2.)) / (2. * pow(Sigma, 2.)));
		vS[i] = vX[i]*vE[i];
		vV[i] = VIS_C * pow(vX[i], (4. / 3.)) * pow(vS[i], (2. / 3.));
		vJ_total += 4 * PI * sqrt(G * M_compact) * vX[i] * vS[i] * delta_X;
	});
	parallel_for(0, N_grids - 2, [=](int i)
	{
		vM_dot[i] = 3. * PI * (vV[i + 1] * vS[i + 1] - vV[i] * vS[i]) / delta_X;
		if (vM_dot[i] >= 0)
			vT_eff[i] = pow(3 * G * M_compact * vM_dot[i] / (8 * PI * a * pow(vX[i], 6)) * (1 - sqrt(X_isco / pow(vX[i], 2))), 0.25);
		else
			vT_eff[i] = 0;
		if (vM_dot[i] >= 0)
			vT_c[i] = pow(3 * OpticalThickness(vE[i]) * pow(vT_eff[i], 4) / 4, 0.25);
			
		else
			vT_c[i] = 0;
		vH[i] = sqrt((vT_c[i] * k * pow(vR[i], 3)) / (mu_p * m_p * G * M_compact));
	});

	L_instant = (vM_dot[0] * G * M_compact) / (2 * R_isco);		// Luminosity in ergs/s

	IrradiationTemperature(vT_irr, N_grids, 2, 0.5, L_instant, vR, vH);

	parallel_for(0, N_grids - 2, [=](int i)
	{
		vT_sur[i + 1] = pow(pow(vT_irr[i + 1], 4) + pow(vT_eff[i + 1], 4), 0.25);
	});

	/*
	* DANGER!!!!!
	*/
	parallel_for(0, N_grids - 2, [=](int i)
	{
		vT_c[i + 1] = pow(3 * OpticalThickness(vE[i + 1]) * pow(vT_eff[i + 1], 4) / 8 + (pow(vT_irr[i + 1], 4)), 0.25);
		vH[i + 1] = sqrt((vT_c[i + 1] * k * pow(vR[i + 1], 3)) / (mu_p * m_p * G * M_compact));
		vV[i + 1] = alpha * sqrt((vT_c[i + 1] * k) / (mu_p * m_p)) * vH[i + 1];
	});
	/*
	*
	*
	*/

	L_instant = (vM_dot[0] * G * M_compact) / (2 * R_isco);		// Luminosity in ergs/s

	ExtractSpectrum(vT_eff, vX, delta_X, N_grids - 2, 1, 15000, 1, false);

	WriteGraphData(vR, vT_irr, N_grids - 1, "TirrvsR.txt", false);
	WriteGraphData(vR, vE, N_grids, "EvsR.txt", false);
	WriteGraphData(vR, vV, N_grids, "VvsR.txt", false);
	WriteGraphData(vR, vM_dot, N_grids-1, "MdotvsR.txt", false);
	WriteGraphData(vR, vT_eff, N_grids-2, "TeffvsR.txt", false);
	WriteGraphData(vR, vT_c, N_grids - 2, "TcvsR.txt", false);
	WriteGraphData(vR, vH, N_grids - 2, "HvsR.txt", false);
	WriteGraphData(vR, vT_sur, N_grids - 2, "TsurvsR.txt", false);



	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	cout << "Initial conditions created in " <<	elapsed.count() << " ms.\n";
	cout << "Total angular momentum is " << vJ_total << " g cm2 s-1.\n";
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
	int l = 0;
	start = chrono::high_resolution_clock::now();

	ofstream file;
	file.open("lightcurve.txt", ios::out);
	if (L_instant > 0)
		file << T << "\t" << L_instant << "\n";						// Write luminosity to file
	file.close();
	while (T < T_max)
	{
		// Determine outer boundary condition*************************************************************
		//************************************************************************************************
		vS[N_grids - 1] = pow(((vM_dot[N_grids - 2] * delta_X / (3. * PI) + vS[N_grids - 2] * vV[N_grids - 2]) 
			/ (VIS_C * pow(vX[N_grids - 1], (4. / 3.)))), 
			(3. / 5.));
		//************************************************************************************************
		//************************************************************************************************

		// Determine time step****************************************************************************
		//************************************************************************************************

		parallel_for(0, N_grids, [=] (int i)
		{
			vdelta_T[i] = (1. / 6.) * (pow(delta_X * vX[i], 2) / vV[i]) / 0.75;
		});
		dT = *min_element(vdelta_T, vdelta_T + N_grids);
		
		//************************************************************************************************
		//************************************************************************************************

		// Iterate in radial coordinate ******************************************************************
		//************************************************************************************************
		parallel_for(0, N_grids - 2, [=](int i)
		{
			
			vS_new[i + 1] = vS[i + 1] + 0.75 * dT / pow((vX[i + 1] * delta_X), 2) *
				(vV[i + 2] * vS[i + 2] + vV[i] * vS[i] - 2. * vV[i + 1] * vS[i + 1]);
			vV_new[i + 1] = VIS_C * pow(vX[i + 1], (4./3.)) * pow(vS_new[i + 1], (2. / 3.));
		});
		//*************************************************************************************************
		//*************************************************************************************************

		// Obtain new values*******************************************************************************
		//*************************************************************************************************
		parallel_for(0, N_grids - 2, [=](int i)
		{
			vS[i + 1] = vS_new[i + 1];
			vV[i + 1] = vV_new[i + 1];
		});
		//*************************************************************************************************
		//*************************************************************************************************
		parallel_for(0, N_grids, [=](int i)
		{
			vE[i] = vS[i] / vX[i];
		});
		parallel_for(0, N_grids - 2, [=](int i)
		{
			vM_dot[i] = 3. * PI * (vV[i + 1] * vS[i + 1] - vV[i] * vS[i]) / delta_X;
			vT_eff[i + 1] = pow(3 * G * M_compact * abs(vM_dot[i + 1]) / (8 * PI * a * pow(vX[i + 1], 6)) * (1 - sqrt(X_isco / pow(vX[i + 1], 2))), 0.25);
			vT_c[i + 1] = pow(3 * OpticalThickness(vE[i + 1]) * pow(vT_eff[i + 1], 4) / 4, 0.25);
			vH[i + 1] = sqrt((vT_c[i + 1] * k * pow(vR[i + 1], 3)) / (mu_p * m_p * G * M_compact));
		});

		L_instant = (vM_dot[0] * G * M_compact) / (2 * R_isco);		// Luminosity in ergs/s

		IrradiationTemperature(vT_irr, N_grids, 2, 0.5, L_instant, vR, vH);

		parallel_for(0, N_grids - 2, [=](int i)
		{
			vT_sur[i + 1] = pow(pow(vT_irr[i + 1], 4) + pow(vT_eff[i + 1], 4), 0.25);
		});

		/*
		* DANGER!!!!!
		*/
		parallel_for(0, N_grids - 2, [=](int i)
		{
			vT_c[i + 1] = pow(3 * OpticalThickness(vE[i + 1]) * pow(vT_eff[i + 1], 4) / 8  + (pow(vT_irr[i + 1], 4)), 0.25);
			vH[i + 1] = sqrt((vT_c[i + 1] * k * pow(vR[i + 1], 3)) / (mu_p * m_p * G * M_compact));
			vV[i + 1] = alpha * sqrt((vT_c[i + 1] * k) / (mu_p * m_p)) * vH[i + 1];
		});
		/*
		*
		*
		*/

		T += dT; // Increase time

		if (T >= (l + 1) * 1800)
		{
			l++;
			file.open("lightcurve.txt", ios::app);
			if (L_instant > 0)
				file << T << "\t" << L_instant << "\n";						// Write luminosity to file
			file.close();
		}
		// Take samples ***********************************************************************************
		//*************************************************************************************************
		if (T >= vT_sample[j])
		{
			j++;
/*
			// Obtain E values for central temperature*********************************************************
			//*************************************************************************************************
			parallel_for(0, N_grids, [=](int i)
			{
				vE[i] = vS[i] / vX[i];
			});
			//*************************************************************************************************
			//*************************************************************************************************

			parallel_for(0, N_grids, [=] (int i)
			{
				// Obtain effective temperature ***********************************************************************************************
				//***************************************************************************************************************************
				if (vM_dot[i] >= 0)
					vT_eff[i] = pow(3 * G * M_compact * vM_dot[i] / (8 * PI * a * pow(vX[i], 6)) * (1 - sqrt(X_isco / pow(vX[i], 2))), 0.25);
				else
					vT_eff[i] = 0;
				//***************************************************************************************************************************
				//***************************************************************************************************************************

				// Obtain central temperature ***********************************************************************************************
				//***************************************************************************************************************************
				if (vM_dot[i] >= 0)
					vT_c[i] = pow(3 * OpticalThickness(vE[i]) * pow(vT_eff[i], 4) / 4, 0.25);
				//vT_c[i] = pow(9 * OpticalThickness(vS[i] / vX[i]) * G * M_compact * vM_dot[i] / (32 * PI * a * pow(vX[i], 6)) * (1 - sqrt(X_isco / pow(vX[i], 2))), 0.25);
				else
					vT_c[i] = 0;
				//***************************************************************************************************************************
				//***************************************************************************************************************************

				vH[i] = sqrt((vT_c[i] * k * pow(vR[i], 3)) / (mu_p * m_p * G * M_compact));
			});

			vT_irr = IrradiationTemperature(N_grids, 0.1, 0.5, L_instant, vR, vH);

			parallel_for(0, N_grids, [=](int i)
			{
				vT_sur[i] = pow(pow(vT_irr[i], 4) + pow(vT_eff[i], 4), 0.25);
			});
			*/
			ExtractSpectrum(vT_eff, vX, delta_X, N_grids - 2, 1, 15000, 1, true);

			WriteGraphData(vR, vE, N_grids, "EvsR.txt", true);
			WriteGraphData(vR, vV, N_grids, "VvsR.txt", true);
			WriteGraphData(vR, vM_dot, N_grids- 1, "MdotvsR.txt", true);
			WriteGraphData(vR, vT_eff, N_grids - 2, "TeffvsR.txt", true);
			WriteGraphData(vR, vT_c, N_grids - 2, "TcvsR.txt", true);
			WriteGraphData(vR, vT_irr, N_grids - 2, "TirrvsR.txt", true);
			WriteGraphData(vR, vT_sur, N_grids - 2, "TsurvsR.txt", true);
			WriteGraphData(vR, vH, N_grids - 2, "HvsR.txt", true);

			end = std::chrono::high_resolution_clock::now();
			elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
			cout << "Current time step is " << dT << " s.\n";
			cout << T/T_max * 100 << fixed << " percent completed! " << elapsed.count() << " ms have elapsed.\n";
		}
		//*************************************************************************************************
		//*************************************************************************************************
	}
	cout << "All done!";
	cin.get();
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

double OpticalThickness(double SurfaceDensity)
{
	const double thompson = 6.6524e-25;                                  // (thompson cross-section for electron) cm-2
	const double m_p = 1.6726231e-24;                                    // (mass of proton) g
	double Opacity = thompson / m_p;
	return SurfaceDensity * Opacity;
}

void ExtractSpectrum(double* T, double* X, double delta_X, int N_grids, double minEnergyEV, double maxEnergyEV, double resolutionEV, bool append)
{
	const double k = 1.380658e-16;                                       // (boltzmann's constant) erg K-1
	const double c = 2.99792458e10;                                      // (speed of light) cm / s
	const double h = 6.6260755e-27;                                      // (planck's constant) erg s
	const double PI = 3.14159265358979323846;
	int numberofchannels = (maxEnergyEV - minEnergyEV) / resolutionEV;
	double* I_bb = new double[numberofchannels];
	double* eV = new double[numberofchannels];
	parallel_for(0, numberofchannels, [=](int i)
	{
		I_bb[i] = 0;
		eV[i] = minEnergyEV + resolutionEV * i;
	});
	parallel_for(0, N_grids - 2, [=](int i)
	{
		parallel_for(0, numberofchannels, [=](int j)
		{
			if(T[i + 1] != 0)
				I_bb[j] += (2 * h * pow(eVtoHz(eV[j]), 3) / pow(c, 2)) / (exp((h * eVtoHz(eV[j])) /
					(k * T[i + 1])) - 1)
					* 4 * PI * pow(X[i], 3) * delta_X; // infinitesimal area
		});
	});
	WriteGraphData(eV, I_bb, numberofchannels, "powerspectrum.txt", append);
}
// Point source assumed
void IrradiationTemperature(double* T_irr, int N_grids, double nu, double epsilon, double L, double* R, double* H)
{
	const double a = 5.67051e-5;                                         // (stefan boltzmann constant) erg cm-2 K-4 s-1
	const double PI = 3.14159265358979323846;

	parallel_for(0, N_grids - 1, [=](int i) 
	{
		parallel_for(0, i + 1, [=](int j)
		{
			if (H[i + 1] < H[j])
				T_irr[i + 1] = 0;
			else
			{
				double C = nu * (1 - epsilon)*((H[i + 1] - H[i]) / (R[i + 1] - R[i]) - H[i + 1] / R[i + 1]);
				if (C > 0 && L > 0)
					T_irr[i + 1] = pow(C * L / (4 * PI * a * R[i + 1]), 0.25);
				else
					T_irr[i + 1] = 0;
			}
		});
	});
}

double* IrradiationTemperaturePower_4(int N_grids, double nu, double epsilon, double L, double* R, double* H)
{
	const double a = 5.67051e-5;                                         // (stefan boltzmann constant) erg cm-2 K-4 s-1
	const double PI = 3.14159265358979323846;

	double* T_irr = new double[N_grids];
	parallel_for(0, N_grids - 1, [=](int i)
	{
		double C = nu * (1 - epsilon)*((H[i + 1] - H[i]) / (R[i + 1] - R[i]) - H[i] / R[i]);
		T_irr[i] = C * L / (4 * PI * a * R[i]);
	});
	return T_irr;
}

double eVtoHz(double eV)
{
	return 2.417990504024e+14 * eV; // not sure...
}