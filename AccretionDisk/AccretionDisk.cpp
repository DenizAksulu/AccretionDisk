#include "stdafx.h"
#include <fstream>
#include <iostream>
//#include <amp_math.h>
#include <algorithm>
//#include <amp.h>
#include <ppl.h>
#include <chrono>
#include <sstream>

using namespace std;
using namespace :: Concurrency;
//using namespace ::Concurrency::fast_math;

// CONSTANTS*******************************************************************
static const double day = 86400e0;                                          // (day) s
static const double month = 2629743.83e0;                                   // (month) s
static const double year = 31556926e0;                                      // (year) s
static const double M_solar = 1.9891e33;                                    // (solar mass) g
static const double AU = 1.496e13;                                          // (astronomical unit) cm
static const double ly = 9.463e17;                                          // (light year) cm
static const double pc = 3.086e18;                                          // (parsec)cm
static const double c = 2.99792458e10;                                      // (speed of light) cm / s
static const double h = 6.6260755e-27;                                      // (planck's constant) erg s
static const double hbar = 1.05457266e-27;                                  // (reduced planck's constant) erg s
static const double G = 6.67259e-8;                                         // (gravitational constant) cm3 g-1 s-2
static const double e = 4.8032068e-10;                                      // (electron charge) esu
static const double m_e = 9.1093897e-28;                                    // (mass of electron) g
static const double m_p = 1.6726231e-24;                                    // (mass of proton) g
static const double m_n = 1.6749286e-24;                                    // (mass of neutron) g
static const double m_H = 1.6733e-24;                                       // (mass of hydrogen) g
static const double k = 1.380658e-16;                                       // (boltzmann's constant) erg K-1
static const double eV = 1.6021772e-12;                                     // (electron volt) erg
static const double a = 5.67051e-5;                                         // (stefan boltzmann constant) erg cm-2 K-4 s-1
static const double thompson = 6.6524e-25;                                  // (thompson cross-section for electron) cm-2
static const double PI = 3.14159265358979323846;
//******************************************************************************************************************

// VARIABLES
static double M_compact;								// Mass of the compact object
static double M_disk;									// Mass of the disk
static double M_dot_boundary;									// Accretion rate
static double R_isco, R_outer, X_isco, X_outer;		// Inner and outer radius
static int N_grids;									// Number of grids
static int N_sample;									// Number of samples
static double delta_X;									// Step size
static double mu_p_hot = 0.63;								// proton ratio (ertan, 2002)
static double mu_p_cold = 0.87;								// proton ratio
static double T_max;
static double T_L = 99999999999;
static double T_corona;
static double L_instant = 0;
static double L_BB = 0;
static double L_optical = 0;
static double L_X = 0;
static double L_previous = 0;
static double alpha_hot = 0.1;
static double alpha_cold = 0.033;
static bool Corona = false;
static double opacity[19][70];
static double logR[19] = { -8, -7.5, -7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1 };
static double logT[70];

double OpticalThickness(double SurfaceDensity);
double NormalizationGaussian(double Mu, double Sigma, double R_isco, double R_outer, double M_disk);
void WriteGraphData(double* X, double* Y, int length, string filename, bool append);
double eVtoHz(double eV);
void ExtractSpectrum(double* T, double* X, double delta_X, int N_grids, double minEnergyEV, double maxEnergyEV, double resolutionEV, bool append);
void IrradiationTemperature(double* T_irr, int N_grids, double nu, double epsilon, double L, double* R, double* H);
double* IrradiationTemperaturePower_4(int N_grids, double nu, double epsilon, double L, double* R, double* H);
double alpha(double T);									// alpha parameter for disk
double mu_p(double T);									// 
double VISC_C(double T_c, double T_irr, double R);
double average(double numbers[], int size); 
double alpha_alternative(double T_c, double T_irr, double R);
double T_c_max(double T_irr, double R);
double T_c_min(double T_irr, double R);
double irr_effect(double T_irr);
double BB_Luminosity(double* T_eff, double* X);
double GetLuminosity(double* T, double* X, double delta_X, int N_grids, double minEnergyEV, double maxEnergyEV, double resolutionEV);

int main()
{
	double T = 0;
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

	// Calculate some initial values
	R_isco = 6 * G * M_compact / pow(c, 2);
	R_outer = 3e11; // cm
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
	double* vO = new double[N_grids];
	double dT = 0;
	double vJ_total = 0;
	vM_dot[N_grids - 1] = M_dot_boundary;
	vM_dot[N_grids - 2] = M_dot_boundary;

	// Read opacity table******************************************************************************************************************************
	//********************************************************************************************************************************************************
	ifstream opal("table#73.txt");
	if (opal.is_open())
	{
		int row = 0, col = 0;
		while (!opal.eof())
		{
			string line;
			getline(opal, line);
			stringstream ss(line);
			ss >> logT[row];
			col = 0;
			while (ss >> opacity[col][row])
			{
				col++;
			}
			row++;
		}
	}
	else
	{
		cout << "Opacity table not found!\n";
	}
	
	// Create initial conditions******************************************************************************************************************************
	//********************************************************************************************************************************************************
	auto start = chrono::high_resolution_clock::now();

	
	double Mu = 5e10; // cm
	double Sigma = 1e10;

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
		vV[i] = VISC_C(0, 0, vR[i]) * pow(vX[i], (4. / 3.)) * pow(vS[i], (2. / 3.));
		vJ_total += 4 * PI * sqrt(G * M_compact) * vX[i] * vS[i] * delta_X;
		vO[i] = thompson / m_p;
		vT_irr[i] = 0;
		vH[i] = 0;
		vdelta_T[i] = 0;
		vT_eff[i] = 0;
		vT_c[i] = 0;
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
		vH[i] = sqrt((vT_c[i] * k * pow(vR[i], 3)) / (mu_p(vT_eff[i]) * m_p * G * M_compact));
		// Find Opacity value*****************************************************
		//************************************************************************
		int a = 0, b = 0;
		for (int m = 0; m < 19; m++)
		{
			if (vH[i] > 0 && vH[i] < 1e13)
			{
				if (log10((vE[i] / vH[i]) / pow(vT_c[i] * 1e-6, 3)) >= logR[18])
				{
					a = 18;
					break;
				}
				if (log10((vE[i] / vH[i]) / pow(vT_c[i] * 1e-6, 3)) < logR[m])
				{
					a = m;
					break;
				}
			}
		}

		for (int n = 0; n < 70; n++)
		{
			if (log10(vT_c[i]) >= logT[69])
			{
				b = 69;
				break;
			}
			if (log10(vT_c[i]) < logT[n])
			{
				b = n;
				break;
			}
		}
		vO[i] = pow(10, opacity[a][b]);
		//************************************************************************
		//************************************************************************
	});

	L_instant = (vM_dot[0] * G * M_compact) / (2 * R_isco);		// Luminosity in ergs/s

	IrradiationTemperature(vT_irr, N_grids, 2, 0.5, L_instant, vR, vH);

	parallel_for(0, N_grids - 2, [=](int i)
	{
		vT_sur[i + 1] = pow(pow(vT_irr[i + 1], 4) + pow(vT_eff[i + 1], 4), 0.25);
	});

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
	WriteGraphData(vR, vO, N_grids, "OvsR.txt", false);



	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	cout << "Initial conditions created in " <<	elapsed.count() << " ms.\n";
	cout << "Total angular momentum is " << vJ_total << " g cm2 s-1.\n";
	//********************************************************************
	//********************************************************************

	cout << "Corona formation time after main outburst? (days)\n";
	cin >> n; T_corona = n*day;

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
	int N_L_sample = 500;
	double* vLT_sample = new double[N_L_sample];
	double deltaLT_sample = log(T_max) / (double)N_L_sample;
	for (int i = 1; i <= N_L_sample; i++)
	{
		vLT_sample[i - 1] = exp(i*deltaLT_sample);
	}
	int j = 0;
	int l = 0;
	start = chrono::high_resolution_clock::now();

	ofstream file;
	file.open("lightcurve.txt", ios::out);
	if (L_instant > 0)
		file << T << "\t" << L_instant << "\n";						// Write luminosity to file
	file.close();

	ofstream file_bolo;
	file_bolo.open("lightcurve_bolo.txt", ios::out);
	if (L_BB > 1e20)
		file_bolo << T << "\t" << L_BB << "\n";						// Write luminosity to file
	file_bolo.close();

	ofstream file_optical;
	file_optical.open("lightcurve_optical.txt", ios::out);
	if (L_optical > 1e20)
		file_optical << T << "\t" << L_optical << "\n";						// Write luminosity to file
	file_optical.close();

	ofstream file_X;
	file_X.open("lightcurve_X.txt", ios::out);
	if (L_X > 1e20)
		file_X << T << "\t" << L_X << "\n";						// Write luminosity to file
	file_X.close();

	L_previous = 0.1;
	T = 0;
	bool message = true;
	while (T < T_max)
	{
		// Determine outer boundary condition*************************************************************
		//************************************************************************************************
		vS[N_grids - 1] = pow(((vM_dot[N_grids - 2] * delta_X / (3. * PI) + vS[N_grids - 2] * vV[N_grids - 2]) 
			/ (VISC_C(vT_c[N_grids - 1], vT_irr[N_grids - 1], vR[N_grids - 1]) * pow(vX[N_grids - 1], (4. / 3.)))),
			(3. / 5.));
		//************************************************************************************************
		//************************************************************************************************

		// Determine time step****************************************************************************
		//************************************************************************************************

		parallel_for(1, N_grids, [=] (int i)
		{
			vdelta_T[i] = (1. / 6.) * (pow(delta_X * vX[i], 2) / vV[i]) / 0.75;
		});
		dT = *min_element(vdelta_T + 1, vdelta_T + N_grids);
		
		//************************************************************************************************
		//************************************************************************************************

		// Iterate in radial coordinate ******************************************************************
		//************************************************************************************************
		parallel_for(0, N_grids - 2, [=](int i)
		{
			
			vS_new[i + 1] = vS[i + 1] + 0.75 * dT / pow((vX[i + 1] * delta_X), 2) *
				(vV[i + 2] * vS[i + 2] + vV[i] * vS[i] - 2. * vV[i + 1] * vS[i + 1]);
			vV_new[i + 1] = VISC_C(vT_c[i + 1], vT_irr[i + 1], vR[i + 1]) * pow(vX[i + 1], (4./3.)) * pow(vS_new[i + 1], (2. / 3.));
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
		parallel_for(0, N_grids - 2, [=](int i)
		{
			vM_dot[i] = 3. * PI * (vV[i + 1] * vS[i + 1] - vV[i] * vS[i]) / delta_X;
		});
		parallel_for(0, N_grids, [=](int i)
		{
			vE[i] = vS[i] / vX[i];
		});
		parallel_for(0, N_grids - 2, [=](int i)
		{
			vT_eff[i + 1] = pow(3 * G * M_compact * abs(vM_dot[i + 1]) / (8 * PI * a * pow(vX[i + 1], 6)) * (1 - sqrt(X_isco / pow(vX[i + 1], 2))), 0.25);
			vT_c[i + 1] = pow(3 * OpticalThickness(vE[i + 1]) * pow(vT_eff[i + 1], 4) / 4, 0.25);
			vH[i + 1] = sqrt((vT_c[i + 1] * k * pow(vR[i + 1], 3)) / (mu_p(vT_eff[i + 1]) * m_p * G * M_compact));


			// Find Opacity value*****************************************************
			//************************************************************************
			int a = 0, b = 0;
			for (int m = 0; m < 19; m++)
			{
				if (log10((vE[i + 1] / vH[i + 1]) / pow(vT_c[i + 1] * 1e-6, 3)) >= logR[18])
				{
					a = 18;
					break;
				}
				if (log10((vE[i + 1] / vH[i + 1]) / pow(vT_c[i + 1] * 1e-6, 3)) < logR[m])
				{
					a = m;
					break;
				}
			}

			for (int n = 0; n < 70; n++)
			{
				if (log10(vT_c[i + 1]) >= logT[69])
				{
					b = 69;
					break;
				}
				if (log10(vT_c[i + 1]) < logT[n])
				{
					b = n;
					break;
				}
			}
			vO[i + 1] = pow(10, opacity[a][b]);
			//************************************************************************
			//************************************************************************
		});

		L_instant = (vM_dot[0] * G * M_compact) / (2 * R_isco);		// Luminosity in ergs/s

		
		if (T > T_L + T_corona)
		{
			if (message)
			{
				cout << "Corona has formed at time T = " << T/day << " days.\n" << elapsed.count() << " ms have elapsed.\n";
				message = false;
			}
			IrradiationTemperature(vT_irr, N_grids, 10, 0.9, L_instant, vR, vH); // Dubus et. al. (2014)
		}
		else
		{
			IrradiationTemperature(vT_irr, N_grids, 0.1, 0.9, L_instant, vR, vH); // Dubus et. al. (2014)
		}
		parallel_for(0, N_grids - 2, [=](int i)
		{
			vT_c[i + 1] = pow(3 * vE[i + 1] * vO[i + 1] * (pow(vT_eff[i + 1], 4)) / 8 + pow(vT_irr[i + 1], 4), 0.25); // Dubus et. al. (2014)
			vH[i + 1] = sqrt((vT_c[i + 1] * k * pow(vR[i + 1], 3)) / (mu_p(vT_eff[i + 1]) * m_p * G * M_compact));
			vV[i + 1] = alpha_alternative(vT_c[i + 1], vT_irr[i + 1], vR[i + 1]) * sqrt((vT_c[i + 1] * k) / (mu_p(vT_eff[i + 1]) * m_p)) * vH[i + 1]; //
		});
		T += dT; // Increase time

		if (l < N_L_sample)
		{
			if (T >= vLT_sample[l])
			{
				if (L_instant < L_previous)
				{
					cout << "Maximum luminosity reached -> L = " << L_instant << " erg/s at time T = " << T << " s.\n" << elapsed.count() << " ms have elapsed.\n";
					T_L = T;
					L_previous = 0;
				}
				else if (L_previous != 0)
				{
					L_previous = L_instant;
				}
				l++;
				file.open("lightcurve.txt", ios::app);
				if (L_instant > 0)
					file << T / day << "\t" << L_instant << "\n";						// Write luminosity to file
				file.close();
				parallel_for(0, N_grids, [=](int i)
				{
					vT_sur[i] = pow(pow(vT_irr[i], 4) + pow(vT_eff[i], 4), 0.25);
				});
				L_BB = BB_Luminosity(vT_sur, vX);
				file_bolo.open("lightcurve_bolo.txt", ios::app);
				if (L_BB > 1e20)
					file_bolo << T / day << "\t" << L_BB << "\n";						// Write luminosity to file
				file_bolo.close();

				L_optical = GetLuminosity(vT_sur, vX, delta_X, N_grids, 1, 4, 0.01);
				file_optical.open("lightcurve_optical.txt", ios::app);
				if (L_optical > 1e20)
					file_optical << T / day << "\t" << L_optical << "\n";						// Write luminosity to file
				file_optical.close();

				L_X = GetLuminosity(vT_sur, vX, delta_X, N_grids, 1000, 10000, 1);
				file_X.open("lightcurve_X.txt", ios::app);
				if (L_X > 1e20)
					file_X << T / day << "\t" << L_X << "\n";						// Write luminosity to file
				file_X.close();
			}
		}
		// Take samples ***********************************************************************************
		//*************************************************************************************************
		if (T >= vT_sample[j])
		{
			j++;

			ExtractSpectrum(vT_eff, vX, delta_X, N_grids - 2, 1, 15000, 1, true);

			WriteGraphData(vR, vE, N_grids, "EvsR.txt", true);
			WriteGraphData(vR, vV, N_grids, "VvsR.txt", true);
			WriteGraphData(vR, vM_dot, N_grids- 1, "MdotvsR.txt", true);
			WriteGraphData(vR, vT_eff, N_grids - 2, "TeffvsR.txt", true);
			WriteGraphData(vR, vT_c, N_grids - 2, "TcvsR.txt", true);
			WriteGraphData(vR, vT_irr, N_grids - 2, "TirrvsR.txt", true);
			WriteGraphData(vR, vT_sur, N_grids - 2, "TsurvsR.txt", true);
			WriteGraphData(vR, vH, N_grids - 2, "HvsR.txt", true);
			WriteGraphData(vR, vO, N_grids, "OvsR.txt", true);

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
	double Opacity = thompson / m_p;
	return SurfaceDensity * Opacity;
}

double BB_Luminosity(double* T_eff, double* X)
{
	double L = 0;
	for(int i = 0; i < N_grids; i++)
	{
		L += pow(X[i], 3) * pow(T_eff[i], 4) * a * 4 * PI * delta_X;
	}
	return L;
}

void ExtractSpectrum(double* T, double* X, double delta_X, int N_grids, double minEnergyEV, double maxEnergyEV, double resolutionEV, bool append)
{
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

double GetLuminosity(double* T, double* X, double delta_X, int N_grids, double minEnergyEV, double maxEnergyEV, double resolutionEV)
{
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
			if (T[i + 1] != 0)
				I_bb[j] += (2 * h * pow(eVtoHz(eV[j]), 3) / pow(c, 2)) / (exp((h * eVtoHz(eV[j])) /
				(k * T[i + 1])) - 1)
				* 4 * PI * pow(X[i], 3) * delta_X; // infinitesimal area
		});
	});
	double L = 0;
	parallel_for(0, numberofchannels, [=, &L](int i)
	{
		L += I_bb[i];
	});
	return L;
}
// Point source assumed
void IrradiationTemperature(double* T_irr, int N_grids, double nu, double epsilon, double L, double* R, double* H)
{
	bool shadow = false;
	parallel_for(0, N_grids - 1, [=, &shadow](int i)
	{
		int n = 0;
		shadow = false;
		parallel_for(0, i + 1, [=, &shadow, &n](int j)
		{
			if (atan((H[i + 1] - 1e6) / R[i + 1]) < atan((H[j] - 1e6) / R[j]))
			{
				n++;
				if (n > 10)
				{
					T_irr[i + 1] = 0;
					shadow = true;
				}
			}
		});
		if (!shadow)
		{
			//double C = nu * (1 - epsilon)*((H[i + 1] - H[i]) / (R[i + 1] - R[i]) - H[i + 1] / R[i + 1]);
			double C = /*5e-3;*/ nu * (1 - epsilon)*(2. / 35.)*(H[i + 1] / R[i + 1]); // lasota (2014)
			if (C > 0 && L > 0)
				T_irr[i + 1] = pow(C * L / (4 * PI * a * R[i + 1] * R[i + 1]), 0.25);
			else
				T_irr[i + 1] = 0;
		}
	});
}

double* IrradiationTemperaturePower_4(int N_grids, double nu, double epsilon, double L, double* R, double* H)
{

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

double alpha(double T)
{
	if (T > 10000)
		return alpha_hot;
	else return alpha_cold;
}

double mu_p(double T)
{
	if (T > 10000)
		return mu_p_hot;
	else return mu_p_cold;
}

double VISC_C(double T_c, double T_irr, double R)
{
	double VIS_C = pow(
		(27. / 32.) * (pow(alpha_alternative(T_c, T_irr, R), 4) * pow(k, 4) * thompson)
		/ (pow(mu_p(T_c), 4) * pow(m_p, 5) * a * G * M_compact)
		, (1. / 3.));
	return VIS_C;
}

double average(double numbers[], int size) {
	double sum = 0;
	parallel_for(0, size, [=, &sum](int i)
	{
		if(numbers[i] < 100)
			sum += numbers[i];
	});
	return sum / (double)size;
}
/*
 *  Dubus et. al. (2001)
 */
double alpha_alternative(double T_c, double T_irr, double R)
{
	//return exp(log(alpha_cold) + (log(alpha_hot) - log(alpha_cold)) * pow(1 + pow(T_critical(T_irr, R) / T_c, 8), -1));
	if (T_c_max(T_irr, R) > T_c_min(T_irr, R))
	{
		return alpha_hot;
	}
	else
		return alpha_cold;
}
double T_critical(double T_irr, double R)
{
	return 0.5 * (T_c_max(T_irr, R) + T_c_min(T_irr, R));
}
double T_c_max(double T_irr, double R)
{
	return 10700 * pow(alpha_cold, -0.1) * pow(R / 1e10, -0.05 * irr_effect(T_irr));
}
double T_c_min(double T_irr, double R)
{
	return (20900 - 11300 * irr_effect(T_irr)) * pow(alpha_hot, -0.22)* pow(M_compact / M_solar, -0.01) * pow(R / 1e10, 0.05 - 0.12 * irr_effect(T_irr));
}
double irr_effect(double T_irr)
{
	return pow(T_irr / 1e4, 2);
}