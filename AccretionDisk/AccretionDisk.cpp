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
static const double rad_const = 7.5646e-15;									// (radiation constant) erg cm-3 K-4
static const double PI = 3.14159265358979323846;
//******************************************************************************************************************

// VARIABLES
static double M_compact;								// Mass of the compact object
static double M_disk;									// Mass of the disk
static double M_dot_boundary;							// Accretion rate
static double R_isco, R_outer, X_isco, X_outer;			// Inner and outer radius
static int N_grids;										// Number of grids
static double delta_X;									// Step size
static double mu_p_hot = 0.63;							// proton ratio (ertan, 2002)
static double mu_p_cold = 0.87;							// proton ratio
static double T_max;
static double T_L = 99999999999;
static double L_instant = 0, L_BB = 0;
static double L_previous = 0;
static double alpha_hot = 0.1;
static double alpha_cold = 0.033;
static double opacity[19][88];
static double logR[19] = { -8, -7.5, -7, -6.5, -6, -5.5, -5, -4.5, -4, -3.5, -3, -2.5, -2, -1.5, -1, -0.5, 0, 0.5, 1 };
static double logT[88];
static double L_edd = 0;
static double M_dot_edd = 0;
static double R_g = 0;
bool CoronaFormed = false;
bool MaximumLuminosityReached = false;
double R_Corona, Theta_Corona;
int trunc_radius_index = 0;



double OpticalThickness(double SurfaceDensity, double opacity);
double NormalizationGaussian(double Mu, double Sigma, double R_isco, double R_outer, double M_disk);
void WriteGraphData(double* X, double* Y, int length, string filename, bool append);
void WriteGraphData(double* X, int* Y, int length, string filename, bool append);
double eVtoHz(double eV);
void ExtractSpectrum(double* T, double* X, double delta_X, int N_grids, double minEnergyEV, double maxEnergyEV, double resolutionEV, bool append);
void IrradiationTemperature(double* T_irr, int N_grids, double nu, double epsilon, double L, double* R, double* H, double offset, bool shadow_enable);
double* IrradiationTemperaturePower_4(int N_grids, double nu, double epsilon, double L, double* R, double* H);
double mu_p(double alpha);
double alpha(double T);
void alpha(double* Alpha, double* T);
double VISC_C(double alpha, double opacity);
double average(double numbers[], int size); 
double alpha_alternative(double T_c, double T_irr, double R);
double alpha_alternative(double* vAlpha, double* T_c, double* T_irr, double* R);
double alpha_alternative_smooth(double* vAlpha, double* T_c, double* T_irr, double* R);
double T_critical(double T_irr, double R);
double T_c_max(double T_irr, double R);
double T_c_min(double T_irr, double R);
double irr_effect(double T_irr);
double BB_Luminosity(double* T_eff, double* X);
double BB_Luminosity(double* T_eff, double* X, int truncation_index);
double GetLuminosity(double* T, double* X, double delta_X, int N_grids, double minEnergyEV, double maxEnergyEV, double resolutionEV);
int FindHotIndex(double* vAlpha);
double EffectiveTemperature_alternative(double M_dot, double R);
double EffectiveTemperature_dubus2014(double E, double R, double V);
double CentralTemperature_dubus2014(double tau, double T_eff, double T_irr);
double E_max(double T_irr, double R);
double E_min(double T_irr, double R);
double alpha_compare(double* vAlpha, double* E, double* T_irr, double* R);
void RadiationPressure(double* P_rad, double* T_c);
void GasPressure(double* P_gas, double* E, double* H, double* T_c, double* Alpha); 
void ScaleHeight_FrankKingRaine(double * H, double * E, double * T_c, double* M_dot, double* R, double* Alpha);
int FindTruncationIndex(double* E, double* H, double* T_c, double* Alpha);
void GetShadows(int* Shadows, double* R, double* H, double R_corona, double Theta_corona); 
void CoronaIrradiationTemperature(double* T_irr, double L, double* R, double* H, double R_corona, double Theta_corona);
void CoronaIrradiationTemperature(double* T_irr, double nu, double epsilon, double L, double* R, double* H, double R_corona);
double TimeDependentCoronaRadius(double T_current, double T_corona, double T_rise, double R_initial, double R_max, bool Formed);
int FindTruncationIndex(double R_corona, double* R);
void GetGAS_RADIATION_ratio(double* GvsR, double* E, double* H, double* Alpha, double* T_c);
void WriteGraphData(double* X, double* Y, double T, int length, string filename, bool append);


int main(int argc, char **argv)
{
	double n;
	double T = 0;
	if (argc == 1)
	{
		cout << "Accretion disk simulation with parallel CPU computing.\n\n";
		cout << "Please enter the mass of the compact object. (M_solar)\n";
		cin >> n;
		M_compact = n * M_solar;
		cout << "Please enter the mass of the disk. (M_solar)\n";
		cin >> n;
		M_disk = n * M_solar;
		cout << "Please enter the accretion rate. (M_solar s-1)\n";
		cin >> n;
		M_dot_boundary = n * M_solar;

		R_g = 2 * G * M_compact / pow(c, 2);
	}
	else
	{
		M_compact = stod(argv[1]) * M_solar;
		M_disk = stod(argv[2]) * M_solar;
		M_dot_boundary = stod(argv[3]) * M_solar;
	}

	// Calculate some initial values
	R_isco = 6 * G * M_compact / pow(c, 2);
	R_outer = 5e10; // cm 
	X_isco = sqrt(R_isco);
	X_outer = sqrt(R_outer);
	L_edd = 1.3e38 * (M_compact / M_solar);
	M_dot_edd = L_edd * 2 * R_isco / (G* M_compact);		// Luminosity in ergs/s
	
	cout << "Mass of the compact object = " << M_compact << " g.\n";
	cout << "Mass of the disk           = " << M_disk << " g.\n";
	cout << "Innermost stable orbit     = " << R_isco << " cm.\n";
	cout << "Outer radius               = " << R_outer << " cm.\n";
	cout << "Accretion rate             = " << M_dot_boundary << " g s-1.\n";
	cout << "Eddington Luminosity       = " << L_edd << " erg s-1.\n";
	if (argc == 1)
	{
		cout << "Please enter the number of grids for the radial coordinate.\n";
		cin >> N_grids;
	}
	else
	{
		N_grids = stod(argv[8]);
	}
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
	double* vdelta_T = new double[N_grids];			// Surface mass density vector
	double* vH = new double[N_grids];
	double* vO = new double[N_grids];
	double* vAlpha = new double[N_grids];
	double* vGvsR = new double[N_grids];
	double dT = 0;
	double vJ_total = 0;
	double vM_total = 0;

	// Read opacity table******************************************************************************************************************************
	//********************************************************************************************************************************************************
	ifstream opal("table_merged.txt");
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

	parallel_for(0, N_grids, [=, &vJ_total, &vM_total](int i)
	{
		vX[i] = X_isco + i * delta_X;
		vR[i] = vX[i] * vX[i];
		vE[i] = E_0 * exp(-1 * (pow((vR[i] - Mu), 2.)) / (2. * pow(Sigma, 2.)));
		vS[i] = vX[i] * vE[i];
		
		vO[i] = thompson / m_p;
		vAlpha[i] = alpha_hot;
		
		vJ_total += 4 * PI * sqrt(G * M_compact) * vX[i] * vS[i] * delta_X;
		vM_total += 4 * PI * vS[i] * vX[i] * vX[i] * delta_X;
		
		vT_irr[i] = 0;
		vH[i] = 0;
		vdelta_T[i] = 0;
		vT_eff[i] = 0;
		vT_c[i] = 10;
		vM_dot[i] = 0;
		vGvsR[i] = 0;
	});


	// Inner and outer boundary condition
	R_Corona = R_isco;
	Theta_Corona = 0;
	vE[0] = 0;
	vS[0] = 0;
	vM_dot[N_grids - 1] = M_dot_boundary;
	vM_dot[N_grids - 2] = M_dot_boundary;
	//*************************

	for (int j = 0; j < 100; j++)
	{
		parallel_for(0, N_grids, [=](int i)
		{
			//************************************************************************
			// Calculate viscosity **************************************
			vV[i] = VISC_C(vAlpha[i], vO[i]) * pow(vX[i], (4. / 3.)) * pow(vS[i], (2. / 3.));
			//************************************************************************
			//************************************************************************

			vT_eff[i] = EffectiveTemperature_dubus2014(vE[i], vR[i], vV[i]);
			vT_c[i] = CentralTemperature_dubus2014(OpticalThickness(vE[i], vO[i]), vT_eff[i], 0);
			vH[i] = sqrt((vT_c[i + 1] * k * pow(vR[i + 1], 3)) / (mu_p(vAlpha[i + 1]) * m_p * G * M_compact));

			// Find Opacity value*****************************************************
			//************************************************************************
			int a = 0, b = 0;
			// 19 R values in the table
			for (int m = 0; m < 19; m++)
			{
				if (log10((vE[i] / (2 * vH[i])) / pow(vT_c[i] * 1e-6, 3)) >= logR[18])
				{
					a = 18;
					break;
				}
				if (log10((vE[i] / (2 * vH[i])) / pow(vT_c[i] * 1e-6, 3)) < logR[m])
				{
					a = m;
					break;
				}
			}
			// 88 T values in the table
			for (int n = 0; n < 88; n++)
			{
				if (log10(vT_c[i]) >= logT[87])
				{
					b = 87;
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
	}

	for (int i = 0; i < N_grids - 2; i++)
	{
		vM_dot[i] = 3. * PI * (vV[i + 1] * vS[i + 1] - vV[i] * vS[i]) / delta_X;
	}

	// Add efficiency constant
	L_instant = 0.1 * (vM_dot[0] * G * M_compact) / (2 * R_isco);		// Luminosity in ergs/s

	GetGAS_RADIATION_ratio(vGvsR, vE, vH, vAlpha, vT_c);				// Calculate P_g/P_r ratio

	WriteGraphData(vR, vGvsR, 0, N_grids, "PGvsR.txt", false);
	WriteGraphData(vR, vAlpha, 0, N_grids, "AvsR.txt", false);
	WriteGraphData(vR, vE, 0, N_grids, "EvsR.txt", false);
	WriteGraphData(vR, vV, 0, N_grids, "VvsR.txt", false);
	WriteGraphData(vR, vM_dot, 0, N_grids, "MdotvsR.txt", false);
	WriteGraphData(vR, vT_eff, 0, N_grids, "TeffvsR.txt", false);
	WriteGraphData(vR, vT_c, 0, N_grids, "TcvsR.txt", false);
	WriteGraphData(vR, vH, 0, N_grids, "HvsR.txt", false);
	WriteGraphData(vR, vO, 0, N_grids, "OvsR.txt", false);

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	cout << "Initial conditions created in " <<	elapsed.count() << " ms.\n";
	cout << "Total angular momentum is " << vJ_total << " g cm2 s-1.\n";
	cout << "Total mass is " << vM_total << " g.\n";
	//********************************************************************
	//********************************************************************
	
	if (argc == 1)
	{
		cout << "Please enter the evolution duration. (months)\n";
		cin >> n; T_max = n*month;
	}
	else
	{
		T_max = stod(argv[9])*month;
	}

	int N_L_sample = 500;			// 1000 yap
	double* vLT_sample = new double[N_L_sample];
	double deltaLT_sample = log(T_max) / (double)N_L_sample;
	for (int i = 1; i <= N_L_sample; i++)
	{
		vLT_sample[i - 1] = i*(T_max/N_L_sample);
	}

	start = chrono::high_resolution_clock::now();

	int j = 0;
	int l = 0;
	L_previous = 0.1;
	T = 0;
	bool message = true;

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

	ofstream file_Rhot;
	file_Rhot.open("Rhot_T.txt", ios::out);
	file_Rhot << T / day << "\t" << vR[FindHotIndex(vAlpha)] << "\n";
	file_Rhot.close(); 

	ofstream file_Rtrunc;
	file_Rtrunc.open("Rtrunc_T.txt", ios::out);
	file_Rtrunc << T / day << "\t" << vR[FindTruncationIndex(vE, vH, vT_c, vAlpha)] << "\n";
	file_Rtrunc.close();

	while (T < T_max)
	{

		// Determine outer boundary condition*************************************************************
		//************************************************************************************************
		int boundary = N_grids;
		vS[boundary - 1] = pow(((vM_dot[boundary - 2] * delta_X / (3. * PI) + vS[boundary - 2] * vV[boundary - 2])
			/ (VISC_C(vAlpha[boundary - 1], vO[boundary - 1]) * pow(vX[boundary - 1], (4. / 3.)))),
			(3. / 5.));
		vV[N_grids - 1] = VISC_C(vAlpha[N_grids - 1], vO[N_grids - 1]) * pow(vX[N_grids - 1], (4. / 3.)) * pow(vS[N_grids - 1], (2. / 3.));
		//************************************************************************************************
		//************************************************************************************************

		// Determine time step****************************************************************************
		//************************************************************************************************
		parallel_for(0, N_grids, [=] (int i)
		{
			vdelta_T[i] = (1. / 3.) * (pow(delta_X * vX[i], 2) / vV[i]) / 0.75; // 1/3 Courant Number better precision at 1/6
		});
		dT = *min_element(vdelta_T, vdelta_T + N_grids);
		//************************************************************************************************
		//************************************************************************************************

		// Iterate in radial coordinate ******************************************************************
		//************************************************************************************************
		parallel_for(1, N_grids - 1, [=](int i)
		{
			vS_new[i] = vS[i] + 0.75 * dT / pow((vX[i] * delta_X), 2) *
				(vV[i + 1] * vS[i + 1] + vV[i - 1] * vS[i - 1] - 2. * vV[i] * vS[i]);
			vV_new[i] = VISC_C(vAlpha[i], vO[i]) * pow(vX[i], (4. / 3.)) * pow(vS_new[i], (2. / 3.));
		});
		//*************************************************************************************************
		//*************************************************************************************************

		// Obtain new values*******************************************************************************
		//*************************************************************************************************
		parallel_for(1, N_grids - 1, [=](int i)
		{
			vS[i] = vS_new[i];
			vV[i] = vV_new[i];
		});
		//*************************************************************************************************
		//*************************************************************************************************

		for (int i = 0; i < N_grids - 2; i++)
		{
			vM_dot[i] = 3. * PI * (vV[i + 1] * vS[i + 1] - vV[i] * vS[i]) / delta_X;
		}

		parallel_for(0, N_grids, [=](int i)
		{
			vE[i] = vS[i] / vX[i];
		});

		// Calculate alpha values*******************
		alpha_compare(vAlpha, vE, vT_irr, vR);
		//******************************************

		for (int j = 0; j < 50; j++)
		{
			parallel_for(0, N_grids, [=](int i)
			{
				vT_eff[i] = EffectiveTemperature_dubus2014(vE[i], vR[i], vV[i]);
				vT_c[i] = CentralTemperature_dubus2014(OpticalThickness(vE[i], vO[i]), vT_eff[i], 0);
				vH[i] = sqrt((vT_c[i] * k * pow(vR[i], 3)) / (mu_p(vAlpha[i]) * m_p * G * M_compact));

				// Find Opacity value*****************************************************
				//************************************************************************
				int a = 0, b = 0;
				for (int m = 0; m < 19; m++)
				{
					if (log10((vE[i] / (2 * vH[i])) / pow(vT_c[i] * 1e-6, 3)) >= logR[18])
					{
						a = 18;
						break;
					}
					if (log10((vE[i] / (2 * vH[i])) / pow(vT_c[i] * 1e-6, 3)) < logR[m])
					{
						a = m;
						break;
					}
				}

				for (int n = 0; n < 88; n++)
				{
					if (log10(vT_c[i]) >= logT[87])
					{
						b = 87;
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
		}

		L_instant = 0.1 * (vM_dot[0] * G * M_compact) / (2 * vR[0]);		// Luminosity in ergs/s

		T += dT; // Increase time

		// Check if maximum luminosity has been reached***************************************************************************************************************
		//************************************************************************************************************************************************************
		if (L_instant < L_previous && L_instant > 1e37 && !MaximumLuminosityReached)
		{
			cout << std::scientific;
			cout << "Maximum luminosity reached -> L = " << L_instant << " erg/s at time T = " << T / day << " days.\n" << elapsed.count() << " ms have elapsed.\n\n";
			MaximumLuminosityReached = true;
		}
		else if(L_instant > L_previous && !MaximumLuminosityReached)
		{
			L_previous = L_instant;
		}
		//************************************************************************************************************************************************************
		//************************************************************************************************************************************************************

		// Record Data************************************************************************************************************************************************
		//************************************************************************************************************************************************************
		if (l < N_L_sample)
		{
			if (T >= vLT_sample[l])
			{
				l++;
				file.open("lightcurve.txt", ios::app);

				if (L_instant > 1e25)
					file << T / day << "\t" << L_instant << "\n";						// Write luminosity to file
				file.close();

				L_BB = BB_Luminosity(vT_eff, vX, trunc_radius_index);
				file_bolo.open("lightcurve_bolo.txt", ios::app);
				if (L_BB > 1e25)
					file_bolo << T / day << "\t" << L_BB << "\n";						// Write luminosity to file
				file_bolo.close();

				file_Rhot.open("Rhot_T.txt", ios::app);
				if (L_BB > 1e25)
					file_Rhot << T / day << "\t" << vR[FindHotIndex(vAlpha)] << "\n";						// Write luminosity to file
				file_Rhot.close();

				file_Rtrunc.open("Rtrunc_T.txt", ios::app);
				if (L_BB > 1e25)
					file_Rtrunc << T / day << "\t" << vR[FindTruncationIndex(vE, vH, vT_c, vAlpha)] << "\n";						// Write luminosity to file
				file_Rtrunc.close();

				if (MaximumLuminosityReached)
				{
					GetGAS_RADIATION_ratio(vGvsR, vE, vH, vAlpha, vT_c);				// Calculate P_g/P_r ratio

					WriteGraphData(vR, vGvsR, 0, N_grids, "PGvsR.txt", true);
					WriteGraphData(vR, vAlpha, T, N_grids, "AvsR.txt", true);
					WriteGraphData(vR, vE, T, N_grids, "EvsR.txt", true);
					WriteGraphData(vR, vV, T, N_grids, "VvsR.txt", true);
					WriteGraphData(vR, vM_dot, T, N_grids, "MdotvsR.txt", true);
					WriteGraphData(vR, vT_eff, T, N_grids, "TeffvsR.txt", true);
					WriteGraphData(vR, vT_c, T, N_grids, "TcvsR.txt", true);
					WriteGraphData(vR, vH, T, N_grids, "HvsR.txt", true);
					WriteGraphData(vR, vO, T, N_grids, "OvsR.txt", true);
				}

				end = std::chrono::high_resolution_clock::now();
				elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
				cout << "Current time step is " << dT << " s.\n";
				cout << std::scientific;
				vJ_total = 0;
				vM_total = 0;
				parallel_for(0, N_grids, [=, &vJ_total, &vM_total](int i)
				{
					vJ_total += 4 * PI * sqrt(G * M_compact) * vX[i] * vS[i] * delta_X;
					vM_total += 4 * PI * vS[i] * vX[i] * vX[i] * delta_X;
				});
				cout << "Total angular momentum is " << vJ_total << " g cm2 s-1.\n";
				cout << "Total mass is             " << vM_total << " g.\n";
				cout << "Hot disk radius is        " << vR[FindHotIndex(vAlpha)] << " cm.\n";
				cout << "Luminosity is             " << (L_instant / L_edd) * 100 << "\% of Eddington Luminosity.\n";
				cout << T / T_max * 100 << fixed << " percent completed! " << elapsed.count() << " ms have elapsed.\n\n";
			}
		}
		//************************************************************************************************************************************************************
		//************************************************************************************************************************************************************
	}
	cout << "All done!";
	if (argc == 1)
	{
		cin.get();
		cin.get();
	}
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
void WriteGraphData(double* X, int* Y, int length, string filename, bool append)
{
	FILE* file;
	if (!append)
		file = fopen(filename.c_str(), "w");
	else
		file = fopen(filename.c_str(), "a");
	parallel_for(0, length, [=, &file](int i)
	{
		fprintf(file, "%lf\t%d\n", X[i], Y[i]);
	});
	fprintf(file, "\n");
	fclose(file);
}
void WriteGraphData(double* X, double* Y, double T, int length, string filename, bool append)
{
	FILE* file;
	if (!append)
		file = fopen(filename.c_str(), "w");
	else
		file = fopen(filename.c_str(), "a");
	for(int i=0; i < length; i++)
	{
		fprintf(file, "%lf\t%lf\t%lf\n", T/day, X[i], Y[i]);
	}
	fprintf(file, "\n");
	fclose(file);
}

double OpticalThickness(double SurfaceDensity, double opacity)
{
	return SurfaceDensity * opacity;
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

double BB_Luminosity(double* T_eff, double* X, int truncation_index)
{
	double L = 0;
	for (int i = truncation_index; i < N_grids; i++)
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
void IrradiationTemperature(double* T_irr, int N_grids, double nu, double epsilon, double L, double* R, double* H, double offset, bool shadow_enable)
{
	bool shadow = false;
	parallel_for(0, N_grids - 1, [=, &shadow](int i)
	{
		if (shadow_enable)
		{
			int n = 0;
			shadow = false;
			parallel_for(0, i + 1, [=, &shadow, &n](int j)
			{
				if (atan((H[i + 1] - offset) / R[i + 1]) < atan((H[j] - offset) / R[j]))
				{
					n++;
				}
			});
			if (n > 1)
			{
				T_irr[i + 1] = 0;
				shadow = true;
			}
		}
		if (!shadow)
		{
			//double C = nu * (1 - epsilon)*((H[i + 1] - H[i]) / (R[i + 1] - R[i]) - H[i + 1] / R[i + 1]);
			double C = nu * 1e-3;/* nu * (1 - epsilon) * (2. / 35.)*(H[i + 1] / R[i + 1]); // lasota (2014)*/
			if (C > 0 && L > 0)
				T_irr[i + 1] = pow(C * L / (4 * PI * a * R[i + 1] * R[i + 1]), 0.25);
			else
				T_irr[i + 1] = 0;
		}
	});
}

void CoronaIrradiationTemperature(double* T_irr, double L, double* R, double* H, double R_corona, double Theta_corona)
{
	int* Shadows = new int[N_grids];
	GetShadows(Shadows, R, H, R_corona, Theta_corona);
	parallel_for(0, N_grids - 1, [=](int i)
	{
		if (Shadows[i + 1] != 1)
		{
			double R_2 = pow(H[i + 1], 2) + pow(R[i + 1], 2);	// square of radial coordinate of disk element
			double Theta = atan(H[i + 1] / R[i + 1]);			// elevation of disk element
			double d_2 = pow(R_corona, 2) + R_2 + 2 * R_corona * sqrt(R_2) * cos(Theta_corona - Theta);	// square of distance between corona element and disk element
			double beta_prime = acos((d_2 + R_2 - pow(R_corona, 2)) / (2 * sqrt(d_2 * R_2)));
			double slope_angle = atan((H[i + 1] - H[i]) / (R[i + 1] / R[i]));
			double beta = beta_prime + slope_angle - Theta;												// viewing angle of corona

			double L_disk = L / (4 * PI * d_2) * sin(beta);
			if (L > 0)
			{
				T_irr[i + 1] = pow(L_disk / a, 0.25);
			}
			else
				T_irr[i + 1] = 0;
		}
	});
}

// Integral
void CoronaIrradiationTemperature(double* T_irr, double nu, double epsilon, double L, double* R, double* H, double R_corona)
{
	int N_corona = 20;
	double* Theta_int = new double[N_corona];
	double Delta_Theta_corona = PI / (2 * (N_corona - 1));
	parallel_for(0, N_corona, [=](int i)
	{
		Theta_int[i] = 0 + i *  Delta_Theta_corona;
	});
	parallel_for(0, N_grids, [=](int i)
	{
		T_irr[i] = 0;
	});
	parallel_for(0, N_grids - 1, [=](int i)
	{
		parallel_for(0, N_corona, [=](int n)
		{
			bool Shadow = false;
			parallel_for(0, i + 1, [=, &Shadow](int j)
			{
				if (atan((H[i + 1] - R_corona * sin(Theta_int[n])) / (R[i + 1] - R_corona * cos(Theta_int[n])))
					<
					atan((H[j] -	 R_corona * sin(Theta_int[n])) / (R[j] -	 R_corona * cos(Theta_int[n]))))
				{
					Shadow = true;
				}
			});

			if (!Shadow && R[i + 1] > R_corona)
			{
				//double C = nu * (1 - epsilon)*((H[i + 1] - H[i]) / (R[i + 1] - R[i]) - H[i + 1] / R[i + 1]);
				//double C = nu * 1e-3;/* nu * (1 - epsilon) * (2. / 35.)*(H[i + 1] / R[i + 1]); // lasota (2014)*/
				
				double R_2 = pow(H[i + 1], 2) + pow(R[i + 1], 2);	// square of radial coordinate of disk element
				double Theta = atan(H[i + 1] / R[i + 1]);			// elevation of disk element
				double d_2 = pow(R_corona, 2) + R_2 + 2 * R_corona * sqrt(R_2) * cos(Theta_int[n] - Theta);	// square of distance between corona element and disk element
				
				double L_p = L * sin(Theta_int[n]) * Delta_Theta_corona * 0.5;								// luminosity of corona element

				double beta_prime = acos((d_2 + R_2 - pow(R_corona, 2)) / (2 * sqrt(d_2 * R_2)));			
				double slope_angle = atan((H[i + 1] - H[i]) / (R[i + 1] / R[i]));
				double beta = beta_prime + slope_angle - Theta;												// viewing angle of corona

				double L_disk = L_p / (4 * PI * d_2) * sin(beta);
				if (L > 0)
				{
						T_irr[i + 1] += pow(L_disk / a, 0.25);
				}
				else
					T_irr[i + 1] += 0;
			}
			else
			{
				T_irr[i + 1] += 0;
			}
		});
	});
}

void GetShadows(int* Shadows, double* R, double* H, double R_corona, double Theta_corona)
{
	Shadows[0] = 0;
	parallel_for(0, N_grids - 1, [=](int i)
	{
		int n = 0;
		parallel_for(0, i + 1, [=, &n](int j)
		{
			if (atan((H[i + 1] - R_corona * sin(Theta_corona)) / abs(R[i + 1] - R_corona * cos(Theta_corona))) <
				atan((H[j] -     R_corona * sin(Theta_corona)) / abs(R[j]     - R_corona * cos(Theta_corona))))
			{
				n++;
			}
		});
		if (n > 1)
		{
			Shadows[i + 1] = 1;
		}
		else Shadows[i + 1] = 0;
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

void alpha(double* Alpha, double* T)
{
	parallel_for(0, N_grids, [=](int i)
	{
		if (T[i] > 10000)
			Alpha[i] =  alpha_hot;
		else Alpha[i] = alpha_cold;
	});
}

double mu_p(double alpha)
{
	if (alpha == alpha_hot)
		return mu_p_hot;
	else return mu_p_cold;
}
int FindHotIndex(double* vAlpha)
{
	int i;
	for (i = 0; i < N_grids - 1; i++)
	{
		if (vAlpha[i] == alpha_cold && vAlpha[i + 1] == alpha_cold)
		{
			return i;
		}
	}
	return N_grids - 1;
}
double VISC_C(double alpha, double opacity)
{
	double VIS_C = pow(
		(27. / 32.) * (pow(alpha, 4) * pow(k, 4) * opacity)
		/ (pow(mu_p(alpha), 4) * pow(m_p, 4) * a * G * M_compact)
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
	double v1 = T_c_max(T_irr, R);
	double v2 = T_c_min(T_irr, R);
	double v3 = irr_effect(T_irr);
	double v4 = pow(R / 1e10, -0.05 * irr_effect(T_irr));
	double v5 = pow(R / 1e10, 0.05 - 0.12 * irr_effect(T_irr));
	double v6 = -0.05 * irr_effect(T_irr);
	double v7 = 0.05 - 0.12 * irr_effect(T_irr);
	double criticaltemp = T_critical(T_irr, R);
	double power = pow(1 + pow(criticaltemp / T_c, 8), -1);
	double calc = alpha_cold * pow((alpha_hot / alpha_cold), power);
	if(calc > 0)
		return calc;
	else
	{
		//cout << "alpha = " << std::scientific << calc;
		return alpha_hot;
	}
	if (T_c_max(T_irr, R) > T_c_min(T_irr, R))
	{
		return alpha_hot;
	}
	else
		return alpha_cold;
}

double alpha_alternative(double* vAlpha, double* T_c, double* T_irr, double* R)
{
	parallel_for(0, N_grids, [=](int i)
	{
		double criticaltemp = T_critical(T_irr[i], R[i]);
		double power = pow(1 + pow(criticaltemp / T_c[i], 8), -1);
		double calc = alpha_cold * pow((alpha_hot / alpha_cold), power);
		if (calc > 0)
			vAlpha[i] = calc;
		else
		{
			//cout << "alpha = " << std::scientific << calc;
			vAlpha[i] = alpha_hot;
		}
	});
}

double alpha_compare(double* vAlpha, double* E, double* T_irr, double* R)
{
	parallel_for(0, N_grids, [=](int i)
	{
		if (E[i] > E_max(T_irr[i], R[i]))
			vAlpha[i] = alpha_hot;
		if (E[i] < E_min(T_irr[i], R[i]))
			vAlpha[i] = alpha_cold;
	});
}

double alpha_alternative_smooth(double* vAlpha, double* T_c, double* T_irr, double* R)
{
	parallel_for(0, N_grids, [=](int i)
	{
		double criticaltemp = T_critical(T_irr[i], R[i]);
		double power = pow(1 + pow(criticaltemp / T_c[i], 8), -1);
		double calc = alpha_cold * pow((alpha_hot / alpha_cold), power);
		if (calc > 0)
			vAlpha[i] = calc;
		else
		{
			vAlpha[i] = alpha_hot;
		}
	});

	int j = FindHotIndex(vAlpha);

	parallel_for(0, N_grids, [=](int i)
	{
		if (i < j)
		{
			vAlpha[i] = alpha_hot;
		}
		else
			vAlpha[i] = alpha_cold;
	});
}
double T_critical(double T_irr, double R)
{
	if (T_c_max(T_irr, R) == std::numeric_limits<double>::infinity())
		return 0;
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
double E_max(double T_irr, double R)
{
	return (10.8 - 10.3 * irr_effect(T_irr)) * pow(alpha_cold, -0.84) * pow(M_compact / M_solar, -0.37 + 0.1 * irr_effect(T_irr))
		* pow(R / 1e10, 1.11 - 0.27 * irr_effect(T_irr));
}
double E_min(double T_irr, double R)
{
	return (8.3 - 7.1 * irr_effect(T_irr)) * pow(alpha_hot, -0.77) * pow(M_compact / M_solar, -0.37)
		* pow(R / 1e10, 1.12 - 0.23 * irr_effect(T_irr));
}

double EffectiveTemperature_alternative(double M_dot, double R)
{
	if (M_dot > 0)
	{
		double T_isco = pow((3 * G * M_compact * M_dot) / (8 * PI * pow(R_isco, 3)) * a, 0.25);
		return T_isco * pow(R / R_isco, -0.75) * pow((1 - pow(R_isco / R, 0.5)), 0.25);
	}
	else return 0;
}

double EffectiveTemperature_dubus2014(double E, double R, double V)
{
	return pow((9. / 8.) * V * E * G * M_compact / (a * pow(R, 3)), 0.25);
}

double CentralTemperature_dubus2014(double tau, double T_eff, double T_irr)
{
	return pow(((3. / 8.) * tau * pow(T_eff, 4) + pow(T_irr, 4)), 0.25);
}

void RadiationPressure(double* P_rad, double* T_c)
{
	parallel_for(0, N_grids, [=](int i)
	{
		P_rad[i] = 4. / 3. * a / c /* rad_const*/ * pow(T_c[i], 4);
	});
}

void GasPressure(double* P_gas, double* E, double* H, double* T_c, double* Alpha)
{
	parallel_for(0, N_grids, [=](int i)
	{
		P_gas[i] = ((E[i] / H[i]) * k * T_c[i]) / (mu_p(Alpha[i]) * m_p);
	});
}

void ScaleHeight_FrankKingRaine(double * H, double * E, double * T_c, double* M_dot, double* R, double* Alpha)
{
	double* P_rad = new double[N_grids];
	double* P_gas = new double[N_grids];
	RadiationPressure(P_rad, T_c);
	GasPressure(P_gas, E, H, T_c, Alpha);

	for (int j = 0; j < 10; j++)
	{
		parallel_for(0, N_grids, [=](int i)
		{
			if (P_gas[i]*10 < P_rad[i])
			{
				H[i] = ((3 * R_isco) / (4 * 0.1)) * (M_dot[i] / M_dot_edd) * (1 - pow(R_isco / R[i], 0.5));
			}
			else
			{
				H[i] = sqrt((T_c[i] * k * pow(R[i], 3)) / (mu_p(Alpha[i]) * m_p * G * M_compact));
			}
		});
		GasPressure(P_gas, E, H, T_c, Alpha);
	}
}

int FindTruncationIndex(double* E, double* H, double* T_c, double* Alpha)
{
	double* P_rad = new double[N_grids];
	double* P_gas = new double[N_grids];

	RadiationPressure(P_rad, T_c);
	GasPressure(P_gas, E, H, T_c, Alpha);

	for (int i = N_grids - 2; i >= 0; i--)
	{
		if (P_rad[i] > P_gas[i] + 0.1 && P_rad[i + 1] > P_gas[i + 1])
		{
			int j = FindHotIndex(Alpha);
			if (i <= j)
			{
				return i;
			}
			else return j;
		}
	}
	return 0;
}

int FindTruncationIndex(double R_corona, double* R)
{
	for (int i = 0; i < N_grids; i++)
	{
		if (R[i] > R_corona)
		{
			return i;
		}
	}
}

double TimeDependentCoronaRadius(double T_current, double T_corona, double T_rise, double R_initial, double R_max, bool Formed)
{
	if (Formed)
	{
		if (T_current < T_corona + T_rise)
		{
			double c_1 = (R_initial - R_max) / pow(T_rise, 2);
			double c_2 = -2 * c_1 * (T_corona + T_rise);
			double c_3 = R_max + pow(T_corona + T_rise, 2) * c_1;
			return c_1 * pow(T_current, 2) + c_2 * T_current + c_3;
		}
		else
			return R_max;
	}
	else return R_isco;
}

void GetGAS_RADIATION_ratio(double* GvsR, double* E, double* H, double* Alpha, double* T_c)
{
	double* P_rad = new double[N_grids];
	double* P_gas = new double[N_grids];

	RadiationPressure(P_rad, T_c);
	GasPressure(P_gas, E, H, T_c, Alpha);

	parallel_for(0, N_grids, [=](int i)
	{
		GvsR[i] = P_rad[i] / P_gas[i];
	});
}


