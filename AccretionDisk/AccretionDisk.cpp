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
static double M_dot_boundary;									// Accretion rate
static double R_isco, R_outer, X_isco, X_outer;		// Inner and outer radius
static int N_grids;									// Number of grids
static int N_sample;									// Number of samples
static double delta_X;									// Step size
static double mu_p_hot = 0.63;								// proton ratio (ertan, 2002)
static double mu_p_cold = 0.87;								// proton ratio
static double T_max;
static double T_L = 99999999999;
static double T_corona = 0;
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
static double L_edd = 0;
static double M_dot_edd = 0;
static double R_point_source = 0;
static double Theta_point_source = 0;

static double R_corona_max = 0;
static double R_corona_initial = 0;
static double T_corona_rise = 0;
bool CoronaFormed = false;



double OpticalThickness(double SurfaceDensity, double opacity);
double NormalizationGaussian(double Mu, double Sigma, double R_isco, double R_outer, double M_disk);
void WriteGraphData(double* X, double* Y, int length, string filename, bool append);
void WriteGraphData(double* X, int* Y, int length, string filename, bool append);
double eVtoHz(double eV);
void ExtractSpectrum(double* T, double* X, double delta_X, int N_grids, double minEnergyEV, double maxEnergyEV, double resolutionEV, bool append);
void IrradiationTemperature(double* T_irr, int N_grids, double nu, double epsilon, double L, double* R, double* H, double offset, bool shadow_enable);
double* IrradiationTemperaturePower_4(int N_grids, double nu, double epsilon, double L, double* R, double* H);
double alpha(double T);									// alpha parameter for disk
double mu_p(double T);									// 
double VISC_C(double alpha, double T_c, double opacity);
double average(double numbers[], int size); 
double alpha_alternative(double T_c, double T_irr, double R);
double alpha_alternative(double* vAlpha, double* T_c, double* T_irr, double* R);
double alpha_alternative_smooth(double* vAlpha, double* T_c, double* T_irr, double* R);
double T_critical(double T_irr, double R);
double T_c_max(double T_irr, double R);
double T_c_min(double T_irr, double R);
double irr_effect(double T_irr);
double BB_Luminosity(double* T_eff, double* X);
double GetLuminosity(double* T, double* X, double delta_X, int N_grids, double minEnergyEV, double maxEnergyEV, double resolutionEV);
int FindHotIndex(double* vAlpha);
double EffectiveTemperature_alternative(double M_dot, double R);
double EffectiveTemperature_dubus2014(double E, double R, double V);
double CentralTemperature_dubus2014(double tau, double T_eff, double T_irr);
double E_max(double T_irr, double R);
double E_min(double T_irr, double R);
double alpha_compare(double* vAlpha, double* E, double* T_irr, double* R);
void RadiationPressure(double* P_rad, double* T_c);
void GasPressure(double* P_gas, double* E, double* H, double* T_c);
void ScaleHeight_FrankKingRaine(double* H, double* E, double* T_c, double* M_dot, double* R);
int FindTruncationIndex(double* E, double* H, double* T_c, double* Alpha);
void GetShadows(int* Shadows, double* R, double* H, double R_corona, double Theta_corona);
void CoronaIrradiationTemperature(double* T_irr, double nu, double epsilon, double L, double* R, double* H, double R_corona, double Theta_corona);
void CoronaIrradiationTemperature(double* T_irr, double nu, double epsilon, double L, double* R, double* H, double R_corona);
double TimeDependentCoronaRadius(double T_current, double T_corona, double T_rise, double R_initial, double R_max);

int main()
{
	double T = 0;
	cout << "Accretion disk simulation with parallel CPU computing.\n\n";

	cout << "Please enter the mass of the compact object. (M_solar)\n";
	double n; cin >> n;
	M_compact = n * M_solar;
	cout << "Please enter the mass of the disk. (M_solar)\n";
	cin >> n;
	M_disk = n * M_solar;
	cout << "Please enter the accretion rate. (M_solar s-1)\n";
	cin >> n;
	M_dot_boundary = n * M_solar;

	cout << "Please enter the initial radius of the corona. (R_g)\n";
	cin >> n;
	R_corona_initial = n * 2 * G * M_compact / pow(c, 2);

	cout << "Please enter the maximum radius of the corona. (R_g)\n";
	cin >> n;
	R_corona_max = n * 2 * G * M_compact / pow(c, 2);

	cout << "Please enter the rise time of the corona. (days)\n";
	cin >> n;
	T_corona_rise = n * day;

	// Calculate some initial values
	R_isco = 6 * G * M_compact / pow(c, 2);
	R_outer = 3e11; // cm
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
	double* vAlpha = new double[N_grids];
	double* vP_rad = new double[N_grids];
	double* vP_gas = new double[N_grids];
	int*    vShadow = new int[N_grids];
	int trunc_radius_index = 0;
	double dT = 0;
	double vJ_total = 0;
	double vM_total = 0;
	vM_dot[N_grids - 1] = M_dot_boundary;
	vM_dot[N_grids - 2] = M_dot_boundary;
	double nu_irr = 0.1, epsilon_irr = 0.9, r_irr = 6 * G * M_compact / pow(c, 2), theta_irr = 0;

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

	parallel_for(0, N_grids, [=, &vJ_total, &vM_total](int i)
	{  
		vR[i] = vX[i]*vX[i];
		vE[i] = E_0 * exp(-1 * (pow((vR[i] - Mu), 2.)) / (2. * pow(Sigma, 2.)));
		vS[i] = vX[i]*vE[i];
		vO[i] = thompson / m_p;
		vAlpha[i] = 0.033;
		vV[i] = VISC_C(0.033, 0, vO[i]) * pow(vX[i], (4. / 3.)) * pow(vS[i], (2. / 3.));
		vJ_total += 4 * PI * sqrt(G * M_compact) * vX[i] * vS[i] * delta_X;
		vM_total += 4 * PI * vS[i] * vX[i] * vX[i] * delta_X;
		vT_irr[i] = 0;
		vH[i] = 0;
		vdelta_T[i] = 0;
		vT_eff[i] = 0;
		vT_c[i] = 0;
		vShadow[i] = 0;
	});
	// 10 iterations are enough
	for (int j = 0; j < 10; j++)
	{
		parallel_for(0, N_grids - 2, [=](int i)
		{
			vM_dot[i] = 3. * PI * (vV[i + 1] * vS[i + 1] - vV[i] * vS[i]) / delta_X;
			if (vM_dot[0] >= 0)
				vT_eff[i + 1] = EffectiveTemperature_dubus2014(vE[i + 1], vR[i + 1], vV[i + 1]);
			else
				vT_eff[i] = 0;
			if (vM_dot[i] >= 0)
				vT_c[i] = CentralTemperature_dubus2014(OpticalThickness(vE[i], vO[i]), vT_eff[i], vT_irr[i]);
			else
				vT_c[i] = 0;
			vH[i] = sqrt((vT_c[i] * k * pow(vR[i], 3)) / (mu_p(vT_c[i]) * m_p * G * M_compact));
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
	}
	L_instant = (vM_dot[0] * G * M_compact) / (2 * R_isco);		// Luminosity in ergs/s

	parallel_for(0, N_grids - 2, [=](int i)
	{
		vT_sur[i + 1] = pow(pow(vT_irr[i + 1], 4) + pow(vT_eff[i + 1], 4), 0.25);
	});

	// Add efficiency constant
	L_instant = 0.1 * (vM_dot[0] * G * M_compact) / (2 * R_isco);		// Luminosity in ergs/s

	ExtractSpectrum(vT_eff, vX, delta_X, N_grids - 2, 1, 15000, 1, false);

	alpha_compare(vAlpha, vT_c, vT_irr, vR);

	RadiationPressure(vP_rad, vT_c);
	GasPressure(vP_gas, vE, vH, vT_c);

	GetShadows(vShadow, vR, vH, R_point_source, Theta_point_source);

	WriteGraphData(vR, vShadow, N_grids, "ShadowvsR.txt", false);
	WriteGraphData(vR, vP_gas, N_grids, "PgasvsR.txt", false);
	WriteGraphData(vR, vP_rad, N_grids, "PradvsR.txt", false);
	WriteGraphData(vR, vAlpha, N_grids, "AvsR.txt", false);
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
	cout << "Total mass is " << vM_total << " g.\n";
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
	int N_L_sample = 10000;
	double* vLT_sample = new double[N_L_sample];
	double deltaLT_sample = log(T_max) / (double)N_L_sample;
	for (int i = 1; i <= N_L_sample; i++)
	{
		vLT_sample[i - 1] = i*(T_max/N_L_sample);
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
	
	ofstream file_R_hot;
	file_R_hot.open("Rhot_T.txt", ios::out);
	file_R_hot << T << "\t" << vR[FindHotIndex(vAlpha)] << "\n";
	file_R_hot.close();

	L_previous = 0.1;
	T = 0;
	bool message = true;
	while (T < T_max)
	{
		// Determine outer boundary condition*************************************************************
		//************************************************************************************************
		vS[N_grids - 1] = pow(((vM_dot[N_grids - 2] * delta_X / (3. * PI) + vS[N_grids - 2] * vV[N_grids - 2]) 
			/ (VISC_C(vAlpha[N_grids - 1], vT_c[N_grids - 1], vO[N_grids - 1]) * pow(vX[N_grids - 1], (4. / 3.)))),
			(3. / 5.));
		
		//************************************************************************************************
		//************************************************************************************************

		// Determine inner boundary condition with truncation*********************************************
		//************************************************************************************************
		/*trunc_radius_index = FindTruncationIndex(vE, vH, vT_c, vAlpha);
		parallel_for(0, trunc_radius_index, [=](int i)
		{
			vT_irr[i] = 0;
			vT_c[i] = 0;
			vV[i] = 0;
			vE[i] = 0;
			vH[i] = 0;
			vS[i] = 0;
		});*/
		//************************************************************************************************
		//************************************************************************************************

		// Determine time step****************************************************************************
		//************************************************************************************************
		parallel_for(trunc_radius_index, N_grids, [=] (int i)
		{
			vdelta_T[i] = (1. / 6.) * (pow(delta_X * vX[i], 2) / vV[i]) / 0.75;
		});
		dT = *min_element(vdelta_T + trunc_radius_index, vdelta_T + N_grids);
		//************************************************************************************************
		//************************************************************************************************

		// Iterate in radial coordinate ******************************************************************
		//************************************************************************************************
		parallel_for(1, N_grids - 1, [=](int i)
		{

			vS_new[i] = vS[i] + 0.75 * dT / pow((vX[i] * delta_X), 2) *
				(vV[i + 1] * vS[i + 1] + vV[i - 1] * vS[i - 1] - 2. * vV[i] * vS[i]);
			vV_new[i] = VISC_C(vAlpha[i], vT_c[i], vO[i]) * pow(vX[i], (4. / 3.)) * pow(vS_new[i], (2. / 3.));
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

		// Calculate alpha values*******************
		alpha_compare(vAlpha, vE, vT_irr, vR);
		//******************************************

		parallel_for(0, N_grids - 2, [=](int i)
		{
			vM_dot[i] = 3. * PI * (vV[i + 1] * vS[i + 1] - vV[i] * vS[i]) / delta_X;
		});
		parallel_for(0, N_grids, [=](int i)
		{
			vE[i] = vS[i] / vX[i];
		});
		/*parallel_for(0, N_grids, [=](int i)
		{
		vT_eff[i] = EffectiveTemperature_dubus2014(vE[i], vR[i], vV[i]);
		vT_c[i] = CentralTemperature_dubus2014(OpticalThickness(vE[i], vO[i]), vT_eff[i], vT_irr[i]);
		vH[i] = sqrt((vT_c[i] * k * pow(vR[i], 3)) / (mu_p(vT_c[i]) * m_p * G * M_compact));
		});
		ScaleHeight_FrankKingRaine(vH, vE, vT_c, vM_dot, vR);*/
		parallel_for(0, N_grids, [=](int i)
		{
			vT_eff[i] = EffectiveTemperature_dubus2014(vE[i], vR[i], vV[i]);
			vT_c[i] = CentralTemperature_dubus2014(OpticalThickness(vE[i], vO[i]), vT_eff[i], vT_irr[i]);
			vH[i] = sqrt((vT_c[i] * k * pow(vR[i], 3)) / (mu_p(vT_c[i]) * m_p * G * M_compact));


			// Find Opacity value*****************************************************
			//************************************************************************
			int a = 0, b = 0;
			for (int m = 0; m < 19; m++)
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

		L_instant = 0.1 * (vM_dot[0] * G * M_compact) / (2 * R_isco);		// Luminosity in ergs/s

		if (T > 10 * day && L_instant < 0.01 * L_edd && !CoronaFormed)
		{
			if (message)
			{
				T_corona = T;
				cout << "Corona has formed at time T = " << T_corona / day << " days.\n" << elapsed.count() << " ms have elapsed.\n\n";
				message = false;
				CoronaFormed = true;
			}
		}
		if (CoronaFormed && T > T_corona + 2 * T_corona_rise)
		{
			cout << "Corona has dissolved T = " << T / day << " days.\n" << elapsed.count() << " ms have elapsed.\n\n";
			CoronaFormed = false;
		}
		if (CoronaFormed)
		{
			CoronaIrradiationTemperature(vT_irr, TimeDependentCoronaRadius(T, T_corona, T_corona_rise, 1, 20), epsilon_irr, L_instant, vR, vH,
				TimeDependentCoronaRadius(T, T_corona, T_corona_rise, R_corona_initial, R_corona_max)); // Dubus et. al. (2014)
		}
		else
		{
			CoronaIrradiationTemperature(vT_irr, 1, epsilon_irr, L_instant, vR, vH, R_corona_initial); // Dubus et. al. (2014)
		}
	
		// Calculate alpha values*******************
		alpha_compare(vAlpha, vE, vT_irr, vR);
		//******************************************

		parallel_for(0, N_grids - 2, [=](int i)
		{
			vM_dot[i] = 3. * PI * (vV[i + 1] * vS[i + 1] - vV[i] * vS[i]) / delta_X;
		});
		parallel_for(0, N_grids, [=](int i)
		{
			vE[i] = vS[i] / vX[i];
		});
		/*parallel_for(0, N_grids, [=](int i)
		{
		vT_eff[i] = EffectiveTemperature_dubus2014(vE[i], vR[i], vV[i]);
		vT_c[i] = CentralTemperature_dubus2014(OpticalThickness(vE[i], vO[i]), vT_eff[i], vT_irr[i]);
		vH[i] = sqrt((vT_c[i] * k * pow(vR[i], 3)) / (mu_p(vT_c[i]) * m_p * G * M_compact));
		});
		ScaleHeight_FrankKingRaine(vH, vE, vT_c, vM_dot, vR);*/
		parallel_for(0, N_grids, [=](int i)
		{
			vT_eff[i] = EffectiveTemperature_dubus2014(vE[i], vR[i], vV[i]);
			vT_c[i] = CentralTemperature_dubus2014(OpticalThickness(vE[i], vO[i]), vT_eff[i], vT_irr[i]);
			vH[i] = sqrt((vT_c[i] * k * pow(vR[i], 3)) / (mu_p(vT_c[i]) * m_p * G * M_compact));


			// Find Opacity value*****************************************************
			//************************************************************************
			int a = 0, b = 0;
			for (int m = 0; m < 19; m++)
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

		T += dT; // Increase time

		if (l < N_L_sample)
		{
			if (T >= vLT_sample[l])
			{
				if (L_instant < L_previous)
				{
					cout << std::scientific;
					cout << "Maximum luminosity reached -> L = " << L_instant << " erg/s at time T = " << T/day << " days.\n" << elapsed.count() << " ms have elapsed.\n\n";
					T_L = T;
					L_previous = 0;
				}
				else if (L_previous != 0)
				{
					L_previous = L_instant;
				}
				l++;
				file.open("lightcurve.txt", ios::app);
				if (L_instant > 1e25)
					file << T / day << "\t" << L_instant << "\n";						// Write luminosity to file
				file.close();
				parallel_for(0, N_grids, [=](int i)
				{
					vT_sur[i] = pow(pow(vT_irr[i], 4) + pow(vT_eff[i], 4), 0.25);
				});
				L_BB = BB_Luminosity(vT_sur, vX);
				file_bolo.open("lightcurve_bolo.txt", ios::app);
				if (L_BB > 1e25)
					file_bolo << T / day << "\t" << L_BB << "\n";						// Write luminosity to file
				file_bolo.close();

				L_optical = GetLuminosity(vT_sur, vX, delta_X, N_grids, 1, 4, 0.001);
				file_optical.open("lightcurve_optical.txt", ios::app);
				//if (L_optical > 1e25)
					file_optical << T / day << "\t" << L_optical << "\n";						// Write luminosity to file
				file_optical.close();

				L_X = GetLuminosity(vT_sur, vX, delta_X, N_grids, 1000, 200000, 100);
				file_X.open("lightcurve_X.txt", ios::app);
				//if (L_X > 1e25)
					file_X << T / day << "\t" << L_X << "\n";						// Write luminosity to file
				file_X.close();

				file_R_hot.open("Rhot_T.txt", ios::app);
				file_R_hot << T / day << "\t" << vR[FindHotIndex(vAlpha)] << "\n";
				file_R_hot.close();
			}
		}
		// Take samples ***********************************************************************************
		//*************************************************************************************************
		if (T >= vT_sample[j])
		{
			j++;

			ExtractSpectrum(vT_eff, vX, delta_X, N_grids - 2, 1, 15000, 1, true);
			
			RadiationPressure(vP_rad, vT_c);
			GasPressure(vP_gas, vE, vH, vT_c);

			GetShadows(vShadow, vR, vH, R_point_source, Theta_point_source);

			WriteGraphData(vR, vShadow, N_grids, "ShadowvsR.txt", true);
			WriteGraphData(vR, vP_gas, N_grids, "PgasvsR.txt", true);
			WriteGraphData(vR, vP_rad, N_grids, "PradvsR.txt", true);
			WriteGraphData(vR, vAlpha, N_grids, "AvsR.txt", true);
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
			cout << "Truncation radius is      " << vR[trunc_radius_index] << " cm.\n";
			cout << "Luminosity is             " << (L_instant / L_edd) * 100 << "\% of Eddington Lumnosity.\n";
			cout << T/T_max * 100 << fixed << " percent completed! " << elapsed.count() << " ms have elapsed.\n\n";
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
			double C = 5e-3;/* nu * (1 - epsilon) * (2. / 35.)*(H[i + 1] / R[i + 1]); // lasota (2014)*/
			if (C > 0 && L > 0)
				T_irr[i + 1] = pow(C * L / (4 * PI * a * R[i + 1] * R[i + 1]), 0.25);
			else
				T_irr[i + 1] = 0;
		}
	});
}

void CoronaIrradiationTemperature(double* T_irr, double nu, double epsilon, double L, double* R, double* H, double R_corona, double Theta_corona)
{
	int* Shadows = new int[N_grids];
	GetShadows(Shadows, R, H, R_corona, Theta_corona);
	parallel_for(0, N_grids - 1, [=](int i)
	{
		if (Shadows[i + 1] != 1)
		{
			//double C = nu * (1 - epsilon)*((H[i + 1] - H[i]) / (R[i + 1] - R[i]) - H[i + 1] / R[i + 1]);
			double C = 5e-3;/* nu * (1 - epsilon) * (2. / 35.)*(H[i + 1] / R[i + 1]); // lasota (2014)*/
			double R_2 = pow(H[i + 1], 2) + pow(R[i + 1], 2);
			double Theta = atan(H[i + 1] / R[i + 1]);
			double d_2 = pow(R_corona, 2) + R_2 + 2 * R_corona * sqrt(R_2) * cos(Theta_corona - Theta);
			if (C > 0 && L > 0)
				T_irr[i + 1] = pow(C * L / (4 * PI * a * d_2), 0.25);
			else
				T_irr[i + 1] = 0;
		}
	});
}

// Integral
void CoronaIrradiationTemperature(double* T_irr, double nu, double epsilon, double L, double* R, double* H, double R_corona)
{
	bool Shadow = false;
	int N_corona = 10;
	double* R_int = new double[N_corona];
	double* Theta_int = new double[N_corona];
	double Delta_R_corona = R_corona / (N_corona - 1);
	double Delta_Theta_corona = PI / (2 * (N_corona - 1));
	parallel_for(0, N_corona, [=](int i)
	{
		R_int[i] = 0 + i * R_corona / Delta_R_corona;
		Theta_int[i] = 0 + i * (PI / 2) / Delta_Theta_corona;
	});
	parallel_for(0, N_grids, [=](int i)
	{
		T_irr[i] = 0;
	});
	L = L / (R_corona * R_corona * PI / 4);
	parallel_for(0, N_grids - 1, [=, &Shadow](int i)
	{
		parallel_for(0, N_corona, [=, &Shadow](int m)
		{
			parallel_for(0, N_corona, [=, &Shadow](int n)
			{
				Shadow = false;
				parallel_for(0, i + 1, [=, &Shadow](int j)
				{
					if (atan((H[i + 1] - R_int[m] * sin(Theta_int[n])) / abs(R[i + 1] - R_int[m] * cos(Theta_int[n]))) <
						atan((H[j] - R_int[m] * sin(Theta_int[n])) / abs(R[j] - R_int[m] * cos(Theta_int[n]))))
					{
						Shadow = true;
					}
				});

				if (!Shadow)
				{
					//double C = nu * (1 - epsilon)*((H[i + 1] - H[i]) / (R[i + 1] - R[i]) - H[i + 1] / R[i + 1]);
					double C = nu * 5e-3;/* nu * (1 - epsilon) * (2. / 35.)*(H[i + 1] / R[i + 1]); // lasota (2014)*/
					double R_2 = pow(H[i + 1], 2) + pow(R[i + 1], 2);
					double Theta = atan(H[i + 1] / R[i + 1]);
					double d_2 = pow(R_int[m], 2) + R_2 + 2 * R_int[m] * sqrt(R_2) * cos(Theta_int[n] - Theta);
					if (C > 0 && L > 0)
						T_irr[i + 1] += pow(C * L / (4 * PI * a * d_2) * R_int[m] * Delta_R_corona * Delta_Theta_corona, 0.25);
					else
						T_irr[i + 1] += 0;

				}
			});
		});
	});
}

void GetShadows(int* Shadows, double* R, double* H, double R_corona, double Theta_corona)
{
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

double mu_p(double T)
{
	if (T > 10000)
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
double VISC_C(double alpha, double T_c, double opacity)
{
	double VIS_C = pow(
		(27. / 32.) * (pow(alpha, 4) * pow(k, 4) * opacity)
		/ (pow(mu_p(T_c), 4) * pow(m_p, 4) * a * G * M_compact)
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
		P_rad[i] = 1. / 3. * rad_const * pow(T_c[i], 4);
	});
}

void GasPressure(double* P_gas, double* E, double* H, double* T_c)
{
	parallel_for(0, N_grids, [=](int i)
	{
		P_gas[i] = ((E[i] / H[i]) * k * T_c[i]) / (mu_p(T_c[i]) * m_p);
	});
}

void ScaleHeight_FrankKingRaine(double * H, double * E, double * T_c, double* M_dot, double* R)
{
	double* P_rad = new double[N_grids];
	double* P_gas = new double[N_grids];
	RadiationPressure(P_rad, T_c);
	GasPressure(P_gas, E, H, T_c);

	for (int j = 0; j < 10; j++)
	{
		parallel_for(0, N_grids, [=](int i)
		{
			if (P_gas[i] > P_rad[i])
			{
				H[i] = ((3 * R_isco) / (4 * 0.1)) * (M_dot[i] / M_dot_edd) * (1 - pow(R_isco / R[i], 0.5));
			}
			else
			{
				H[i] = sqrt((T_c[i + 1] * k * pow(R[i + 1], 3)) / (mu_p(T_c[i + 1]) * m_p * G * M_compact));
			}
		});
		RadiationPressure(P_rad, T_c);
		GasPressure(P_gas, E, H, T_c);
	}
}

int FindTruncationIndex(double* E, double* H, double* T_c, double* Alpha)
{
	double* P_rad = new double[N_grids];
	double* P_gas = new double[N_grids];

	RadiationPressure(P_rad, T_c);
	GasPressure(P_gas, E, H, T_c);

	for (int i = 0; i < N_grids - 1; i++)
	{
		if (P_rad[i] < P_gas[i] + 0.1 && P_rad[i + 1] < P_gas[i + 1])
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

double TimeDependentCoronaRadius(double T_current, double T_corona, double T_rise, double R_initial, double R_max)
{
	double c_1 = (R_max - R_initial) / (pow(T_rise, 2) - 4 * T_corona * T_rise);
	double c_2 = 2 * (T_corona + T_rise) * (R_max - R_initial) / (4 * T_corona * T_rise - pow(T_rise, 2));
	double c_3 = R_initial - pow(T_corona, 2) * (R_max - R_initial) / (4 * T_corona * T_rise - pow(T_rise, 2));

	return c_1 * pow(T_current, 2) + c_2 * T_current + c_3;
}


