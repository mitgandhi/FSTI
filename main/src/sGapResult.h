#pragma once

using namespace std;

struct sGapResult
{
	//General results
	double time,phi;
	//Revcounter
	int revcounter;
	//Pressures
	double pDC;						//Displacement chamber pressure [bar]
	double pHP;						//High pressure port pressure [bar]
	double PLP;						//Low pressure port pressure [bar]
	//Positions and velocities
	vector<double> gappositions;		//parts positions vector [m]
	vector<double> gapvelocities;	//parts velocities vector [m/s]
	vector<double> QSKTot;			//Total leakage from pistons [l/min]
	vector<double> MSKTot;			//Total torque loss from pistons [Nm]
	vector<double> MFT;				//Total torque loss due to pistons spin in chamber [Nm]
	vector<double> PhiD_mech;			//Total energy dissipated [W]
	vector<double> PhiD_vol;			//Total energy dissipated [W]


};
