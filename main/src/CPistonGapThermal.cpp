#include "CPistonGap.h"
#include "..\..\caspar_input\input.h"
#include "logger.h"
#pragma once


extern class CGapInput myGapInput;
extern struct sGapResult myGapResult;
extern class CMesh myMesh;
extern class CThermal myThermal;
extern class CFEMThermal myFEMThermal;
extern class CGapUtils myGapUtils;
extern class input myinput;


//Calculate istantaneous thermal flux generated from gap to piston and cylinder
void CPistonGap::PistonCylinderCalcGapThermalFlux(void)
{	
	double lambda;
	Range all = Range::all();


	//Heat fluxes clear
	QgapK = 0.0;
	QgapB = 0.0;

	//Oil conduction coefficient
	lambda = my_oil ->get_lambda();//oilpistongap.oillambda;

	//CONDUCTIVE HEAT FLUX TO SOLID PARTS FROM FLUID [W/m2]
	Array<double,1> To(N*M);
	//Cylinder
	To = T(Range(0,N*M-1));
	QgapB = lambda/dz2 * (To - TB_surf_gap);
	//Piston
	To = T(Range(N*M*(Q-1),N*M*Q-1));
	QgapK = lambda/dz2 * (To - TK_surf_gap);


};
//Calculate thermal flux to solid surfaces considering boundary flux over shaft revolution
void CPistonGap::PistonCylinderCalcBodyThermalFlux(void)
{
	double phi_deg,phi_rad,speedK,timeold,TDC,TCase,
		AlphaCase,AlphaDC,dt,lvar,lA,lch,lB;
	int phioffset;
	Range all = Range::all();
	Array<double,1> Ebody;
	Array<double,1> Egap;
	Array<double,1> qDC;
	Array<double,1> qCase;
	

	//Parameters
	speedK = operatingpistongap.speedK;
	phi_rad = operatingpistongap.phi_rad;
	phi_deg = operatingpistongap.phi_deg;
	phioffset = (int) (floor(phi_rad/dphi)*speedK);
	timeold = myGapResult.time;


	//Geometry
	lvar = geometrypistongap.lvar;
	lA = geometrypistongap.lA;
	lB = geometrypistongap.lB;
	lch = geometrypistongap.lch;


	//Temperatures
	if(phi_deg < 180.0)
	{
		TDC = temperaturepistongap.THP;
	}
	else
	{ 
		TDC = temperaturepistongap.TLP;
	}
	TCase = temperaturepistongap.TCase;	
	//Covection Case side
	AlphaDC = oilpistongap.AlphaDC;								
	//Covection DC side
	AlphaCase = oilpistongap.AlphaCase;
	//Delta Time
	dt = timenew-timeold;


	//----------------------PISTON INTERPOLATION-------------------------//
	Egap.resize(QgapK.extent(0));					Egap = 0.0;
	Ebody.resize(myGapInput.xyzfK_th.extent(0));	Ebody = 0.0;
	qDC.resize(myGapInput.xyzfK_th.extent(0));		qDC = 0.0;
	qCase.resize(myGapInput.xyzfK_th.extent(0));	qCase = 0.0;
	//Gap istantaneous energy [J/m2] = [W/m2] * [s]
	Egap = QgapK * dt;
	//Energy flux rotation due to piston relative motion in gap
	if(phioffset>0)
	{
		Array<double,1> EgapCopy;
		EgapCopy.resize(N*M);
		EgapCopy = Egap;
		Egap(Range(phioffset*M,N*M-1)) = EgapCopy(Range(0,(N-phioffset)*M-1));
		Egap(Range(0,phioffset*M-1)) = EgapCopy(Range((N-phioffset)*M,N*M-1));
	};
	//Energy flux to solid body interpolation
	Ebody = myGapUtils.InterpolateFields(FaceIdK_f2s_th,FaceDistK_f2s_th,Egap);
	//Boundary heat fluxes [W/m2]
	qDC = AlphaDC * (TDC - TK_surf);
	qCase = AlphaCase * (TCase - TK_surf);
	//Boundary energy flux [J/m2]
	qDC *= dt;	qCase *= dt;
	//Boundaries
	Ebody = where( zfK_th<(lch+lA), qDC , Ebody );
	Ebody = where( zfK_th>(lch+lA+lvar), qCase , Ebody );
	//Final
	EbodyK += Ebody;


	//----------------------CYLINDER INTERPOLATION-------------------------//
	Egap.resize(QgapB.extent(0));					Egap = 0.0;
	Ebody.resize(myGapInput.xyzfB_th.extent(0));	Ebody = 0.0;
	qDC.resize(myGapInput.xyzfB_th.extent(0));		qDC = 0.0;
	qCase.resize(myGapInput.xyzfB_th.extent(0));	qCase = 0.0;
	//Gap istantaneous energy [J/m2] = [W/m2] * [s]
	Egap = QgapB * dt;
	//Energy flux to solid body interpolation
	Ebody = myGapUtils.InterpolateFields(FaceIdB_f2s_th,FaceDistB_f2s_th,Egap);
	//Boundary heat fluxes [W/m2]
	qDC = AlphaDC * (TDC - TB_surf);
	qCase = AlphaCase * (TCase - TB_surf);
	//Boundary energy flux [J/m2]
	qDC *= dt;	qCase *= dt;
	//Boundaries
	Ebody = where( zfB_th<lB, qDC , Ebody );
	Ebody = where( zfB_th>(lB+lvar), qCase , Ebody );
	//Final
	EbodyB += Ebody;

}

//Calculate piston and cylinder body temperatures and thermal deformation
void CPistonGap::PistonCylinderSolveBodyThermal(void)
{

	/*Log << "EbodyK = " << "\n";
	for( int i = 0; i < EbodyK.size();i++){
		Log << EbodyK(i) << "\n";
	}

		Log << "EbodyB = " << "\n";
	for( int i = 0; i < EbodyB.size();i++){
		Log << EbodyB(i) << "\n";
	}*/

	//Dump IM's dwm
	/*myGapInput.xyzfK_p.free(); myGapInput.xyznK_p.free(); myGapInput.xyzfB_p.free();
	myGapInput.xyznB_p.free();*/ myGapInput.IM_piston.free(); myGapInput.IM_cylinder.free();

	//shaft speed [rev/s]
	double speed = operatingpistongap.speed/60.0;

	//---------------Calculate temperature piston------------//
	qbi_piston[0] = 0.0;
	//[J/(m2 rev) * rev/s] = [W/m2] 
	EbodyK *= speed;
	//limit flux
	//EbodyK = where(EbodyK>=1.0e5,1.0e5,EbodyK);
	EbodyK = where(EbodyK<=-1.0e3,-1.0e3,EbodyK);
	EbodyK = where(EbodyK>=5.0e4,5.0e4,EbodyK);
	//EbodyK = where(EbodyK<=-5.0e4,-5.0e4,EbodyK);
	//damp flux
	EbodyK = EbodyK_old + AlphaTh * (EbodyK - EbodyK_old);
	//assign flux to flux vector - from faces heat flux to surface nodes
	qbi_piston[0] = EbodyK;
	/*Log << "Boundary EbodyK " << "\n";
	for(int b = 0; b < EbodyK.size(); b++){
		Log << EbodyK(b) << "\n";
	}*/
	for(int i = 0;i<myinput.data.thermal.piston.neumann_bc.size();i++){
		Array<double,1> temp ;
		temp.resize(myMesh.nCells,myinput.data.thermal.piston.neumann_bc[i].q[0]);
		temp[0]= myinput.data.thermal.piston.neumann_bc[i].q[0];
		qbi_piston.push_back(temp);
	};
	qbi_piston[0] = where(qbi_piston[0]==0.0,1.0e-9,qbi_piston[0]);
	//read mesh and solve thermal
	myThermal.readMeshThermal("piston");
	for(int a = 0; a < qbi_piston.size(); a++){
		/*Log << "Boundary qbi_piston " << a << "\n";
		for(int b = 0; b < qbi_piston[a].size(); b++){
			Log << qbi_piston[a](b) << "\n";
		}*/
	}
	myThermal.ThermalSolve("piston",qbi_piston,TK_body,TK_surf);
	//assign flux to old array
	EbodyK_old = EbodyK;
	//reset flux
	EbodyK = 0.0;

	//---------------Solve piston thermal deformation------------//
	if(ThermalDeformation)
	{
		myFEMThermal.FEMThermalSolve("piston",TK_body,defK_th,myMesh.piston_name);
	}
	//free variables
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.cxyz_th.free(); myMesh.E.free();
	myMesh.v.free(); myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear();



	//---------------Calculate temperature cylinder------------//
	qbi_cylinder[0] = 0.0;
	//[J/(m2 rev) * rev/s] = [W/m2] 
	EbodyB *= speed;
	//limit flux
	//EbodyB = where(EbodyB>=1.0e5,1.0e5,EbodyB);
	EbodyB = where(EbodyB<=-1.0e4,-1.0e4,EbodyB);
	EbodyB = where(EbodyB>=5.0e4,5.0e4,EbodyB);
	//EbodyB = where(EbodyB<=-5.0e4,-5.0e4,EbodyB);
	//damp flux
	EbodyB = EbodyB_old + AlphaTh * (EbodyB - EbodyB_old);
	//assign flux to flux vector
	qbi_cylinder[0] = EbodyB;
	/*Log << "Boundary EbodyB " << "\n";
	for(int b = 0; b < EbodyB.size(); b++){
		Log << EbodyB(b) << "\n";
	}*/
	qbi_cylinder[0] = where(qbi_cylinder[0]==0.0,1.0e-9,qbi_cylinder[0]);
	//read mesh
	myThermal.readMeshThermal("cylinder");
	//calculate other bores heat flux
	if(myinput.data.options_piston.general.EHDTestRig || myinput.data.options_piston.general.TriboTestRig)
		qbi_cylinder.resize(1);
	else
		qbi_cylinder.resize(myinput.data.operating_conditions.npistons);
	PistonCylinderBlockFlux(qbi_cylinder);
	//add other defined flux boundaries
	for(int i = 0;i<myinput.data.thermal.block.neumann_bc.size();i++){
		Array<double,1> temp; 
		temp.resize(myMesh.nNodes);//should be more than plenty...
		temp = myinput.data.thermal.block.neumann_bc[i].q[0];
		//Log << "defining qb " << i << "\n";
		qbi_cylinder.push_back(temp);
	};
	//solve thermal
	for(int a = 0; a < qbi_cylinder.size(); a++){
		//Log << "Boundary qbi_cylinder " << a << "\n";
		for(int b = 0; b < qbi_cylinder[a].size(); b++){
			//Log << qbi_cylinder[a](b) << "\n";
		}
	}
	myThermal.ThermalSolve("cylinder",qbi_cylinder,TB_body,TB_surf);
	//assign flux to old array
	EbodyB_old = EbodyB;
	//reset flux
	EbodyB = 0.0;
	
	//---------------Solve cylinder thermal deformation------------//
	if(ThermalDeformation)
	{
		myFEMThermal.FEMThermalSolve("cylinder",TB_body,defB_th,myMesh.cylinder_name);
	}
	//free variables
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.cxyz_th.free(); myMesh.E.free();
	myMesh.v.free(); myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear(); 

	//Load IM's dwm
	if(myinput.data.options_piston.general.PressureDeformation)
	{
		//read surface coordinates
		//Log << "Reading pressure mesh coordinates file... " << "\t";
		//myGapInput.readBodySurfacexyzPressure();
		//Log << "done!" << "\n";
		//Log << " " << "\n";
		//read matrices
		Log << "Reading influence matrices... " << "\t";
		myGapInput.readInfluenceMatricesPistonCylinder(myinput.data.options_piston.general.IM_piston_path,myinput.data.options_piston.general.IM_bushing_path);
		Log << "done!" << "\n";
		Log << " " << "\n";
	};

};
//Interpolate bore flux to other cylinder bores
void CPistonGap::PistonCylinderBlockFlux(vector<Array<double,1>> &qbi)
{

	Range all = Range::all();

	//Assign cylinder gap coordinates for reference bore
	int nFaces_gap = (int) myThermal.faceid_qb[0].size();
	Array<double,2> xyz_0;	xyz_0.resize(nFaces_gap,3);	xyz_0 = 0.0;
	for(int i=0;i<nFaces_gap;i++)
	{
		xyz_0(i,0) = myThermal.xyzf(myThermal.faceid_qb[0][i],0);
		xyz_0(i,1) = myThermal.xyzf(myThermal.faceid_qb[0][i],1);
		xyz_0(i,2) = myThermal.xyzf(myThermal.faceid_qb[0][i],2);
	}
	//Array of coordinates rotation of reference bore to other bores
	Array<double,2> xyz_r;	xyz_r.resize(nFaces_gap,3);	xyz_r = 0.0;
	//Array of coordinates nodes actual i-th bore
	Array<double,2> xyz_i;
	//Angle between cylinder bores
	double dtheta = 2.0*PI/(double) qbi.size();
	double theta = dtheta;
	//Interpolate reference heat flux to other bores coordinates
	int nqb = (int) qbi.size();
	//Log << "nqb=" << nqb << "\n";
	for(int i=0;i<nqb-1;i++)
	{
		//Rotate reference bore according to actual i-th bore position
		xyz_r(all,0) = xyz_0(all,0)*cos(theta) - 1.0 * xyz_0(all,1)*sin(theta);
		xyz_r(all,1) = xyz_0(all,0)*sin(theta) + xyz_0(all,1)*cos(theta);
		xyz_r(all,2) = xyz_0(all,2);

		//Increment theta
		theta -= dtheta;
		//Assign actual i-th cylinder gap coordinates
		nFaces_gap = (int) myThermal.faceid_qb[i+1].size();
		xyz_i.resize(nFaces_gap,3);	xyz_i = 0.0;
		for(int j=0;j<nFaces_gap;j++)
		{
			xyz_i(j,0) = myThermal.xyzf(myThermal.faceid_qb[i+1][j],0);
			xyz_i(j,1) = myThermal.xyzf(myThermal.faceid_qb[i+1][j],1);
			xyz_i(j,2) = myThermal.xyzf(myThermal.faceid_qb[i+1][j],2);
		}
		
		//Search neighbours nodes from reference coordinates to actual i-th cylinder
		Array<double,2> FaceDist; Array<int,2> FaceId;
		myGapUtils.SearchNeighbours(xyz_r,xyz_i,FaceId,FaceDist,10);
		//Interpolate reference heat flux to i-th heat flux
		qbi[i+1] = myGapUtils.InterpolateFields(FaceId,FaceDist,qbi[0]);
			
	}

};