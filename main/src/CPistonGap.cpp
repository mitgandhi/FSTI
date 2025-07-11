#include "CPistonGap.h"
#include "logger.h"
#include <iostream>
#include <iomanip>
#include "../../caspar_input\input.h"
#pragma once


extern class CGapInput myGapInput;
extern class CGapUtils myGapUtils;
extern struct sGapResult myGapResult;
extern class CMesh myMesh;
extern class CThermal myThermal;
extern class CFEMThermal myFEMThermal;
extern class input myinput;
extern bool force_balance_iterative;
extern double mixed_friction;

CPistonGap::CPistonGap()
{
	//Initialize Processor Load
	GetSystemTimes(&IdleTime, &KernelTime, &UserTime);

	//Initialize counters and angles
	myGapResult.phi = 0.0 ;
	myGapResult.revcounter = 0 ;


	//Initialize variables
	timenew = 0.0;
	lvarold = 0.0;
	PhiD_mech = 0.0;
	PhiD_vol = 0.0;


	//Initialize leakages
	operatingpistongap.QSK = 0.0;
	operatingpistongap.QSK_p = 0.0;
	operatingpistongap.QSK_c = 0.0;


	//Initialize forces
	forcespistongap.FTKx = 0.0;		forcespistongap.FTKy = 0.0;
	forcespistongap.FTKx_p = 0.0;	forcespistongap.FTKy_p = 0.0;
	forcespistongap.FTKx_c = 0.0;	forcespistongap.FTKy_c = 0.0;
	forcespistongap.FTbx = 0.0;		forcespistongap.FTby = 0.0;
	forcespistongap.FTbx_p = 0.0;	forcespistongap.FTby_p = 0.0;
	forcespistongap.FTbx_c = 0.0;	forcespistongap.FTby_c = 0.0;
	

	//Parameters
	R_h = 1.0;
	R_p = 1.0;
	Rold_h = 1.0;
	Rold_p = 1.0;
	nmax = myinput.data.options_piston.numeric.nmax;
	Rmin_p = myinput.data.options_piston.numeric.Rmin_p;
	Rmin_h = myinput.data.options_piston.numeric.Rmin_h;
	Rmin_R = myinput.data.options_piston.numeric.Rmin_R;
	Rmin_E = myinput.data.options_piston.numeric.Rmin_E;
	epsilonK = myinput.data.options_piston.numeric.epsilonK;
	AlphaP = myinput.data.options_piston.numeric.AlphaP;
	AlphaMu = myinput.data.options_piston.numeric.AlphaMu;
	AlphaDef = myinput.data.options_piston.numeric.AlphaDef;
	AlphaTh = myinput.data.options_piston.numeric.AlphaTh;
	EmodK = myinput.data.options_piston.numeric.EmodK;
	EmodB = myinput.data.options_piston.numeric.EmodB;
	vK = myinput.data.options_piston.numeric.vK;
	vB = myinput.data.options_piston.numeric.vB;
	Eprime = 2.0 / ( ( 1.0 - pow(vK,2.0) )/EmodK + ( 1.0 - pow(vB,2.0) )/EmodB ) ; 


	//Copy piston simulation options
	EHDTestRig = myinput.data.options_piston.general.EHDTestRig;
	TriboTestRig = myinput.data.options_piston.general.TriboTestRig;
	PistonMacroGeometry = myinput.data.options_piston.general.McrK;
	CylinderMacroGeometry = myinput.data.options_piston.general.McrB;
	PressureDeformation = myinput.data.options_piston.general.PressureDeformation;
	PressureDeformationOMP = myinput.data.options_piston.general.PressureDeformationOMP;
	ThermalDeformation = myinput.data.options_piston.general.ThermalDeformation;
	ReynoldsMultiGrid = myinput.data.options_piston.general.ReynoldsMultiGrid;
	ReadpFile = myinput.data.options_piston.general.ReadpFile;
	EnergyEquation = myinput.data.options_piston.general.EnergyEquation;
	HeatTransfer = myinput.data.options_piston.general.HeatTransfer;


	//Gap computing grid structure
	gridpistongap.N = myinput.data.options_piston.fluid_grid.GS.N;
	gridpistongap.M = myinput.data.options_piston.fluid_grid.GS.M;
	gridpistongap.Q = myinput.data.options_piston.fluid_grid.GS.Q;
	N = gridpistongap.N;							//Number of fluid elements x direction (circumferential direction)
	M = gridpistongap.M;							//Number of fluid elements y direction (axial direction)
	Q = gridpistongap.Q;							//Number of fluid elements z direction (film thickness)
	

	//If multigrid overwrite inputs from grids file
	if(ReynoldsMultiGrid)
	{
		gridpistongap.N = myinput.data.options_piston.fluid_grid.MG.MG_N[0];
		gridpistongap.M = myinput.data.options_piston.fluid_grid.MG.MG_M[0];
		gridpistongap.Q = myinput.data.options_piston.fluid_grid.MG.Q;
		N = gridpistongap.N;			
		M = gridpistongap.M;
		Q = gridpistongap.Q;
		//myinput.data.options_piston.fluid_grid.GS.N = N; //fixed dwm
		//myinput.data.options_piston.fluid_grid.GS.M = M;
	}


	//Geometric parameters
	geometrypistongap.dK = myinput.data.geometry.dK;
	geometrypistongap.dDK = myinput.data.geometry.dDK;
	geometrypistongap.dZ = myinput.data.geometry.dZ;
	geometrypistongap.dB = myinput.data.geometry.dB;
	geometrypistongap.le = myinput.data.geometry.le;
	geometrypistongap.lF = myinput.data.geometry.lF;
	geometrypistongap.lK = myinput.data.geometry.lK;
	geometrypistongap.lch = myinput.data.geometry.lch;
	geometrypistongap.lSK = myinput.data.geometry.lSK;
	geometrypistongap.lKG = myinput.data.geometry.lKG;
	geometrypistongap.lZ0 = myinput.data.geometry.lZ0;
	geometrypistongap.lengthB = myinput.data.geometry.lengthB;
	geometrypistongap.lengthcanalB = myinput.data.geometry.lengthcanalB;
	geometrypistongap.mK = myinput.data.geometry.mK;
	geometrypistongap.hmin = myinput.data.options_piston.numeric.hmin;
	geometrypistongap.doutG = myinput.data.geometry.doutG;
	geometrypistongap.dinG = myinput.data.geometry.dinG;
	geometrypistongap.rK_red = myinput.data.geometry.rK_red;
	geometrypistongap.lK_hs = myinput.data.geometry.lK_hs;
	geometrypistongap.npistons = myinput.data.operating_conditions.npistons;


	//Cylinder block pitch radius
	rB = geometrypistongap.dB * 0.5;
	//Bushing radius
	rZ = geometrypistongap.dZ * 0.5;
	//Piston radius
	rK = geometrypistongap.dK * 0.5;
	//Slipper surface area
	AreaG = PI / 4.0 * ( pow(geometrypistongap.doutG,2) - pow(geometrypistongap.dinG,2) );	
	//Piston displacement chamber area
	AreaK = PI / 4.0 * ( pow(geometrypistongap.dK,2) - pow(geometrypistongap.dDK,2) );	
	//Delta x in gap circumferential direction
	dx = geometrypistongap.dK * PI / N;	
	//Delta phi in gap circumference
	dphi = 2.0*PI / N;
	//Moment arm piston axis
	zKj.resize(N*M);				zKj=0;
	//Angular field gap
	phi.resize(N*M);				phi=0;
	for(int i=0;i<N;i++)
	{
		phi(Range(i*M,(i+1)*M-1)) = (i+0.5)*dphi;
	};


	//Multigrid levels initialization
	if(ReynoldsMultiGrid)
	{
		//Size fields
		int nL = myinput.data.options_piston.fluid_grid.MG.nL;
		levels.resize(nL);
		//Sizing
		int l=0;
		while(l<nL)
		{
			levels[l].VW = myinput.data.options_piston.fluid_grid.MG.VW;
			levels[l].MGInt = myinput.data.options_piston.fluid_grid.MG.MGInt;
			levels[l].v1 = myinput.data.options_piston.fluid_grid.MG.v1;
			levels[l].v2 = myinput.data.options_piston.fluid_grid.MG.v2;
			levels[l].nL = nL;
			levels[l].N = myinput.data.options_piston.fluid_grid.MG.MG_N[l];
			levels[l].M = myinput.data.options_piston.fluid_grid.MG.MG_M[l];

			levels[l].ap.resize(levels[l].N*levels[l].M);			levels[l].ap=0.0;
			levels[l].an.resize(levels[l].N*levels[l].M);			levels[l].an=0.0;
			levels[l].as.resize(levels[l].N*levels[l].M);			levels[l].as=0.0;
			levels[l].ae.resize(levels[l].N*levels[l].M);			levels[l].ae=0.0;
			levels[l].aw.resize(levels[l].N*levels[l].M);			levels[l].aw=0.0;
			levels[l].b.resize(levels[l].N*levels[l].M);			levels[l].b=0.0;
			levels[l].r.resize(levels[l].N*levels[l].M);			levels[l].r=0.0;
			levels[l].e.resize(levels[l].N*levels[l].M);			levels[l].e=0.0;
			levels[l].x.resize(levels[l].N*levels[l].M);			levels[l].x=0.0;
			levels[l].h.resize(levels[l].N*levels[l].M);			levels[l].h=0.0;
			levels[l].mu.resize(levels[l].N*levels[l].M);			levels[l].mu=0.0;
			levels[l].phi.resize(levels[l].N*levels[l].M);			levels[l].phi=0.0;	
			levels[l].xyz.resize(levels[l].N*levels[l].M,3);		levels[l].xyz=0.0;			

			l++;
		}
	};


	//Operating parameters
	operatingpistongap.phi_deg = 0.0;
	operatingpistongap.phi_rad = 0.0;
	operatingpistongap.omega = myinput.data.operating_conditions.speed;
	operatingpistongap.speed = operatingpistongap.omega*(30/PI);
	operatingpistongap.beta_rad = myinput.data.operating_conditions.beta;
	operatingpistongap.beta_deg = operatingpistongap.beta_rad*(180/PI);
	operatingpistongap.gamma_rad = myinput.data.geometry.gamma;
	operatingpistongap.gamma_deg = operatingpistongap.gamma_rad*(180/PI);
	operatingpistongap.betamax_rad = myinput.data.operating_conditions.betamax;
	operatingpistongap.betamax_deg = operatingpistongap.betamax_rad*(180/PI);
	operatingpistongap.mode = myinput.data.operating_conditions.mode;
	operatingpistongap.pHP = myinput.data.operating_conditions.HP;
	operatingpistongap.pLP = myinput.data.operating_conditions.LP;
	operatingpistongap.pCase = myinput.data.operating_conditions.pCase;
	

	//If EHD Test Rig is simulated, maximum displacmenet angle is hardcoded
	if(EHDTestRig)
	{
		operatingpistongap.beta_deg = 17.0;
		operatingpistongap.beta_rad = 17.0 * (PI/180.0);
		operatingpistongap.betamax_deg = 17.0;
		operatingpistongap.betamax_rad = 17.0 * (PI/180.0);
		geometrypistongap.lF = 28.660 * 1e-3;
	}
	operatingpistongap.speedK = myinput.data.geometry.speedK;


	//Switch beta if motoring mode
	if(operatingpistongap.mode==2)
	{
		operatingpistongap.beta_deg *= -1.0;
		operatingpistongap.beta_rad *= -1.0;
	}


	//If EHD Test Rig is simulated, no relative piston rotation is considered
	if(EHDTestRig)
	{
		operatingpistongap.speedK = 0.0;
	}
	

	//Fluid parameters
	//oilpistongap.oiltype = myinput.data.oil.general.oiltype;
	//oilpistongap.oildensity = myinput.data.oil.oildensity;
	//oilpistongap.oilviscosity = myinput.data.oil.oilviscosity;
	//oilpistongap.oilbetaP = myinput.data.oil.oilbetaP;
	//oilpistongap.oilbetaT = myinput.data.oil.oilbetaT;
	//oilpistongap.oilPc1 = myinput.data.oil.oilPc1;
	///oilpistongap.oilPc2 = myinput.data.oil.oilPc2;
	//oilpistongap.oilTc1 = myinput.data.oil.oilTc1;
	//oilpistongap.oilTc2 = myinput.data.oil.oilTc2;
	//oilpistongap.oilW = myinput.data.oil.oilW;
	//oilpistongap.alpha1 = myinput.data.oil.alpha1;
	//oilpistongap.alpha2 = myinput.data.oil.alpha2;
	//oilpistongap.alpha3 = myinput.data.oil.alpha3;
	//oilpistongap.oilC = myinput.data.oil.oilC;
	//oilpistongap.oillambda = myinput.data.oil.oillambda;
	oilpistongap.AlphaDC = myinput.data.options_piston.numeric.AlphaDC;
	oilpistongap.AlphaCase = myinput.data.options_piston.numeric.AlphaCase;
	if(myinput.data.oil.general.oiltype == 0){
		Log << "Using Constant Oil Properties" << "\n";
		my_oil = new constant_oil(myinput);
	}
	else if(myinput.data.oil.general.oiltype == 1){
		Log << "Using User Defined Oil Model" << "\n";
		my_oil=new user_defined_oil(myinput);
	}
	/*else if(myinput.data.oil.general.oiltype == 2){
		Log << "Using HLP32 Oil Model" << "\n";
		my_oil=new HLP32_oil(myinput);
	}
	else if(myinput.data.oil.general.oiltype == 3){
		Log << "Using Skydrol Oil Model" << "\n";
		my_oil=new skydrol_oil(myinput);
	}	
	else if(myinput.data.oil.general.oiltype == 4){
		Log << "Using Red Oil Model" << "\n";
		my_oil=new red_oil(myinput);
	}
	else if(myinput.data.oil.general.oiltype == 5){
		Log << "Using SAE 10 Oil Model" << "\n";
		my_oil=new SAE10W(myinput);
	}
	else if(myinput.data.oil.general.oiltype == 6){
		Log << "Using Mil-H-5606 Oil Model" << "\n";
		my_oil=new milh5606_oil(myinput);
	}
	else if(myinput.data.oil.general.oiltype == 7){
		Log << "Using Exxon DTE Excel 32 Oil Model" << "\n";
		my_oil=new ExxonDTE10Excel32(myinput);
	}*/
	else if(myinput.data.oil.general.oiltype == 8){
		Log << "Using ISO 46 Oil Model" << "\n";
		my_oil=new ISO46_oil(myinput);
	}
	else if(myinput.data.oil.general.oiltype == 9){
		Log << "Using Skydrol_LD4 Model" << "\n";
		my_oil=new Parker_skydrol_oil(myinput);
	}
	/*else if(myinput.data.oil.general.oiltype == 50){
		Log << "Using Water Model" << "\n";
		my_oil=new water_oil(myinput);
	}*/
	else{
		Log << "Oil type " << myinput.data.oil.general.oiltype << " not supported!" << endl;
		exit(1);
	}
	

	//Temepratures parameters
	temperaturepistongap.Tmax = myinput.data.options_piston.numeric.Tmax;
	temperaturepistongap.THP = myinput.data.operating_conditions.T_HP;
	temperaturepistongap.TLP = myinput.data.operating_conditions.T_LP;
	temperaturepistongap.TCase = myinput.data.operating_conditions.T_Leak;

	
	//Surface face centers coordinates
	xyzf_gap.resize(N*M,3);				xyzf_gap = 0;
	

	//Deformation variables
	defK_gap.resize(N*M);				defK_gap = 0.0;
	defB_gap.resize(N*M);				defB_gap = 0.0;
	defK_p_gap.resize(N*M);				defK_p_gap = 0.0;
	defK_p_gap_old.resize(N*M);			defK_p_gap_old = 0.0;
	defK_p_gap_squeeze.resize(N*M);		defK_p_gap_squeeze = 0.0;
	defK_th_gap.resize(N*M);			defK_th_gap = 0.0;
	defB_p_gap.resize(N*M);				defB_p_gap = 0.0;
	defB_p_gap_old.resize(N*M);			defB_p_gap_old = 0.0;
	defB_p_gap_squeeze.resize(N*M);		defB_p_gap_squeeze = 0.0;
	defB_th_gap.resize(N*M);			defB_th_gap = 0.0;


	//Macro-geometry variables
	McrK.resize(N*M);					McrK=0.0;
	McrB.resize(N*M);					McrB=0.0;


	//Gap flow variables
	dht.resize(N*M);					dht=0;			
	sigma.resize(N*M);					sigma=0;
	T.resize(N*M*Q);					T=0;
	r_T.resize(N*M*Q);					r_T=0;
	TK_surf_gap.resize(N*M);			TK_surf_gap=0;
	TB_surf_gap.resize(N*M);			TB_surf_gap=0;
	Tfluid.resize(N*M*Q);				Tfluid=0;			
	vx.resize(N*M*Q);					vx=0;			
	vy.resize(N*M*Q);					vy=0;
	vx_p.resize(N*M*Q);					vx_p=0;			
	vy_p.resize(N*M*Q);					vy_p=0;
	vx_c.resize(N*M*Q);					vx_c=0;			
	vy_c.resize(N*M*Q);					vy_c=0;
	dz3.resize(N*M*Q);					dz3=0;			
	dz2.resize(N*M*Q);					dz2=0;
	dAy.resize(N*M*Q);					dAy=0;	
	oilviscosity.resize(N*M*Q);			oilviscosity=0;	
	oilviscosity_old.resize(N*M*Q);		oilviscosity_old=0;
	oildensity.resize(N*M*Q);			oildensity=0;	


	//Sizing piston fields
	p.resize(N*M);					p=0;
	pnew.resize(N*M);				pnew=0;	
	r_p.resize(N*M);				r_p=0;
	pold.resize(N*M);				pold=0;
	pfluid.resize(N*M);				pfluid=0;
	ploop.resize(N*M);				ploop=0;
	if(EHDTestRig)
	{
		pEHD.resize(1620);			pEHD=0.0;
	}
	h.resize(N*M);					h=0;
	hT.resize(N*M);					hT=0;
	xm.resize(N*M);					xm=0;
	ym.resize(N*M);					ym=0;
	hold.resize(N*M);				hold=0;
	h1.resize(N*M);					h1=0;
	hK.resize(N*M);					hK=0;
	mu.resize(N*M);					mu=0;
	muT.resize(N*M);				muT=0.0;
	dpx.resize(N*M);				dpx=0;			
	dpy.resize(N*M);				dpy=0;			
	ap.resize(N*M);					ap=0;			
	an.resize(N*M);					an=0;			
	as.resize(N*M);					as=0;			
	ae.resize(N*M);					ae=0;			
	aw.resize(N*M);					aw=0;			
	b.resize(N*M);					b=0;			
	

	//PistonEnergyGS vectors definition
	Tnew.resize(N*M*Q);				Tnew=0;			
	apE.resize(N*M*Q);				apE=0;		
	anE.resize(N*M*Q);				anE=0;			
	asE.resize(N*M*Q);				asE=0;		
	aeE.resize(N*M*Q);				aeE=0;			
	awE.resize(N*M*Q);				awE=0;			
	atE.resize(N*M*Q);				atE=1;			
	abE.resize(N*M*Q);				abE=1;			
	bE.resize(N*M*Q);				bE=0;			
	dvxz.resize(N*M*Q);				dvxz=0;
	dvyz.resize(N*M*Q);				dvyz=0;
	Dx.resize(N*M*Q);				Dx=0;
	Dy.resize(N*M*Q);				Dy=0;
	Dz.resize(N*M*Q);				Dz=0;
	Fx.resize(N*M*Q);				Fx=0;
	Fy.resize(N*M*Q);				Fy=0;
	Px.resize(N*M*Q);				Px=0;
	Py.resize(N*M*Q);				Py=0;
	Ax.resize(N*M*Q);				Ax=0;
	Ay.resize(N*M*Q);				Ay=0;
	QgapK.resize(N*M);				QgapK=0;
	QgapB.resize(N*M);				QgapB=0;
	dvxT_p.resize(N*M);				dvxT_p=0;
	dvyT_p.resize(N*M);				dvyT_p=0;
	dvxT_c.resize(N*M);				dvxT_c=0;
	dvyT_c.resize(N*M);				dvyT_c=0;
	taux.resize(N*M);				taux=0;
	tauy.resize(N*M);				tauy=0;
	FfK.resize(N*M);				FfK=0;
	FfKx.resize(N*M);				FfKx=0;
	FfKy.resize(N*M);				FfKy=0;
	MfKx.resize(N*M);				MfKx=0;
	MfKy.resize(N*M);				MfKy=0;
	FcK.resize(N*M);				FcK=0;
	FcKx.resize(N*M);				FcKx=0;
	FcKy.resize(N*M);				FcKy=0;
	McKx.resize(N*M);				McKx=0;
	McKy.resize(N*M);				McKy=0;


	//Piston forces vectors sizing
	forcespistongap.F_fluid.resize(4);
	forcespistongap.F_external.resize(4);
	forcespistongap.F_contact.resize(4);
	forcespistongap.dF.resize(4);

	//Check to see if there is a resume file
	/*WIN32_FIND_DATA ResumeData;
	HANDLE hResume;

	ResumeFile = "0";
	Log << "Looking for resume file." << "\n";
	hResume = FindFirstFile("./output/piston/resume.*.bin",&ResumeData);
	if(hResume == INVALID_HANDLE_VALUE){
		Log << "File Not Found: ";
		unsigned int temp = GetLastError();
		Log << temp << "\n";
	}
	else
	{
		Log << "Resume File Found: " << ResumeData.cFileName << "\n";
		FindClose(hResume);
		ResumeFile = ResumeData.cFileName;
	}*/


	//-------------------Read piston mesh and size fields-----------------//
	if(HeatTransfer || ThermalDeformation)
	{
		//read mesh 
		myThermal.readMeshThermal("piston");
		//read surface coordinates piston body
		myGapInput.readBodySurfacexyzThermal("piston",myinput.data.thermal.piston.meshfile);
	}
		//Load IM's dwm
	if(myinput.data.options_piston.general.PressureDeformation)
	{
		//read surface coordinates
		Log << "Reading pressure mesh coordinates file... " << "\t";
		myGapInput.readBodySurfacexyzPressure();
		Log << "done!" << "\n";
		Log << " " << "\n";
	};
	//size all piston fields
	SizeFieldsPiston();
	//if(!ResumeFile.string::compare("0")){
		//---------------------Calculate temperature piston---------------------//
		if(HeatTransfer)
		{
			myThermal.ThermalSolve("piston",qbi_piston,TK_body,TK_surf);
		};
		//-------------------Solve piston thermal deformation------------------//
		if(ThermalDeformation)
		{
			myFEMThermal.FEMThermalSolve("piston",TK_body,defK_th,myMesh.piston_name);
		};
	//}
	//free variables
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.E.free(); myMesh.v.free();
	myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear();
	

	//-------------------Read cylinder mesh and size fields-----------------//
	if(HeatTransfer || ThermalDeformation)
	{
		//read mesh 
		myThermal.readMeshThermal("cylinder");
		//read surface coordinates cylinder body
		myGapInput.readBodySurfacexyzThermal("cylinder",myinput.data.thermal.block.meshfile);
	}
	//size all cylinder fields based on mesh
	SizeFieldsCylinder();
	//if(!ResumeFile.string::compare("0")){
		//-------------------Calculate temperature cylinder------------------//
		if(HeatTransfer)
		{
			myThermal.ThermalSolve("cylinder",qbi_cylinder,TB_body,TB_surf);
		}
		//----------------Solve cylinder thermal deformation------------------//
		if(ThermalDeformation)
		{
			myFEMThermal.FEMThermalSolve("cylinder",TB_body,defB_th,myMesh.cylinder_name);
		}
	//}
	//free variables
	myMesh.xyz.free(); myMesh.conn.free(); myMesh.E.free(); myMesh.v.free();
	myMesh.alpha.free(); myThermal.nodeid_qb.clear(); myThermal.phi.clear(); 


	//Load IM's dwm
	if(myinput.data.options_piston.general.PressureDeformation)
	{
		//read matrices
		Log << "Reading influence matrices... " << "\t";
		myGapInput.readInfluenceMatricesPistonCylinder(myinput.data.options_piston.general.IM_piston_path,myinput.data.options_piston.general.IM_bushing_path);
		Log << "done!" << "\n";
		Log << " " << "\n";
	};
	
	//Set displacement chamber pressure
	PistonSetPressure();
	//Initialize fields
	PistonInitializePressure();
	PistonInitializeTemperature();
	PistonInitializeViscosity(1.0e8,120.0);
	//Initialize density --- Lizhi 03/19/2015
	PistonCalcDensity();

	//First revolution
	myGapResult.revcounter = 1 ;

	//cout.unsetf(ios::scientific); dwm
	
};
CPistonGap::~CPistonGap()	//Destructor
{


}; 
//Solve piston gap
void CPistonGap::PistonGap(vector<double> &xp,vector<double> &vp)
{
	Log << "Piston/Cylinder Convergence Loop..." << "\n";

	//Initialize parameters
	int FLAG_P = 1;
	counter = 0;
	AlphaDef = myinput.data.options_piston.numeric.AlphaDef;
	AlphaMu = myinput.data.options_piston.numeric.AlphaMu;
	AlphaP = myinput.data.options_piston.numeric.AlphaP;
	R_p = 1.0;	R_h = 1.0;	R_mu = 1.0;
	pnew = ploop;	pold = ploop;

	//Initialize groove
	numgvlizhi=0;
	if(myinput.data.options_piston.numeric.cgv.size() > 0)
	{
		Pistongroove();
		//pgvold.resize(numgvlizhi);
		//for(int j;j<numgvlizhi;j++)
		//	pgvold(j) = ploop(nlimitlizhi(j)+1);
	}
	//cout<<"Pistongroove()"<<"\n";
	//FSI fixed-point iteration convergence loop
	do
	{
		//calculate density --- Lizhi 03/17/2015
		//PistonCalcDensity();

		//Calculate film thickness
		PistonCalch(xp,vp);
		//cout<<"PistonCalch()"<<"\n";

		//Calculate piston squeeze velocity
		PistonCalcdht(vp);
		//cout<<"PistonCalcdht()"<<"\n";
		//Calculate fluid pressure - Reynolds equation
		if(ReynoldsMultiGrid)
		{

			if(myinput.data.options_piston.numeric.cgv.size() > 0)
			{ 
				Log << "In order to solve Reynolds equation for groove profiled piston, turn Mutigrid switch off!!" << endl;
				exit(1);
			};
			//Reynolds equation coefficients
			PistonReynoldsCalcCoefficientsMG();
			//Solve Reynolds equation
			PistonReynoldsMG();
		}
		else
		{
			//Reynolds equation coefficients
			PistonReynoldsCalcCoefficientsGS();
			//Solve Reynolds equation
			PistonReynoldsGS();
		};
		//Relax fluid pressure
		p = pold + AlphaP * (p - pold);
		//for(int j=0;j<numgvlizhi;j++)
		//	pgvlizhi(j) = pgvold(j) + AlphaP * (pgvlizhi(j) - pgvold(j));

		//Calculate fluid viscosity
		PistonCalcViscosity(1.0e8,120.0);
		//Relax fluid viscosity
		oilviscosity = oilviscosity_old + AlphaMu * (oilviscosity - oilviscosity_old);

		//Assign old pressure deformations
		defK_p_gap_old = defK_p_gap;
		defB_p_gap_old = defB_p_gap;

		//Update contact forces
		PistonCalcContactForces();

		//EHD
		if(PressureDeformation)
		{

			//Calculate pressure deformations
			PistonCylinderInterpolatePressureDeformations();
			//Relax pressure deformations
			defK_p_gap = defK_p_gap_old + AlphaDef * (defK_p_gap - defK_p_gap_old);
			defB_p_gap = defB_p_gap_old + AlphaDef * (defB_p_gap - defB_p_gap_old);
		}
		
		//Loop residual
		FLAG_P = PistonCalcLoopResidual();

		//Counter
		counter++;

	}while(FLAG_P == 1 && counter < nmax);
	//Counter output
	fout.open("./output/piston/matlab/PistonLoopCounter.dat",ios::app);
	fout << counter-1 << "\n";
	fout.close();
	fout.clear();
	//Residual spacing
	fout.open("./output/piston/matlab/Rfsi.dat",ios::app);
	fout << "\n";
	fout.close();
	fout.clear();

	int tempint = counter - 1;
	Log << "Loop Iterations: " << tempint << "\n";
	Log << "Done!" << "\n";


	//Assign pressure array to result
	ploop = p;
	//pgvold = pgvlizhi;

	defK_p_gap_squeeze = defK_p_gap;
	defB_p_gap_squeeze = defB_p_gap;

	//Calculate fluid velocities
	//cout<<"beforePistonCalcV"<<"\n";
	PistonCalcV();
	//cout<<"PistonCalcV"<<"\n";
	//PistonEnergyGS();
	//Calculate fluid temperature - Energy equation
	if(EnergyEquation)
	{
		//Solve Energy equation
		PistonEnergyGS();
		//cout<<"PistonEnergyGS()"<<"\n";
		//Calculate fluid density
		PistonCalcDensity();
		//cout<<"PistonCalcDensity"<<"\n";
	}


	//Calculate ehd test rig pressure field
	if(EHDTestRig)
	{
		PistonCalcEHDTestRigPressureField( );
	}


	//Calculate fluid film leakage in axial direction
	PistonCalcLeakage( );


	//Calculate fluid film viscous friction forces
	PistonCalcFrictionForces( );
	//PistonGuessFTG(); dwm
	PistonReadFTG();


	//Calculate heat fluxes coming from the fluid film
	PistonCylinderCalcGapThermalFlux( );


};
//Solve piston gap for Jacobian calculations
void CPistonGap::PistonGapJacobian(vector<double> &xp,vector<double> &vp,vector<double> &dF)
{

	//Calculate film thickness
	PistonCalch(xp,vp);

	//Calculate piston squeeze velocity
	PistonCalcdht(vp);

	//Reynolds equation
	if(ReynoldsMultiGrid)
	{

			if(myinput.data.options_piston.numeric.cgv.size() > 0)
			{ Log << "In order to solve Reynolds equation for groove profiled piston, turn Mutigrid switch off!!" << endl;
			exit(1);
			}
		//Coefficients
		PistonReynoldsCalcCoefficientsMG();
		//Solve
		PistonReynoldsMG();		
	}
	else
	{
		//Coefficients
		PistonReynoldsCalcCoefficientsGS();
		//Solve
		PistonReynoldsGS();
	}

	//Calculate the fluid forces due to gap pressure field
	PistonCalcFluidForces( );

	//Calculate the external forces
	if(myinput.data.options_piston.general.inclined_piston)
	{
		PistonCalcExternalForces_inclined_piston( );
	}
	else
	{
		
		PistonCalcExternalForces( );
	}

	//Calculate the force balance dF considering all the forces
	PistonCalcdF(dF);

};
//Calculate actual shaft angle
void CPistonGap::PistonCalcPhi(double time)
{
	//Phi angle [rad] between 0 and 2PI
	operatingpistongap.phi_rad = fmod(operatingpistongap.omega*time,2.0*PI);
	//Current simualation time
	timenew = time;
	//Return angle phi to variable in myNewtonIteration class 
	operatingpistongap.phi_deg = operatingpistongap.phi_rad*180/PI;
	Log << "PHI: " << operatingpistongap.phi_deg << " [deg]" << "\n";

};
//Determine displacement chamber pressure from file or ideal
void CPistonGap::PistonSetPressure(void)
{

	vector<double> pressure;
	if(ReadpFile)
	{ 
		pressure = PistonGetFilePressure( 0.0 );
	}
	else 
	{
		pressure = PistonGetIdealPressure( 0.0 );
	}
	//Assign calculated DC pressure to PistonGap variable and myNewtonIteration variable 
	operatingpistongap.pDC = pressure[0];
	operatingpistongap.pHP = pressure[1];
	operatingpistongap.pLP = pressure[2];

};
//Calculate ideal displacment chamber pressure profile
vector<double> CPistonGap::PistonGetIdealPressure(double deltaangle)
{
	vector<double> pressure;
	pressure.resize(3);
	double pHP = operatingpistongap.pHP;
	double pLP = operatingpistongap.pLP;
	double phi = fmod(operatingpistongap.phi_deg + deltaangle,360.0);
	

    if( (fabs(phi) < 3.0) && ( fabs(phi)>=0) )
	{         
		pressure[0] = pLP + (pHP - pLP)/3.0*phi;      
	}
    else if( (fabs(phi) < 179) && ( fabs(phi)>=3) )
	{        
		pressure[0] = pHP;
	}
    else if(  (fabs(phi) < 182) && ( fabs(phi) >=179) )
	{
		pressure[0] = pHP - (pHP-pLP)/3.0*( fabs(phi) - 179.00 );
	}
    else
	{
		pressure[0] = pLP;      
	}
	pressure[1] = pHP;
	pressure[2] = pLP;

    return pressure;

};
//Read simulated displacement chamber pressure profile
vector<double> CPistonGap::PistonGetFilePressure(double deltaangle)
{
	double phi = fmod(operatingpistongap.phi_deg + deltaangle,360.0);
	double newTime = 0;
	double newAngle = 0;
	double speed = 0;
	int size = (int) myGapInput.pFile.time.size();
	vector<double> pressure;
	speed = operatingpistongap.speed;
	for(int i=0; i < size ; i++)
	{
		newTime = myGapInput.pFile.time[i];
		newAngle = newTime * 360 * (speed / 60 );
		if(newAngle > 360.0)
		{
			newAngle -= 360.0;
		}
		if(phi<=newAngle)
		{
			pressure.push_back(myGapInput.pFile.pDC[i]);
			pressure.push_back(myGapInput.pFile.pHP[i]);
			pressure.push_back(myGapInput.pFile.pLP[i]);
			break;
		}
	}
	//In case simulated angle is bigger then last time in file get last pressure
	if(phi>newAngle)
	{
		pressure.push_back(myGapInput.pFile.pDC[size-1]);
		pressure.push_back(myGapInput.pFile.pHP[size-1]);
		pressure.push_back(myGapInput.pFile.pLP[size-1]);
	}
	return pressure;
		
};
//Initialize variables in between time steps
void CPistonGap::PistonInitializeVariables(double time)
{

	//Calculate shaft angle
	PistonCalcPhi(time);

	//Set displacement chamber pressure
	PistonSetPressure();
	double tempdouble = operatingpistongap.pDC*1.0e-5;
	Log << "pDC: " << tempdouble << " [bar]" << "\n";

	if(myinput.data.options_piston.general.inclined_piston)
	{
		//Calculate the piston stroke
		PistonCalcsK_inclined();

		//Calculate the piston velocity
		PistonCalcvK_inclined();
	}
	else
	{
		//Calculate the piston stroke
		PistonCalcsK();

		//Calculate the piston velocity
		PistonCalcvK();
	}
	//Calculate the variable gap length
	PistonCalclvar();
	
	//Update geometrical parameters only if lvar changes
	if( geometrypistongap.lvar != lvarold )
	{
		//Calculate fluid element main dimensions
		dy = geometrypistongap.lvar / (1.0*M);	
		dAz = dx * dy;
		//Calculate fluid film coordinates
		PistonCylinderGapSurfaceCoordinates();
		//Set levels coordinates and faces for MG
		if(ReynoldsMultiGrid)
		{
			PistonSetMG();
		};
	}
	//Assign old gap length value
	lvarold = geometrypistongap.lvar;

	//Calculate nodes and face centers in case of FSI
	if(HeatTransfer || PressureDeformation || ThermalDeformation)
	{
		//Search fluid-structure neighbors
		int nb = 10;
		PistonCylinderSearchNodesNeighbours(nb);
	}
	
	//Interpolate surface temperatures for current position
	if(HeatTransfer)
	{
		PistonCylinderInterpolateSurfaceTemperatures();
	}

	//Interpolate thermal surface deformations for current position
	if(ThermalDeformation)
	{
		PistonCylinderInterpolateThermalDeformations( );
	}

	//Results
	myGapResult.pDC = operatingpistongap.pDC*1.0e-5;	//Displacement chamber pressure
	myGapResult.pHP = operatingpistongap.pHP*1.0e-5;	//High pressure port pressure
	myGapResult.PLP = operatingpistongap.pLP*1.0e-5;	//Low pressure port pressure


};
//Initialize fluid fields: pressure, temperature and viscosity
void CPistonGap::PistonInitializePressure(void)
{

	for(int i=0;i<N*M;i++)
	{
		if(i%M==(M-1))
		{
			p(i) = operatingpistongap.pCase;
		}
		else if(i%M==0)
		{
			p(i) = operatingpistongap.pDC;
		}
		else
		{
			p(i) = (operatingpistongap.pDC + operatingpistongap.pCase)*0.5;
		};
	};
	//Other field
	pnew = p;
	pold = p;
	ploop = p;
}
void CPistonGap::PistonInitializeTemperature(void)
{
	double TDC,TCase,Tmax,phi;
	double skt,M1,M2,b,dt_1,dt_2;
	int st;


	//Fixed surface temperature profile for piston and bushing
	phi = operatingpistongap.phi_rad * 180/PI;
	if(phi < 180)
	{
		TDC = temperaturepistongap.THP;
	}
	else
	{
		TDC = temperaturepistongap.TLP;
	}
	TCase = temperaturepistongap.TCase;
	Tmax = temperaturepistongap.Tmax;
  

	//Numerical field with temperature gradient
	M1 = M*3.0/4.0;
	M2 = M/4.0;
	st = (int) ceil(M1);
	skt = floor(M2);
	dt_1= (Tmax - TDC)/st;
	dt_2= (TCase - Tmax)/skt;
	b = Tmax - dt_2*st;
	//Temperature
	T(Range(0,st)) = TDC + tensor::i*dt_1;
	T(Range(st+1,M-1)) = b + (tensor::i+(st+1))*dt_2;
	for(int i=0;i<N*Q;i++)
	{
		T(Range(i*M,(i+1)*M-1)) = T(Range(0,M-1));
	};
	Tnew = T;
	//Gap surface temperatures
	TK_surf_gap = T(Range(N*M*(Q-1),N*M*Q-1));
	TB_surf_gap = T(Range(0,N*M-1));

};
void CPistonGap::PistonInitializeViscosity(double plim,double Tlim)
{
	/*double oilPc1,oilPc2,oilTc1,oilTc2,oilbetaP,oilbetaT,oilW;
	double mu_0 = 1.702e-2;
	double A0 = 2.278e-8;
	double A1 = 1.477e-10;
	double A2 = 5.127e-13;
	double B = 5.428;
	double C = 94.72;

	oilPc1 = oilpistongap.oilPc1;
	oilPc2 = oilpistongap.oilPc2;
	oilTc1 = oilpistongap.oilTc1;
	oilTc2 = oilpistongap.oilTc2;
	oilbetaP = oilpistongap.oilbetaP;
	oilbetaT = oilpistongap.oilbetaT;
	oilW = oilpistongap.oilW;

	//Limit p and T
	Tfluid = where(T>Tlim,Tlim,T);
	pfluid = where(p>plim,plim,p);

	//Type 0
	if (oilpistongap.oiltype == 0)
	{
		oilviscosity = oilpistongap.oilviscosity;
	}

	//Type 1
	if (oilpistongap.oiltype == 1)
	{
		Array<double,1> a(N*M*Q);
		Array<double,1> b(N*M*Q);
		Array<double,1> mu_kin(N*M*Q);
		a=0.0; b=0.0; mu_kin=0.0;
		//Fluid average viscosity
		for(int i=0;i<Q;i++)
		{
			a(Range(i*N*M,(i+1)*N*M-1)) = oilPc1 * pfluid * pow( (Tfluid(Range(i*N*M,(i+1)*N*M-1))+273.15),oilPc2 );        //change for viscosity second option 10.11.07 Daniel
		};
		b = pow( 10.0, oilTc1 - oilTc2*log10(Tfluid+273.15) );
		mu_kin = oilW * exp( a ) * ( pow(10.0,b) - 0.5 );
		for(int i=0;i<Q;i++)
		{
			oilviscosity(Range(i*N*M,(i+1)*N*M-1)) = 1.0e-6 * mu_kin(Range(i*N*M,(i+1)*N*M-1)) * ( oilpistongap.oildensity*(1.0 + oilbetaP * pfluid - oilbetaT * (Tfluid(Range(i*N*M,(i+1)*N*M-1))- 20.0)) );
		};
	}

	//Type 2
	if (oilpistongap.oiltype == 2)
	{
		Array<double,1> alpha(N*M*Q);			alpha=0.0; 
		Array<double,1> Temp_coeff(N*M*Q);		Temp_coeff=0.0; 
		Array<double,1> Temp_coeff2(N*M*Q);		Temp_coeff2=0.0;
		//Temperature coefficient
		alpha = A0 - A1*Tfluid + A2*( pow(Tfluid,2.0) );
		for(int i=0;i<Q;i++)
		{
			Temp_coeff(Range(i*N*M,(i+1)*N*M-1)) =  B*( 50.0 - Tfluid(Range(i*N*M,(i+1)*N*M-1)) )/( C + Tfluid(Range(i*N*M,(i+1)*N*M-1)) ) + alpha(Range(i*N*M,(i+1)*N*M-1)) * pfluid ;
		};
		//Limit Temp_coeff to 4 for pressure calculations
		Temp_coeff = where( Temp_coeff>=4.5, 4.5, Temp_coeff );
		Temp_coeff2 = exp(Temp_coeff);
		oilviscosity = mu_0 * Temp_coeff2;
	}

	//Type 3
	if (oilpistongap.oiltype == 3)
	{
		//Skydrol - model developed using AMEsim & SAE specs by Marco Zecchi
		double pref = 1.00000000000000e005;
		double Tref = 1.00000000000000e001;
		double rhoref = 1.06690035585407e003;
		double muref = 4.78435154853258e-002;
		double cpref = 1.67717744042304e003;
		double lambdaref = 1.36331614988375e-001;

		double P1 = 1.02059515692251e-008;
		double T1 = -1.72947118956637e-002;
		double T2 = 6.53539084574889e-005;
		double PT = -1.40663906488398e-011;
		double PT2 = -3.12090395313133e-013;
		double T3 = -9.56687972187690e-008;
		double PT3 = 1.01110563823919e-015;

		Array<double,1> dp (p - pref);
		Array<double,1> dT (T - Tref);
		double Psi;
		Array<double,1> temp(N*M*Q);
		
		for(int j=0;j<Q;j++)
		{
			for(int i=0;i<(N*M);i++)
			{
				Psi = P1*dp(i) + T1*dT(j*N*M+i) + T2*pow(dT(j*N*M+i), 2) + PT*dp(i)*dT(j*N*M+i) + PT2*dp(i)*pow(T1,2) 
				+ T3*pow(dT(j*N*M+i), 3) + PT3*dp(i)*pow(T1, 3)
				;
				temp(j*N*M+i) = muref * pow(10.0,Psi);
			};
		};
		oilviscosity = temp;
	}*/

	/*for(int a=0;a<M*N;a++){
		oilviscosity(a) = my_oil -> get_mu(p(a),T(a));
	}*/

	//Lizhi 03/23/2015
	for(int b = 0;b<Q;b++)
	{
		for(int a = 0;a<M*N;a++)
		{
			//myp = mixed_friction * sigma(a) + p(a);
			oilviscosity(b*M*N+a) = my_oil -> get_mu(p(a),T(b*M*N+a));
		}
	}


	//Old value
	oilviscosity_old = oilviscosity;
};
//Interpolate structure fields from solid to fluid
void CPistonGap::PistonCylinderInterpolateSurfaceTemperatures(void)
{
	double phi_rad,speedK;
	int phioffset;
	Range all = Range::all();
	Array<double,1> Tsurf;


	//Parameters
	speedK = operatingpistongap.speedK;
	phi_rad = operatingpistongap.phi_rad;
	phioffset = (int) (floor(phi_rad/dphi)*speedK);



	//-------------Piston temperature interpolation-------------//
	TK_surf_gap = 0.0;
	TK_surf_gap = myGapUtils.InterpolateFields(FaceIdK_s2f_th,FaceDistK_s2f_th,TK_surf);


	//------------Cylinder temperature interpolation------------//
	TB_surf_gap = 0.0;
	TB_surf_gap = myGapUtils.InterpolateFields(FaceIdB_s2f_th,FaceDistB_s2f_th,TB_surf);


	//Smooth interpolated fields
	PistonCylinderSmoothField(TK_surf_gap,10);
	PistonCylinderSmoothField(TB_surf_gap,10);
	

	//Rotate piston temperature as seen by fluid
	if(phioffset>0)
	{
		Array<double,1> T_copy;
		T_copy.resize(N*M);
		T_copy = TK_surf_gap;
		TK_surf_gap(Range(phioffset*M,N*M-1)) = T_copy(Range(0,(N-phioffset)*M-1));
		TK_surf_gap(Range(0,phioffset*M-1)) = T_copy(Range((N-phioffset)*M,N*M-1));
	};

};
void CPistonGap::PistonCylinderInterpolatePressureDeformations(void)
{

	//int phioffset;
	//double phi_rad;
	//PistonCalcContactForces();

	//Calculate piston pressure deformation from IM
	PistonFEMPressureSurfaceDeformation( );
	
	//Calculate cylinder pressure deformation from IM
	PistonCylinderFEMPressureSurfaceDeformation( );

	//Interpolation piston
	defK_p_gap = 0.0;
	defK_p_gap = myGapUtils.InterpolateFields(NodeIdK_s2f_p,NodeDistK_s2f_p,defK_p);//added a function to rotate piston below
		
	//Interpolation cylinder
	defB_p_gap = 0.0;
	defB_p_gap = myGapUtils.InterpolateFields(NodeIdB_s2f_p,NodeDistB_s2f_p,defB_p);

	//Smooth interpolated fields
	PistonCylinderSmoothField(defK_p_gap,10);
	PistonCylinderSmoothField(defB_p_gap,10);

	//Rotate piston pressure deformation as seen by fluid
	/*phi_rad = operatingpistongap.phi_rad;
	phioffset = (int) (floor(phi_rad/dphi));
	if(phioffset>0)
	{
		Array<double,1> defK_copy;
		defK_copy.resize(N*M);
		defK_copy = defK_p_gap;
		defK_p_gap(Range(phioffset*M,N*M-1)) = defK_copy(Range(0,(N-phioffset)*M-1));
		defK_p_gap(Range(0,phioffset*M-1)) = defK_copy(Range((N-phioffset)*M,N*M-1));
	};*/

};
void CPistonGap::PistonCylinderInterpolateThermalDeformations(void)
{
	int phioffset;
	double phi_rad,speedK;
	
	speedK = operatingpistongap.speedK;
	phi_rad = operatingpistongap.phi_rad;
	phioffset = (int) (floor(phi_rad/dphi)*speedK);


	//----------------------PISTON INTERPOLATION----------------------//
	defK_th_gap = 0.0;
	//Interpolation
	defK_th_gap = myGapUtils.InterpolateFields(NodeIdK_s2f_th,NodeDistK_s2f_th,defK_th);


	//-----------------------CYLINDER INTERPOLATION----------------------//
	defB_th_gap = 0.0;
	//Interpolation
	defB_th_gap = myGapUtils.InterpolateFields(NodeIdB_s2f_th,NodeDistB_s2f_th,defB_th);


	//Smooth interpolated fields
	PistonCylinderSmoothField(defK_th_gap,10);
	PistonCylinderSmoothField(defB_th_gap,10);


	//Rotate piston thermal deformation as seen by fluid
	if(phioffset>0)
	{
		Array<double,1> defK_copy;
		defK_copy.resize(N*M);
		defK_copy = defK_th_gap;
		defK_th_gap(Range(phioffset*M,N*M-1)) = defK_copy(Range(0,(N-phioffset)*M-1));
		defK_th_gap(Range(0,phioffset*M-1)) = defK_copy(Range((N-phioffset)*M,N*M-1));
	};

};
void CPistonGap::PistonCylinderSmoothField(Array<double,1> &field,int iter)
{
	double fn,fs,fe,fw,fp;
	int counter;

	counter=0;
	do
	{
		for(int i=0;i<N*M;i++)
		{
			//North
			if(i%M==(M-1))
			{
				fn = 0.0;
			}
			else
			{
				fn = field(i+1);
			}
			//South
			if(i%M==0)
			{
				fs = 0.0;
			}
			else
			{
				fs = field(i-1);
			}
			//East
			if(i>=(N-1)*M)
			{
				fe = field(i-(N-1)*M);
			}
			else
			{
				fe = field(i+M);
			}
			//West
			if(i<M)
			{
				fw = field(i+(N-1)*M);
			}
			else
			{
				fw = field(i-M);
			}
			//Point
			fp = field(i);

			//Average smoothing
			double fsum = fn + fs + fe + fw;
			if(i%M==(M-1) || i%M==0)
			{
				fsum /= 3.0;
			}
			else
			{
				fsum /= 4.0;
			}
			field(i) = 0.5 * ( fp + fsum );
		}

		counter++;

	}while(counter<=iter);


};


void CPistonGap::PistonCalcvK_inclined(void){
	
	double omega, beta, phi,zeta;
	omega = operatingpistongap.omega;
	beta = operatingpistongap.beta_rad;
	phi = operatingpistongap.phi_rad;
	zeta = myinput.data.options_piston.general.zeta * PI / 180;
	double r_temp= rB - (2*rB*tan(beta)*tan(zeta)*sin(phi/2)*sin(phi/2));
	//Center of mass inertia force z directions
	

	operatingpistongap.vK = -1.0 * operatingpistongap.omega * r_temp *tan(beta)/cos(zeta) * sin(phi);
};
//Calculate piston stroke	
void CPistonGap::PistonCalcsK_inclined(void){
	double beta, betamax, phi, sK_beta, sK_gamma, delta_psi, phi_odp, zeta;

	beta = operatingpistongap.beta_rad;
	betamax = operatingpistongap.betamax_rad;
	
	phi = operatingpistongap.phi_rad;
	zeta = myinput.data.options_piston.general.zeta * PI / 180;

	sK_beta = - rB * tan(beta)/cos(zeta)*( 1 - cos(phi));
	geometrypistongap.sK = sK_beta - rB * (tan(betamax) - tan(beta));

};
//Calculate piston velocity
void CPistonGap::PistonCalcvK(void)
{

	double omega, beta, phi, gamma;
	omega = operatingpistongap.omega;
	beta = operatingpistongap.beta_rad;
	phi = operatingpistongap.phi_rad;
	gamma = operatingpistongap.gamma_rad;

	operatingpistongap.vK = -1.0 * operatingpistongap.omega * rB *( tan(beta) * sin(phi) + tan(gamma) * cos(phi) / cos(beta));
	
	/*//Test MG
	operatingpistongap.vK = -10.0;*/

};
//Calculate piston stroke
void CPistonGap::PistonCalcsK(void)
{
	double beta, betamax, phi, sK_beta, sK_gamma, gamma, delta_psi, phi_odp;

	beta = operatingpistongap.beta_rad;
	betamax = operatingpistongap.betamax_rad;
	gamma = operatingpistongap.gamma_rad;
	phi = operatingpistongap.phi_rad;

	phi_odp = (beta == 0) ? 0 : -atan(tan(gamma)/sin(beta));
	// offset in sk(0) due to gamma angle
	delta_psi = rB*tan(beta)*(1 - cos(phi_odp)) + rB*tan(gamma)*sin(phi_odp)/cos(beta);

	sK_beta = - rB * tan(beta) * ( 1 - cos(phi) );
	sK_gamma = - rB * tan(gamma) * sin(phi) / cos(beta);

	
	geometrypistongap.sK = sK_beta + sK_gamma - rB * ( tan(betamax) - tan(beta) )+delta_psi;
	//geometrypistongap.sK = sK1 - rB * ( tan(betamax) - tan(beta) );

};
//Calculate variable gap length
void CPistonGap::PistonCalclvar(void)
{
	double Tbeta = operatingpistongap.beta_rad;
	double Tbetamax = operatingpistongap.betamax_rad;
	double pos_A = 0.0;			//Position of the front flat piston edge from the bottom of the cylinder [m] (pos_R old CASPAR)
	double pos_B = 0.0;			//Position of the back flat piston edge from beginning of bushing [m] (pos_L old CASPAR)
	double pos_B_mesh = 0.0;;	//Position of the back flat piston edge from beginning of bushing [m] for meshing purposes (pos_L old CASPAR)
	double temp_B = 0.0;;		//Variable distance from piston edge B to cylinder block beginning edge [m] (pistonB old CASPAR)
	double delta_A = 0.0;;		//Distance from front flat piston surface to bushing end when piston is at ODC [m] (deltlmin old CASPAR)
	double delta_B = 0.0;;		//Distance from piston head center to cylinder block beginning edge when piston is at ODC [m] (deltlminB old CASPAR)
	double lvar = 0.0;										//Variable gap length [m]
	double lengthB = geometrypistongap.lengthB;				//Cylinder block length [m]
	double lengthcanalB = geometrypistongap.lengthcanalB;	//Cylinder block canal lenght [m]
	double lZ0 = geometrypistongap.lZ0;						//Displacment chamber lenght at ODC [m]
	double lF = geometrypistongap.lF;						//Bushing lenght [m]
	double le = geometrypistongap.le;						//Bushing beginning position [m]
	double lK = geometrypistongap.lK;						//Piston length [m]
	double lKG = geometrypistongap.lKG;						//Piston surface length [m]
	double sK = geometrypistongap.sK;						//Piston stroke [m]
	double T;					//Translation of piston due to J and K offsets.
	double K = myinput.data.geometry.offset_K;
	double J = myinput.data.geometry.offset_J;

	//Translation of piston due to J and K offsets.
	T = tan(Tbeta)*(J-K*tan(0.5*Tbeta)) - tan(Tbetamax)*(J-K*tan(0.5*Tbetamax));

	//Distance from front flat piston surface to bushing end when piston is at ODC [m]
	geometrypistongap.delta_A = lengthB - lengthcanalB - lZ0 - lF - le - T;
	delta_A = geometrypistongap.delta_A;
	

	//Position of the front flat piston edge from the bottom of the cylinder
	geometrypistongap.pos_A = lZ0 + sK + T;
	pos_A = geometrypistongap.pos_A;


	//Distance from piston head center to cylinder block beginning edge when piston is at ODC [m]
	geometrypistongap.delta_B = lK - delta_A - lF - le;
	delta_B = geometrypistongap.delta_B;


	//Position of the back flat piston edge from beginning of bushing [m]
	temp_B = sK + delta_B - (lK - lKG);
	pos_B = - temp_B - le;
	pos_B_mesh = pos_B;
	if (pos_B < 0.0)
	{
		pos_B = 0.0;
	}
	geometrypistongap.pos_B = pos_B;
	

	//Variable gap length [m]
	lvar = lengthB - lengthcanalB - pos_A - pos_B - le ;      
	if ( lvar >= lF - pos_B )
	{
		lvar = lF  - pos_B ;
	}
	//cout << "lvar is: "<< lvar << "\n";
	geometrypistongap.lvar = lvar ;
	//cout << "lvar is: "<< geometrypistongap.lvar << "\n";


	//Variable distance between piston head and beginning edge of cylinder block [m]
	geometrypistongap.lout = delta_B + sK ;


	//Variable distance between piston head and far end of the gap [m]
	geometrypistongap.zRK = geometrypistongap.lout + lvar + le + pos_B ; 


	//Determination of eventual portion of piston out of gap in DC [m]
	if(pos_B_mesh < 0.0)
	{
		geometrypistongap.lA = lKG - ( lvar + fabs(pos_B_mesh) ) ;
	}
	else
	{
		geometrypistongap.lA = lKG - lvar;
	}
	if(geometrypistongap.lA<=1.0e-6)
	{
		geometrypistongap.lA = 0.0;
	}
	

	//Determination of eventual portion of cylinder out of gap in DC [m]
	geometrypistongap.lB = lF - geometrypistongap.lvar - pos_B ;
	if(geometrypistongap.lB <= 1.0e-6)
	{
		geometrypistongap.lB = 0.0;
	}


	//Distance between cylinder block origin and cylinder block beginning edge [m]
	double beta = operatingpistongap.beta_rad ;
	double sKmax = 2.0 * rB * tan(beta) ;
	geometrypistongap.delta_0 = delta_B - 0.5 * sKmax ;


	//Distance from piston DC coordinate system origin to cylinder block origin [m]
	geometrypistongap.z_A = ( lF - geometrypistongap.lB ) + le + geometrypistongap.delta_0 ;
	
	
	//Distance from piston Case coordinate system origin to cylinder block origin [m]
	geometrypistongap.z_B = ( lF - geometrypistongap.lB - geometrypistongap.lvar ) + le + geometrypistongap.delta_0 ;


	//Remove extra decimals
	geometrypistongap.lA = floor(geometrypistongap.lA*1.0e6) / 1.0e6 ;
	geometrypistongap.lB = floor(geometrypistongap.lB*1.0e6) / 1.0e6 ;
	geometrypistongap.lvar = floor(geometrypistongap.lvar*1.0e6 )/ 1.0e6 ;
	geometrypistongap.zRK = floor(geometrypistongap.zRK*1.0e6) / 1.0e6 ; 

}
//Calculate fluid density
void CPistonGap::PistonCalcDensity(void)
{
	/*double oilRhozero = oilpistongap.oildensity;	 //[kg/m3]
	double oilbetaT = oilpistongap.oilbetaT;		//[1/K]
	double oilbetaP = oilpistongap.oilbetaP;		//[1/Pa]
	Array<double,1> oilRhoT(N*M*Q);

	double rs = 1047.03;		//[kg/m3]
	double als = 5.761668e-4;	//[1/K]
	double a1 = 0.0732965;
	double a2 = 1965.02;		//[bar]
	double a3 = -2.96813;		//[bar/K]

	//Type 0
	if (oilpistongap.oiltype == 0)
	{
		oildensity = oilRhozero;
	}

	//Type 1
	if (oilpistongap.oiltype == 1)
	{
		for(int i=0;i<Q;i++)
		{
			oildensity(Range(i*N*M,(i+1)*N*M-1)) = oilRhozero * (1 + oilbetaP * pfluid - oilbetaT * (Tfluid(Range(i*N*M,(i+1)*N*M-1)) - 20));
		};
	}

	//Type 2
	if (oilpistongap.oiltype == 2)
	{
		oilRhoT = rs*(1 - als*(Tfluid + 273.15));
		for(int i=0;i<Q;i++)
		{
			oildensity(Range(i*N*M,(i+1)*N*M-1)) = oilRhoT(Range(i*N*M,(i+1)*N*M-1)) / ( 1-a1*log( (a2+a3*(Tfluid(Range(i*N*M,(i+1)*N*M-1)) + 273.15) + 1e-5 * pfluid )/(a2+a3*(Tfluid(Range(i*N*M,(i+1)*N*M-1)) + 273.15)) ) );
		};
	}
	
	//Type 3
	if (oilpistongap.oiltype == 3)
	{
		//Skydrol - model developed using AMEsim & SAE specs by Marco Zecchi
		double pref = 1.00000000000000e005;
		double Tref = 1.00000000000000e001;
		double rhoref = 1.06690035585407e003;
		double muref = 4.78435154853258e-002;
		double cpref = 1.67717744042304e003;
		double lambdaref = 1.36331614988375e-001;

		double T1 = 7.80474256715724e-004;
		double T2 = 8.33673780452984e-007;
		double P1 = -8.76030682143087e-010;
		double P2 = 3.76207238974272e-018;
		double PT = -5.44707978002713e-012;
		double PT2 = 1.17343267760062e-014;
		double P2T = 2.78987659880380e-020;
		double P2T2 = -7.01924728092232e-023;
		
		Array<double,1> dp(pfluid - pref);
		Array<double,1> dT(Tfluid - Tref);

		double vs;
		for(int j=0;j<Q;j++)
		{
			for(int i=0;i<(N*M);i++)
			{
			vs = ( 
				(1.0/rhoref)*( 
				1 + P1*dp(i) + P2*pow(dp(i), 2.0) + T1*dT(N*M*j+i) + T2*pow(dT(N*M*j+i), 2.0) +  
				PT*dp(i)*dT(N*M*j+i) + PT2*dp(i)*pow(dT(N*M*j+i),2) + P2T*dT(N*M*j+i)*pow(dp(i),2) +
				P2T2*pow(dp(i), 2)*pow(dT(N*M*j+i), 2)
				));
			oildensity(i) = 1.0/vs;
			};
		};
		


	}*/
	for(int b=0;b< Q ;b++){
		for(int a=0;a<M*N;a++){
			oildensity(b*M*N+a) = my_oil -> get_rho(p(a),T(b*M*N+a));//maybe pfluid and Tfluid...
		}
	}
};
//Calculate fluid viscosity
void CPistonGap::PistonCalcViscosity(double plim,double Tlim)
{
	/*double oilPc1,oilPc2,oilTc1,oilTc2,oilbetaP,oilbetaT,oilW;
	double mu_0 = 1.702e-2;
	double A0 = 2.278e-8;
	double A1 = 1.477e-10;
	double A2 = 5.127e-13;
	double B = 5.428;
	double C = 94.72;

	oilPc1 = oilpistongap.oilPc1;
	oilPc2 = oilpistongap.oilPc2;
	oilTc1 = oilpistongap.oilTc1;
	oilTc2 = oilpistongap.oilTc2;
	oilbetaP = oilpistongap.oilbetaP;
	oilbetaT = oilpistongap.oilbetaT;
	oilW = oilpistongap.oilW;

	//Limit p and T
	Tfluid = where(T>Tlim,Tlim,T);
	pfluid = where(p>plim,plim,p);

	//Type 0
	if (oilpistongap.oiltype == 0)
	{
		oilviscosity = oilpistongap.oilviscosity;
	}

	//Type 1
	if (oilpistongap.oiltype == 1)
	{
		Array<double,1> a(N*M*Q);
		Array<double,1> b(N*M*Q);
		Array<double,1> mu_kin(N*M*Q);
		a=0.0; b=0.0; mu_kin=0.0;
		//Fluid average viscosity
		for(int i=0;i<Q;i++)
		{
			a(Range(i*N*M,(i+1)*N*M-1)) = oilPc1 * pfluid * pow( (Tfluid(Range(i*N*M,(i+1)*N*M-1))+273.15),oilPc2 );        //change for viscosity second option 10.11.07 Daniel
		};
		b = pow( 10.0, oilTc1 - oilTc2*log10(Tfluid+273.15) );
		mu_kin = oilW * exp( a ) * ( pow(10.0,b) - 0.5 );
		for(int i=0;i<Q;i++)
		{
			oilviscosity(Range(i*N*M,(i+1)*N*M-1)) = 1.0e-6 * mu_kin(Range(i*N*M,(i+1)*N*M-1)) * ( oilpistongap.oildensity * (1 + oilbetaP * pfluid - oilbetaT * (Tfluid(Range(i*N*M,(i+1)*N*M-1))- 20)) );
		};
	}

	//Type 2
	if (oilpistongap.oiltype == 2)
	{
		Array<double,1> alpha(N*M*Q);			alpha=0.0; 
		Array<double,1> Temp_coeff(N*M*Q);		Temp_coeff=0.0; 
		Array<double,1> Temp_coeff2(N*M*Q);		Temp_coeff2=0.0;
		//Temperature coefficient
		alpha = A0 - A1*Tfluid + A2*( pow(Tfluid,2.0) );
		for(int i=0;i<Q;i++)
		{
			Temp_coeff(Range(i*N*M,(i+1)*N*M-1)) =  B*( 50.0 - Tfluid(Range(i*N*M,(i+1)*N*M-1)) )/( C + Tfluid(Range(i*N*M,(i+1)*N*M-1)) ) + alpha(Range(i*N*M,(i+1)*N*M-1)) * pfluid ;
		};
		//Limit Temp_coeff to 4 for pressure calculations
		Temp_coeff = where( Temp_coeff>=4.5, 4.5, Temp_coeff );
		Temp_coeff2 = exp(Temp_coeff);
		oilviscosity = mu_0 * Temp_coeff2;
	}
	
	//Type 3
	if (oilpistongap.oiltype == 3)
	{
		//Skydrol - model developed using AMEsim & SAE specs by Marco Zecchi
		double pref = 1.00000000000000e005;
		double Tref = 1.00000000000000e001;
		double rhoref = 1.06690035585407e003;
		double muref = 4.78435154853258e-002;
		double cpref = 1.67717744042304e003;
		double lambdaref = 1.36331614988375e-001;

		double P1 = 1.02059515692251e-008;
		double T1 = -1.72947118956637e-002;
		double T2 = 6.53539084574889e-005;
		double PT = -1.40663906488398e-011;
		double PT2 = -3.12090395313133e-013;
		double T3 = -9.56687972187690e-008;
		double PT3 = 1.01110563823919e-015;

		Array<double,1> dp(pfluid - pref);
		Array<double,1> dT(Tfluid - Tref);
		
		double Psi;
		Array<double,1> temp(N*M*Q);
		/*for(int i=0;i<(N*M*Q-1);i++)
		{
			Psi = ( 
			P1*dp(i) + T1*dT(i) + T2*pow(dT(i), 2) + PT*dp(i)*dT(i) + PT2*dp(i)*pow(T1,2) 
			+ T3*pow(dT(i), 3) + PT3*dp(i)*pow(T1, 3)
									);
			temp(i) = muref * pow(10.0,Psi);
		};
		for(int j=0;j<Q;j++)
		{
			for(int i=0;i<(N*M);i++)
			{
				Psi = P1*dp(i) + T1*dT(j*N*M+i) + T2*pow(dT(j*N*M+i), 2) + PT*dp(i)*dT(j*N*M+i) + PT2*dp(i)*pow(T1,2) 
				+ T3*pow(dT(j*N*M+i), 3) + PT3*dp(i)*pow(T1, 3)
				;
				temp(j*N*M+i) = muref * pow(10.0,Psi);
			};
		};
		oilviscosity = temp;
	}*/

	double myp;

	for(int b = 0;b<Q;b++){

		for(int a = 0;a<M*N;a++)
		{
			myp = mixed_friction * sigma(a) + p(a);
			oilviscosity(b*M*N+a) = my_oil -> get_mu(myp,T(b*M*N+a));
		}
	}
};
void CPistonGap::PistonGuessViscosity(double &oilviscosityguess)
{
	/*double Temperature,Pressure,a,b,oilPc1,oilPc2,oilTc1,oilTc2,oilbetaP,oilbetaT,oilW,oilkinematicviscosity,oildynamicviscosity;
	
	double mu_0 = 1.702e-02;
	double A0 = 2.278e-8;
	double A1 = 1.477e-10;
	double A2 = 5.127e-13;
	double B = 5.428;
	double C = 94.72;
	double alpha,Temp_coeff,Temp_coeff2;

	oilPc1 = oilpistongap.oilPc1;
	oilPc2 = oilpistongap.oilPc2;
	oilTc1 = oilpistongap.oilTc1;
	oilTc2 = oilpistongap.oilTc2;
	oilbetaP = oilpistongap.oilbetaP;
	oilbetaT = oilpistongap.oilbetaT;
	oilW = oilpistongap.oilW;

	//Setting temperature and pressure
	Temperature = temperaturepistongap.TCase;
	Pressure = operatingpistongap.pDC;
				
	//---------VISCOSITY CALCULATION-----------//
	//Oil Type 0: Constant kinematic viscosity insert by the user
	if(oilpistongap.oiltype == 0)
	{
		oilviscosityguess = oilpistongap.oilviscosity;
	}
	//Oil Type 1: Linear formula defined by the user
	if(oilpistongap.oiltype == 1)
	{
		a = oilPc1 * Pressure * pow( Temperature+273.4,oilPc2 );        //change for viscosity second option 10.11.07 Daniel
		b = pow( 10.0, oilTc1 - oilTc2*log10(Temperature+273.4) );
		oilkinematicviscosity = oilW * exp(a) * (pow(10.0,b) - 0.5);
		oildynamicviscosity = 1.0e-6 * oilkinematicviscosity*( oilpistongap.oildensity*(1 + oilbetaP * Pressure - oilbetaT * (Temperature - 20)) );
		oilviscosityguess = oildynamicviscosity; 
	}
	//Oil Type 2: HLP32 oil
	if(oilpistongap.oiltype == 2)
	{
		alpha = A0 - A1*Temperature + A2*(Temperature*Temperature);
		Temp_coeff =  B*(50-Temperature)/(C+Temperature) + alpha*Pressure;
		if(Temp_coeff >= 4.5) 
		{
			Temp_coeff = 4.5;
		}
		Temp_coeff2 = exp(Temp_coeff);
		oilviscosityguess = mu_0 * Temp_coeff2;
	}
	
	//Type 3
	if (oilpistongap.oiltype == 3)
	{
		//Skydrol - model developed using AMEsim & SAE specs by Marco Zecchi
		double pref = 1.00000000000000e005;
		double Tref = 1.00000000000000e001;
		double rhoref = 1.06690035585407e003;
		double muref = 4.78435154853258e-002;
		double cpref = 1.67717744042304e003;
		double lambdaref = 1.36331614988375e-001;

		double P1 = 1.02059515692251e-008;
		double T1 = -1.72947118956637e-002;
		double T2 = 6.53539084574889e-005;
		double PT = -1.40663906488398e-011;
		double PT2 = -3.12090395313133e-013;
		double T3 = -9.56687972187690e-008;
		double PT3 = 1.01110563823919e-015;

		double dp = (Pressure - pref);
		double dT = (Temperature - Tref);
		
		double Psi = 
			P1*dp + T1*dT + T2*pow(dT, 2) + PT*dp*dT + PT2*dp*pow(T1,2) 
			+ T3*pow(dT, 3) + PT3*dp*pow(T1, 3)
									;
		oilviscosityguess = muref * pow(10.0,Psi);

	}*/

	oilviscosityguess = my_oil ->get_mu(operatingpistongap.pDC,temperaturepistongap.TCase);

};
//Calculate rigid film thickness from control points
void CPistonGap::PistonCalcRigidh(vector<double> &xp)
{
	//Calculate rigid film thickness linear slopes
	double dxm,dym,xp0,xp1,xp2,xp3,lvar;

	//Variable gap length
	lvar = geometrypistongap.lvar;
	
	//Piston eccentricities
	xp0 = xp[0];
	xp1 = xp[1];
	xp2 = xp[2];
	xp3 = xp[3];

	//Gap linear gradient
	dxm = (xp2 - xp0)/lvar;
	dym = (xp3 - xp1)/lvar;
	
	//Rigid linear gap variation
	xm=0.0; ym=0.0;
	for(int i=0;i<N;i++)
	{
		xm(Range(i*M,(i+1)*M-1)) = xp0 + (0.5 + tensor::i)*dy*dxm;
		ym(Range(i*M,(i+1)*M-1)) = xp1 + (0.5 + tensor::i)*dy*dym;
	};

	//Rigid film thickness relative to piston sliding surface
	hK = 0.0;
	hK = sqrt( pow((rZ * cos(phi) - xm),2.0) + pow((rZ * sin(phi) - ym),2.0) ) - rK;

};
//Calculate total film thickness considering macrogeometry and/or deformations
void CPistonGap::PistonCalch(vector<double> &xp,vector<double> &vp)
{
	vector<double> myxp; myxp.resize(4);
	double dt = myinput.data.options_piston.numeric.Simalphastep * (PI/180.0)/myinput.data.operating_conditions.speed;

	//Implicit ODE
	if(force_balance_iterative)
		myxp = xp;
	else
		for(short i = 0;i<4;i++)
			myxp[i] = xp[i] + vp[i] * dt;

	//Calculate rigid film thickness
	PistonCalcRigidh(myxp);

	//Calculate piston surface total gap deformation
	defK_gap = defK_p_gap + defK_th_gap;

	//Calculate cylinder surface total gap deformation
	defB_gap = defB_p_gap + defB_th_gap;

	//Calculate piston macro-geometry
	PistonCalcMcrK( );

	//Calculate cylinder macro-geometry
	PistonCalcMcrB( );

	//Relative deformed & macro-geometric film thickness (sliding surface)
	hK += - defK_gap - McrK;

	//Total deformed & macro-geometric film thickness
	h = hK + defB_gap - McrB;
	
	//Uncut gap
	h1 = 0.0;
	h1 = h;

	//Add step profile on gap height - Lizhi
	if(myinput.data.options_piston.numeric.stploc > 0)
	{
		h += stepfield;
	}

	//Minimum gap
	double hmin = geometrypistongap.hmin;	
	h = where(h<hmin, hmin, h);
	hK = where(hK<hmin, hmin, hK);

};
//Relax viscosity and thickness fields to smooth friction
void CPistonGap::PistonRelaxFields(void)
{
	Array<double,1> muT_old(N*M);
	Array<double,1> hT_old(N*M);

	//Recalculate viscosity with wider pressure boundary
	PistonCalcViscosity(5.0e8,120.0);

	//Previous time step values
	muT_old = muT;
	hT_old = hT;

	//Average fluid viscosity over fluid film
	muT = 0.0;
	for(int i=0;i<Q;i++)
	{
		muT(Range(0,N*M-1)) += oilviscosity(Range(i*N*M,(i+1)*N*M-1)) ;
	};
	muT /= (1.0*Q) ;
	//Under-relax viscosity to avoid force oscillation
	muT = muT_old + 0.1 * (muT - muT_old) ;

	//Fluid film thickness current time step
	hT = h ;
	//Under-relax fluid film to avoid force oscillation
	hT = hT_old + 0.1 * (hT - hT_old) ;

	//Gap delta z
	dz2=0.0;
	dz3=0.0;
	for(int i=0;i<Q;i++)
	{
		dz2(Range(i*N*M,(i+1)*N*M-1)) = hT / (1.0*Q);
		dz3(Range(i*N*M,(i+1)*N*M-1)) += (i+0.5)*dz2(Range(i*N*M,(i+1)*N*M-1));
	};

	//Fluid film area y direction
	dAy =  dx * dz2;

	//Reassign loop viscosity
	oilviscosity = oilviscosity_old ;

}
//Calculate piston macrogeometry
void CPistonGap::PistonCalcMcrK(void)
{

	//Gap geometrical variables
	double lKG = geometrypistongap.lKG;
	double rK_red = geometrypistongap.rK_red;
	double lK_hs = geometrypistongap.lK_hs;
	double lA = geometrypistongap.lA;

	//Cylindrical piston gap surface
	if(PistonMacroGeometry==0)
	{
		McrK = 0.0;
	}
	//Step-wise piston gap surface
	else if(PistonMacroGeometry==1)
	{
		//Diameters difference
		double d_difference = 0;
		//Diameter gradient over length	
		double d_tan_alpha = 0;
		//Gap axial cell coordinate
		Array<double,1> r1(N*M); r1 = 0.0;
		for(int i=0;i<N;i++)
		{
			r1(Range(i*M,(i+1)*M-1)) = (0.5 + tensor::i) * dy;
		};
		r1 += lA;
		//Total steps to count
		int n_max = (int) myGapInput.Geometry.stepwisegap_d_K.size() - 1 ;
		//If the first value of length is bigger than 0.0
		//Log << "Size of stepwisegap_d_k = " << n_max << "\n";
		McrK = where(r1 < myGapInput.Geometry.stepwisegap_l_K[0], myGapInput.Geometry.stepwisegap_d_K[0], McrK );
		//Loop over piston length
		int n = 0;
		Array<double,1> d_x(N*M); d_x = 0.0;
	    while(n < n_max)
		{
			//Case 1
			McrK = where(fabs(myGapInput.Geometry.stepwisegap_l_K[n]-r1)<=1e-12, myGapInput.Geometry.stepwisegap_d_K[n], McrK );
			//Case 2
			McrK = where(fabs(myGapInput.Geometry.stepwisegap_l_K[n+1]-r1)<=1e-12, myGapInput.Geometry.stepwisegap_d_K[n+1], McrK );
			//Case 3
			d_difference = myGapInput.Geometry.stepwisegap_d_K[n+1] - myGapInput.Geometry.stepwisegap_d_K[n];
			d_tan_alpha = d_difference/(myGapInput.Geometry.stepwisegap_l_K[n+1] - myGapInput.Geometry.stepwisegap_l_K[n]);
			d_x = d_tan_alpha * (r1 - myGapInput.Geometry.stepwisegap_l_K[n]);
			McrK = where(r1 > myGapInput.Geometry.stepwisegap_l_K[n] && r1 < myGapInput.Geometry.stepwisegap_l_K[n+1], myGapInput.Geometry.stepwisegap_d_K[n] + d_x, McrK );
			//Counter
			n++;
		}
	}
	//2D cylinder gap surface profile
	else if(PistonMacroGeometry == 2){
		//Gap geometrical variable
		//double lB = geometrypistongap.lB;
		//Diameters difference
		//double d_difference = 0;
		//Diameter gradient over length	
		//double d_tan_alpha = 0;
		//Gap axial cell coordinate
		Array<double,1> r1(N*M); r1=0.0;
		for(int i=0;i<N;i++)
		{
			r1(Range(i*M,(i+1)*M-1)) = (0.5 + tensor::i) * dy;
		};
		r1 += geometrypistongap.lA;
		//Gap circumferential cell coordinate
		Array<double,1> r2(N*M); r2 = 0.0;
		//Array<double,1> d_x(N*M); d_x=0.0;
		double pi = 4.0 * atan(1.0);
		double shift = -1.0 * myinput.data.geometry.speedK * operatingpistongap.phi_rad;
		while(shift < 0)
			shift += 2 * pi;
		for(int i = 0;i<N;i++)
		{
			r2(Range(i*M,(i+1)*M-1)) = fmod(i * (2 * pi / N) + shift,2*pi);
		};

		//Bilinear Interpolation Loop
		for(int i=0;i<M*N;i++){
			double x = r2(i);
			//Log << "x=" << x << "\n";
			double y = r1(i);
			//Log << "y=" << y << "\n";
			int y1M = floor(y/(myinput.data.geometry.lKG/(myGapInput.Geometry.PistonAx-1))); //Cell coordinate in axial direction
			//Log << "Y-index 1=" << y1M << "\n";
			int y2M = ceil(y/(myinput.data.geometry.lKG/(myGapInput.Geometry.PistonAx-1)));
			//Log << "Y-index 2=" << y2M << "\n";
			int x1N = floor(x*(myGapInput.Geometry.PistonCirc-1)/(2*pi));//Cell coordinate in circumferential direction.
			//Log << "X-index 1=" << x1N << "\n";
			int x2N = ceil(x*(myGapInput.Geometry.PistonCirc-1)/(2*pi));
			//if(x2N >= myGapInput.Geometry.BushingCirc)
			//	x2N = 0;
				
			//Log << "X-index 2=" << x2N << "\n";
			
			double q11 = myGapInput.Geometry.pistonsurface[x1N*(myGapInput.Geometry.PistonAx+1)+y1M];
			//Log << "q11=" << q11 << "\n";
			double q12 = myGapInput.Geometry.pistonsurface[x1N*(myGapInput.Geometry.PistonAx+1)+y2M];
			//Log << "q12=" << q12 << "\n";
			double q21 = myGapInput.Geometry.pistonsurface[x2N*(myGapInput.Geometry.PistonAx+1)+y1M];
			//Log << "q21=" << q21 << "\n";
			double q22 = myGapInput.Geometry.pistonsurface[x2N*(myGapInput.Geometry.PistonAx+1)+y2M];
			//Log << "q22=" << q22 << "\n";
			double x1 = x1N * (2*pi/(myGapInput.Geometry.PistonCirc-1));
			//Log << "x1=" << x1 << "\n";
			double x2 = x2N * (2*pi/(myGapInput.Geometry.PistonCirc-1));
			//Log << "x2=" << x2 << "\n";
			double y1 = y1M * (myinput.data.geometry.lKG/(myGapInput.Geometry.PistonAx-1));
			//Log << "y1=" << y1 << "\n";
			double y2 = y2M * (myinput.data.geometry.lKG/(myGapInput.Geometry.PistonAx-1));
			//Log << "y2=" << y2 << "\n";
			
			if((y1M == y2M) && (x1N == x2N)){
				McrK(i) = myGapInput.Geometry.pistonsurface[x1N*myGapInput.Geometry.PistonAx+y1M];
				//Log << "McrB=" << McrB(i) << "\n";
			}
			else if(y1M == y2M){
				McrK(i) = q11 + (x-x1)*(q21-q11)/(x2-x1);
				//Log << "McrB=" << McrB(i) << "\n";
			}
			else if(x1N == x2N){
				McrK(i) = q11 + (y-y1)*(q12-q11)/(y2-y1);
				//Log << "McrB=" << McrB(i) << "\n";
			}
			else{
				McrK(i) = (q11*(x2-x)*(y2-y)/((x2-x1)*(y2-y1))+q21*(x-x1)*(y2-y)/((x2-x1)*(y2-y1))+q12*(x2-x)*(y-y1)/((x2-x1)*(y2-y1))+q22*(x-x1)*(y-y1)/((x2-x1)*(y2-y1)));
				//Log << "McrB=" << McrB(i) << "\n";
			};
		};
		
	}
	//Half-spherical piston gap surface
	else if(PistonMacroGeometry==3)
	{
		Array<double,1> r1(N*M); r1 = 0.0;
		for(int i=0;i<N;i++)
		{
			r1(Range(i*M,(i+1)*M-1)) = (0.5 + tensor::i) * dy;
		};
		double rS = ( pow(lK_hs,2.0) + pow(rK_red,2.0) ) / ( 2.0*rK_red );
		double ym = lKG - lK_hs - lA;
		McrK = where( r1 < ym, 0.0 , rS - sqrt(pow(rS,2.0) - pow((r1-ym),2.0)) );
	}
	//Spherical piston gap surface
	else if(PistonMacroGeometry==4)
	{
		Array<double,1> r1(N*M); r1 = 0.0;
		for(int i=0;i<N;i++)
		{
			r1(Range(i*M,(i+1)*M-1)) = (0.5 + tensor::i) * dy;
		};
		double rS = ( pow(lKG,2.0)/4.0 + pow(rK_red,2.0) ) / ( 2.0*rK_red );
		Array<double,1> s(N*M);	s = 0.0;
		s = lKG/2.0 - (r1 + lA);
		McrK = rS - sqrt(pow(rS,2.0) - pow(s,2.0));
	}
	//Polynomial piston gap surface
	else
	{
		Array<double,1> r1(N*M); r1 = 0.0;
		for(int i=0;i<N;i++)
		{
			r1(Range(i*M,(i+1)*M-1)) = (0.5 + tensor::i) * dy;
		};
		int n = (int) myGapInput.Geometry.polygap_coeff.size() - 1 ;
		Array<double,1> x(N*M); x = 0.0;
		x = (r1 + lA)*1.0e3;  //x axis for polynomial is in mm
		double a0 = myGapInput.Geometry.polygap_coeff[n];
		Array<double,1> temp(N*M); temp = 0.0;
		for(int i=0; i<n ;i++)
		{
			 temp += myGapInput.Geometry.polygap_coeff[i]*pow(x,n-i);
		}
		temp += a0;
		McrK = temp / 1.0e3;
	};

};
//Calculate cylinder macrogeometry
void CPistonGap::PistonCalcMcrB(void)
{

	//Step-wise cylinder gap surface
	if(CylinderMacroGeometry == 1)
	{
		//Gap geometrical variable
		double lB = geometrypistongap.lB;
		//Diameters difference
		double d_difference = 0;
		//Diameter gradient over length	
		double d_tan_alpha = 0;
		//Gap axial cell coordinate
		Array<double,1> r1(N*M); r1=0.0;
		for(int i=0;i<N;i++)
		{
			r1(Range(i*M,(i+1)*M-1)) = (0.5 + tensor::i) * dy;
		};
		r1 += lB;
		//Total steps to count
		int n_max = (int) myGapInput.Geometry.stepwisegap_d_B.size()-1;
		//If the first value of length is bigger than 0.0
		McrB = where(r1 < myGapInput.Geometry.stepwisegap_l_B[0], myGapInput.Geometry.stepwisegap_d_B[0], McrB );
		//Loop over cylinder length
		int n = 0;
		Array<double,1> d_x(N*M); d_x=0.0;
		while(n < n_max)
		{
			//Case 1
			McrB = where(fabs(myGapInput.Geometry.stepwisegap_l_B[n]-r1)<=1e-12, myGapInput.Geometry.stepwisegap_d_B[n], McrB );
			//Case 2
			McrB = where(fabs(myGapInput.Geometry.stepwisegap_l_B[n+1]-r1)<=1e-12, myGapInput.Geometry.stepwisegap_d_B[n+1], McrB );
			//Case 3
			d_difference = myGapInput.Geometry.stepwisegap_d_B[n+1] - myGapInput.Geometry.stepwisegap_d_B[n];
			d_tan_alpha = d_difference/(myGapInput.Geometry.stepwisegap_l_B[n+1] - myGapInput.Geometry.stepwisegap_l_B[n]);
			d_x = d_tan_alpha * (r1 - myGapInput.Geometry.stepwisegap_l_B[n]);
			McrB = where(r1 > myGapInput.Geometry.stepwisegap_l_B[n] && r1 < myGapInput.Geometry.stepwisegap_l_B[n+1], myGapInput.Geometry.stepwisegap_d_B[n] + d_x, McrB );
			//Counter
			n++;
		}
	};

	//2D cylinder gap surface profile
	if(CylinderMacroGeometry == 2){
		//Gap geometrical variable
		//double lB = geometrypistongap.lB;
		//Diameters difference
		//double d_difference = 0;
		//Diameter gradient over length	
		//double d_tan_alpha = 0;
		//Gap axial cell coordinate
		Array<double,1> r1(N*M); r1=0.0;
		for(int i=0;i<N;i++)
		{
			r1(Range(i*M,(i+1)*M-1)) = (0.5 + tensor::i) * dy;
		};
		r1 += geometrypistongap.lB;
		//Gap circumferential cell coordinate
		Array<double,1> r2(N*M); r2 = 0.0;
		//Array<double,1> d_x(N*M); d_x=0.0;
		double pi = 4.0 * atan(1.0);
		for(int i = 0;i<N;i++)
		{
			r2(Range(i*M,(i+1)*M-1)) = i * (2 * pi / N);
		};

		//Bilinear Interpolation Loop
		for(int i=0;i<M*N;i++){
			double x = r2(i);
			//Log << "x=" << x << "\n";
			double y = r1(i);
			//Log << "y=" << y << "\n";
			int y1M = floor(y/(myinput.data.geometry.lF/(myGapInput.Geometry.BushingAx-1))); //Cell coordinate in axial direction
			//Log << "Y-index 1=" << y1M << "\n";
			int y2M = ceil(y/(myinput.data.geometry.lF/(myGapInput.Geometry.BushingAx-1)));
			//Log << "Y-index 2=" << y2M << "\n";
			int x1N = floor(x*(myGapInput.Geometry.BushingCirc-1)/(2*pi));//Cell coordinate in circumferential direction.
			//Log << "X-index 1=" << x1N << "\n";
			int x2N = ceil(x*(myGapInput.Geometry.BushingCirc-1)/(2*pi));
			//if(x2N >= myGapInput.Geometry.BushingCirc)
			//	x2N = 0;
				
			//Log << "X-index 2=" << x2N << "\n";
			
			double q11 = myGapInput.Geometry.bushingsurface[x1N*(myGapInput.Geometry.BushingAx+1)+y1M];
			//Log << "q11=" << q11 << "\n";
			double q12 = myGapInput.Geometry.bushingsurface[x1N*(myGapInput.Geometry.BushingAx+1)+y2M];
			//Log << "q12=" << q12 << "\n";
			double q21 = myGapInput.Geometry.bushingsurface[x2N*(myGapInput.Geometry.BushingAx+1)+y1M];
			//Log << "q21=" << q21 << "\n";
			double q22 = myGapInput.Geometry.bushingsurface[x2N*(myGapInput.Geometry.BushingAx+1)+y2M];
			//Log << "q22=" << q22 << "\n";
			double x1 = x1N * (2*pi/(myGapInput.Geometry.BushingCirc-1));
			//Log << "x1=" << x1 << "\n";
			double x2 = x2N * (2*pi/(myGapInput.Geometry.BushingCirc-1));
			//Log << "x2=" << x2 << "\n";
			double y1 = y1M * (myinput.data.geometry.lF/(myGapInput.Geometry.BushingAx-1));
			//Log << "y1=" << y1 << "\n";
			double y2 = y2M * (myinput.data.geometry.lF/(myGapInput.Geometry.BushingAx-1));
			//Log << "y2=" << y2 << "\n";
			
			if((y1M == y2M) && (x1N == x2N)){
				McrB(i) = myGapInput.Geometry.bushingsurface[x1N*myGapInput.Geometry.BushingAx+y1M];
				//Log << "McrB=" << McrB(i) << "\n";
			}
			else if(y1M == y2M){
				McrB(i) = q11 + (x-x1)*(q21-q11)/(x2-x1);
				//Log << "McrB=" << McrB(i) << "\n";
			}
			else if(x1N == x2N){
				McrB(i) = q11 + (y-y1)*(q12-q11)/(y2-y1);
				//Log << "McrB=" << McrB(i) << "\n";
			}
			else{
				McrB(i) = (q11*(x2-x)*(y2-y)/((x2-x1)*(y2-y1))+q21*(x-x1)*(y2-y)/((x2-x1)*(y2-y1))+q12*(x2-x)*(y-y1)/((x2-x1)*(y2-y1))+q22*(x-x1)*(y-y1)/((x2-x1)*(y2-y1)));
				//Log << "McrB=" << McrB(i) << "\n";
			};
		};
		
	};

};
//Calculate squeeze motion from control points
void CPistonGap::PistonCalcdht(vector<double> &vp)
{
	double dvxm,dvym,vp0,vp1,vp2,vp3,lvar,hmin;
	lvar = geometrypistongap.lvar;
	hmin = geometrypistongap.hmin;
	
	//Piston velocities
	vp0 = vp[0];
	vp1 = vp[1];
	vp2 = vp[2];
	vp3 = vp[3];		
	
	//Velocity gradients in gap
	dvxm = (vp2 - vp0)/lvar;
	dvym = (vp3 - vp1)/lvar;

	//Velocity linear profile
	Array<double,1> vxm(N*M);
	Array<double,1> vym(N*M);
	for(int i=0;i<N;i++)
	{
		vxm(Range(i*M,(i+1)*M-1)) = vp0 + (0.5 + tensor::i)*dy*dvxm;
		vym(Range(i*M,(i+1)*M-1)) = vp1 + (0.5 + tensor::i)*dy*dvym;
	};

	//Squeeze velocity field
	dht = vxm * cos(phi) + vym * sin(phi);

};
//Calculate gap axial and circumferential fluid velocity
void CPistonGap::PistonCalcV(void)
{
	double vK,speedK,omega,pDC,pCase,hmin;
	Range all = Range::all( );


	//Variables
	vK = operatingpistongap.vK;				
	omega = operatingpistongap.omega;		
	speedK = operatingpistongap.speedK;
	pDC = operatingpistongap.pDC;
	pCase = operatingpistongap.pCase;
	hmin = geometrypistongap.hmin;
	zerolizhi.resize(N*M);


	//Pressure gradients - central difference
	for(int i=0;i<N*M;i++)
	{
		//dpy backward - 3rd order
		if(i%M==(M-1))
		{
			dpy(i) = ( 11.0/6.0 * p(i) - 3.0 * p(i-1) + 3.0/2.0 * p(i-2) - 1.0/3.0 * p(i-3) ) / dy ;
		}
		//dpy forward - 3rd order
		else if(i%M==0)
		{
			dpy(i) = ( -11.0/6.0 * p(i) + 3.0 * p(i+1) - 3.0/2.0 * p(i+2) + 1.0/3.0 * p(i+3) ) / dy ;
		}
		//dpy central - 2nd order
		else
		{
			dpy(i) = ( p(i+1) -  p(i-1) ) / (2.0*dy) ;
		};
		//dpx central - 2nd order
		if(i>=(N-1)*M)
		{
			dpx(i) = ( p(i-(N-1)*M) - p(i-M) ) / (2.0*dx) ;
		}
		else if(i<M)
		{
			dpx(i) = ( p(i+M) - p(i+(N-1)*M) ) / (2.0*dx) ;
		}
		else
		{
			dpx(i) = ( p(i+M) - p(i-M) ) / (2.0*dx) ;
		};
	};



	//Relax fields to smooth results
	PistonRelaxFields( );


	zerolizhi = 1;//initialized 

	//lizhi try to artificially assign velocities in this HUGE and COMPLICATIED ARRAY!
	if(numgvlizhi > 0)
	{	
		//cout<<numgvlizhi<<"\n";
		for(int j=0;j<numgvlizhi;j++)
		{
			for(int i=0;i<N;i++)
			{
				zerolizhi(Range(nlimitlizhi(j) + i * M + 1, slimitlizhi(j) + i * M - 1) ) = 0;
			};
		};
	
	};

	//cout<<"check 1"<<"\n";

	//Fluid velocities
	for(int i=0;i<Q;i++)
	{
		//Poiseuille
		vx_p(Range(i*N*M,(i+1)*N*M-1)) = zerolizhi * 0.5 / muT * dpx * ( dz3(Range(i*N*M,(i+1)*N*M-1)) * dz3(Range(i*N*M,(i+1)*N*M-1)) - hT * dz3(Range(i*N*M,(i+1)*N*M-1)) ) ;
		vy_p(Range(i*N*M,(i+1)*N*M-1)) = zerolizhi * 0.5 / muT * dpy * ( dz3(Range(i*N*M,(i+1)*N*M-1)) * dz3(Range(i*N*M,(i+1)*N*M-1)) - hT * dz3(Range(i*N*M,(i+1)*N*M-1)) ) ;
		//Couette
		vx_c(Range(i*N*M,(i+1)*N*M-1)) = zerolizhi * omega * rK * speedK * dz3(Range(i*N*M,(i+1)*N*M-1)) / hT ;
		vy_c(Range(i*N*M,(i+1)*N*M-1)) = zerolizhi * vK * dz3(Range(i*N*M,(i+1)*N*M-1)) / hT ;
	};

	//Total
	vx = vx_p + vx_c ;
	vy = vy_p + vy_c ;

	//cout<<"check 2"<<"\n";

	//Velocity gradients in gap direction
	dvxz(Range(0,N*M-1)) = zerolizhi * ( vx(Range(N*M,2*N*M-1)) - vx(Range(0,N*M-1)) )  / dz2(Range(0,N*M-1));
	dvyz(Range(0,N*M-1)) = zerolizhi * ( vy(Range(N*M,2*N*M-1)) - vy(Range(0,N*M-1)) )  / dz2(Range(0,N*M-1));

	//cout<<"check 3"<<"\n";

	dvxz(Range((Q-1)*N*M,Q*N*M-1)) = zerolizhi * ( vx(Range((Q-1)*N*M,Q*N*M-1)) - vx(Range((Q-2)*N*M,(Q-1)*N*M-1)) ) / dz2(Range((Q-1)*N*M,Q*N*M-1));
	dvyz(Range((Q-1)*N*M,Q*N*M-1)) = zerolizhi * ( vy(Range((Q-1)*N*M,Q*N*M-1)) - vy(Range((Q-2)*N*M,(Q-1)*N*M-1)) ) / dz2(Range((Q-1)*N*M,Q*N*M-1));

	//cout<<"check 4"<<"\n";

	for(int i=1;i<Q-1;i++)
	{
		dvxz(Range(i*N*M,(i+1)*N*M-1)) = zerolizhi * ( vx(Range((i+1)*N*M,(i+2)*N*M-1)) - vx(Range((i-1)*N*M,(i)*N*M-1)) ) / (2.0*dz2(Range(i*N*M,(i+1)*N*M-1)));
		dvyz(Range(i*N*M,(i+1)*N*M-1)) = zerolizhi * ( vy(Range((i+1)*N*M,(i+2)*N*M-1)) - vy(Range((i-1)*N*M,(i)*N*M-1)) ) / (2.0*dz2(Range(i*N*M,(i+1)*N*M-1)));
	}

	//cout<<"check 5"<<"\n";
}
//Calculate gap axial and circumferential leakage
void CPistonGap::PistonCalcLeakage(void)
{
	Array<int,1> zerotemp;
	zerotemp.resize(N*M*Q);

	for(int k=0; k<Q; k++)
	{
		zerotemp(Range(k*N*M,(k+1)*N*M-1)) = zerolizhi;

		//cout<<"check1"<<"\n";

		for(int j=0; j<cgvlizhi.size(); j++)
		{
			if(cgvlizhi(j) == 3)
			{
				for(int i=0;i<N;i++)
				{
					zerotemp(Range(k * N * M + i * M, k * N * M + slimitlizhi(j) + i * M - 1) ) = 0;
					//cout<<"check2"<<"\n";
				};
			}
		}
	}

	//Poiseuille
	operatingpistongap.QSK_p = sum(vy_p * dAy * zerotemp) / (1.0*M);

	//Couette
	operatingpistongap.QSK_c = sum(vy_c * dAy * zerotemp) / (1.0*M);

	//Total
	operatingpistongap.QSK = operatingpistongap.QSK_p + operatingpistongap.QSK_c;
	//cout<<"QSK1: "<<operatingpistongap.QSK<<"\n";

	//lizhi correct leakage calculation for groove profile here!
	if(myinput.data.options_piston.numeric.cgv.size() > 0)
	{		
		operatingpistongap.QSK = operatingpistongap.QSK * (M*N*Q) / sum(zerotemp); 
		//cout<<"Ratio: "<<M / ( M - sum(slimitlizhi) + sum(nlimitlizhi) + numgvlizhi )<<"\n";
		//cout<<"M: "<<M<<"\n";
		//cout<<"QSK2: "<<operatingpistongap.QSK<<"\n";
		//for(int j=0;j<numgvlizhi;j++)
		//{
		//	cout<<"slimit: "<<slimitlizhi(j)<<"\n";
		//	cout<<"nlimit: "<<nlimitlizhi(j)<<"\n";
		//}
		//cout<<"sumslimit: "<<sum(slimitlizhi)<<"\n";
		//cout<<"numnlimit: "<<sum(nlimitlizhi)<<"\n";

	};

	//Volumetric power loss [W]
	double delta_p = operatingpistongap.pDC - operatingpistongap.pCase;
	PhiD_vol = operatingpistongap.QSK * delta_p;
	/*PhiD_vol = 0.0;
	int k;
	double dPcell;
	
	for(unsigned int i = 0;i<M*N;i++){
			
		if(i%M == 0){
			dPcell = (p(i)-p(i+1));
		}
		else if(i%M == (M-1)){
			dPcell = (p(i-1)-p(i));
		}
		else{
			dPcell = (p(i-1)-p(i+1))/2;
		}
		for(unsigned short j = 0;j<Q;j++){
			k = j*M*N+i;
			PhiD_vol += (vy_p(k) + vy_c(k))*dAy(k)*dPcell;
		}
	}*/
}
//Calculate external loop fluid-structure residual
int CPistonGap::PistonCalcLoopResidual(void)
{

	//Check the residual
	int FLAG = 1;
	Array<double,1> R2(N*M);	R2=0.0;
	Array<double,1> R3(N*M*Q);	R3=0.0;

	//Calculate pressure residual
	Rold_p = R_p;
	R2 = fabs( p - pold ) / p;
	R_p = sum(R2)/(N*M);

	//Log << R_p << "\n";

	//Calculate viscosity residual
	Rold_mu = R_mu;
	R3 = fabs( oilviscosity - oilviscosity_old ) / oilviscosity;
	R_mu = sum(R3)/(N*M*Q);

	//Calculate structure residual
	if(PressureDeformation)
	{
		Rold_h = R_h;
		R2 = fabs( h - hold ) / h;
		R_h = sum(R2)/(N*M);
	}

	//Control under-relaxation up
	if(R_p<Rold_p)
	{
		double iSecret = (rand() % 10 + 1)*1.0e-3;
		AlphaP += AlphaP*iSecret;
	}
	if(R_mu<Rold_mu)
	{
		double iSecret = (rand() % 10 + 1)*1.0e-3;
		AlphaMu += AlphaMu*iSecret;
	}
	if(PressureDeformation)
	{
		if(R_h<Rold_h)
		{
			double iSecret = (rand() % 10 + 1)*1.0e-3;
			AlphaDef += AlphaDef*iSecret;
		}
	}
	//Control under-relaxation down
	if(R_p>Rold_p)
	{
		double iSecret = (rand() % 100 + 1)*1.0e-3;
		AlphaP -= AlphaP*iSecret;
	}
	if(R_mu>Rold_mu)
	{
		double iSecret = (rand() % 100 + 1)*1.0e-3;
		AlphaMu -= AlphaMu*iSecret;
	}
	if(PressureDeformation)
	{
		if(R_h>Rold_h)
		{
			double iSecret = (rand() % 100 + 1)*1.0e-3;
			AlphaDef -= AlphaDef*iSecret;
		}
	}
	//Upper limit
	if( AlphaP > 0.5 )
	{
		AlphaP = 0.5;
	}
	if( AlphaMu > 0.5 )
	{
		AlphaP = 0.5;
	}
	if( AlphaDef > 0.5 )
	{
		AlphaP = 0.5;
	}
	//Lower limit - force convergence
	if( AlphaP < 1.0e-3 || AlphaMu < 1.0e-3 || AlphaDef < 1.0e-3 )
	{
		FLAG = 0;
		fout.open("./output/piston/matlab/Rfsi.dat",ios::app);
		fout << "Forced Convergence!" << "\n";
		fout.close();
		fout.clear();
	}
	
	//Output log
	fout.open("./output/piston/matlab/Rfsi.dat",ios::app);
	fout << scientific << R_p << "\t" <<  R_mu << "\t" << R_h << "\t" << AlphaP << "\t" <<  AlphaMu << "\t" << AlphaDef << "\n";
	fout.close();
	fout.clear();

	//Convergence
	if( R_p < Rmin_p || R_h < Rmin_h )
	{
		FLAG = 0;
	}

	//Assign old fluid pressure
	pold = p;
	//Assign old fluid viscosity
	oilviscosity_old = oilviscosity;

	//Assign old structure
	hold = h;

	return FLAG;

};
//Size piston fields
void CPistonGap::SizeFieldsPiston(void)
{
	Range all = Range::all();

	int n;

	//Assign & initialize fields pressure FEM
	if(PressureDeformation)
	{
		//piston face centers z coordinate
		n = (int) myGapInput.xyzfK_p.extent(0);
		zfK_p.resize(n);
		zfK_p = myGapInput.xyzfK_p(all,2);
		//piston nodes
		n = (int) myGapInput.xyznK_p.extent(0);
		//piston surface deformation
		defK_p.resize(n); defK_p = 0.0;
	};

	//Assign & initialize common fields thermal and FEM thermal
	if(HeatTransfer || ThermalDeformation)
	{
		//piston face centers z coordinate
		n = (int) myGapInput.xyzfK_th.extent(0);
		zfK_th.resize(n);
		zfK_th = myGapInput.xyzfK_th(all,2);
		//temperatures
		double TDC = 0.5*(myinput.data.operating_conditions.T_HP+myinput.data.operating_conditions.T_LP);
		double TCase = myinput.data.operating_conditions.T_Leak;
		TK_body.resize(myThermal.nCells);
		TK_body = 0.5*(TDC+TCase);
	}

	//Assign & initialize fields thermal
	if(HeatTransfer)
	{
		n = (int) myGapInput.xyzfK_th.extent(0);
		//gap surface temperature
		TK_surf.resize(n);
		//gap heat flux
		EbodyK.resize(n); EbodyK=0.0;
		EbodyK_old.resize(n); EbodyK_old=0.0;
		//heat flux vector
		int nqb = (int) myThermal.faceid_qb.size();
		for(int i=0;i<nqb;i++)
		{
			Array<double,1> qbi;
			qbi.resize((int) myThermal.faceid_qb[i].size());	qbi = 1.0e-9;
			qbi_piston.push_back(qbi);
		}
	}

	//Assign & initialize fields FEM thermal
	if(ThermalDeformation)
	{
		n = (int) myGapInput.xyznK_th.extent(0);
		defK_th.resize(n); defK_th = 0.0;
	};

}
//Size cylinder fields
void CPistonGap::SizeFieldsCylinder(void)
{
	Range all = Range::all();

	int n;

	//Assign & initialize fields pressure FEM
	if(PressureDeformation)
	{
		//cylinder face centers z direction
		n = (int) myGapInput.xyzfB_p.extent(0);
		zfB_p.resize(n);
		zfB_p = myGapInput.xyzfB_p(all,2);
		//cylinder nodes
		n = (int) myGapInput.xyznB_p.extent(0);
		//cylinder surface deformation
		defB_p.resize(n); defB_p = 0.0;
	};

	//Assign & initialize common fields thermal and FEM thermal
	if(HeatTransfer || ThermalDeformation)
	{
		//cylinder face centers z coordinate
		n = (int) myGapInput.xyzfB_th.extent(0);
		zfB_th.resize(n);
		zfB_th = myGapInput.xyzfB_th(all,2);
		double TDC = 0.5*(myinput.data.operating_conditions.T_HP+myinput.data.operating_conditions.T_LP);
		double TCase = myinput.data.operating_conditions.T_Leak;
		TB_body.resize(myThermal.nCells);
		TB_body = 0.5*(TDC+TCase);
	};

	//Assign & initialize fields thermal
	if(HeatTransfer)
	{
		n = (int) myGapInput.xyzfB_th.extent(0);
		//gap surface temperature
		TB_surf.resize(n);
		//heat flux
		EbodyB.resize(n); EbodyB=0.0;
		EbodyB_old.resize(n); EbodyB_old=0.0;
		//gap heat flux vector
		int nqb = (int) myThermal.faceid_qb.size();
		for(int i=0;i<nqb;i++)
		{
			Array<double,1> qbi;
			qbi.resize((int) myThermal.faceid_qb[i].size());	qbi = 1.0e-9;
			qbi_cylinder.push_back(qbi);
		}
	};

	//Assign & initialize fields FEM thermal
	if(ThermalDeformation)
	{
		n = (int) myGapInput.xyznB_th.extent(0);
		defB_th.resize(n); defB_th = 0.0;
	};


}
//Define gap cell centers coordinates
void CPistonGap::PistonCylinderGapSurfaceCoordinates(void)
{
	Range all=Range::all();

	//Fluid surface coordinates for structure to fluid interpolation
	xyzf_gap(all,0) = rK * cos(phi);
	xyzf_gap(all,1) = rK * sin(phi);
	for(int i=0;i<N;i++)
	{
		//Centroids coordinates coarse2fine
		xyzf_gap(Range(i*M,(i+1)*M-1),2) = (0.5 + tensor::i)*dy;
		//Fluid and contact forces moment arm
		zKj(Range(i*M,(i+1)*M-1)) = (0.5 + tensor::i)*dy;
	};
};
//Search fluid-structure neighbours
void CPistonGap::PistonCylinderSearchNodesNeighbours(int nb)
{

	double lA,lB,lch;
	Range all=Range::all();

	//Lenghts
	lA = geometrypistongap.lA;
	lB = geometrypistongap.lB;
	lch = geometrypistongap.lch;

	
	//---------------------------------PISTON PRESSURE MESH-----------------------------------//
	if(PressureDeformation)
	{
		//Shift gap coordinates according to piston position
		xyzf_gap(all,2) += (lA+lch);
		//Search structure to fluid nodes
		myGapUtils.SearchNeighbours(myGapInput.xyznK_p,xyzf_gap,NodeIdK_s2f_p,NodeDistK_s2f_p,nb); 
		//Search fluid to strucutre face centers
		myGapUtils.SearchNeighbours(xyzf_gap,myGapInput.xyzfK_p,FaceIdK_f2s_p,FaceDistK_f2s_p,nb);
		//Reset gap coordinates
		xyzf_gap(all,2) -= (lA+lch);
	}
	//---------------------------------PISTON THERMAL MESH-----------------------------------//
	if(HeatTransfer || ThermalDeformation)
	{
		//Shift gap coordinates according to piston position
		xyzf_gap(all,2) += (lA+lch);
		//Search structure to fluid face centers
		myGapUtils.SearchNeighbours(myGapInput.xyzfK_th,xyzf_gap,FaceIdK_s2f_th,FaceDistK_s2f_th,nb); 
		//Search structure to fluid nodes
		myGapUtils.SearchNeighbours(myGapInput.xyznK_th,xyzf_gap,NodeIdK_s2f_th,NodeDistK_s2f_th,nb); 	
		//Search fluid to structure face centers
		myGapUtils.SearchNeighbours(xyzf_gap,myGapInput.xyzfK_th,FaceIdK_f2s_th,FaceDistK_f2s_th,nb);
		//Reset gap coordinates 
		xyzf_gap(all,2) -= (lA+lch);
	}


	//------------------------------CYLINDER PRESSURE MESH-----------------------------------//
	if(PressureDeformation)
	{
		//Shift gap coordinates according to piston position
		xyzf_gap(all,2) += lB;
		//Search structure to fluid nodes
		myGapUtils.SearchNeighbours(myGapInput.xyznB_p,xyzf_gap,NodeIdB_s2f_p,NodeDistB_s2f_p,nb); 
		//Search fluid to strucutre face centers
		myGapUtils.SearchNeighbours(xyzf_gap,myGapInput.xyzfB_p,FaceIdB_f2s_p,FaceDistB_f2s_p,nb); 
		//Reset gap coordinates
		xyzf_gap(all,2) -= lB;
	}
	//---------------------------------CYLINDER THERMAL MESH-----------------------------------//
	if(HeatTransfer || ThermalDeformation)
	{
		//Shift gap coordinates according to piston position
		xyzf_gap(all,2) += lB;
		//Shift reference ccordinate system from block center to cylinder bore center
		if(EHDTestRig==0 && TriboTestRig==0)
		{
			myGapInput.xyzfB_th(all,1) -= rB;
			myGapInput.xyznB_th(all,1) -= rB;
		}
		//Search structure to fluid face centers
		myGapUtils.SearchNeighbours(myGapInput.xyzfB_th,xyzf_gap,FaceIdB_s2f_th,FaceDistB_s2f_th,nb);
		//Search structure to fluid nodes
		myGapUtils.SearchNeighbours(myGapInput.xyznB_th,xyzf_gap,NodeIdB_s2f_th,NodeDistB_s2f_th,nb);
		//Search fluid to strucutre face centers
		myGapUtils.SearchNeighbours(xyzf_gap,myGapInput.xyzfB_th,FaceIdB_f2s_th,FaceDistB_f2s_th,nb);
		//Reset gap coordinates
		xyzf_gap(all,2) -= lB;
		//Reset structure coordinates
		if(EHDTestRig==0 && TriboTestRig==0)
		{
			myGapInput.xyzfB_th(all,1) += rB;
			myGapInput.xyznB_th(all,1) += rB;
		}
	}
	

};
//Calculate EHD test rig pressure field in sensors position
void CPistonGap::PistonCalcEHDTestRigPressureField(void)
{
	Range all=Range::all();

	//Calculate EHD pressure sensors coordinates positions
	Array<double,1> phi(1620);
	Array<double,1> z(9);
	Array<double,2> xyzEHD(1620,3);

	z(0) = 2.5;		z(3) = 8.0;		z(6) = 23.66;
	z(1) = 2.5;		z(4) = 14.33;	z(7) = 26.16;
	z(2) = 5.0;		z(5) = 20.66;	z(8) = 26.16;	

	//Angular field EHD surface
	double dphi = 2.0*PI/180.0;
	for(int i=0;i<180;i++)
	{
		phi(Range(i*9,(i+1)*9-1)) = i*dphi;
	};

	//Cylindrical field EHD surface
	xyzEHD(all,0) = rK*cos(phi);
	xyzEHD(all,1) = rK*sin(phi);
	for(int i=0;i<180;i++)
	{
		for(int j=0;j<9;j++)
		{
			xyzEHD(i*9+j,2) = z(j)*1.0e-3;
		}
	}

	//Search neighbours
	Array<double,2> FaceDist; Array<int,2> FaceId;
	myGapUtils.SearchNeighbours(xyzf_gap,xyzEHD,FaceId,FaceDist,1);
	//Interpolate
	pEHD = 0.0;
	pEHD = myGapUtils.InterpolateFields(FaceId,FaceDist,p);
	//Correct pressure if sensor falls outside gap length
	double pCase = operatingpistongap.pCase;
	double lvar = geometrypistongap.lvar;
	pEHD = where(xyzEHD(all,2)>=lvar,pCase,pEHD);

};
//Calculate groove stuff outside the loop - Lizhi
void CPistonGap::Pistongroove(void)
{
	double pDC = operatingpistongap.pDC;
	double pCase = operatingpistongap.pCase;
	//lizhi initial the parameter of the groove here, should be able to read in the input file
	//strposilizhi = myinput.data.options_piston.numeric.stgv;//[m]distance between the first groove to the end of the piston near DC 0.0055 for EATON
	//wgvlizhi = myinput.data.options_piston.numeric.wgv;//[m]Radius of the groove
	//stpgvlizhi = myinput.data.options_piston.numeric.spgv;//[m]distance between each groove
	//numgvlizhi = myinput.data.options_piston.numeric.ngv;//number of the grooves
	//poblizhi = myinput.data.options_piston.numeric.pobgv;//groove on piston or bushing? 1: piston linear; 2: bushing linear; 3:piston constant; 4:bushing constant

	//lgvlizhi.resize(numgvlizhi);
	//pgvlizhi.resize(numgvlizhi);
	//nlimitlizhi.resize(numgvlizhi);
	//slimitlizhi.resize(numgvlizhi);

	//temparary hard code here Lizhi
	vector <double> pgv;
	vector <double> wgv;
	vector <int> cgv;
	/*pgv.push_back(10e-3);
	pgv.push_back(5e-3);
	pgv.push_back(3e-3);
	wgv.push_back(0.866e-3);
	wgv.push_back(0.866e-3);
	wgv.push_back(0.866e-3);
	cgv.push_back(1);
	cgv.push_back(3);
	cgv.push_back(4);*/

	numgvlizhi = myinput.data.options_piston.numeric.cgv.size();
	pgv.resize(numgvlizhi);
	wgv.resize(numgvlizhi);
	cgv.resize(numgvlizhi);

	pgv = myinput.data.options_piston.numeric.pgv;
	wgv = myinput.data.options_piston.numeric.wgv;
	cgv = myinput.data.options_piston.numeric.cgv;

	

	//cout<<"check 1 pgv.size: "<<pgv.size()<<"\n";
	//cout<<"pgv[0]: "<<pgv[0]<<" pgv[1]: "<<pgv[1]<<" pgv[2]: "<<pgv[2]<<"\n";
	//cout<<"wgv[0]: "<<wgv[0]<<" wgv[1]: "<<wgv[1]<<" wgv[2]: "<<wgv[2]<<"\n";
	//cout<<"cgv[0]: "<<cgv[0]<<" pgv[1]: "<<cgv[1]<<" pgv[2]: "<<cgv[2]<<"\n";
	//calculate the position of the groove relative to the DC end of the gap
	for (int igv = 0; igv < pgv.size(); igv++)
	{
		if (cgv[igv] == 1 || cgv[igv] == 3)
		{
			pgv[igv] -= geometrypistongap.lA;
		}
		else
		{
			pgv[igv] -= geometrypistongap.lB;
			cgv[igv] -= 1;//cgv = 1 or 3, 1 = linear, 3 = DC
		}
	}
	//cout<<"check 2"<<"\n";
	//sort the groove by position: [0] is DC end
	for (int i_sort = 0; i_sort < pgv.size() - 1; i_sort++)
	{
		for (int j_sort = 0; j_sort < pgv.size() - i_sort - 1; j_sort++)
		{
			if (pgv[j_sort] > pgv[j_sort + 1])
			{
				double templizhi = pgv[j_sort];
				pgv[j_sort] = pgv[j_sort + 1];
				pgv[j_sort + 1] = templizhi;
				templizhi = wgv[j_sort];
				wgv[j_sort] = wgv[j_sort + 1];
				wgv[j_sort + 1] = templizhi;
				templizhi = cgv[j_sort];
				cgv[j_sort] = cgv[j_sort + 1];
				cgv[j_sort + 1] = templizhi;
			}
		}
	}
	//cout<<"check 3"<<"\n";
	//merge grooves if overlap
	vector <double> pgv_new;
	vector <double> wgv_new;
	vector <double> cgv_new;
	for (int i = 0; i < pgv.size() - 1; i++)
	{
		if (pgv[i+1] - pgv[i] > 0.5 * (wgv[i+1] + wgv[i]) + dy)
		{
			pgv_new.push_back(pgv[i]);
			wgv_new.push_back(wgv[i]);

			cgv_new.push_back(cgv[i]);
		}
		else
		{
			double lb = min((pgv[i] - 0.5 * wgv[i]), (pgv[i+1] - 0.5 * wgv[i+1]));
			double ub = max((pgv[i] + 0.5 * wgv[i]), (pgv[i+1] + 0.5 * wgv[i+1]));
			pgv[i+1] = 0.5 * (lb + ub);
			wgv[i+1] = ub - lb;
			cgv[i+1] = max(cgv[i],cgv[i+1]);
		}
	}
	pgv_new.push_back(pgv[pgv.size()-1]);
	wgv_new.push_back(wgv[wgv.size()-1]);
	cgv_new.push_back(cgv[cgv.size()-1]);
	cout<<"check 4  pgv_new.size: "<<pgv_new.size()<<"\n";
	cout<<pgv_new[0]<<"\t"<<wgv_new[0]<<"\t"<<cgv_new[0]<<"\n";
	//check groove position
	if (pgv_new[0] - 0.5 * wgv_new[0] < 0)
	{
		pgv_new.erase(pgv_new.begin());
		wgv_new.erase(wgv_new.begin());
		cgv_new.erase(cgv_new.begin());
	}
	//cout<<"check 4.5  pgv_new.size: "<<pgv_new.size()<<"\n";
	//cout<<"pgv_new[0]: "<<pgv_new[0]<<" pgv_new[1]: "<<pgv_new[1]<<" pgv_new[2]: "<<pgv_new[2]<<"\n";
	while (pgv_new[pgv_new.size()-1] + 0.5 * wgv_new[wgv_new.size()-1] > geometrypistongap.lvar)
	{
		pgv_new.pop_back();
		wgv_new.pop_back();
		cgv_new.pop_back();
	}
	//cout<<"check 5  pgv_new.size: "<<pgv_new.size()<<"\n";
	//draw pressure profile
	vector <double> pgv_collapse;
	vector <double> px;
	vector <double> py;

	double collapse_temp = 0;
	for (int i = 0; i < pgv_new.size(); i++)
	{
		pgv_collapse.push_back(pgv_new[i] - 0.5 * wgv_new[i] - collapse_temp);
		collapse_temp += wgv_new[i];
	}


	px.push_back(0);
	py.push_back(pDC);
	for (int i = 0; i < pgv_collapse.size(); i++)
	{
		if (cgv_new[i] == 3)
		{
			px.push_back(pgv_collapse[i]);
			py.push_back(pDC);
		}
	}
	px.push_back(geometrypistongap.lvar - collapse_temp);
	py.push_back(pCase);

	//for (int i = 0; i < py.size(); i++)
	//{
	//	cout<<py[i]<<"	";
	//}
	//cout<<"\n";
	//for (int i = 0; i < py.size(); i++)
	//{
	//	cout<<px[i]<<"	";
	//}
	//cout<<"\n"<<"\n";
	//cout<<"check 6"<<"\n";

	
	numgvlizhi = pgv_new.size();
	cgvlizhi.resize(numgvlizhi);
	//lgvlizhi.resize(numgvlizhi);
	pgvlizhi.resize(numgvlizhi);
	nlimitlizhi.resize(numgvlizhi);
	slimitlizhi.resize(numgvlizhi);

	//assign pressure into grooves
	for (int i = 0; i < pgv_collapse.size(); i++)
	{
		for (int j = 0; j < px.size() - 1; j++)
		{
			if (pgv_collapse[i] <= px[j+1])
			{
				//pgvlizhi(i) = py[j] + (py[j+1] - py[j])/(px[j+1] - px[j]) * (pgv_new[i] - px[j]);
				pgvlizhi(i) = py[j] + (py[j+1] - py[j])/(px[j+1] - px[j]) * (pgv_collapse[i] - px[j]);
				//cout<<pgvlizhi(i)<<"	";
				break;
			}
		}
	}
	//cout<<"\n";

	//numgvlizhi = pgv_new.size();
	//lgvlizhi.resize(numgvlizhi);
	//pgvlizhi.resize(numgvlizhi);
	//nlimitlizhi.resize(numgvlizhi);
	//slimitlizhi.resize(numgvlizhi);

	//jlizhi=0;
	//cout<<"check 7"<<"\n";
	for(int jlizhi = 0; jlizhi < pgv_new.size(); jlizhi++)
	{
		cgvlizhi(jlizhi) = cgv_new[jlizhi];
		for(int i=0;i<M;i++)
		{
			if(fabs(dy*i+dy-(pgv_new[jlizhi]-wgv_new[jlizhi]/2))<dy/2)
			{
				nlimitlizhi(jlizhi)=i;
				//cout<<jlizhi<<"\n";
				//cout<<"n = "<<i<<"	";
			}
			if(fabs(dy*i-(pgv_new[jlizhi]+wgv_new[jlizhi]/2))<dy/2)
			{
				slimitlizhi(jlizhi)=i;
				//cout<<jlizhi<<"\n";
				//cout<<"s = "<<i<<"\n";
			}
		}
	}
	//cout<<"check 8"<<"\n";


	/*if (poblizhi==0)
	{
		jlizhi=0;
		for(int i=0;i<numgvlizhi;i++)
		{
			if ((strposilizhi+stpgvlizhi*i-geometrypistongap.lA) > wgvlizhi)
			{
				lgvlizhi(jlizhi)=strposilizhi+stpgvlizhi*i-geometrypistongap.lA;
				if (groove_condition[jlizhi] == 1)
				{
					pgvlizhi(jlizhi)=pDC-(pDC-pCase)*lgvlizhi(jlizhi)/(dy*M);
				}
				else if (groove_condition[jlizhi] == 0)
				{
				}
				pgvlizhi(jlizhi)=pDC-(pDC-pCase)*lgvlizhi(jlizhi)/(dy*M);
				jlizhi++;
			}
		
		};
	}
	
	else if (poblizhi==1)//linear on piston Lizhi
	{
		jlizhi=0;
		for(int i=0;i<numgvlizhi;i++)
		{
			if ((strposilizhi+stpgvlizhi*i-geometrypistongap.lA) > wgvlizhi)
			{
				lgvlizhi(jlizhi)=strposilizhi+stpgvlizhi*i-geometrypistongap.lA;
				pgvlizhi(jlizhi)=pDC-(pDC-pCase)*lgvlizhi(jlizhi)/(dy*M);
				jlizhi++;
			}
		
		};
	}
	else if (poblizhi==2)//linear on bushing Lizhi
	{
		for(int i=0;i<numgvlizhi;i++)
		{
			lgvlizhi(i)=strposilizhi+stpgvlizhi*i-geometrypistongap.lA;
			pgvlizhi(i)=pDC-(pDC-pCase)*lgvlizhi(i)/(dy*M);		
		};
	}
	else if (poblizhi==3)//constant dc pressure on piston Lizhi
	{
		jlizhi=0;
		for(int i=0;i<numgvlizhi;i++)
		{
			if ((strposilizhi+stpgvlizhi*i-geometrypistongap.lA) > wgvlizhi)
			{
				lgvlizhi(jlizhi)=strposilizhi+stpgvlizhi*i-geometrypistongap.lA;
				pgvlizhi(jlizhi)=pDC;
				jlizhi++;
			}
		
		};
	}
	else if (poblizhi==4)//constant dc pressure on bushing Lizhi
	{
		for(int i=0;i<numgvlizhi;i++)
		{
			lgvlizhi(i)=strposilizhi+stpgvlizhi*i-geometrypistongap.lA;
			pgvlizhi(i)=pDC;		
		};
	}


	for(int i=0;i<numgvlizhi;i++)
	{
		lgvlizhi(i)=strposilizhi+stpgvlizhi*i-geometrypistongap.lA;
		pgvlizhi(i)=pDC-(pDC-pCase)*lgvlizhi(i)/(dy*M);		
	};

	//Calculateing the locating of the groove boundary LIZHI
	jlizhi=0;
	for(int i=0;i<M;i++)
	{
		if(fabs(dy*i+dy-(lgvlizhi(jlizhi)-wgvlizhi/2))<dy/2 && jlizhi<numgvlizhi)
		{
			nlimitlizhi(jlizhi)=i;
			//cout<<jlizhi<<"\n";
			//cout<<i<<"\n";
		}
		if(fabs(dy*i-(lgvlizhi(jlizhi)+wgvlizhi/2))<dy/2 && jlizhi<numgvlizhi)
		{
			slimitlizhi(jlizhi)=i;
			//cout<<jlizhi<<"\n";
			//cout<<i<<"\n";
			jlizhi++;
		}
	};*/

	};
//Calculate step location here - Lizhi
void CPistonGap::Pistonsteplocation(void)
{
	stepfield.resize(M*N);
	stepfield = 0;

	for(int i=0;i<M;i++)
	{
		if(fabs(i*dy-myinput.data.options_piston.numeric.stploc) < dy/2)
		{
			stepboundary = i;
		}
	};
	
	//cout<<stepboundary<<"\n";

	for(int i=0;i<N;i++)
	{
		stepfield(Range( (i+1)*M-1-stepboundary , (i+1)*M-1) ) += myinput.data.options_piston.numeric.stpdep;
	};
};
//Calculate average number of cores used in IM calculations.
void CPistonGap::IMParallel(int i,int load,int n){
	if(i==0){
		processorload.resize(0);
		nprocessors.resize(0);
	}
	else if(i==1){
		processorload.push_back(load);
		nprocessors.push_back(n);
	}
	else if(i==2){
		float loadi, proc;
		loadi = 0;
		proc = 0;
		for(int i = 0;i<processorload.size();i++){
			loadi += processorload[i];
			proc += nprocessors[i];
		}
		loadi /= processorload.size();
		proc /= processorload.size();
		Log << "Processor Load: " << loadi << "\n";
		Log << "Parallel Processes: " << proc << "\n";
	}
};













//Calculate piston mesh motion for .vtk file sequence
void CPistonGap::PistonMeshMotion(void)
{

	Range all=Range::all();

	//Read .pos file
	ifstream fin;
	string line;
	vector<Array<double,1>> pos;
	fin.open("./outputs_piston/Parts_Positions.txt");
	int i=0;
	while(!fin.eof())
	{
		getline(fin,line);
		istringstream iss(line,istringstream::in);
		Array<double,1> m; m.resize(22); m=0.0;
		pos.push_back(m);
		for(int j=0;j<(int)m.extent(0);j++)
		{
			double c=0.0;
			iss >> c;
			pos[i](j)= c;
		};
		i++;
	}

	//angle and piston positions arrays
	int sz = (int) pos.size();
	Array<double,1> angle(sz);	angle=0.0;
	Array<double,2> xp(sz,4);	xp=0.0;
	for(int i=0;i<sz;i++)
	{
		angle(i) = pos[i](1);
		xp(i,0) = pos[i](2);
		xp(i,1) = pos[i](3);
		xp(i,2) = pos[i](4);
		xp(i,3) = pos[i](5);
	};

	//read mesh 
	myThermal.readMeshThermal("piston");


	//generate moving mesh .vtk
	double phi = 0.0;
	double dphi = 5.0;
	do
	{
		//deg to rad
		operatingpistongap.phi_rad = phi*PI/180.0;

		//calculate piston stroke
		PistonCalcsK( );
		
		//calculate variable gap lenght
		PistonCalclvar();
		double lch = geometrypistongap.lch;
		double lA = geometrypistongap.lA+lch;
		if(phi==0.0)
		{
			myMesh.xyz(all,2) -= lA;
		}
	
		//shaft angle
		int i=0;
		double phi_i=0.0;
		do
		{
			phi_i = angle(i);
			i++;
		}while(phi_i<phi);
		//positions
		double xp0 = xp(i-1,0);
		double xp1 = xp(i-1,1);
		double xp2 = xp(i-1,2);
		double xp3 = xp(i-1,3);

		fout.open("./xp.dat",ios::app);
		fout << scientific << xp0 << "\t" << xp1 << "\t" << xp2 << "\t" << xp3 << "\n";
		fout.close();
		fout.clear();

		//fluid film linear gradient
		double lvar = geometrypistongap.lvar;
		double sK = geometrypistongap.sK;
		double dxm = (xp2 - xp0)/lvar;
		double dym = (xp3 - xp1)/lvar;

		//file name based on angle
		string angle_i;
		char buffer[500];
		//angle number
		_itoa_s((int) phi,buffer,10);
		angle_i = buffer;


		double scale = 4.0e1;
		//output .vtk
		string output = "./phi_" + angle_i + ".vtk";
		fout.open(output.c_str());
		fout << "# vtk DataFile Version 2.0" << "\n";
		fout << "vtk output" << "\n";
		fout << "ASCII" << "\n";
		fout << "DATASET UNSTRUCTURED_GRID" << "\n";
		fout << "POINTS " << myMesh.nNodes << " double" << "\n";
		fout.precision(8);
		for(int i=0;i<myMesh.nNodes;i++)
		{
			fout << scientific <<  myMesh.xyz(i,0) + scale *  ( xp0 + dxm*myMesh.xyz(i,2) ) << "\t" <<  myMesh.xyz(i,1) + scale * ( xp1 + dym*myMesh.xyz(i,2) ) << "\t" <<  myMesh.xyz(i,2) + sK << "\n";
		}
		fout << "CELLS " << "\t" << myMesh.nCells << "\t" << 5*myMesh.nCells << "\n";
		for(int i = 0; i < myMesh.nCells; i++)
		{
			fout << 4 << "\t" << myMesh.conn(i,0) << "\t" << myMesh.conn(i,1)  << "\t" << myMesh.conn(i,2) << "\t" << myMesh.conn(i,3) << "\n";
		}
		fout << "CELL_TYPES " << myMesh.nCells << "\n";
		for(int i = 0; i < myMesh.nCells; i++)
		{
			fout << 10 << "\n";
		}
		fout.close();
		fout.clear();

		//increment angle
		phi+=dphi;

	}while(phi<360.0);


	//system("PAUSE");
};