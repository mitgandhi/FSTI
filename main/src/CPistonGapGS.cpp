#include "CPistonGap.h"
#include "../../caspar_input/input.h"
#include "logger.h"
#pragma once



extern class CGapInput myGapInput;
extern struct sGapResult myGapResult;
extern class CGapUtils myGapUtils;
extern class input myinput;

//Solve Reynolds equation using GS SOR method
void CPistonGap::PistonReynoldsGS(void)
{
	//cout<<"check before GS"<<"\n";
	//SOR loop
	double pn,ps,pe,pw;
	double alpha = 1.8;
	int iterations = 0;
	int iter_max = 10000;
	double R = 1.0;
	pn = 0.0; ps = 0.0; pe = 0.0; pw = 0.0;
	pcon.resize(pnew.size()); pcon = pnew;

	//coefficient of mass conservation to calculate pressure in grooves: p = (C-A)/(B-D) --- Lizhi 03/17/2015
	/*Array<double,1> A_groove;
	Array<double,1> B_groove;
	Array<double,1> C_groove;
	Array<double,1> D_groove;
	A_groove.resize(numgvlizhi);
	B_groove.resize(numgvlizhi);
	C_groove.resize(numgvlizhi);
	D_groove.resize(numgvlizhi);

	A_groove = 0;
	B_groove = 0;
	C_groove = 0;
	D_groove = 0;*/


	/*Array<double,1> ps_test;
	Array<double,1> pn_test;
	ps_test.resize(numgvlizhi);
	pn_test.resize(numgvlizhi);*/

	//Initialize the viscosity and density for mass conservation calculation --- Lizhi 03/17/2015
	/*mu = 0;
	rho = 0;
	for(int i=0;i<Q;i++)
	{
		mu(Range(0,N*M-1)) += oilviscosity(Range(i*N*M,(i+1)*N*M-1));
		rho(Range(0,N*M-1)) += oildensity(Range(i*N*M,(i+1)*N*M-1));
	};
	mu/=Q;
	rho/=Q;*/

	/*//Test MG
	clock_t start,end,ticks;
	start = clock();*/
	do
	{
		FLAGlizhi=0; //FlAG > 0 means groove
		for(int i=0;i<N*M;i++)
		{
			
			//North
			if(i%M==(M-1))
			{
				pn=0.0;
			}
			else
			{
				pn = pnew(i+1);
			}
			//South
			if(i%M==0)
			{
				ps=0.0;
			}
			else
			{
				ps = pnew(i-1);
			}
			//East
			if(i>=(N-1)*M)
			{
				pe = pnew(i-(N-1)*M);
			}
			else
			{
				pe = pnew(i+M);
			}
			//West
			if(i<M)
			{
				pw = pnew(i+(N-1)*M);
			}
			else
			{
				pw = pnew(i-M);
			}
			
			//SOR
			//cout<<i<<"  "<<FLAGlizhi<<"\n";
			if(FLAGlizhi==0)
			{
			pnew(i) = pnew(i) + alpha * ( ( ( an(i) * pn + as(i) * ps + ae(i) * pe + aw(i) * pw + b(i) ) / ap(i) ) - pnew(i) );
			pnew(i) = pnew(i) < 1e4? 1e4:pnew(i);
			}
			else
			{
				pnew(i)=pgvlizhi(FLAGlizhi-1); //Groove pressure Lizhi
			};

			
			for(int j=0;j<numgvlizhi;j++)
			{
				//Groove boundary recognition #1: next cell is south bound, next cell is out groove --- Lizhi 03/17/2015
				if((i+1-slimitlizhi(j))%M==0 && i>0)
				{
					FLAGlizhi=0;
					//cout<<"#1 "<<i<<"\n";
				};
				//Groove boundary recognition #2: this cell is north bound, next cell is in groove --- Lizhi 03/17/2015
				if((i-nlimitlizhi(j))%M==0 && i>0)
				{
					FLAGlizhi=j+1;
					pgvlizhi(j) -= 1.0 * Bs_g[j][(i-nlimitlizhi(j))/M] * (pnew(i) - pcon(i)) / (Cs_g(j) - Cn_g(j));
					//cout<<"#2 "<<i<<"\n";
					//cout<<i<<"---"<<ceil((i-nlimitlizhi(j))/M-0.5)<<"\n";
					//cout<<"check pgvlizhi at north "<<pgvlizhi(0)<<"  "<<pnew(i)<<"\n";
				};
				//Groove boundary recognition #3: this cell is south bound --- Lizhi 03/17/2015
				if((i-slimitlizhi(j))%M==0 && i>0)
				{
					pgvlizhi(j) += 1.0 * Bn_g[j][(i-slimitlizhi(j))/M] * (pnew(i) - pcon(i)) / (Cs_g(j) - Cn_g(j));
					//cout<<"#3 "<<i<<"\n";
					//cout<<"check pgvlizhi at south "<<pgvlizhi(0)<<"  "<<pnew(i)<<"\n";
				};
				//If groove is connect to the DC, overwrite the groove pressure to DC pressure
				if(cgvlizhi(j) == 3)
				{
					pgvlizhi(j) = operatingpistongap.pDC;
				}
			};
		}
		//Lizhi do test here
		/*Array<double,1> pgvlizhi_test;
		pgvlizhi_test.resize(numgvlizhi);
		for(int j=0;j<numgvlizhi;j++)
		{
			pgvlizhi_test(j) = (C_groove(j) - A_groove(j)) / (B_groove(j) - D_groove(j));
			cout<<"pgv("<<j<<") = "<<pgvlizhi(j)<<"; pgv_test("<<j<<") = "<<pgvlizhi_test(j)<<"\n";
		}
		//This is the real thing
		for(int j=0;j<numgvlizhi;j++)
		{
			pnew(nlimitlizhi(j)+1) = pnew(nlimitlizhi(j)+1) + alpha * ((C_groove(j) - A_groove(j)) / (B_groove(j) - D_groove(j)) - pnew(nlimitlizhi(j)+1));
			//cout<<"pgv("<<j<<") = "<<pgvlizhi(j)<<"\n";
			for(int i=0;i<N;i++)
				for(int ii=nlimitlizhi(j)+1;ii<slimitlizhi(j);ii++)
					pnew(i*M+ii) = pnew(nlimitlizhi(j)+1);
		}*/




		//Counter
		iterations++;

		//Residual GS
		PistonReynoldsCalcResidualGS(R);
		//Log<<R<<"\n";
		//cout<<"R = "<<R<<"\n";
		//for(int j=0;j<numgvlizhi;j++)
		//	cout<<"trouble("<<j<<") = "<<(C_groove(j) - A_groove(j)) / (B_groove(j) - D_groove(j)) - (pgvlizhi(j) - alpha * ((C_groove(0) - A_groove(0)) / (B_groove(0) - D_groove(0))))/(1 - alpha)<<"\n";
		//Residual in MPa
	}while(R > (Rmin_R * 1e6) && iterations < iter_max);

	//if(numgvlizhi > 0)
	//	cout<<"check pgvlizhi "<<pgvlizhi(0)<<"\n";


	/*//Test MG
	end = clock();
	ticks = end - start;
	double time = ticks / (double) CLOCKS_PER_SEC;

	fout.open("./outputs_piston/nGS.dat");
	fout << iterations << "\n";
	fout.close();

	fout.open("./outputs_piston/tGS.dat");
	fout << time << "\n";
	fout.close();*/


	//Log
	if(iterations>=iter_max)
	{
		fout.open("./output/piston/matlab/PistonGapLoopLog.txt",ios::app);
		fout << "Max Number of Iterations Reached in Reynolds Equation Loop" << "\n";
		fout << "Shaft Angle: " << "\t" << operatingpistongap.phi_rad*180/PI << "\n";
		fout << "Residual: " << "\t" << R << "\n";
		fout << " " << "\n";
		fout.close();
		fout.clear();
	}

	//Limit pressure
	if(PressureDeformation)
	{
		p = where(pnew>1.0e9,1.0e9,pnew);
	}
	else
	{
		p = where(pnew>1.0e9,1.0e9,pnew);
	};
	//p = where(p<1.0e4,1.0e4,p);


	/*//Test MG
	fout.open("./outputs_piston/pGS.dat");
	for(int i = 0; i < N ; i++)
	{
		for(int j = 0; j < M ; j++)
        {
			fout << scientific << pnew(j+i*M)*1.0e-5 << "\t" ;
        }
		fout << "\n";
	}
	fout.close();
	fout.clear();
	system("PAUSE");*/
	//cout<<"check after GS"<<"\n";
};

//Calculate Reynolds diffusive coeffcients and source term
void CPistonGap::PistonReynoldsCalcCoefficientsGS(void)
{
	//cout<<"check before GS coefficient"<<"\n";
	double vK,omega,speedK,pDC,pCase;
	Range all = Range::all( );

	//Geometry
	vK = operatingpistongap.vK;				
	omega = operatingpistongap.omega;		
	speedK = operatingpistongap.speedK;
	
	//Initializing the pressure 2D vector 
	pDC = operatingpistongap.pDC;
	pCase = operatingpistongap.pCase;
	mu = 0;
	rho2d.resize(N*M);
	rho2d = 0;
	for(int i=0;i<Q;i++)
	{
		mu(Range(0,N*M-1)) += oilviscosity(Range(i*N*M,(i+1)*N*M-1));
		rho2d(Range(0,N*M-1)) += oildensity(Range(i*N*M,(i+1)*N*M-1));
	};
	mu/=Q;
	rho2d/=Q;
	//Face diffusivity values
	double mun,mus,mue,muw,hKn,hKs,hKe,hKw,hn,hs,he,hw,dhKx,dhKy,defsqueezeK,defsqueezeB;
	if(numgvlizhi > 0)
	{
	As_g.resize(numgvlizhi);
	An_g.resize(numgvlizhi);
	Cs_g.resize(numgvlizhi);
	Cn_g.resize(numgvlizhi);
	As_g = 0;
	An_g = 0;
	Cs_g = 0;
	Cn_g = 0;
	Bs_g.resize(numgvlizhi, vector<double>(N,0));
	Bn_g.resize(numgvlizhi, vector<double>(N,0));
	Bsps.resize(numgvlizhi);
	Bnpn.resize(numgvlizhi);	
	Bsps = 0;
	Bnpn = 0;
	}
	//double hplus = 0;//lizhi study step profile here
	b = 0.0;
	for(int i=0;i<N*M;i++)
	{
		//North
		if(i%M==(M-1))
		{
			mun = mu(i);
			hn = h(i);
			hKn = hK(i);
		}
		else
		{
			mun = 0.5*( mu(i) + mu(i+1) );
			hn = 0.5*( h(i) + h(i+1) );
			hKn = 0.5*( hK(i) + hK(i+1) );
		}
		//South
		if(i%M==0)
		{
			mus = mu(i);
			hs = h(i);
			hKs = hK(i);
		}
		else
		{
			mus = 0.5*( mu(i) + mu(i-1) );
			hs = 0.5*( h(i) + h(i-1) );
			hKs = 0.5*( hK(i) + hK(i-1) );
		}
		//East
		if(i>=(N-1)*M)
		{
			mue = 0.5*( mu(i) + mu(i-(N-1)*M) );
			he = 0.5*( h(i) + h(i-(N-1)*M) );
			hKe = 0.5*( hK(i) + hK(i-(N-1)*M) );
		}
		else
		{
			mue = 0.5*( mu(i) + mu(i+M) );
			he = 0.5*( h(i) + h(i+M) );
			hKe = 0.5*( hK(i) + hK(i+M) );
		}
		//West
		if(i<M)
		{
			muw = 0.5*( mu(i) + mu(i+(N-1)*M) );
			hw = 0.5*( h(i) + h(i+(N-1)*M) );
			hKw = 0.5*( hK(i) + hK(i+(N-1)*M) );
		}
		else
		{
			muw = 0.5*( mu(i) + mu(i-M) );
			hw = 0.5*( h(i) + h(i-M) );
			hKw = 0.5*( hK(i) + hK(i-M) );
		}

		//Sharp the step boundary - Lizhi
		/*if(myinput.data.options_piston.numeric.stploc > 0)
		{
			if(i%M==(M-1-stepboundary))
			{
				hs -= myinput.data.options_piston.numeric.stpdep/2;
			};

			if(i%M==(M-1-stepboundary-1))
			{
				hn += myinput.data.options_piston.numeric.stpdep/2;
			};
		}*/

		//North
		an(i) = pow(hn,3.0)*dx / (dy*6.0*mun);
		//Dirichlet boundary
		if(i%M==(M-1) && i>0)
		{
			an(i) *= 2.0;
			b(i) += an(i)*pCase;
		};

		//South
		as(i) = pow(hs,3.0)*dx / (dy*6.0*mus);
		//Dirichlet boundary
		if(i%M==0)
		{
			as(i) *= 2.0;
			b(i) += as(i)*pDC;
		};

		//East
		ae(i) = pow(he,3.0)*dy / (dx*6.0*mue);
		//West
		aw(i) = pow(hw,3.0)*dy / (dx*6.0*muw);
		//Point
		ap(i) = an(i) + as(i) + ae(i) + aw(i);
		//Gradients sliding part
		dhKx = (hKe - hKw)/dx;
		dhKy = (hKn - hKs)/dy;
		//Change in deformation
		defsqueezeK = (CPistonGap::defK_p_gap(i) - CPistonGap::defK_p_gap_squeeze(i))/CPistonGap::dx;
		defsqueezeB = (CPistonGap::defB_p_gap(i) - CPistonGap::defB_p_gap_squeeze(i))/CPistonGap::dx;
		//Source
		b(i) += ( - speedK * omega * rK * (he - hw)/dx - vK * (hn - hs)/dy 
			+ 2.0 * (dht(i) + speedK * omega * rK * dhKx + vK * dhKy + defsqueezeK - defsqueezeB ) ) * dx * dy;
		//Groove coefficient
		if(numgvlizhi > 0)
		{
			for(int j=0;j<numgvlizhi;j++)
			{
				//Groove boundary recognition #2: this cell is north bound, next cell is in groove --- Lizhi 03/17/2015
				if((i-nlimitlizhi(j))%M==0 && i>0)
				{
					//FLAGlizhi=j+1;
					As_g(j) += 6 * vK * rho2d(i) * h(i);
					//cout<<"i = "<<i<<"\n";
					//cout<<"rho "<<rho2d(i)<<"; mu "<<mu(i)<<"\n";
					//cout<<"h = "<<h(i)<<"\n";
					//cout<<"As = "<<As_g(0)<<"\n";
					Cs_g(j) += -1/dy * rho2d(i) * h(i) * h(i) * h(i) / mu(i);
					//cout<<"Cs = "<<Cs_g(0)<<"\n";
					Bs_g[j][(i-nlimitlizhi(j))/M] = 1/dy * rho2d(i) * h(i) * h(i) * h(i) / mu(i);
					//cout<<"Bs_g["<<j<<"]["<<(i-nlimitlizhi(j))/M<<"] "<<Bs_g[j][(i-nlimitlizhi(j))/M]<<"\n";
					Bsps(j) += 1/dy * rho2d(i) * pnew(i) * h(i) * h(i) * h(i) / mu(i);
					//cout<<"Bsps = "<<Bsps(0)<<"\n";
				};
				//Groove boundary recognition #3: this cell is south bound --- Lizhi 03/17/2015
				if((i-slimitlizhi(j))%M==0 && i>0)
				{
					An_g(j) += 6 * vK * rho2d(i) * h(i);
					//cout<<"i = "<<i<<"\n";
					//cout<<"rho "<<rho2d(i)<<"; mu "<<mu(i)<<"\n";
					//cout<<"h = "<<h(i)<<"\n";
					//cout<<"An = "<<An_g(0)<<"\n";
					Cn_g(j) += 1/dy * rho2d(i) * h(i) * h(i) * h(i) / mu(i);
					//cout<<"Cn = "<<Cn_g(0)<<"\n";
					Bn_g[j][(i-slimitlizhi(j))/M] = -1/dy * rho2d(i) * h(i) * h(i) * h(i) / mu(i);
					//cout<<"Bn_g["<<j<<"]["<<(i-slimitlizhi(j))/M<<"] "<<Bn_g[j][(i-slimitlizhi(j))/M]<<"\n";
					Bnpn(j) += -1/dy * rho2d(i) * pnew(i) * h(i) * h(i) * h(i) / mu(i);
					//cout<<"Bnpn = "<<Bnpn(0)<<"\n";
				};
			};
		};
	};
	for(int j=0;j<numgvlizhi;j++)
	{
		pgvlizhi(j) = (An_g(j) - As_g(j) + Bnpn(j) - Bsps(j))/(Cs_g(j) - Cn_g(j));
		/*cout<<"An_g("<<j<<") = "<<An_g(j)<<"\n";
		cout<<"As_g("<<j<<") = "<<As_g(j)<<"\n";
		cout<<"Cn_g("<<j<<") = "<<Cn_g(j)<<"\n";
		cout<<"Cs_g("<<j<<") = "<<Cs_g(j)<<"\n";
		cout<<"Bnpn("<<j<<") = "<<Bnpn(j)<<"\n";
		cout<<"Bsps("<<j<<") = "<<Bsps(j)<<"\n";*/
		//cout<<"check pgvlizhi("<<j<<") = "<<pgvlizhi(j)/1e5<<" [bar]"<<"\n";
	}

};
//Calculate Reynolds residual for GS SOR method
void CPistonGap::PistonReynoldsCalcResidualGS(double &R)
{

	R = max(fabs(pcon - pnew));

	pcon = pnew;

	/*double pn,ps,pe,pw;

	FLAGlizhi=0;

	//Calculate residual field
	r_p = 0.0;
	for(int i=0;i<N*M;i++)
	{
		//North
		if(i%M==(M-1) && i>0)
		{
			pn = 0.0;
		}
		else
		{
			pn = pnew(i+1);
		}
		//South
		if(i%M==0)
		{
			ps = 0.0;
		}
		else
		{
			ps = pnew(i-1);
		}
		//East
		if(i>=(N-1)*M)
		{
			pe = pnew(i-(N-1)*M);
		}
		else
		{
			pe = pnew(i+M);
		}
		//West
		if(i<M)
		{
			pw = pnew(i+(N-1)*M);
		}
		else
		{
			pw = pnew(i-M);
		}

		//Groove boundary recognition #1 Lizhi
		if(numgvlizhi > 0)
		{
		for(int j=0;j<numgvlizhi;j++)
		{
			if((i-slimitlizhi(j))%M==0)
			{
				FLAGlizhi=0;
			};
		};
		};
		
		//Residual LIZHI
		r_p(i) = b(i)  - ( ap(i) * pnew(i) - an(i) * pn - as(i) * ps - ae(i) * pe - aw(i) * pw  );

		if(FLAGlizhi)
			{
				r_p(i) = 0;//residual is 0 in groove
			}
		
		//Groove boundary recognition #2 Lizhi		
		if(numgvlizhi > 0)
		{
		for(int j=0;j<numgvlizhi;j++)
			{
				if((i-nlimitlizhi(j))%M==0)
				{
					FLAGlizhi=j+1;
				};
			};
		};

	};

	//Calculate scaled residual
	R = 0.0;
	double n = sum(fabs(r_p));
	double d = sum(fabs(ap*pnew));
	//Test MG
	double n = sqrt(sum(pow(r_p,2.0)));
	double d = sqrt(sum(pow(b,2.0)));
	R = n/d;
	*/


}
//Solve Energy equation using GS SOR method
void CPistonGap::PistonEnergyGS(void) 
{
	double cp,alpha,lambda,TDC,TCase,phi;

	//Fluid properties
	lambda = my_oil -> get_lambda();//oilpistongap.oillambda;
	cp = my_oil -> get_C();//oilpistongap.oilC;

	//Boundaries
	phi = operatingpistongap.phi_rad * 180.0/PI;
	if(phi < 180.0)
	{
		TDC = temperaturepistongap.THP;
	}
	else
	{
		TDC = temperaturepistongap.TLP;
	}
	TCase = temperaturepistongap.TCase;

	//Source term
	bE = oilviscosity * ( pow(dvxz,2.0) + pow(dvyz,2.0) ) * dx * dy * dz2;
	//Diffusive coefficients
	Dx = lambda*dy*dz2 / dx;
	Dy = lambda*dx*dz2 / dy;
	Dz = lambda*dx*dy / dz2;
	//Convective coefficients
	Fx = oildensity * cp * vx * dy * dz2;
	Fy = oildensity * cp * vy * dx * dz2;
	//Peclet number
	Px = Fx / Dx;
	Py = Fy / Dy;
	//Power Law Scheme (Patankar Book)
	Ax = pow((1 - 0.1*fabs(Px)),5.0);
	Ax = max(Ax,0.0);
	Ay = pow((1 - 0.1*fabs(Py)),5.0);
	Ay = max(Ay,0.0);
	//Spatial coefficients
	anE = Ay * Dy + max(-1.0 * Fy,0.0); 
    asE = Ay * Dy + max(Fy,0.0);
    aeE = Ax * Dx + max(-1.0 * Fx,0.0); 
    awE = Ax * Dx + max(Fx,0.0); 
	atE = Dz;
	abE = Dz;
	apE = anE + asE + aeE + awE + atE + abE;
	
	//SOR Loop
	Tnew = T;
	alpha = 1.3;
	int iterations = 0;
	int iter_max = 10000;
	double R = 1.0;
	double Tn,Ts,Te,Tw,Tt,Tb;
	do
	{

		int q=0;
		//SOR
		for(int i=0;i<N*M*Q;i++)
		{
			if(i%(N*M)==0)
			{
				q++;
			}
			//North
			if(i%M==(M-1) && i>0)
			{
				Tn = TCase;
			}
			else
			{
				Tn = Tnew(i+1);
			}
			//South
			if(i%M==0)
			{
				Ts = TDC;
			}
			else
			{
				Ts = Tnew(i-1);
			}

			//East
			if(i>=q*(N-1)*M)
			{
				Te = Tnew(i-(N-1)*M);
			}
			else
			{
				Te = Tnew(i+M);
			}
			//West
			if(i<M+(q-1)*N*M)
			{
				Tw = Tnew(i+(N-1)*M);
			}
			else
			{
				Tw = Tnew(i-M);
			}

			//Top
			if(i>=N*M*(Q-1))
			{
				Tt = TK_surf_gap(i-N*M*(Q-1));
			}
			else
			{
				Tt = Tnew(i+N*M);
			}
			//Bottom
			if(i<N*M)
			{
				Tb = TB_surf_gap(i);
			}
			else
			{
				Tb = Tnew(i-N*M);
			}
			//SOR
			Tnew(i) += alpha * ( ( ( anE(i) * Tn + asE(i) * Ts + aeE(i) * Te + awE(i) * Tw + atE(i) * Tt + abE(i) * Tb + bE(i) ) / apE(i) ) - Tnew(i) );
			//cout<<"bE("<<i<<") "<<bE(i)<<"   Tnew("<<i<<") "<<Tnew(i)<<"\n";
		};

		//Counter
		iterations++;

		//Residual GS
		PistonEnergyCalcResidualGS(R);

	}while(R > Rmin_E && iterations < iter_max);


	//Log
	if(iterations>=iter_max)
	{
		fout.open("./output/piston/matlab/PistonGapLoopLog.txt",ios::app);
		fout << "Max Number of Iterations Reached in Energy Equation Loop" << "\n";
		fout << "Shaft Angle: " << "\t" << operatingpistongap.phi_rad*180/PI << "\n";
		fout << "Residual: " << "\t" << R << "\n";
		fout << " " << "\n";
		fout.close();
		fout.clear();
	};

	//Assign temperature
	T = Tnew;

};
//Calculate Energy residual for GS SOR method
void CPistonGap::PistonEnergyCalcResidualGS(double &R)
{

	double Tn,Ts,Te,Tw,Tt,Tb,TDC,TCase,phi;

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

	//Calculate residual field
	r_T = 0.0;
	int q=0;
	for(int i=0;i<N*M*Q;i++)
	{
		if(i%(N*M)==0)
		{
			q++;
		}
		//North
		if(i%M==(M-1) && i>0)
		{
			Tn = TCase;
		}
		else
		{
			Tn = Tnew(i+1);
		}
		//South
		if(i%M==0)
		{
			Ts = TDC;
		}
		else
		{
			Ts = Tnew(i-1);
		}

		//East
		if(i>=q*(N-1)*M)
		{
			Te = Tnew(i-(N-1)*M);
		}
		else
		{
			Te = Tnew(i+M);
		}
		//West
		if(i<M+(q-1)*N*M)
		{
			Tw = Tnew(i+(N-1)*M);
		}
		else
		{
			Tw = Tnew(i-M);
		}

		//Top
		if(i>=N*M*(Q-1))
		{
			Tt = TK_surf_gap(i-N*M*(Q-1));
		}
		else
		{
			Tt = Tnew(i+N*M);
		}
		//Bottom
		if(i<N*M)
		{
			Tb = TB_surf_gap(i);
		}
		else
		{
			Tb = Tnew(i-N*M);
		}
		//SOR
		r_T(i) = bE(i) - ( apE(i) * Tnew(i) - anE(i) * Tn - asE(i) * Ts - aeE(i) * Te - awE(i) * Tw - atE(i) * Tt - abE(i) * Tb );
	};

	//Calculate scaled residual
	R = 0.0;
	double n = sum(r_T);
	double d = sum(apE*Tnew);
	//double n = sqrt(sum(pow(r_T,2.0)));
	//double d = sqrt(sum(pow(bE,2.0)));
	R = n/d;

};
