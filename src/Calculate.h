#ifndef calculate_H_
#define calculate_H_


using namespace std;


/////////////////////////                 method1 /////////////////////////////////
int cal_RR_MA ( vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  double &  CalResult, statementVar & Var )
{
	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;
	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++ ;
	}
	
	Var.tmpAA=(Var.DDE[1][1])+(Var.DDE[1][0]);
	if (Var.tmpAA==0)
	{
		return 0 ;
	}
	if ( (Var.DDE[1][1]+Var.DDE[0][1])==0)
	{
		return 0 ;
	}

	Var.ALL_count=Var.DDE[0][0]+Var.DDE[0][1]+Var.tmpAA;
	Var.probHaps[0]=((Var.DDE[0][0])/Var.ALL_count);
	Var.probHaps[1]=((Var.DDE[0][1])/Var.ALL_count);
	Var.probHaps[2]=((Var.DDE[1][0])/Var.ALL_count);

	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];

	Var.Cal_B=(Var.pA1)*(Var.pA2);
	Var.Cal_A  = 1.0-(Var.pA1+Var.pA2)+Var.Cal_B ;

	if  (Var.Cal_A==0  || Var.Cal_B==0 )
	{
		if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
		if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
		if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}

		Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
		Var.pA2 = Var.probHaps[0]+Var.probHaps[2];

		Var.Cal_B = (Var.pA1)*(Var.pA2);
	    Var.Cal_A  = 1.0-(Var.pA1+Var.pA2)+Var.Cal_B ;

	}

	Var.D_A = Var.probHaps[0]-Var.Cal_B ;
	CalResult = (Var.D_A*Var.D_A)/(Var.Cal_A*Var.Cal_B);

	return 1;
}










 int cal_RR_D_MA(vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  PairInfoV1 &  CalResult, statementVar & Var)
{

	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;

	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++;
	}

	Var.tmpAA=Var.DDE[1][1]+Var.DDE[1][0];

	if (Var.tmpAA==0)
	{
		return 0 ;
	}
	if ( (Var.DDE[1][1]+Var.DDE[0][1])==0)
	{
		return 0 ;
	}

	Var.ALL_count=Var.DDE[0][0]+Var.DDE[0][1]+Var.tmpAA;
	Var.probHaps[0]=((Var.DDE[0][0])/Var.ALL_count);
	Var.probHaps[1]=((Var.DDE[0][1])/Var.ALL_count);
	Var.probHaps[2]=((Var.DDE[1][0])/Var.ALL_count);

	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pB1 = 1.0-Var.pA1;
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
	Var.pB2 = 1.0-Var.pA2;

	Var.XpA1_pA2=Var.pA1*Var.pA2;
	Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

	if  (Var.D_A>0)
	{
		Var.Cal_A=Var.pB1*Var.pA2;
		Var.Cal_B=Var.pA1*Var.pB2; 
	}
	else 
	{
		Var.D_A = 0.0-Var.D_A;
		Var.Cal_A = (Var.pB1)*(Var.pB2);
		Var.Cal_B = Var.XpA1_pA2 ;
	}

	Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{		Var.D_max=Var.Cal_B;	}

	if  (Var.D_max==0)
	{
		if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
		if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
		if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}

		Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
		Var.pB1 = 1.0-Var.pA1;
		Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
		Var.pB2 = 1.0-Var.pA2;

		Var.XpA1_pA2=Var.pA1*Var.pA2;
		Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

		if  (Var.D_A>0)
		{
			Var.Cal_A=Var.pB1*Var.pA2;
			Var.Cal_B=Var.pA1*Var.pB2; 
		}
		else 
		{
			Var.D_A = 0.0-Var.D_A;
			Var.Cal_A = (Var.pB1)*(Var.pB2);
			Var.Cal_B = Var.XpA1_pA2 ;
		}
		Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{		Var.D_max=Var.Cal_B;	}
	}

	CalResult.D = Var.D_A/Var.D_max;
	CalResult.RR = (Var.D_A/Var.Cal_A)*(Var.D_A/Var.Cal_B);

	return 1;
}



int cal_RR_D2_MA( vector<BaseType> & Base1 , vector<BaseType>  & Base2 , PairInfoV2 &  CalResult , statementVar & Var )
{

	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;
	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++;
	}
	Var.known[0]=Var.DDE[0][0];
	Var.known[1]=Var.DDE[0][1];
	Var.known[2]=Var.DDE[1][0];
	Var.known[3]=Var.DDE[1][1];


	Var.tmpAA=Var.known[3]+Var.known[2];
	if (  Var.tmpAA==0 )
	{
		return 0 ;
	}
	else if ( (Var.known[3]+Var.known[1])==0)
	{
		return 0 ;
	}

	Var.ALL_count=Var.known[0]+Var.known[1]+Var.tmpAA;

	Var.probHaps[0]=(Var.known[0])/Var.ALL_count;
	Var.probHaps[1]=(Var.known[1])/Var.ALL_count;
	Var.probHaps[2]=(Var.known[2])/Var.ALL_count;
	Var.probHaps[3]=1-Var.probHaps[0]-Var.probHaps[1]-Var.probHaps[2];

	if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
	if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
	if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}
	if (Var.probHaps[3] < 1e-10) { Var.probHaps[3]=1e-10;}


	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pB1 = 1.0-Var.pA1;
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
	Var.pB2 = 1.0-Var.pA2;



	Var.XpA1_pA2=Var.pA1*Var.pA2;
	Var.XpA1_pB2=Var.pA1*Var.pB2;
	Var.XpB1_pA2=Var.pB1*Var.pA2;
	Var.XpB1_pB2=Var.pB1*Var.pB2;

	Var.loglike1 =( Var.known[0]*log(Var.probHaps[0]) + Var.known[1]*log(Var.probHaps[1]) + Var.known[2]*log(Var.probHaps[2]) + Var.known[3]*log(Var.probHaps[3]))/Var.LN10;
	Var.loglike0 =( Var.known[0]*log(Var.XpA1_pA2) + Var.known[1]*log(Var.XpA1_pB2) + Var.known[2]*log(Var.XpB1_pA2) + Var.known[3]*log(Var.XpB1_pB2))/Var.LN10;

	Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

	if (Var.D_A < 0 ) 
	{
		Var.D_A = 0.0-Var.D_A;
		Var.Cal_A = Var.XpB1_pB2;
		Var.Cal_B = Var.XpA1_pA2 ;
	}
	else
	{
		Var.Cal_A = Var.XpB1_pA2;
		Var.Cal_B = Var.XpA1_pB2;
	}

	//Var.D_max=min(Var.Cal_A,Var.Cal_B);
	Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{		Var.D_max=Var.Cal_B;	}
	CalResult.D = Var.D_A/Var.D_max;
	CalResult.RR = (Var.D_A/Var.Cal_A)*(Var.D_A/Var.Cal_B);
	CalResult.LOD=(Var.loglike1-Var.loglike0);

	return 1;
}










 int cal_RR_D3_MA( vector<BaseType> & Base1  , vector<BaseType>  & Base2   ,  PairInfoV3 &  CalResult , statementVar & Var  )
{

	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;
	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++;
	}
	Var.known[0]=Var.DDE[0][0];
	Var.known[1]=Var.DDE[0][1];
	Var.known[2]=Var.DDE[1][0];
	Var.known[3]=Var.DDE[1][1];

	Var.tmpAA=Var.known[3]+Var.known[2];
	if (  Var.tmpAA==0 )
	{
		return 0 ;
	}
	else if ( (Var.known[3]+Var.known[1])==0)
	{
		return 0 ;
	}

	Var.ALL_count=Var.known[0]+Var.known[1]+Var.tmpAA;

	Var.probHaps[0]=(Var.known[0])/Var.ALL_count;
	Var.probHaps[1]=(Var.known[1])/Var.ALL_count;
	Var.probHaps[2]=(Var.known[2])/Var.ALL_count;
	Var.probHaps[3]=1-Var.probHaps[0]-Var.probHaps[1]-Var.probHaps[2];



	if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
	if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
	if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}
	if (Var.probHaps[3] < 1e-10) { Var.probHaps[3]=1e-10;}


	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pB1 = 1.0-Var.pA1;
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
	Var.pB2 = 1.0-Var.pA2;


	Var.XpA1_pA2=Var.pA1*Var.pA2;
	Var.XpA1_pB2=Var.pA1*Var.pB2;
	Var.XpB1_pA2=Var.pB1*Var.pA2;
	Var.XpB1_pB2=Var.pB1*Var.pB2;
	//	if (XpA1_pA2<1e-10 ){XpA1_pA2=1e-10;} 	if (XpA1_pB2<1e-10 ){XpA1_pB2=1e-10;}	if (XpB1_pA2<1e-10 ){XpB1_pA2=1e-10;} 	if (XpB1_pB2<1e-10 ){XpB1_pB2=1e-10;}

	Var.loglike1 =( Var.known[0]*log(Var.probHaps[0]) + Var.known[1]*log(Var.probHaps[1]) + Var.known[2]*log(Var.probHaps[2]) + Var.known[3]*log(Var.probHaps[3]))/Var.LN10;
	Var.loglike0 =( Var.known[0]*log(Var.XpA1_pA2) + Var.known[1]*log(Var.XpA1_pB2) + Var.known[2]*log(Var.XpB1_pA2) + Var.known[3]*log(Var.XpB1_pB2))/Var.LN10;

	Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

	if (Var.D_A< 0 ) 
	{ //  ALL_count  is tmp
		Var.ALL_count=Var.probHaps[0]; Var.probHaps[0]=Var.probHaps[1]; Var.probHaps[1]=Var.ALL_count;
		Var.ALL_count=Var.probHaps[3]; Var.probHaps[3]=Var.probHaps[2]; Var.probHaps[2]=Var.ALL_count;
		Var.pA2 = Var.pA2 + Var.pB2;		Var.pB2 = Var.pA2 - Var.pB2;		Var.pA2 = Var.pA2 - Var.pB2;
		Var.D_A = 0.0-Var.D_A;
		Var.ALL_count=Var.known[0]; Var.known[0]=Var.known[1]; Var.known[1]=Var.ALL_count;
		Var.ALL_count=Var.known[3]; Var.known[3]=Var.known[2]; Var.known[2]=Var.ALL_count;
		Var.Cal_A = (Var.pA2)*(Var.pB1);
		Var.Cal_B = (Var.pA1)*(Var.pB2);
	}
	else
	{
		Var.Cal_A = Var.XpB1_pA2;
		Var.Cal_B = Var.XpA1_pB2;
	}

	//Var.D_max=min(Var.Cal_A,Var.Cal_B);
	Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{		Var.D_max=Var.Cal_B;	}

	CalResult.D = Var.D_A/Var.D_max;
	CalResult.RR = (Var.D_A/Var.Cal_A)*(Var.D_A/Var.Cal_B);
	CalResult.LOD=(Var.loglike1-Var.loglike0);


	Var.XpA1_pA2=Var.pA1*Var.pA2;
	for (Var.i=0; Var.i<100; Var.i++)
	{
		Var.dpr = (double)Var.i*0.01;
		Var.tmpAA = Var.dpr*Var.D_max + Var.XpA1_pA2 ;
		Var.tmpAB = Var.pA1-Var.tmpAA;
		Var.tmpBA = Var.pA2-Var.tmpAA;
		Var.tmpBB = Var.pB1-Var.tmpBA;
		Var.lsurface[Var.i] = (Var.known[0]*log(Var.tmpAA) + Var.known[1]*log(Var.tmpAB) + Var.known[2]*log(Var.tmpBA) + Var.known[3]*log(Var.tmpBB))/Var.LN10;
	}
	//   i=100;
	Var.dpr = (double)100*0.01;
	Var.tmpAA = Var.dpr*Var.D_max + Var.XpA1_pA2;
	Var.tmpAB = Var.pA1-Var.tmpAA;
	Var.tmpBA = Var.pA2-Var.tmpAA;
	Var.tmpBB = Var.pB1-Var.tmpBA;
	/* one value will be 0 */
	if (Var.tmpAA < 1e-10) { Var.tmpAA=1e-10;}
	if (Var.tmpAB < 1e-10) { Var.tmpAB=1e-10;}
	if (Var.tmpBA < 1e-10) { Var.tmpBA=1e-10;}
	if (Var.tmpBB < 1e-10) { Var.tmpBB=1e-10;}
	Var.lsurface[100] = (Var.known[0]*log(Var.tmpAA) + Var.known[1]*log(Var.tmpAB) + Var.known[2]*log(Var.tmpBA) + Var.known[3]*log(Var.tmpBB))/Var.LN10;

	Var.total_prob=0.0;
	Var.sum_prob=0.0;

	for (Var.i=0; Var.i<=100; Var.i++) {
		Var.lsurface[Var.i] -= Var.loglike1;
		Var.lsurface[Var.i] = pow(10.0,Var.lsurface[Var.i]);
		Var.total_prob += Var.lsurface[Var.i];
	}

	Var.cut5off=Var.total_prob*0.05;

	for (Var.i=0; Var.i<=100; Var.i++) 
	{
		Var.sum_prob += Var.lsurface[Var.i];
		if (Var.sum_prob > Var.cut5off &&	Var.sum_prob-Var.lsurface[Var.i] < Var.cut5off ) 
		{
			Var.low_i = Var.i-1;
			break;
		}
	}

	Var.sum_prob=0.0;
	for (Var.i=100; Var.i>=0; (Var.i)--)
	{
		Var.sum_prob += Var.lsurface[Var.i];
		if (Var.sum_prob > Var.cut5off &&	Var.sum_prob-Var.lsurface[Var.i] < Var.cut5off ) 
		{
			Var.high_i = Var.i+1;
			break;
		}
	}


	if (Var.high_i > 100){ Var.high_i = 100; }

	CalResult.low_i=Var.low_i;
	CalResult.high_i=Var.high_i;

	return 1;

}

/////////////////////////////////////////////////////////////////////////////////////////////////





/////////////////////////                    ethod1 /////////////////////////////////
int cal_RR_MB(vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  PairInfoV1 &  CalResult, statementVar & Var )
{
	
	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;
	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++ ;
	}

	Var.tmpAA=(Var.DDE[1][1])+(Var.DDE[1][0]);

	if (Var.tmpAA==0)
	{
		CalResult.D=-1;
		return 0 ;
	}
	if ( (Var.DDE[1][1]+Var.DDE[0][1])==0)
	{
		CalResult.D=-1;
		return 0 ;
	}

	Var.ALL_count=Var.DDE[0][0]+Var.DDE[0][1]+Var.tmpAA;
	Var.probHaps[0]=((Var.DDE[0][0])/Var.ALL_count);
	Var.probHaps[1]=((Var.DDE[0][1])/Var.ALL_count);
	Var.probHaps[2]=((Var.DDE[1][0])/Var.ALL_count);

	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];

	Var.Cal_B = (Var.pA1)*(Var.pA2);
	Var.Cal_A = 1- Var.pA1-Var.pA2+Var.Cal_B;

	if  (Var.Cal_A==0  || Var.Cal_B==0 )
	{
		if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
		if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
		if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}

		Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
		Var.pA2 = Var.probHaps[0]+Var.probHaps[2];

		Var.Cal_B = (Var.pA1)*(Var.pA2);
		Var.Cal_A = 1- Var.pA1-Var.pA2+Var.Cal_B;
	}
	CalResult.D=1;
	Var.D_A = Var.probHaps[0]-Var.Cal_B ;
	CalResult.RR = (Var.D_A/Var.Cal_A)*(Var.D_A/Var.Cal_B);

	return 1;
}












int cal_RR_D_MB(vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  PairInfoV1 &  CalResult, statementVar & Var)
{

	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;

	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++;
	}

	Var.tmpAA=Var.DDE[1][1]+Var.DDE[1][0];

	if (Var.tmpAA==0)
	{
		CalResult.D=-1;
		return 0 ;
	}
	else if ( (Var.DDE[1][1]+Var.DDE[0][1])==0)
	{
		CalResult.D=-1;
		return 0 ;
	}

	Var.ALL_count=Var.DDE[0][0]+Var.DDE[0][1]+Var.tmpAA;
	Var.probHaps[0]=((Var.DDE[0][0])/Var.ALL_count);
	Var.probHaps[1]=((Var.DDE[0][1])/Var.ALL_count);
	Var.probHaps[2]=((Var.DDE[1][0])/Var.ALL_count);

	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pB1 = 1.0-Var.pA1;
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
	Var.pB2 = 1.0-Var.pA2;

	Var.XpA1_pA2=Var.pA1*Var.pA2;
	Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

	if  (Var.D_A>0)
	{
		Var.Cal_A=Var.pB1*Var.pA2;
		Var.Cal_B=Var.pA1*Var.pB2; 
	}
	else 
	{
		Var.D_A = 0.0-Var.D_A;
		Var.Cal_A = (Var.pB1)*(Var.pB2);
		Var.Cal_B = Var.XpA1_pA2 ;
	}

	Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{		Var.D_max=Var.Cal_B;	}

	if  (Var.D_max==0)
	{
		if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
		if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
		if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}

		Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
		Var.pB1 = 1.0-Var.pA1;
		Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
		Var.pB2 = 1.0-Var.pA2;

		Var.XpA1_pA2=Var.pA1*Var.pA2;
		Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

		if  (Var.D_A>0)
		{
			Var.Cal_A=Var.pB1*Var.pA2;
			Var.Cal_B=Var.pA1*Var.pB2; 
		}
		else 
		{
			Var.D_A = 0.0-Var.D_A;
			Var.Cal_A = (Var.pB1)*(Var.pB2);
			Var.Cal_B = Var.XpA1_pA2 ;
		}
		Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{		Var.D_max=Var.Cal_B;	}
	}

	CalResult.D = Var.D_A/Var.D_max;
	CalResult.RR = (Var.D_A/Var.Cal_A)*(Var.D_A/Var.Cal_B);

	return 1;
}



int cal_RR_D2_MB( vector<BaseType> & Base1 , vector<BaseType>  & Base2 , PairInfoV2 &  CalResult , statementVar & Var )
{

	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;
	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++;
	}
	Var.known[0]=Var.DDE[0][0];
	Var.known[1]=Var.DDE[0][1];
	Var.known[2]=Var.DDE[1][0];
	Var.known[3]=Var.DDE[1][1];


	Var.tmpAA=Var.known[3]+Var.known[2];
	if (  Var.tmpAA==0 )
	{
		CalResult.D=-1;
		return 0 ;
	}
	else if ( (Var.known[3]+Var.known[1])==0)
	{
		CalResult.D=-1;
		return 0 ;
	}

	Var.ALL_count=Var.known[0]+Var.known[1]+Var.tmpAA;

	Var.probHaps[0]=(Var.known[0])/Var.ALL_count;
	Var.probHaps[1]=(Var.known[1])/Var.ALL_count;
	Var.probHaps[2]=(Var.known[2])/Var.ALL_count;
	Var.probHaps[3]=1-Var.probHaps[0]-Var.probHaps[1]-Var.probHaps[2];

	if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
	if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
	if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}
	if (Var.probHaps[3] < 1e-10) { Var.probHaps[3]=1e-10;}


	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pB1 = 1.0-Var.pA1;
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
	Var.pB2 = 1.0-Var.pA2;



	Var.XpA1_pA2=Var.pA1*Var.pA2;
	Var.XpA1_pB2=Var.pA1*Var.pB2;
	Var.XpB1_pA2=Var.pB1*Var.pA2;
	Var.XpB1_pB2=Var.pB1*Var.pB2;

	Var.loglike1 =( Var.known[0]*log(Var.probHaps[0]) + Var.known[1]*log(Var.probHaps[1]) + Var.known[2]*log(Var.probHaps[2]) + Var.known[3]*log(Var.probHaps[3]))/Var.LN10;
	Var.loglike0 =( Var.known[0]*log(Var.XpA1_pA2) + Var.known[1]*log(Var.XpA1_pB2) + Var.known[2]*log(Var.XpB1_pA2) + Var.known[3]*log(Var.XpB1_pB2))/Var.LN10;

	Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

	if (Var.D_A < 0 ) 
	{
		Var.D_A = 0.0-Var.D_A;
		Var.Cal_A = Var.XpB1_pB2;
		Var.Cal_B = Var.XpA1_pA2 ;
	}
	else
	{
		Var.Cal_A = Var.XpB1_pA2;
		Var.Cal_B = Var.XpA1_pB2;
	}

	//Var.D_max=min(Var.Cal_A,Var.Cal_B);
	Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{		Var.D_max=Var.Cal_B;	}
	CalResult.D = Var.D_A/Var.D_max;
	CalResult.RR = (Var.D_A/Var.Cal_A)*(Var.D_A/Var.Cal_B);
	CalResult.LOD=(Var.loglike1-Var.loglike0);

	return 1;
}










int cal_RR_D3_MB( vector<BaseType> & Base1  , vector<BaseType>  & Base2   ,  PairInfoV3 &  CalResult , statementVar & Var  )
{

	Var.DDE[0][0]=0;	Var.DDE[0][1]=0;	Var.DDE[1][0]=0;	Var.DDE[1][1]=0;
	for (Var.i=0 ;Var.i<Var.Asize ; (Var.i)++)
	{
		Var.DDE[(Base1[Var.i].Value)][(Base2[Var.i].Value)]++;
	}
	Var.known[0]=Var.DDE[0][0];
	Var.known[1]=Var.DDE[0][1];
	Var.known[2]=Var.DDE[1][0];
	Var.known[3]=Var.DDE[1][1];

	Var.tmpAA=Var.known[3]+Var.known[2];
	if (  Var.tmpAA==0 )
	{
		CalResult.D=-1;
		return 0 ;
	}
	else if ( (Var.known[3]+Var.known[1])==0)
	{
		CalResult.D=-1;
		return 0 ;
	}

	Var.ALL_count=Var.known[0]+Var.known[1]+Var.tmpAA;

	Var.probHaps[0]=(Var.known[0])/Var.ALL_count;
	Var.probHaps[1]=(Var.known[1])/Var.ALL_count;
	Var.probHaps[2]=(Var.known[2])/Var.ALL_count;
	Var.probHaps[3]=1-Var.probHaps[0]-Var.probHaps[1]-Var.probHaps[2];



	if (Var.probHaps[0] < 1e-10) { Var.probHaps[0]=1e-10;}
	if (Var.probHaps[1] < 1e-10) { Var.probHaps[1]=1e-10;}
	if (Var.probHaps[2] < 1e-10) { Var.probHaps[2]=1e-10;}
	if (Var.probHaps[3] < 1e-10) { Var.probHaps[3]=1e-10;}


	Var.pA1 = Var.probHaps[0]+Var.probHaps[1];
	Var.pB1 = 1.0-Var.pA1;
	Var.pA2 = Var.probHaps[0]+Var.probHaps[2];
	Var.pB2 = 1.0-Var.pA2;


	Var.XpA1_pA2=Var.pA1*Var.pA2;
	Var.XpA1_pB2=Var.pA1*Var.pB2;
	Var.XpB1_pA2=Var.pB1*Var.pA2;
	Var.XpB1_pB2=Var.pB1*Var.pB2;
	//	if (XpA1_pA2<1e-10 ){XpA1_pA2=1e-10;} 	if (XpA1_pB2<1e-10 ){XpA1_pB2=1e-10;}	if (XpB1_pA2<1e-10 ){XpB1_pA2=1e-10;} 	if (XpB1_pB2<1e-10 ){XpB1_pB2=1e-10;}

	Var.loglike1 =( Var.known[0]*log(Var.probHaps[0]) + Var.known[1]*log(Var.probHaps[1]) + Var.known[2]*log(Var.probHaps[2]) + Var.known[3]*log(Var.probHaps[3]))/Var.LN10;
	Var.loglike0 =( Var.known[0]*log(Var.XpA1_pA2) + Var.known[1]*log(Var.XpA1_pB2) + Var.known[2]*log(Var.XpB1_pA2) + Var.known[3]*log(Var.XpB1_pB2))/Var.LN10;

	Var.D_A = Var.probHaps[0]-Var.XpA1_pA2 ;

	if (Var.D_A< 0 ) 
	{ //  ALL_count  is tmp
		Var.ALL_count=Var.probHaps[0]; Var.probHaps[0]=Var.probHaps[1]; Var.probHaps[1]=Var.ALL_count;
		Var.ALL_count=Var.probHaps[3]; Var.probHaps[3]=Var.probHaps[2]; Var.probHaps[2]=Var.ALL_count;
		Var.pA2 = Var.pA2 + Var.pB2;		Var.pB2 = Var.pA2 - Var.pB2;		Var.pA2 = Var.pA2 - Var.pB2;
		Var.D_A = 0.0-Var.D_A;
		Var.ALL_count=Var.known[0]; Var.known[0]=Var.known[1]; Var.known[1]=Var.ALL_count;
		Var.ALL_count=Var.known[3]; Var.known[3]=Var.known[2]; Var.known[2]=Var.ALL_count;
		Var.Cal_A = (Var.pA2)*(Var.pB1);
		Var.Cal_B = (Var.pA1)*(Var.pB2);
	}
	else
	{
		Var.Cal_A = Var.XpB1_pA2;
		Var.Cal_B = Var.XpA1_pB2;
	}

	//Var.D_max=min(Var.Cal_A,Var.Cal_B);
	Var.D_max=Var.Cal_A;	if  (Var.Cal_A>Var.Cal_B)	{		Var.D_max=Var.Cal_B;	}

	CalResult.D = Var.D_A/Var.D_max;
	CalResult.RR = (Var.D_A/Var.Cal_A)*(Var.D_A/Var.Cal_B);
	CalResult.LOD=(Var.loglike1-Var.loglike0);


	Var.XpA1_pA2=Var.pA1*Var.pA2;
	for (Var.i=0; Var.i<100; Var.i++)
	{
		Var.dpr = (double)Var.i*0.01;
		Var.tmpAA = Var.dpr*Var.D_max + Var.XpA1_pA2 ;
		Var.tmpAB = Var.pA1-Var.tmpAA;
		Var.tmpBA = Var.pA2-Var.tmpAA;
		Var.tmpBB = Var.pB1-Var.tmpBA;
		Var.lsurface[Var.i] = (Var.known[0]*log(Var.tmpAA) + Var.known[1]*log(Var.tmpAB) + Var.known[2]*log(Var.tmpBA) + Var.known[3]*log(Var.tmpBB))/Var.LN10;
	}
	//   i=100;
	Var.dpr = (double)100*0.01;
	Var.tmpAA = Var.dpr*Var.D_max + Var.XpA1_pA2;
	Var.tmpAB = Var.pA1-Var.tmpAA;
	Var.tmpBA = Var.pA2-Var.tmpAA;
	Var.tmpBB = Var.pB1-Var.tmpBA;
	/* one value will be 0 */
	if (Var.tmpAA < 1e-10) { Var.tmpAA=1e-10;}
	if (Var.tmpAB < 1e-10) { Var.tmpAB=1e-10;}
	if (Var.tmpBA < 1e-10) { Var.tmpBA=1e-10;}
	if (Var.tmpBB < 1e-10) { Var.tmpBB=1e-10;}
	Var.lsurface[100] = (Var.known[0]*log(Var.tmpAA) + Var.known[1]*log(Var.tmpAB) + Var.known[2]*log(Var.tmpBA) + Var.known[3]*log(Var.tmpBB))/Var.LN10;

	Var.total_prob=0.0;
	Var.sum_prob=0.0;

	for (Var.i=0; Var.i<=100; Var.i++) {
		Var.lsurface[Var.i] -= Var.loglike1;
		Var.lsurface[Var.i] = pow(10.0,Var.lsurface[Var.i]);
		Var.total_prob += Var.lsurface[Var.i];
	}

	Var.cut5off=Var.total_prob*0.05;

	for (Var.i=0; Var.i<=100; Var.i++) 
	{
		Var.sum_prob += Var.lsurface[Var.i];
		if (Var.sum_prob > Var.cut5off &&	Var.sum_prob-Var.lsurface[Var.i] < Var.cut5off ) 
		{
			Var.low_i = Var.i-1;
			break;
		}
	}

	Var.sum_prob=0.0;
	for (Var.i=100; Var.i>=0; (Var.i)--)
	{
		Var.sum_prob += Var.lsurface[Var.i];
		if (Var.sum_prob > Var.cut5off &&	Var.sum_prob-Var.lsurface[Var.i] < Var.cut5off ) 
		{
			Var.high_i = Var.i+1;
			break;
		}
	}


	if (Var.high_i > 100){ Var.high_i = 100; }

	CalResult.low_i=Var.low_i;
	CalResult.high_i=Var.high_i;

	return 1;

}

/////////////////////////////////////////////////////////////////////////////////////////////////













#endif // calculate_H_  ;


