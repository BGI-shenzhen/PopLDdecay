#ifndef calculate_H_
#define calculate_H_

using namespace std;

int cal_RR_MA ( vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  double &  CalResult, statementVar & Var )
{
	memset(Var.DDE, 0, sizeof(Var.DDE));
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

	memset(Var.DDE, 0, sizeof(Var.DDE));

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

	memset(Var.DDE, 0, sizeof(Var.DDE));
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

	memset(Var.DDE, 0, sizeof(Var.DDE));
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

int cal_RR_MB(vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  PairInfoV1 &  CalResult, statementVar & Var )
{
	unsigned short int DDE[3][3] = {{0}};
	int Asize = Var.Asize;
	for (int _i = 0; _i < Asize; _i++) {
		DDE[(Base1[_i].Value)][(Base2[_i].Value)]++;
	}
	double tmpAA = DDE[1][1] + DDE[1][0];
	if (tmpAA == 0) { CalResult.D = -1; return 0; }
	if ((DDE[1][1] + DDE[0][1]) == 0) { CalResult.D = -1; return 0; }
	double ALL_count = DDE[0][0] + DDE[0][1] + tmpAA;
	double prob0 = (DDE[0][0]) / ALL_count;
	double prob1 = (DDE[0][1]) / ALL_count;
	double prob2 = (DDE[1][0]) / ALL_count;
	double pA1 = prob0 + prob1;
	double pA2 = prob0 + prob2;
	double Cal_B = pA1 * pA2;
	double Cal_A = 1.0 - pA1 - pA2 + Cal_B;
	if (Cal_A == 0 || Cal_B == 0) {
		if (prob0 < 1e-10) prob0 = 1e-10;
		if (prob1 < 1e-10) prob1 = 1e-10;
		if (prob2 < 1e-10) prob2 = 1e-10;
		pA1 = prob0 + prob1; pA2 = prob0 + prob2;
		Cal_B = pA1 * pA2; Cal_A = 1.0 - pA1 - pA2 + Cal_B;
	}
	CalResult.D = 1;
	double D_A = prob0 - Cal_B;
	CalResult.RR = (D_A / Cal_A) * (D_A / Cal_B);
	return 1;
}

int cal_RR_D_MB(vector<BaseType>  & Base1  , vector<BaseType> &  Base2 ,  PairInfoV1 &  CalResult, statementVar & Var)
{
	unsigned short int DDE[3][3] = {{0}};
	int Asize = Var.Asize;
	for (int _i = 0; _i < Asize; _i++) {
		DDE[(Base1[_i].Value)][(Base2[_i].Value)]++;
	}
	double tmpAA = DDE[1][1] + DDE[1][0];
	if (tmpAA == 0) { CalResult.D = -1; return 0; }
	if ((DDE[1][1] + DDE[0][1]) == 0) { CalResult.D = -1; return 0; }
	double ALL_count = DDE[0][0] + DDE[0][1] + tmpAA;
	double prob0 = (DDE[0][0]) / ALL_count;
	double prob1 = (DDE[0][1]) / ALL_count;
	double prob2 = (DDE[1][0]) / ALL_count;
	double pA1 = prob0 + prob1;
	double pB1 = 1.0 - pA1;
	double pA2 = prob0 + prob2;
	double pB2 = 1.0 - pA2;
	double XpA1_pA2 = pA1 * pA2;
	double D_A = prob0 - XpA1_pA2;
	double Cal_A, Cal_B;
	if (D_A > 0) {
		Cal_A = pB1 * pA2; Cal_B = pA1 * pB2;
	} else {
		D_A = -D_A; Cal_A = pB1 * pB2; Cal_B = XpA1_pA2;
	}
	double D_max = Cal_A;
	if (Cal_A > Cal_B) D_max = Cal_B;
	if (D_max == 0) {
		if (prob0 < 1e-10) prob0 = 1e-10;
		if (prob1 < 1e-10) prob1 = 1e-10;
		if (prob2 < 1e-10) prob2 = 1e-10;
		pA1 = prob0 + prob1; pB1 = 1.0 - pA1;
		pA2 = prob0 + prob2; pB2 = 1.0 - pA2;
		XpA1_pA2 = pA1 * pA2;
		D_A = prob0 - XpA1_pA2;
		if (D_A > 0) {
			Cal_A = pB1 * pA2; Cal_B = pA1 * pB2;
		} else {
			D_A = -D_A; Cal_A = pB1 * pB2; Cal_B = XpA1_pA2;
		}
		D_max = Cal_A;
		if (Cal_A > Cal_B) D_max = Cal_B;
	}
	CalResult.D = D_A / D_max;
	CalResult.RR = (D_A / Cal_A) * (D_A / Cal_B);
	return 1;
}

int cal_RR_D2_MB( vector<BaseType> & Base1 , vector<BaseType>  & Base2 , PairInfoV2 &  CalResult , statementVar & Var )
{
	const double LN10 = 2.3025850929940456840;
	unsigned short int DDE[3][3] = {{0}};
	int Asize = Var.Asize;
	for (int _i = 0; _i < Asize; _i++) {
		DDE[(Base1[_i].Value)][(Base2[_i].Value)]++;
	}
	int known0 = DDE[0][0], known1 = DDE[0][1];
	int known2 = DDE[1][0], known3 = DDE[1][1];
	double tmpAA = known3 + known2;
	if (tmpAA == 0) { CalResult.D = -1; return 0; }
	if ((known3 + known1) == 0) { CalResult.D = -1; return 0; }
	double ALL_count = known0 + known1 + tmpAA;
	double prob0 = (known0) / ALL_count;
	double prob1 = (known1) / ALL_count;
	double prob2 = (known2) / ALL_count;
	double prob3 = 1.0 - prob0 - prob1 - prob2;
	if (prob0 < 1e-10) prob0 = 1e-10;
	if (prob1 < 1e-10) prob1 = 1e-10;
	if (prob2 < 1e-10) prob2 = 1e-10;
	if (prob3 < 1e-10) prob3 = 1e-10;
	double pA1 = prob0 + prob1;
	double pB1 = 1.0 - pA1;
	double pA2 = prob0 + prob2;
	double pB2 = 1.0 - pA2;
	double XpA1_pA2 = pA1 * pA2;
	double XpA1_pB2 = pA1 * pB2;
	double XpB1_pA2 = pB1 * pA2;
	double XpB1_pB2 = pB1 * pB2;
	double loglike1 = (known0 * log(prob0) + known1 * log(prob1) + known2 * log(prob2) + known3 * log(prob3)) / LN10;
	double loglike0 = (known0 * log(XpA1_pA2) + known1 * log(XpA1_pB2) + known2 * log(XpB1_pA2) + known3 * log(XpB1_pB2)) / LN10;
	double D_A = prob0 - XpA1_pA2;
	double Cal_A, Cal_B;
	if (D_A < 0) {
		D_A = -D_A; Cal_A = XpB1_pB2; Cal_B = XpA1_pA2;
	} else {
		Cal_A = XpB1_pA2; Cal_B = XpA1_pB2;
	}
	double D_max = Cal_A;
	if (Cal_A > Cal_B) D_max = Cal_B;
	CalResult.D = D_A / D_max;
	CalResult.RR = (D_A / Cal_A) * (D_A / Cal_B);
	CalResult.LOD = loglike1 - loglike0;
	return 1;
}

int cal_RR_D3_MB( vector<BaseType> & Base1  , vector<BaseType>  & Base2   ,  PairInfoV3 &  CalResult , statementVar & Var  )
{
	const double LN10 = 2.3025850929940456840;
	unsigned short int DDE[3][3] = {{0}};
	int Asize = Var.Asize;
	for (int _i = 0; _i < Asize; _i++) {
		DDE[(Base1[_i].Value)][(Base2[_i].Value)]++;
	}
	int k0 = DDE[0][0], k1 = DDE[0][1], k2 = DDE[1][0], k3 = DDE[1][1];
	double tmpAA = k3 + k2;
	if (tmpAA == 0) { CalResult.D = -1; return 0; }
	if ((k3 + k1) == 0) { CalResult.D = -1; return 0; }
	double ALL_count = k0 + k1 + tmpAA;
	double pHap0 = (k0) / ALL_count;
	double pHap1 = (k1) / ALL_count;
	double pHap2 = (k2) / ALL_count;
	double pHap3 = 1.0 - pHap0 - pHap1 - pHap2;
	if (pHap0 < 1e-10) pHap0 = 1e-10;
	if (pHap1 < 1e-10) pHap1 = 1e-10;
	if (pHap2 < 1e-10) pHap2 = 1e-10;
	if (pHap3 < 1e-10) pHap3 = 1e-10;
	double pA1 = pHap0 + pHap1;
	double pB1 = 1.0 - pA1;
	double pA2 = pHap0 + pHap2;
	double pB2 = 1.0 - pA2;
	double XpA1_pA2 = pA1 * pA2;
	double XpA1_pB2 = pA1 * pB2;
	double XpB1_pA2 = pB1 * pA2;
	double XpB1_pB2 = pB1 * pB2;
	double loglike1 = (k0 * log(pHap0) + k1 * log(pHap1) + k2 * log(pHap2) + k3 * log(pHap3)) / LN10;
	double loglike0 = (k0 * log(XpA1_pA2) + k1 * log(XpA1_pB2) + k2 * log(XpB1_pA2) + k3 * log(XpB1_pB2)) / LN10;
	double D_A = pHap0 - XpA1_pA2;
	double Cal_A, Cal_B;
	if (D_A < 0) {
		double t = pHap0; pHap0 = pHap1; pHap1 = t;
		t = pHap3; pHap3 = pHap2; pHap2 = t;
		t = pA2; pA2 = pB2; pB2 = t;
		D_A = -D_A;
		int ti = k0; k0 = k1; k1 = ti;
		ti = k3; k3 = k2; k2 = ti;
		Cal_A = pA2 * pB1;
		Cal_B = pA1 * pB2;
	} else {
		Cal_A = XpB1_pA2;
		Cal_B = XpA1_pB2;
	}
	double D_max = Cal_A;
	if (Cal_A > Cal_B) D_max = Cal_B;
	CalResult.D = D_A / D_max;
	CalResult.RR = (D_A / Cal_A) * (D_A / Cal_B);
	CalResult.LOD = loglike1 - loglike0;
	XpA1_pA2 = pA1 * pA2;
	double lsurface[101];
	for (int i = 0; i < 100; i++) {
		double dpr = (double)i * 0.01;
		double tAA = dpr * D_max + XpA1_pA2;
		double tAB = pA1 - tAA;
		double tBA = pA2 - tAA;
		double tBB = pB1 - tBA;
		lsurface[i] = (k0 * log(tAA) + k1 * log(tAB) + k2 * log(tBA) + k3 * log(tBB)) / LN10;
	}
	{
		double tAA = 1.0 * D_max + XpA1_pA2;
		double tAB = pA1 - tAA;
		double tBA = pA2 - tAA;
		double tBB = pB1 - tBA;
		if (tAA < 1e-10) tAA = 1e-10;
		if (tAB < 1e-10) tAB = 1e-10;
		if (tBA < 1e-10) tBA = 1e-10;
		if (tBB < 1e-10) tBB = 1e-10;
		lsurface[100] = (k0 * log(tAA) + k1 * log(tAB) + k2 * log(tBA) + k3 * log(tBB)) / LN10;
	}
	double total_prob = 0.0, sum_prob = 0.0;
	for (int i = 0; i <= 100; i++) {
		lsurface[i] -= loglike1;
		lsurface[i] = pow(10.0, lsurface[i]);
		total_prob += lsurface[i];
	}
	double cut5off = total_prob * 0.05;
	int low_i = 0;
	sum_prob = 0.0;
	for (int i = 0; i <= 100; i++) {
		sum_prob += lsurface[i];
		if (sum_prob > cut5off && sum_prob - lsurface[i] < cut5off) {
			low_i = i - 1;
			break;
		}
	}
	int high_i = 100;
	sum_prob = 0.0;
	for (int i = 100; i >= 0; i--) {
		sum_prob += lsurface[i];
		if (sum_prob > cut5off && sum_prob - lsurface[i] < cut5off) {
			high_i = i + 1;
			break;
		}
	}
	if (high_i > 100) high_i = 100;
	CalResult.low_i = low_i;
	CalResult.high_i = high_i;
	return 1;
}

#endif // calculate_H_  ;
