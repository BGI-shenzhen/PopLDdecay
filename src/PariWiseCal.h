#ifndef PariWiseCal_H_
#define PariWiseCal_H_

using namespace std;

int PairWiseComNewOUT_A( In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPList, int & Flag_for_pro )
{


	map <string,map<llong, vector <BaseType>  > > :: iterator key1_it ;

	map<llong, vector <BaseType> >  :: iterator key2 ;
	map<llong, vector <BaseType> >  :: iterator key2_se  ;
	llong Dis=0;
	statementVar  Var ;	

	key1_it=SNPList.begin(); 
	key2=(key1_it->second).begin();
	Var.Asize= (key2->second).size();

	cerr<<"##begin pair-wise R^2 cal/D'.\nAfter filter Remain SNP Number :\t"<<Flag_for_pro<<"\n##SampleSize*2 =\t"<<Var.Asize<<endl;

	string Stat=(paraFA04->InStr2)+".stat.gz";
	ogzstream OUT ((Stat).c_str());
	if((!OUT.good()))
	{
		cerr << "open OUT File error: "<<(paraFA04->InStr2)<<endl;
		return  0;
	}



	if ((paraFA04->TF) ==5) 
	{
		int MaxDis=(paraFA04->InInt)+1;
		int **RRNum = new int*[MaxDis];

		for(int i = 0; i <MaxDis ; i++)
		{
			RRNum[i] = new int [101];
			for (int j =0 ; j< 101 ; j++)
			{
				RRNum[i][j]=0;
			}
		}

		OUT <<"#Dist\tR^2\tR^2_count\n";
		OUT <<"#2SampleSize\t"<<Var.Asize<<"\n";


		double  CalResult ;
		int RRint ;
		for (key1_it=SNPList.begin(); key1_it!=SNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();

			for(   ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ;  key2_se++ ; 
				for(   ; key2_se!=(key1_it->second).end(); key2_se++)
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					if(cal_RR_MA( key2->second , key2_se->second  ,CalResult, Var) ==1)
					{
						RRint=int(CalResult*100);
						RRNum[Dis][RRint]++;
					}
				}
			}
		}





		for(int i = 0; i <(paraFA04->InInt) ; i++)
		{

			for (int j =0 ; j< 101 ; j++)
			{
				if (RRNum[i][j]==0)
				{
					continue;
				}
				double RR_double=j*1.0/100;
				OUT<<i+1<<"\t"<<setprecision(2)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<RR_double<<"\t"<<RRNum[i][j]<<"\n";
			}

			delete[]  RRNum[i] ;
		}
		delete[]  RRNum;
	}



	////////////////////////////////////////////////////////////
	else if ((paraFA04->TF) ==4) 
	{


		int MaxDis=(paraFA04->InInt)+1;
		int **RRNum = new int*[MaxDis];
		int **DNum = new int*[MaxDis];

		for(int i = 0; i <MaxDis ; i++)
		{
			RRNum[i] = new int [101];
			DNum[i] = new int [101];

			for (int j =0 ; j< 101 ; j++)
			{
				RRNum[i][j]=0;
				DNum[i][j]=0;
			}
		}

		OUT <<"#Dist\tR^2\tR^2_count\tD'\tD_count\n";
		OUT <<"#2SampleSize\t"<<Var.Asize<<"\n";


		PairInfoV1  CalResult ;
		int Dint=0;
		int RRint=0;

		for (key1_it=SNPList.begin(); key1_it!=SNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();

			for(  ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ;  key2_se++ ; 
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					else
					{
						if(	cal_RR_D_MA( key2->second , key2_se->second  ,CalResult, Var) ==1)
						{
							Dint=int(CalResult.D*100);
							RRint=int(CalResult.RR*100);
							RRNum[Dis][RRint]++;
							DNum[Dis][Dint]++;
						}
					}
				}
			}
		}





		for(int i = 0; i <(paraFA04->InInt) ; i++)
		{

			for (int j =0 ; j< 101 ; j++)
			{
				if (RRNum[i][j]==0 && DNum[i][j]==0 )
				{
					continue;
				}
				double RR_double=j*1.0/100;
				OUT<<i+1<<"\t"<<setprecision(2)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<RR_double<<"\t"<<RRNum[i][j]<<"\t"<<RR_double<<"\t"<<DNum[i][j]<<"\n";
			}

			delete[]  RRNum[i] ;
			delete[]  DNum[i] ;
		}
		delete[]  RRNum;
		delete[]  DNum;
	}



	OUT.close();
	return 0;
}

#endif

