#ifndef processV3_H_
#define processV3_H_

using namespace std;


int PairWiseComV1( In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPList ,  StarRsult *  All_Stat , int & Flag_for_pro )
{

	map <string,map<llong, vector <BaseType>  > > :: iterator key1_it ;
	int bin_count=int(Flag_for_pro/100)+1;
	int Flag=1;
	int Flag_Now=0;

	map<llong, vector <BaseType> >  :: iterator key2 ;
	map<llong, vector <BaseType> >  :: iterator key2_se  ;
	llong Dis=0;
	statementVar  Var ;	

	key1_it=SNPList.begin(); 
	key2=(key1_it->second).begin();
	Var.Asize= (key2->second).size();

	cerr<<"##begin pair-wise R^2 cal. after filter Remain SNP Number :\t"<<Flag_for_pro<<"\n#\% number bin is\t"<<bin_count<<endl;



//#####################


	if ((paraFA04->TF) ==1)
	{
		for (key1_it=SNPList.begin(); key1_it!=SNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();

			for(   ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ;  key2_se++ ; 
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					{
						vector<BaseType> & g1 = key2->second;
						vector<BaseType> & g2 = key2_se->second;
						unsigned short int dde[3][3] = {{0}};
						int asize = Var.Asize;
						for (int _i = 0; _i < asize; _i++) {
							dde[g1[_i].Value][g2[_i].Value]++;
						}
						int tmpAA = dde[1][1] + dde[1][0];
						if (tmpAA == 0) continue;
						if ((dde[1][1] + dde[0][1]) == 0) continue;
						int ALL_count = dde[0][0] + dde[0][1] + tmpAA;
						double prob0 = (double)dde[0][0] / ALL_count;
						double prob1 = (double)dde[0][1] / ALL_count;
						double prob2 = (double)dde[1][0] / ALL_count;
						double pA1 = prob0 + prob1;
						double pA2 = prob0 + prob2;
						double Cal_B = pA1 * pA2;
						double Cal_A = 1.0 - pA1 - pA2 + Cal_B;
						if (Cal_A == 0 || Cal_B == 0) {
							if (prob0 < 1e-10) prob0 = 1e-10;
							if (prob1 < 1e-10) prob1 = 1e-10;
							if (prob2 < 1e-10) prob2 = 1e-10;
							pA1 = prob0 + prob1;
							pA2 = prob0 + prob2;
							Cal_B = pA1 * pA2;
							Cal_A = 1.0 - pA1 - pA2 + Cal_B;
						}
						double D_A = prob0 - Cal_B;
						double rr = (D_A * D_A) / (Cal_A * Cal_B);
						All_Stat[Dis].Count++;
						All_Stat[Dis].sumRR += rr;
					}
				}
			}
		}
	}


	else if  ((paraFA04->TF) ==2  )
	{

		PairInfoV1  CalResult ;

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
							All_Stat[Dis].Count++;
							All_Stat[Dis].sumD+=CalResult.D;
							All_Stat[Dis].sumRR+=CalResult.RR;

						}
					}
				}
			}
		}
	}


	else if  ((paraFA04->TF) ==3 )
	{
		double  CalResult ;
		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tr^2\tDist\n";
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
					if (cal_RR_MA( key2->second , key2_se->second  ,CalResult, Var) ==1)
					{
						All_Stat[Dis].Count++;
						All_Stat[Dis].sumRR+=CalResult;
						LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<CalResult<<"\t"<<Dis<<"\n";
					}
				}
			}
		}
		LDOUT.close();
	}




	if ((paraFA04->TF) ==6 )
	{
		PairInfoV1  CalResult ;
		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tD'\tr^2\tDist\n";
		for (key1_it=SNPList.begin(); key1_it!=SNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();
			for( 	; key2!=(key1_it->second).end(); key2++ )
			{
				key2_se=key2 ;  key2_se++ ; 
				Flag_Now++;
				if ((Flag_Now/bin_count)>Flag)
				{
					cerr<<Flag<<"%...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";
					Flag++;
				}
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					else
					{
						if(	cal_RR_D_MA( key2->second , key2_se->second  , CalResult , Var) ==1)
						{
							LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<CalResult.D<<"\t"<<CalResult.RR<<"\t"<<Dis<<"\n";
							All_Stat[Dis].Count++;
							All_Stat[Dis].sumD+=CalResult.D;
							All_Stat[Dis].sumRR+=CalResult.RR;
						}
					}
				}
			}	
		}
		LDOUT.close();
	}


	else if  ((paraFA04->TF) ==7  )
	{
		PairInfoV2  CalResult ;
		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tD'\tLOD\tr^2\tDist\n";
		for (key1_it=SNPList.begin(); key1_it!=SNPList.end() ;key1_it++ )
		{
			key2=(key1_it->second).begin();

			for( ; key2!=(key1_it->second).end(); key2++ )
			{
				key2_se=key2 ;  key2_se++ ; 
				Flag_Now++;
				if ((Flag_Now/bin_count)>Flag)
				{
					cerr<<Flag<<"%...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";
					Flag++;
				}
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					else
					{
						if(	cal_RR_D2_MA( key2->second , key2_se->second  , CalResult, Var ) ==1)
						{
							LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<CalResult.D<<"\t"<<CalResult.LOD<<"\t"<<CalResult.RR<<"\t"<<Dis<<"\n";
							All_Stat[Dis].Count++;
							All_Stat[Dis].sumD+=CalResult.D;
							All_Stat[Dis].sumRR+=CalResult.RR;
						}
					}
				}
			}
		}
		LDOUT.close();


	}


	else if  ((paraFA04->TF) ==8   )
	{

		PairInfoV3  CalResult ;
		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tD'\tLOD\tr^2\tCIlow\tCIhi\tDist\n";
		for (key1_it=SNPList.begin(); key1_it!=SNPList.end() ;key1_it++)		
		{
			key2=(key1_it->second).begin();

			for( ; key2!=(key1_it->second).end(); key2++ )
			{
				key2_se=key2 ;  key2_se++ ; 
				Flag_Now++;
				if ((Flag_Now/bin_count)>Flag)
				{
					cerr<<Flag<<"%...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";
					Flag++;
				}
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					else
					{
						if(	cal_RR_D3_MA( key2->second , key2_se->second  , CalResult, Var ) ==1)
						{
							LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<CalResult.D<<"\t"<<CalResult.LOD<<"\t"<<CalResult.RR<<"\t"<<setprecision(2)<<CalResult.low_i/100.0<<"\t"<<CalResult.high_i/100.0<<"\t"<<Dis<<"\n";
							All_Stat[Dis].Count++;
							All_Stat[Dis].sumD+=CalResult.D;
							All_Stat[Dis].sumRR+=CalResult.RR;

						}
					}
				}
			}
		}
		LDOUT.close();


	}




	








	cerr<<Flag<<"%......-->100%.......\t\tALL done"<<endl;

	return 1;
}



#endif // processV3_H_ //
///////// swimming in the sky and flying in the sea ////////////


