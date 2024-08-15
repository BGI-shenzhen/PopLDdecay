#ifndef processV1_H_
#define processV1_H_

using namespace std;


inline int PairWiseComV3 ( In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPList , StarRsult *All_Stat  , int & Flag_for_pro )
{
	map < string,map <llong, int > >   NewSNPList ;
	vector < vector <BaseType> >  NewBase ;
	changSNPList (SNPList, NewSNPList ,  NewBase ) ;
	SNPList.clear();

	map <string,map<llong, int  > > :: iterator key1_it ;
	int bin_count=int(Flag_for_pro/100)+1;
	int Flag=1;
	int Flag_Now=0;

	map<llong, int >  :: iterator key2 ;
	map<llong, int >  :: iterator key2_se  ;
	llong Dis=0 ;
	statementVar  Var ;	

	Var.Asize=NewBase[1].size();
	int FirstVeSite=NewBase.size();


	cerr<<"##begin pair-wise R^2 cal. after filter Remain SNP Number :\t"<<Flag_for_pro<<"\tSites In Memory Numbers :"<<FirstVeSite<<"\n#\% number bin is\t"<<bin_count<<endl;




	if ((paraFA04->TF) ==1)
	{





		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tD'\tr^2\tDist\n";

		PairInfoV1 ** ArryCalResult = new PairInfoV1 *[FirstVeSite];
		ArryCalResult[0]= new PairInfoV1 [FirstVeSite*FirstVeSite] ;
		for(int i = 1; i < FirstVeSite; i++)
		{
			ArryCalResult[i] = ArryCalResult[i-1]+FirstVeSite;
		}


		for (key1_it=NewSNPList.begin(); key1_it!=NewSNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();
			for(  ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ;  key2_se++ ; 
				Flag_Now++;
				if ((Flag_Now/bin_count)>Flag)
				{
					Flag++;
					cerr<<Flag<<"% done...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";
				}
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}

					if(	(ArryCalResult[key2->second][key2_se->second]).D > -1)
					{
						LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).D<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).RR<<"\t"<<Dis<<"\n";
						All_Stat[Dis].Count++;
						All_Stat[Dis].sumD+=(ArryCalResult[key2->second][key2_se->second]).D;
						All_Stat[Dis].sumRR+=(ArryCalResult[key2->second][key2_se->second]).RR;

					}
					else if  ((ArryCalResult[key2->second][key2_se->second]).D==-2)
					{

						cal_RR_D( NewBase[key2->second], NewBase[key2_se->second],ArryCalResult[key2->second][key2_se->second], Var) ;
						if(	(ArryCalResult[key2->second][key2_se->second]).D > -1)
						{
							LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).D<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).RR<<"\t"<<Dis<<"\n";
							All_Stat[Dis].Count++;
							All_Stat[Dis].sumD+=(ArryCalResult[key2->second][key2_se->second]).D;
							All_Stat[Dis].sumRR+=(ArryCalResult[key2->second][key2_se->second]).RR;

						}
					}
				}
			}
		}



		delete [] ArryCalResult ;
		LDOUT.close();


	}



	else if ((paraFA04->TF) ==2)
	{





		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tD'\tLOD\tr^2\tDist\n";

		PairInfoV2 ** ArryCalResult = new PairInfoV2 *[FirstVeSite];
		ArryCalResult[0]= new PairInfoV2 [FirstVeSite*FirstVeSite] ;
		for(int i = 1; i < FirstVeSite; i++)
		{
			ArryCalResult[i] = ArryCalResult[i-1]+FirstVeSite;
		}


		for (key1_it=NewSNPList.begin(); key1_it!=NewSNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();
			for(  ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ;  key2_se++ ; 
				Flag_Now++;
				if ((Flag_Now/bin_count)>Flag)
				{
					Flag++;
					cerr<<Flag<<"% done...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";
				}
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}

					if(	(ArryCalResult[key2->second][key2_se->second]).D > -1)
					{
						LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).D<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).LOD<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).RR<<"\t"<<Dis<<"\n";
						All_Stat[Dis].Count++;
						All_Stat[Dis].sumD+=(ArryCalResult[key2->second][key2_se->second]).D;
						All_Stat[Dis].sumRR+=(ArryCalResult[key2->second][key2_se->second]).RR;
					}
					else if ((ArryCalResult[key2->second][key2_se->second]).D==-2)
					{

						cal_RR_D2( NewBase[key2->second], NewBase[key2_se->second],ArryCalResult[key2->second][key2_se->second], Var) ;
						if(	(ArryCalResult[key2->second][key2_se->second]).D > -1)
						{
							LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).D<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).LOD<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).RR<<"\t"<<Dis<<"\n";

							All_Stat[Dis].Count++;
							All_Stat[Dis].sumD+=(ArryCalResult[key2->second][key2_se->second]).D;
							All_Stat[Dis].sumRR+=(ArryCalResult[key2->second][key2_se->second]).RR;

						}

					}
				}
			}
		}



		delete [] ArryCalResult ;
		LDOUT.close();


	}




	else if ((paraFA04->TF) ==3)
	{





		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tD'\tLOD\tr^2\tCIlow\tCIhi\tDist\n";

		PairInfoV3 ** ArryCalResult = new PairInfoV3 *[FirstVeSite];
		ArryCalResult[0]= new PairInfoV3 [FirstVeSite*FirstVeSite] ;
		for(int i = 1; i < FirstVeSite; i++)
		{
			ArryCalResult[i] = ArryCalResult[i-1]+FirstVeSite;
		}


		for (key1_it=NewSNPList.begin(); key1_it!=NewSNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();
			for(  ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ;  key2_se++ ; 
				Flag_Now++;
				if ((Flag_Now/bin_count)>Flag)
				{
					Flag++;
					cerr<<Flag<<"% done...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";
				}
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}

					if(	(ArryCalResult[key2->second][key2_se->second]).D > -1)
					{
						LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).D<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).LOD<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).RR<<"\t"<<setprecision(2)<<(ArryCalResult[key2->second][key2_se->second]).low_i/100.0<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).high_i/100.0<<"\t"<<Dis<<"\n";

						All_Stat[Dis].Count++;
						All_Stat[Dis].sumD+=(ArryCalResult[key2->second][key2_se->second]).D;
						All_Stat[Dis].sumRR+=(ArryCalResult[key2->second][key2_se->second]).RR;

					}
					else if ((ArryCalResult[key2->second][key2_se->second]).D==-2)
					{

						cal_RR_D3( NewBase[key2->second], NewBase[key2_se->second],ArryCalResult[key2->second][key2_se->second], Var) ;
						if(	(ArryCalResult[key2->second][key2_se->second]).D > -1)
						{
							LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).D<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).LOD<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).RR<<"\t"<<setprecision(2)<<(ArryCalResult[key2->second][key2_se->second]).low_i/100.0<<"\t"<<(ArryCalResult[key2->second][key2_se->second]).high_i/100.0<<"\t"<<Dis<<"\n";

							All_Stat[Dis].Count++;
							All_Stat[Dis].sumD+=(ArryCalResult[key2->second][key2_se->second]).D;
							All_Stat[Dis].sumRR+=(ArryCalResult[key2->second][key2_se->second]).RR;

						}

					}
				}
			}
		}


		delete [] ArryCalResult ;
		LDOUT.close();


	}






	else
	{



		PairInfoV1 ** ArryCalResult = new PairInfoV1 *[FirstVeSite];
		ArryCalResult[0]= new PairInfoV1 [FirstVeSite*FirstVeSite] ;
		for(int i = 1; i < FirstVeSite; i++)
		{
			ArryCalResult[i] = ArryCalResult[i-1]+FirstVeSite;
		}


		for (key1_it=NewSNPList.begin(); key1_it!=NewSNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();
			for(  ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ;  key2_se++ ; 
				Flag_Now++;
				if ((Flag_Now/bin_count)>Flag)
				{
					Flag++;
					cerr<<Flag<<"% done...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";
				}
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}

					if(	(ArryCalResult[key2->second][key2_se->second]).D > -1)
					{

						All_Stat[Dis].Count++;
						All_Stat[Dis].sumD+=(ArryCalResult[key2->second][key2_se->second]).D;
						All_Stat[Dis].sumRR+=(ArryCalResult[key2->second][key2_se->second]).RR;

					}
					else if  ((ArryCalResult[key2->second][key2_se->second]).D==-2)
					{
						cal_RR_D( NewBase[key2->second], NewBase[key2_se->second],ArryCalResult[key2->second][key2_se->second], Var) ;
						if(	(ArryCalResult[key2->second][key2_se->second]).D > -1)
						{

							All_Stat[Dis].Count++;
							All_Stat[Dis].sumD+=(ArryCalResult[key2->second][key2_se->second]).D;
							All_Stat[Dis].sumRR+=(ArryCalResult[key2->second][key2_se->second]).RR;

						}
					}
				}
			}
		}



		delete [] ArryCalResult ;

	}




	cerr<<Flag<<"%......-->100%.......\t\tALL done"<<endl;


	return 1;
}



#endif // processV1_H_ //
///////// swimming in the sky and flying in the sea ////////////


