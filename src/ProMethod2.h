#ifndef processV2_H_
#define processV2_H_

using namespace std;


int ComPairDEst( vector <BaseType> &  AA,  vector <BaseType> & BB ,int  & Asize )
{
	for (int RR=0; RR<Asize  ; RR++)
	{
		if  (   AA[RR].Value <   BB[RR].Value ) 
		{
			return 1 ;
		}
		else if  (   AA[RR].Value >   BB[RR].Value )
		{
			return -1;
		}
	}
	return 0;
}

int  CompairedTwoVec( vector < vector <BaseType> > & NewBase  , vector <BaseType> & NowVector , int & Asize,                         vector < vector < pair <string, llong> > >  & InfoSite, pair <string, llong> & Site   )
{
	int mid;
	int start=0;
	int end=NewBase.size()-1;

	int Key1=ComPairDEst(NewBase[0],NowVector, Asize);

	vector < vector <BaseType> > :: iterator  theIterator ;
	vector < vector < pair <string, llong> > >  :: iterator The_it_Site ;
	vector < pair <string, llong> > ::  iterator   It_second ;

	while (start <= end)
	{
		mid = start + (end - start) / 2; //直接平均可能會溢位，所以用此算法
		Key1=ComPairDEst(NewBase[mid],NowVector, Asize);
		if (Key1 < 0)
		{
			end = mid - 1;
		}
		else if (Key1 > 0 )
		{
			start = mid + 1 ;
		}
		else
		{
			InfoSite[mid].push_back(Site);
			return 1 ;			
		}
	}

	if (end<1) {end=0;}

	for (int ii=end ; ii<=start  ;ii++)
	{
		Key1=ComPairDEst(NewBase[ii],NowVector, Asize);
		if (Key1<0)
		{			
			continue;
		}
		else if  (Key1>0)
		{
			theIterator= NewBase.begin()+ii+1;
			NewBase.insert( theIterator, 1, NowVector );
			The_it_Site=InfoSite.begin()+ii+1;
			vector < pair <string, llong> >  Atemp;
			Atemp.push_back(Site);
			InfoSite.insert( The_it_Site, 1, Atemp  );
			return 1;
		}
		else 
		{
			InfoSite[ii].push_back(Site);
			return 1 ;
		}
	}

	theIterator= NewBase.begin();
	NewBase.insert( theIterator, 1, NowVector );
	The_it_Site=InfoSite.begin();
	vector < pair <string, llong> >  Atemp;
	Atemp.push_back(Site);
	InfoSite.insert( The_it_Site, 1, Atemp  );
	return 1 ;

}

int changSNPList ( map <string,map <llong, vector <BaseType> > >  & SNPList  ,   map < string,map <llong, int > >  & NewSNPList , vector < vector <BaseType> > & NewBase )
{
	map <string,map<llong, vector <BaseType>  > > :: iterator key1_it ;
	key1_it=SNPList.begin();
	map<llong, vector <BaseType> >  :: iterator key2 ;
	key2=(key1_it->second).begin();
	int Asize= (key2->second).size();

	NewBase.push_back(key2->second);

	vector < vector < pair <string, llong> > >  InfoSite ;
	pair <string, llong> Site ;


	Site=make_pair(key1_it->first,key2->first);
	vector < pair <string, llong> > tempA;
	tempA.push_back(Site);
	InfoSite.push_back(tempA);

	for (key1_it=SNPList.begin(); key1_it!=SNPList.end() ;key1_it++ )
	{
		for( key2=(key1_it->second).begin() ; key2!=(key1_it->second).end(); key2++ )
		{
			Site=make_pair(key1_it->first,key2->first);
			CompairedTwoVec( NewBase , (key2->second), Asize, InfoSite , Site );
		}
	}


	map < string,map <llong, int > >  :: iterator itSNP ;
	int FirstVeSite=InfoSite.size();
	for (int ii=0 ; ii<FirstVeSite  ; ii++)
	{
		int SeconDVe=InfoSite[ii].size();
		for (int jj=0 ; jj<SeconDVe ; jj++)
		{
			Site=InfoSite[ii][jj] ;
			itSNP=NewSNPList.find(Site.first);
			if (itSNP == NewSNPList.end())
			{
				map <llong, int > DD;
				DD.insert(map <llong,int >::value_type(Site.second,ii));
				NewSNPList.insert(map <string,map <llong,int> > ::value_type(Site.first,DD));
			}
			else
			{
				(itSNP->second).insert(map <llong, int >  :: value_type(Site.second,ii)) ;
			}
		}
	}

	return 1;
	//  double middel sort vector
}


int PairWiseComV2 ( In3str1v * paraFA04, Para_18 * para_18 , map <string,map <llong, vector <BaseType> > >  & SNPList , StarRsult * All_Stat , int & Flag_for_pro )
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
		PairInfoV1  CalResult ;
		map < int ,map <int, PairInfoV1 > >   Hash_CalResult ;
		map < int ,map <int, PairInfoV1 > >  :: iterator  OUT_it_key1 ;
		map < int, PairInfoV1 >  :: iterator  OUT_it_key2;


		for (key1_it=NewSNPList.begin(); key1_it!=NewSNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();
			for(  ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ;  key2_se++ ; 
				//Flag_Now++;
				//if ((Flag_Now/bin_count)>Flag) {	Flag++;		cerr<<Flag<<"% done...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";			}

				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					OUT_it_key1=Hash_CalResult.find(key2->second);
					if  (OUT_it_key1!=Hash_CalResult.end())
					{
						OUT_it_key2=(OUT_it_key1->second).find(key2_se->second);
						if  (OUT_it_key2!=(OUT_it_key1->second).end())
						{
							if ((OUT_it_key2->second).D>-1)
							{
								All_Stat[Dis].Count++;
								All_Stat[Dis].sumRR+=(OUT_it_key2->second).RR;
							}
						}
						else
						{
							cal_RR_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
							(OUT_it_key1->second).insert(map <int,PairInfoV1 > :: value_type(key2_se->second,CalResult));
							if(	CalResult.D> -1)
							{

								All_Stat[Dis].Count++;
								All_Stat[Dis].sumRR+=CalResult.RR;
							}
						}
					}
					else
					{
						cal_RR_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
						map <int,PairInfoV1 > DDtemp;
						DDtemp.insert(map <int,PairInfoV1 > :: value_type(key2_se->second,CalResult));
						Hash_CalResult.insert(map < int ,map <int, PairInfoV1 > > :: value_type (key2->second,DDtemp));
						if(	CalResult.D> -1)
						{
							All_Stat[Dis].Count++;
							All_Stat[Dis].sumRR+=CalResult.RR;
						}
					}
				}
			}
		}

	}



	else if ((paraFA04->TF) ==2)
	{
		PairInfoV1  CalResult ;
		map < int ,map <int, PairInfoV1 > >   Hash_CalResult ;
		map < int ,map <int, PairInfoV1 > >  :: iterator  OUT_it_key1 ;
		map < int, PairInfoV1 >  :: iterator  OUT_it_key2;

		for (key1_it=NewSNPList.begin(); key1_it!=NewSNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();
			for(  ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ;  key2_se++ ; 
//				Flag_Now++;				if ((Flag_Now/bin_count)>Flag)				{					Flag++;					cerr<<Flag<<"% done...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";				}
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					OUT_it_key1=Hash_CalResult.find(key2->second);
					if  (OUT_it_key1!=Hash_CalResult.end())
					{
						OUT_it_key2=(OUT_it_key1->second).find(key2_se->second);
						if  (OUT_it_key2!=(OUT_it_key1->second).end())
						{
							if ((OUT_it_key2->second).D>-1)
							{
								All_Stat[Dis].Count++;
								All_Stat[Dis].sumD+=(OUT_it_key2->second).D;
								All_Stat[Dis].sumRR+=(OUT_it_key2->second).RR;
							}
						}
						else
						{
							cal_RR_D_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
							(OUT_it_key1->second).insert(map <int,PairInfoV1 > :: value_type(key2_se->second,CalResult));
							if(	CalResult.D> -1)
							{

								All_Stat[Dis].Count++;
								All_Stat[Dis].sumD+=CalResult.D;
								All_Stat[Dis].sumRR+=CalResult.RR;


							}

						}
					}
					else
					{
						cal_RR_D_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
						map <int,PairInfoV1 > DDtemp;
						DDtemp.insert(map <int,PairInfoV1 > :: value_type(key2_se->second,CalResult));
						Hash_CalResult.insert(map < int ,map <int, PairInfoV1 > > :: value_type (key2->second,DDtemp));
						if(	CalResult.D> -1)
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




	else if ((paraFA04->TF) ==3)
	{


		PairInfoV1  CalResult ;
		map < int ,map <int, PairInfoV1 > >   Hash_CalResult ;
		map < int ,map <int, PairInfoV1 > >  :: iterator  OUT_it_key1 ;
		map < int, PairInfoV1 >  :: iterator  OUT_it_key2;


		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tr^2\tDist\n";




		for (key1_it=NewSNPList.begin(); key1_it!=NewSNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();
			for(  ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ;  key2_se++ ; 
//				Flag_Now++;				if ((Flag_Now/bin_count)>Flag)				{					Flag++;					cerr<<Flag<<"% done...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";				}
				for(   ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					OUT_it_key1=Hash_CalResult.find(key2->second);
					if  (OUT_it_key1!=Hash_CalResult.end())
					{
						OUT_it_key2=(OUT_it_key1->second).find(key2_se->second);
						if  (OUT_it_key2!=(OUT_it_key1->second).end())
						{
							if ((OUT_it_key2->second).D>-1)
							{
								All_Stat[Dis].Count++;
								All_Stat[Dis].sumRR+=(OUT_it_key2->second).RR;
								LDOUT<<(key1_it->first)<<"\t"<<(key2->first)<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<(OUT_it_key2->second).RR<<"\t"<<Dis<<"\n";
							}
						}
						else
						{
							cal_RR_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
							(OUT_it_key1->second).insert(map <int,PairInfoV1 > :: value_type(key2_se->second,CalResult));
							if(	CalResult.D> -1)
							{
								LDOUT<<(key1_it->first)<<"\t"<<(key2->first)<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<CalResult.RR<<"\t"<<Dis<<"\n";
								All_Stat[Dis].Count++;
								All_Stat[Dis].sumRR+=CalResult.RR;
							}
						}
					}
					else
					{
						cal_RR_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
						map <int,PairInfoV1 > DDtemp;
						DDtemp.insert(map <int,PairInfoV1 > :: value_type(key2_se->second,CalResult));
						Hash_CalResult.insert(map < int ,map <int, PairInfoV1 > > :: value_type (key2->second,DDtemp));
						if(	CalResult.D> -1)
						{
							LDOUT<<(key1_it->first)<<"\t"<<(key2->first)<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<CalResult.RR<<"\t"<<Dis<<"\n";
							All_Stat[Dis].Count++;
							All_Stat[Dis].sumRR+=CalResult.RR;
						}
					}
				}
			}
		}





		LDOUT.close();

	}

	else if ((paraFA04->TF) ==6 )
	{
		PairInfoV1  CalResult ;
		map < int ,map <int, PairInfoV1 > >   Hash_CalResult ;
		map < int ,map <int, PairInfoV1 > >  :: iterator  OUT_it_key1 ;
		map < int, PairInfoV1 >  :: iterator  OUT_it_key2;

		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tD'\tr^2\tDist\n";

		for (key1_it=NewSNPList.begin(); key1_it!=NewSNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();
			for(  ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ; 
				Flag_Now++;
				if ((Flag_Now/bin_count)>Flag)
				{
					Flag++;
					cerr<<Flag<<"% done...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";
				}
				for(  key2_se++  ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					OUT_it_key1=Hash_CalResult.find(key2->second);
					if  (OUT_it_key1!=Hash_CalResult.end())
					{
						OUT_it_key2=(OUT_it_key1->second).find(key2_se->second);
						if  (OUT_it_key2!=(OUT_it_key1->second).end())
						{
							if ((OUT_it_key2->second).D>-1)
							{
								LDOUT<<(key1_it->first)<<"\t"<<(key2->first)<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<(OUT_it_key2->second).D<<"\t"<<(OUT_it_key2->second).RR<<"\t"<<Dis<<"\n";
								All_Stat[Dis].Count++;
								All_Stat[Dis].sumD+=(OUT_it_key2->second).D;
								All_Stat[Dis].sumRR+=(OUT_it_key2->second).RR;
							}
						}
						else
						{
							cal_RR_D_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
							(OUT_it_key1->second).insert(map <int,PairInfoV1 > :: value_type(key2_se->second,CalResult));
							if(	CalResult.D> -1)
							{
								LDOUT<<(key1_it->first)<<"\t"<<(key2->first)<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<CalResult.D<<"\t"<<CalResult.RR<<"\t"<<Dis<<"\n";
								All_Stat[Dis].Count++;
								All_Stat[Dis].sumD+=CalResult.D;
								All_Stat[Dis].sumRR+=CalResult.RR;
							}
						}
					}
					else
					{
						cal_RR_D_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
						map <int,PairInfoV1 > DDtemp;
						DDtemp.insert(map <int,PairInfoV1 > :: value_type(key2_se->second,CalResult));
						Hash_CalResult.insert(map < int ,map <int, PairInfoV1 > > :: value_type (key2->second,DDtemp));
						if(	CalResult.D> -1)
						{
							LDOUT<<(key1_it->first)<<"\t"<<(key2->first)<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<CalResult.D<<"\t"<<CalResult.RR<<"\t"<<Dis<<"\n";
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
	else if ((paraFA04->TF) ==7)
	{



		PairInfoV2  CalResult ;
		map < int ,map <int, PairInfoV2 > >   Hash_CalResult ;
		map < int ,map <int, PairInfoV2 > >  :: iterator  OUT_it_key1 ;
		map < int, PairInfoV2 >  :: iterator  OUT_it_key2;




		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tD'\tLOD\tr^2\tDist\n";


		for (key1_it=NewSNPList.begin(); key1_it!=NewSNPList.end() ;key1_it++)
		{
			key2=(key1_it->second).begin();
			for(  ; key2!=(key1_it->second).end(); key2++)
			{
				key2_se=key2 ; 
				Flag_Now++;
				if ((Flag_Now/bin_count)>Flag)
				{
					Flag++;
					cerr<<Flag<<"% done...\t"<<key1_it->first<<"\t"<<(key2->first)<<"\r";
				}
				for(  key2_se++ ; key2_se!=(key1_it->second).end(); key2_se++ )
				{
					Dis=(key2_se->first)-(key2->first);
					if ( Dis> (paraFA04->InInt) )
					{
						break ;
					}
					OUT_it_key1=Hash_CalResult.find(key2->second);
					if (OUT_it_key1!=Hash_CalResult.end())
					{
						OUT_it_key2=(OUT_it_key1->second).find(key2_se->second);
						if (OUT_it_key2!=(OUT_it_key1->second).end())
						{
							if ((OUT_it_key2->second).D>-1)
							{
								LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<(OUT_it_key2->second).D<<"\t"<<(OUT_it_key2->second).LOD<<"\t"<<(OUT_it_key2->second).RR<<"\t"<<Dis<<"\n";
								All_Stat[Dis].Count++;
								All_Stat[Dis].sumD+=(OUT_it_key2->second).D;
								All_Stat[Dis].sumRR+=(OUT_it_key2->second).RR;
							}
						}
						else
						{
							cal_RR_D2_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
							(OUT_it_key1->second).insert(map <int,PairInfoV2 > :: value_type(key2_se->second,CalResult));
							if(	CalResult.D> -1)
							{
								LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<CalResult.D<<"\t"<<CalResult.LOD<<"\t"<<CalResult.RR<<"\t"<<Dis<<"\n";

								All_Stat[Dis].Count++;
								All_Stat[Dis].sumD+=CalResult.D;
								All_Stat[Dis].sumRR+=CalResult.RR;

							}
						}
					}
					else
					{
						cal_RR_D2_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
						map <int,PairInfoV2 > DDtemp;
						DDtemp.insert(map <int,PairInfoV2 > :: value_type(key2_se->second,CalResult));
						Hash_CalResult.insert(map < int ,map <int, PairInfoV2 > > :: value_type (key2->second,DDtemp));
						if(	CalResult.D> -1)
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
	else if ((paraFA04->TF) ==8 )
	{




		PairInfoV3  CalResult ;
		map < int ,map <int, PairInfoV3 > >   Hash_CalResult ;
		map < int ,map <int, PairInfoV3 > >  :: iterator  OUT_it_key1 ;
		map < int, PairInfoV3 >  :: iterator  OUT_it_key2;



		string LDOUT_file=(paraFA04->InStr2)+".LD.gz";
		ogzstream LDOUT (LDOUT_file.c_str());
		LDOUT<<"#chr\tSite1\tSite2\tD'\tLOD\tr^2\tCIlow\tCIhi\tDist\n";



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
					OUT_it_key1=Hash_CalResult.find(key2->second);
					if (OUT_it_key1!=Hash_CalResult.end())
					{
						OUT_it_key2=(OUT_it_key1->second).find(key2_se->second);
						if (OUT_it_key2!=(OUT_it_key1->second).end())
						{
							if ((OUT_it_key2->second).D>-1)
							{
								LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<(OUT_it_key2->second).D<<"\t"<<(OUT_it_key2->second).LOD<<"\t"<<(OUT_it_key2->second).RR<<"\t"<<setprecision(2)<<(OUT_it_key2->second).low_i/100.0<<"\t"<<(OUT_it_key2->second).high_i/100.0<<"\t"<<Dis<<"\n";
								All_Stat[Dis].Count++;
								All_Stat[Dis].sumD+=(OUT_it_key2->second).D;
								All_Stat[Dis].sumRR+=(OUT_it_key2->second).RR;

							}
						}
						else
						{
							cal_RR_D3_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
							(OUT_it_key1->second).insert(map <int,PairInfoV3 > :: value_type(key2_se->second,CalResult));
							if(	 CalResult.D> -1 )
							{
								LDOUT<<key1_it->first<<"\t"<<key2->first<<"\t"<<(key2_se->first)<<setprecision(4)<<setiosflags(ios::right)<<setiosflags(ios::fixed)<<"\t"<<CalResult.D<<"\t"<<CalResult.LOD<<"\t"<<CalResult.RR<<"\t"<<setprecision(2)<<CalResult.low_i/100.0<<"\t"<<CalResult.high_i/100.0<<"\t"<<Dis<<"\n";
								All_Stat[Dis].Count++;
								All_Stat[Dis].sumD+=CalResult.D;
								All_Stat[Dis].sumRR+=CalResult.RR;
							}
						}
					}
					else
					{
						cal_RR_D3_MB( NewBase[key2->second], NewBase[key2_se->second],CalResult, Var) ;
						map <int,PairInfoV3 > DDtemp;
						DDtemp.insert(map <int,PairInfoV3 > :: value_type(key2_se->second,CalResult));
						Hash_CalResult.insert(map < int ,map <int, PairInfoV3 > > :: value_type (key2->second,DDtemp));
						if(	CalResult.D> -1)
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



#endif // processV2_H_ //
///////// swimming in the sky and flying in the sea ////////////


