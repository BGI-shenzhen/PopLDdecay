#ifndef comm_H_
#define comm_H_


using namespace std;
typedef int  llong ;


//////////////////////////////// some small funtion /////////////////////////////////////////////


inline void  LogLackArg( string  flag )
{
	cerr << "\t\tLack Argument for [ -"<<flag<<" ]"<<endl;
}

/*
inline string add_Asuffix ( string path )
{
	string ext =path.substr(path.rfind('.') ==string::npos ? path.length() : path.rfind('.') + 1);
	if (ext != "gz")
	{
		path=path+".gz" ;
	}
	return path ;
}
*/

inline string &  replace_all(string &  str,const  string &  old_Avalue,const string &  new_Avalue)
{
	while(true)   {
		string::size_type  pos(0);
		if(   (pos=str.find(old_Avalue))!=string::npos   )
			str.replace(pos,old_Avalue.length(),new_Avalue);
		else   break;
	}
	return   str;
}


inline void split(const string& str,vector<string>& tokens,  const string& delimiters = " ")
{
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	while (string::npos != pos || string::npos != lastPos)
	{
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		lastPos = str.find_first_not_of(delimiters, pos);
		pos = str.find_first_of(delimiters, lastPos);
	}
}





///////////////////swimming in the sky & flying in the sea/////////////////////////////





#endif // comm_H_  ;


