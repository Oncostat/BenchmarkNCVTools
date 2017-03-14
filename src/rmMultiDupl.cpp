#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h> 
#include <time.h>  
#include <vector>


using namespace std;




int main(int argc, char* argv[]){
	ifstream data(argv[1], std::ios::in);
	if (!data) {cerr << "The input file does not exist." << std::endl;return 1;}

	string ldata,elts;
	vector<string> chr,pos,ref,alt,AF,tmpchr(1),tmppos(1),tmpref(1),tmpalt(1),tmpAF(1);
	int test=0;

	getline(data,ldata,'\n');
	istringstream line(ldata);
	int count=0;
	while(getline(line,elts,'\t')){
		if(count==0){chr.push_back(elts);}
		if(count==1){pos.push_back(elts);}
		if(count==2){ref.push_back(elts);}
		if(count==3){alt.push_back(elts);}
		if(count==4){AF.push_back(elts);}
		count++;
	}


	while(test==0){
		getline(data,ldata,'\n');
		istringstream line(ldata);
		count=0;
		while(getline(line,elts,'\t')){
			if(count==0){tmpchr.at(0)=elts;}
			if(count==1){tmppos.at(0)=elts;}
			if(count==2){tmpref.at(0)=elts;}
			if(count==3){tmpalt.at(0)=elts;}
			if(count==4){tmpAF.at(0)=elts;}
			count++;
		}
		if(tmpchr.at(0)!=chr.at(0) || (tmpchr.at(0)==chr.at(0) && tmppos.at(0)!=pos.at(0))){
			for(int i=0; i<chr.size(); i++){
				cout << chr.at(i) << '\t' << pos.at(i) << '\t' << ref.at(i) << '\t' << alt.at(i) << '\t' << AF.at(i) << endl;
			}
			chr.clear();
			pos.clear();
			ref.clear();
			alt.clear();
			AF.clear();
			chr = tmpchr;
			pos = tmppos;
			ref = tmpref;
			alt = tmpalt;
			AF = tmpAF;
		}else{
			chr.push_back(tmpchr.at(0));
			pos.push_back(tmppos.at(0));
			ref.push_back(tmpref.at(0));
			alt.push_back(tmpalt.at(0));
			AF.push_back(tmpAF.at(0));
			if(tmpalt.at(0)==alt.at(0)){
				chr.erase(chr.begin());
				pos.erase(pos.begin());
				ref.erase(ref.begin());
				alt.erase(alt.begin());
				AF.erase(AF.begin());
			}
			if(data.eof()){
				for(int i=0; i<chr.size(); i++){
					cout << chr.at(i) << '\t' << pos.at(i) << '\t' << ref.at(i) << '\t' << alt.at(i) << '\t' << AF.at(i) << endl;
				}
				test++;
			}
		}
	}
	return 0;
}



