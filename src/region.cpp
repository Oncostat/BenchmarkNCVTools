#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <algorithm>
#include <vector>



using namespace std;
int main(int argc, char* argv[]){
	ifstream data(argv[1], std::ios::in);
	ifstream patho(argv[2], std::ios::in);

	if(!data){cerr << "The 1000 genomes input file does not exist." << endl;return 1;}
	if(!patho){cerr << "The pathologic variant file does not exist." << endl;return 1;}

	string ldata,lpatho,elts;
	int count,stop,tmp;
	vector<int> vidata(2);
	vector<int> viscore(2);

	getline(data,ldata,'\n');
	streampos posref = data.tellg();
	streampos pos = posref;

	while(getline(patho,lpatho,'\n')){
		istringstream w(lpatho);
		count=0;
		while(getline(w,elts,'\t') && count<2){
			istringstream(elts) >> tmp;
			viscore.at(count)=tmp;
			count++;
		}
		stop=0;
		while(stop==0 && !data.eof()){
			getline(data,ldata,'\n');
			istringstream line(ldata);
			count=0;
			while(getline(line,elts,'\t') && count <2){
				if(count==0 || count==1){
					istringstream(elts) >> tmp;
					vidata.at(count)=tmp;
				}
				count++;
			}
			if(vidata.at(0)>viscore.at(0) || (vidata.at(0)==viscore.at(0) && vidata.at(1)>(viscore.at(1)+1000))){
				stop++;
				data.clear();data.seekg(pos);
			}else{
				if(vidata.at(0)==viscore.at(0) && vidata.at(1)>=(viscore.at(1)-1000) && vidata.at(1)<=(viscore.at(1)+1000)){
					cout << ldata << endl;
					stop++;
				}else{
					pos = data.tellg();
				}
			}
		}
	}
	return 0;
}




