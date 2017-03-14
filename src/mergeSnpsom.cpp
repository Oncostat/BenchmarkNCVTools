#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>
#include <vector>

using namespace std;

int main(int argc, char* argv[]){
	ifstream data(argv[1], std::ios::in);
        ifstream scores("Scores/SNPSOM.bed", std::ios::in);
	ofstream output(argv[2], ios::out | ios::trunc);


	if(!data){cerr << "The input file does not exist." << endl;return 1;}
	if(!scores){cerr << "The file of SNP and SOM scores does not exist." << endl;return 1;}
	if(!output){cerr << "The output file cannot be created." << std::endl;return 1;}

	streampos pos = 0,spos = 0;
	string ldata,lscores,elts,sc,na;
	int count,tmp,stop=0, test=0;
	vector<int> vscores(3);
	vector<int> savcoord(3);
	vector<int> vdata(2);
	na= "\tNA\tNA\tNA\tNA\tNA";

	while(test==0){
		getline(scores,lscores,'\n');
		spos = scores.tellg();
		istringstream linescores(lscores);
		count=0;
		sc="";
		while(getline(linescores,elts,'\t')){
			if(count<3){
				istringstream(elts) >> tmp;
				vscores.at(count)=tmp;
			}else{
				sc += '\t' + elts;
			}
			count++;
		}

		stop=0;
		while(stop==0){
			pos = data.tellg();
			getline(data,ldata,'\n');
			istringstream line(ldata);
			count=0;
			while(getline(line,elts,'\t') && count<2){
				istringstream(elts) >> tmp;
				vdata.at(count)=tmp;
				count++;
			}
			if((vdata.at(0)==vscores.at(0) && vdata.at(1)<vscores.at(1) && vdata.at(1)>savcoord.at(2)) ||
					(vdata.at(0)<vscores.at(0) && vdata.at(0)==savcoord.at(0) && vdata.at(1)>savcoord.at(1))){
				if(ldata.length()>0){output << ldata << na;}
				if(!data.eof()){output << endl;}else{stop++;}
			}else{
				if(vdata.at(0)==vscores.at(0) && vdata.at(1)>vscores.at(1) && vdata.at(1)<=vscores.at(2)){
					if(ldata.length()>0){output << ldata << sc;}
					if(!data.eof()){output << endl;}else{stop++;}
				}else{
					if(vdata.at(0)>vscores.at(0) || (vdata.at(0)==vscores.at(0) && vdata.at(1)>vscores.at(2))){
						data.clear();data.seekg(pos);
						stop++;
					}
				}
			}
		}
		for(int i=0;i<3;i++){
			savcoord.at(i)=vscores.at(i);
		}
		if(data.eof() || scores.eof()){test++;}
	}

	data.close();
	scores.close();

	return 0;
}

