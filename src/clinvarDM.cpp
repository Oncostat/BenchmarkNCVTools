/*
 * clinvar_rm_useless.cpp
 *
 *  Created on: 25 janv. 2016
 *      Author: drubay
 */



#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <stdio.h>
#include <cstring>
#include <vector>


using namespace std;



int main(int argc, char* argv[]){
	ifstream data(argv[1], std::ios::in);
	ofstream output(argv[2], ios::out | ios::trunc);

	if (!data) {cerr << "The input does not exist." << std::endl;return 1;}
	if (!output) {cerr << "The output cannot be created." << std::endl;return 1;}

	string ldata;
	int stop=0,stopint;
	int count=0, chr, pos, ref, alt,inf,countLabel=0;
	string hnames;
	string elts,tmp,tmp2,tmp3;
	int test;
	vector<string> vdata;
	vector<string> valt;
	vector<string> vinf;
	while(getline(data,ldata,'\n')){
		istringstream line(ldata);
		int count=0;
		while(getline(line,elts,'\t')){
			vdata.push_back(elts);
		}
		istringstream vec(vdata.at(4));
		stopint=0;
		while(getline(vec,tmp,';') && stopint==0){
			if(tmp.length()>10){
				int found = tmp.find("=");
				if(tmp.substr(0,found)=="CLNDSDBID"){
					tmp2 = tmp.substr(found+1,tmp.size());
					int init =tmp2.find("C");
					if(init>=0){
						int size=tmp2.size()-init;
						vdata.at(4)=tmp2.substr(init,min(size,8));
					}else{
						vdata.at(4)="NF"+to_string(countLabel);
						countLabel++;
					}
					stopint++;
				}
			}
		}
		if(stopint==0){vdata.at(4)="";}

		istringstream linalt(vdata.at(3)); 
		while(getline(linalt,elts,',')){
			valt.push_back(elts);
		}
		istringstream lininf(vdata.at(4)); 
		test=0;
		for(int i=0;i<(valt.size());i++){
			if(vdata.at(2).length()==1 && valt.at(i).length()==1){
				output << vdata.at(0) << '\t' << vdata.at(1) << '\t' << vdata.at(2) << '\t' << valt.at(i) << '\t' << vdata.at(4) << endl;
			}
		}
		vdata.clear();
		valt.clear();
		vinf.clear();
	}
	return 0;
}


