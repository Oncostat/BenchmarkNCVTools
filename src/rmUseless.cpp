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

	double AF;
	string ldata;
	int stop=0;
	int count=0, chr, pos, ref, alt,inf;
	string hnames;
	while(stop==0){
		getline(data,ldata,'\n');
		if(ldata[1]!='#'){ 
			istringstream header(ldata);
			stop=1;
			while(getline(header,hnames,'\t')){
				if(hnames=="CHROM") chr=count;
				if(hnames=="POS") pos=count;
				if(hnames=="REF") ref=count;
				if(hnames=="ALT") alt=count;
				if(hnames=="INFO") inf=count;
				count++;
				}
		}
	}
	string elts,tmp;
	int test;
	vector<string> vdata;
	vector<string> valt;
	vector<string> vinf;
	while(getline(data,ldata,'\n')){ 
		istringstream line(ldata);
		int count=0;
		while(getline(line,elts,'\t')){
			if(count==chr || count==pos || count==ref || count==alt || count==inf){
					vdata.push_back(elts);
				}
			count++;
		}
		istringstream linalt(vdata.at(3));
		while(getline(linalt,elts,',')){
			valt.push_back(elts);
		}
		istringstream lininf(vdata.at(4));
		test=0;
		while(getline(lininf,elts,';') && test==0){
			if(elts.length()>3){
				if(elts.substr(0,3)=="AF="){
					elts = elts.substr(3,elts.length()); 
					istringstream lineinf(elts);
					while(getline(lineinf,elts,',')){
						vinf.push_back(elts); 
					}
					test++;
				}
			}
		}
		for(int i=0;i<(valt.size());i++){
			if(vdata.at(2).length()==1 && valt.at(i).length()==1){
				cout << vdata.at(0) << '\t' << vdata.at(1) << '\t' << vdata.at(2) << '\t' << valt.at(i) << '\t' << vinf.at(i) << endl;
			}
		}
		vdata.clear();
		valt.clear();
		vinf.clear();
	}
	return 0;
}
