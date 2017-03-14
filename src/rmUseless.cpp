/*
 * first_step.cpp
 *
 *  Created on: 8 déc. 2015
 *      Author: drubay
 */


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
#include <vector>


using namespace std;



int main(int argc, char* argv[])
//int main() //add argc argv for the choice of the frequency with a default value >1 (>100% to be sure to include all variants)
{
	ifstream data(argv[1], std::ios::in);

//	ifstream data("../Databases/1KGp1.vcf", std::ios::in);
//	ofstream output("../Databases/data1KGp1.txt", ios::out | ios::trunc);

	/*checks*/
	if (!data) {cerr << "The input file does not exist." << std::endl;return 1;}

	/*Data management*/
	double AF;
	string ldata;
	int stop=0;
	int count=0, chr, pos, ref, alt,inf;
	string hnames;
	while(stop==0){
		getline(data,ldata,'\n');
		if(ldata[1]!='#'){ //skip comments
			istringstream header(ldata);
			stop=1;
			while(getline(header,hnames,'\t')){ //find column pos, ref & alt
				if(hnames=="CHROM") chr=count;
				if(hnames=="POS") pos=count;
				if(hnames=="REF") ref=count;
				if(hnames=="ALT") alt=count;
				if(hnames=="INFO") inf=count;
				count++;
				}
		}
	}// Include only the colnames, thus, next getline begin at the first informative line
	string elts,tmp;
	int test;
	vector<string> vdata;
	vector<string> valt;
	vector<string> vinf;
	while(getline(data,ldata,'\n')){ // && testext<10
		istringstream line(ldata);
		int count=0;
		while(getline(line,elts,'\t')){//elements of chain
			if(count==chr || count==pos || count==ref || count==alt || count==inf){
					vdata.push_back(elts);
				}
			count++;
		}
		istringstream linalt(vdata.at(3));   //extract the different ALTs
		while(getline(linalt,elts,',')){
			valt.push_back(elts);
		}
		istringstream lininf(vdata.at(4));   //extract the different AFs
		test=0;
		while(getline(lininf,elts,';') && test==0){
			if(elts.length()>3){
				if(elts.substr(0,3)=="AF="){
					elts = elts.substr(3,elts.length()); //remove the first 3 characters ("AF=")
					istringstream lineinf(elts);
					while(getline(lineinf,elts,',')){
						vinf.push_back(elts); 		//stock the different values of AF separated by ','
					}
					test++;
				}
			}
		}
		//for each values of the vector, re-write the line if it is not an indel
		for(int i=0;i<(valt.size());i++){
			if(vdata.at(2).length()==1 && valt.at(i).length()==1){ //if not an indel (1 character corresponding to one nucleotide)
				cout << vdata.at(0) << '\t' << vdata.at(1) << '\t' << vdata.at(2) << '\t' << valt.at(i) << '\t' << vinf.at(i) << endl;
			}
		}
		vdata.clear();
		valt.clear();
		vinf.clear();
	}
	return 0;
}
