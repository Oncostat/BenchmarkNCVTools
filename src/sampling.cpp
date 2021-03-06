/*
 * variant_sampling.cpp
 *
 *  Created on: 13 nov. 2015
 *      Author: drubay
 */


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <stdio.h>
#include <cstring>
#include <stdlib.h>     
#include <time.h>       
#include <algorithm>	
#include <vector>



using namespace std;

int main(int argc, char* argv[])
{
	/* import */
	ifstream data(argv[1], std::ios::in); //allows to change the data (for example, if several according the AF)

	int nbs, seed; //number of SNPs to sample and seed
	/*checks*/
	if(!data){cerr << "The input file does not exist." << std::endl;return 1;}
	if(argc < 3){
		cerr << "Please indicate the number of variants to sample" << std::endl;
	}else{
		istringstream ss1(argv[2]); //nb samples required
		if (!(ss1 >> nbs)) cerr << "Invalid number " << argv[2] << '\n';
	}
	if(argc < 4){ //Notice that if there is only 2 argument, seed = default. In the other cases, the provided number is check and used as seed
		time_t tim=time(NULL);
		srand (tim); //default seed
		cout << "The seed was initialized at " << tim << endl;
	}else{
		istringstream ss2(argv[3]); //seed
		if (!(ss2 >> seed)){
			cerr << "Invalid number " << argv[3] << '\n';
		}else{
			srand(seed);
		}
	}

	/*Count the number of lines to attribute rdm number to lines
	 * add argv for number of lines to skip this step (possibly time consuming)*/
	string ldata;
	int count=0;
	while(getline(data,ldata,'\n')){
		count++;
	}
	data.clear();  // Clears all error flags.
	data.seekg(0);
	/*Generate random number without replacement*/
	std::vector <int> vec(count); //created a vector from 1 to number of line
	for(int i=1; i<=count; i++){
		vec[i]=i;
	}
	random_shuffle(vec.begin(),vec.end()); //shuffle line numbers
	std::vector <int> rdm(nbs);
	for(int i=0; i<nbs; i++){    //select the first nbs shuffled lines
		rdm[i]=vec[i];
	}
	sort(rdm.begin(),rdm.end());

	/*Sampling of corresponding rows*/
	getline(data,ldata,'\n'); //skip header
	count=0;
	int pos=0;
	while(pos<nbs && !data.eof()){
		count++;
		getline(data,ldata,'\n');
		if(count==rdm[pos]){
			cout << ldata << endl;
		pos++;}
	}
	data.close();
	return 0;
}



