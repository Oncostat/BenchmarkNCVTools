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

int main(int argc, char* argv[]){
	ifstream data(argv[1], std::ios::in);
	int nbs, seed;
	if(!data){cerr << "The input file does not exist." << std::endl;return 1;}
	if(argc < 3){
		cerr << "Please indicate the number of variants to sample" << std::endl;
	}else{
		istringstream ss1(argv[2]);
		if (!(ss1 >> nbs)) cerr << "Invalid number " << argv[2] << '\n';
	}
	if(argc < 4){
		time_t tim=time(NULL);
		srand (tim);
		cout << "The seed was initialized at " << tim << endl;
	}else{
		istringstream ss2(argv[3]);
		if (!(ss2 >> seed)){
			cerr << "Invalid number " << argv[3] << '\n';
		}else{
			srand(seed);
		}
	}
	string ldata;
	int count=0;
	while(getline(data,ldata,'\n')){
		count++;
	}
	data.clear();
	data.seekg(0);
	std::vector <int> vec(count);
	for(int i=1; i<=count; i++){
		vec[i]=i;
	}
	random_shuffle(vec.begin(),vec.end());
	std::vector <int> rdm(nbs);
	for(int i=0; i<nbs; i++){
		rdm[i]=vec[i];
	}
	sort(rdm.begin(),rdm.end());

	getline(data,ldata,'\n');
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


