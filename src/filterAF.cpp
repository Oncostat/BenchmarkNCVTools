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
	double minAF,AFcheck;
	double nbs;
	if(argc < 3){
		cerr << "Please indicate the minimal expected AF or go to next step if you don't want to filter according the AF" << endl;
		return 1;
	}
	stringstream(argv[2]) >> minAF; 

	int count,stop;
	string ldata,elts,AF;
	while(getline(data,ldata,'\n')){
			istringstream line(ldata);
			count=0,stop=0;
			while(getline(line,elts,'\t') && stop==0){
				if(count==4){
					stringstream(elts) >> AFcheck;
					stop++;
				}
				count++;
			}
			if(AFcheck>=minAF){
				cout << ldata << endl;
			}
	}
	data.close();
}



