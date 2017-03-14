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

	if (!data){cerr << "The file " << data << " does not exist." << std::endl;}

	string ldata, elt, str_nbline;
	int nbline=0;
	string tmp, checkelt;

	while(getline(data,ldata,'\n')){
		nbline++;
		ostringstream convertnbl;
		convertnbl << nbline;
		str_nbline = convertnbl.str();
		istringstream line(ldata);
		int count=0, testloop=0;
		tmp="";
		while(getline(line,elt,'\t') && count<2){
			if(count==0){
				if(elt=="23"){elt="X";}
				if(elt=="24"){elt="Y";}
				if(elt=="25"){elt="M";}
				tmp = elt + '\t';
				}else{
				if(elt!=checkelt){
					checkelt = elt;
					istringstream tmpelt(elt);
					int int_elt;
					tmpelt >> int_elt;
					int_elt--;
					string str_elt;
					ostringstream convert;
					convert << int_elt;
					str_elt = convert.str();
					tmp = tmp + str_elt  + '\t' + elt;
					tmp = "chr" + tmp + '\t' + str_nbline + '\n'; 		
					cout << tmp;
				}
			}
			count++;
		}
	}
	data.close();
	return 0;
}


