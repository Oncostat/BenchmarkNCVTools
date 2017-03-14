

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <boost/algorithm/string.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_io.hpp>
#include <boost/unordered_map.hpp>
#include <boost/timer.hpp>



using namespace std;

int main(int argc, char* argv[]){
	ifstream data(argv[1], std::ios::in);
	ifstream gwava(argv[2], std::ios::in);

	if(!data){cerr << "The input file does not exist." << std::endl;return 1;}
	if(!gwava){cerr << "The input GWAVA scored file does not exist." << std::endl;return 1;}


        string ldata,lscore,tmp,elts,tmpsave;

        int count;
	streampos pos=0;
	vector<string> strsc(7);
        vector<string> strsave(7);

        getline(gwava,lscore);
        boost::split(strsave, lscore, boost::is_any_of("\t"));
        tmpsave=strsave.at(0).substr(3,strsave.at(0).length()-3);

        while(getline(data,ldata)){
            vector<string> strsd(3,"");
            istringstream w(ldata);
            count=0;
            while(getline(w,elts,'\t') && count<3){strsd.at(count)=elts;count++;}
            if(strsd[0]==tmpsave && strsd[1]==strsave[2]){
                    cout << ldata << '\t' << strsave[4] << '\t' << strsave[5] << '\t' << strsave[6] << endl;
            }else{
                pos = gwava.tellg();
                getline(gwava,lscore);
                boost::split(strsc, lscore, boost::is_any_of("\t"));
                tmp=strsc.at(0).substr(3,strsc.at(0).length()-3);
                if(strsd[0]==tmp && strsd[1]==strsc[2]){
                    cout << ldata << '\t' << strsc[4] << '\t' << strsc[5] << '\t' << strsc[6] << endl;
                    strsave=strsc;
                }else{
                    cout << ldata << "\tNA\tNA\tNA" << endl;
                    gwava.clear();
                    gwava.seekg(pos);
                }
            }
        }

	return 0;
}
