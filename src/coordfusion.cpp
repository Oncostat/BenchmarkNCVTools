#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <istream>
#include <stdio.h>
#include <cstring>
#include <vector>


using namespace std;

int main(int argc, char* argv[])
{
	ifstream data(argv[1], std::ios::in);

	if(!data){cerr << "The file does not exist." << endl;return 1;}

	string ldata,elts,chr,savechr,rec="";
	int count,next,stop,begin,end,savebegin,saveend,test;
	vector<string> score,savescore;

	getline(data,ldata,'\n'); 
	istringstream line(ldata);
	count=0;
	while(getline(line,elts,'\t')){
		if(count==0){savechr=elts;}
		if(count==1){istringstream(elts) >> savebegin;}
		if(count==2){istringstream(elts) >> saveend;}
		if(count>2){savescore.push_back(elts);}
		count++;
	}
	stop=0;
	while(stop==0){
		getline(data,ldata,'\n');
		istringstream line(ldata);
		count=0;
		next=0;
		while(getline(line,elts,'\t')){
			if(count==0){chr=elts;}
			if(count==1){istringstream(elts) >> begin;}
			if(count==2){istringstream(elts) >> end;}
			if(count>2){score.push_back(elts);}
			count++;
		}
		if(data.eof() && score.size()==0){
			stop++;
			for(int i=0;i<savescore.size();i++){rec += '\t' + savescore.at(i);}
			cout << savechr << '\t' << savebegin << '\t' << saveend << rec << endl;
		}else{
			if(chr!=savechr){
				for(int i=0;i<savescore.size();i++){rec += '\t' + savescore.at(i);}
				cout << savechr << '\t' << savebegin << '\t' << saveend << rec << endl;
				rec="";
				savechr=chr;
				savebegin=begin;
				saveend=end;
				for(int i=0;i<savescore.size();i++){savescore.at(i)=score.at(i);}
			}else{
				if(begin!=saveend){
					for(int i=0;i<savescore.size();i++){rec += '\t' + savescore.at(i);}
					cout << savechr << '\t' << savebegin << '\t' << saveend << rec << endl;
					rec="";
					savechr=chr;
					savebegin=begin;
					saveend=end;
					for(int i=0;i<savescore.size();i++){savescore.at(i)=score.at(i);}
				}else{
					test=0;
					for(int i=0;i<savescore.size();i++){if(savescore.at(i)!=score.at(i)){test++;}}
					if(test!=0){
						for(int i=0;i<savescore.size();i++){rec += '\t' + savescore.at(i);}
						cout << savechr << '\t' << savebegin << '\t' << saveend << rec << endl;
						rec="";
						savechr=chr;
						savebegin=begin;
						saveend=end;
						for(int i=0;i<savescore.size();i++){savescore.at(i)=score.at(i);}
					}else{
						saveend=end;
					}
				}
			}
			if(data.eof()){
				stop++;
				for(int i=0;i<savescore.size();i++){rec += '\t' + savescore.at(i);}
				cout << savechr << '\t' << savebegin << '\t' << saveend << rec;
			}
			score.clear();
		}
	}
	data.close();
	return 0;
}
