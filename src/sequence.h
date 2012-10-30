/*
	Code DNA sequence. 
	A/a : 001
	G/g : 010
	C/c : 100
	T/t : 000
*/
#ifndef SEQUENCE_H
#define SEQUENCE_H
#include <stdio.h>
#include <string>
#include <vector>
using namespace std;

string reverse_dna_seq(string seq)
{
	string seq_rev(seq.size(),'N');
	for (unsigned int i=0;i<seq.size();i++){
		switch(seq[seq.size() - 1 -i ]){
			case 'A':
				seq_rev[i] = 'T';
				break;
                        case 'C':
                                seq_rev[i] = 'G';
                                break;
                        case 'G':
                                seq_rev[i] = 'C';
                                break;
                        case 'T':
                                seq_rev[i] = 'A';
                                break;
			case 'N':
				seq_rev[i] = 'N';
				break;
		};
	}
	return seq_rev;
}
bool code_dna_seq(const string &seq, vector<int> &code)
{
	bool isOK = true;
	for (unsigned int i=0;i<seq.size();i++){
		switch (seq[i]){
			case 'A': 
				code.push_back(0);
				code.push_back(0);
				code.push_back(1);	
				break;
			case 'a':
                                code.push_back(0);
                                code.push_back(0);
                                code.push_back(1);
                                break;
			case 'G':
                                code.push_back(0);
                                code.push_back(1);
                                code.push_back(0);
                                break;
			case 'g':
                                code.push_back(0);
                                code.push_back(1);
                                code.push_back(0);
                                break;
			case 'C':
                                code.push_back(1);
                                code.push_back(0);
                                code.push_back(0);
                                break;
			case 'c':
                                code.push_back(1);
                                code.push_back(0);
                                code.push_back(0);
                                break;
			case 'T':
                                code.push_back(0);
                                code.push_back(0);
                                code.push_back(0);
                                break;
			case 't':
                                code.push_back(0);
                                code.push_back(0);
                                code.push_back(0);
                                break;
			default:
				isOK = false;
				break;
		}
	}
	if(isOK==false){
		code.clear();
	}
	return isOK;
}

#endif //SEQUENCE_H



