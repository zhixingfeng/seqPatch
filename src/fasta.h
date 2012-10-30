/*
fasta.h - Header file for fasta file manipulation
written by JIANG Hui, 
Institute for Computational and Mathematical Engineering, Stanford University
March, 2008 -
*/

#ifndef FASTA_H
///Define this macro to prevent from including this header file more than once.
#define FASTA_H

#include "stl.h"
#include "file_operation.h"
#include "string_operation.h"
class fasta{
public:
	bool verbose;
	vector<string> tags;
	vector<string> sequences;

	fasta(){verbose = true; clear();};
	void clear(){tags.clear(); sequences.clear();};
	string get_sequence(const string tag){
		for (int i = 0; i < (int)tags.size(); i++) {
			if (tag == tags[i]) return sequences[i];
		}
		return "";
	}
	string get_sequence(const string tag, int start, int end){
		for (int i = 0; i < (int)tags.size(); i++) {
			if (tag == tags[i]) return sequences[i].substr(start, end - start);
		}
		return "";
	}
	string get_sequence(const string tag, vector<pair<int, int> > &regions, int separate_N = 0){
		string result = "";
		string seperator = "";
		seperator.insert(seperator.begin(), separate_N, 'N');
		for (int i = 0; i < (int)tags.size(); i++) {
			if (tag == tags[i]) {
				for (int j = 0; j < (int)regions.size(); j++) {
					result += sequences[i].substr(regions[j].first, regions[j].second - regions[j].first);
					if (j < (int)regions.size() - 1 && separate_N > 0) result += seperator;
				}
				break;
			}
		}
		return result;
	}
	bool read_from_file(const string filename, bool check_sequence = true);
	bool write_to_file(FILE *file, const string tag, const string &sequence);
	bool write_to_file(const string filename);
	bool enumerate(const string filename, int read_length = 25);
};

inline bool is_fasta_file(const string file_name){
	if (!file_exists(file_name)) return false;
	ifstream ifs(file_name.c_str());
	string readline;
	bool result = true;

	while(getline(ifs, readline)) {
		if (readline == "") continue;
		if (readline[0] != '>') result = false;
		break;
	}

	ifs.close();
	return result;
}

inline bool fasta::read_from_file(const string filename, bool check_sequence){
	if (!is_fasta_file(filename)) {
		if (verbose) cout << filename << " is not a vaid fasta file.\n";
		return false;
	}

	if (verbose) printf("Loading fasta file %s.\n", filename.c_str());

	char dem = '\n';
	string tag = "";
	string sequence = "";
	string readline;
	ifstream ifs(filename.c_str());
	// check delimiter, \r or \n
	getline(ifs, readline);
	if (readline[readline.size()-1] == '\r'){
		dem = '\r';
		printf("delimiter is \\r \n");
	}
	else{
		dem = '\n';
		printf("delimiter is \\n \n");
	}
	ifs.seekg(0);
	
	int total_len = 0;
	while(getline(ifs, readline)) {
		if (readline == "") continue;
		if (readline[readline.size()-1]=='\r')
			readline.resize(readline.size()-1);
		if (readline[0] == '>') {
			if (sequence != "") {
				tags.push_back(tag);
				sequences.push_back(sequence);
				total_len += (int)sequence.length();
			}
			tag = readline.substr(1);
			sequence = "";
		} else {
			if (check_sequence) {
				for (int i = 0; i < (int)readline.length(); i++) {
					if (!isalpha(readline[i])) {
						printf("%d\n",i);
						if (verbose) printf("failed: error letter %c occurred when reading fasta file %s.\n", readline[i], filename.c_str());
						ifs.close();
						return false;
					}
				}
			}
			sequence += readline;
		}
	}

	if (sequence != "") {
		tags.push_back(tag);
		sequences.push_back(sequence);
		total_len += (int)sequence.length();
	}
	ifs.close();
	if (verbose) printf("%d sequences read, total length is %d.\n",(int) sequences.size(), total_len);
	return true;
}

inline bool fasta::write_to_file(FILE *file, const string tag, const string &sequence) {
	fprintf(file, ">%s\n", tag.c_str());
	for (int j = 0; j < (int)sequence.length(); j+=50) {
		fprintf(file, "%s\n", sequence.substr(j,50).c_str());
	}
	return true;
}

inline bool fasta::write_to_file(const string filename) {
	FILE *file=fopen(filename.c_str(), "wt");
	if (!file){
		cout << "error opening file for writing:" << filename << endl;
		return false;
	}
	if (verbose) printf("writing fasta file %s.\n", filename.c_str());
	int total_len = 0;
	for (int i = 0; i < (int)tags.size(); i++) {
		if (!write_to_file(file, tags[i], sequences[i])) return false;
		total_len += (int)sequences[i].length();
	}
	if (verbose) printf("%d sequences wrote, total length is %d.\n", (int)sequences.size(), total_len);
	fclose(file);
	return true;
}

inline bool is_dna_file(const string file_name){
	if (!file_exists(file_name)) return false;
	ifstream ifs(file_name.c_str());
	string readline;
	bool result = true;

	while(getline(ifs, readline)) {
		readline = trim_space(toupper(readline));
		if (readline == "") continue;
		if (readline[0] != 'A' && readline[0] != 'C' && readline[0] != 'G' && readline[0] != 'T' && readline[0] != 'N' && readline[0] != '.') result = false;
		break;
	}

	ifs.close();
	return result;
}

inline bool read_dna_from_file(const string filename, vector<string> &sequences) {
	if (!is_dna_file(filename)) {
		cout << filename << " is not a vaid DNA sequence file.\n";
		return false;
	}

	printf("Loading DNA sequence file %s.\n", filename.c_str());

	string readline;
	ifstream ifs(filename.c_str());
	int total_len = 0;
	while(getline(ifs, readline)) {
		if (readline == "") continue;
		sequences.push_back(readline);
		total_len += (int)readline.length();
	}

	ifs.close();
	printf("%d sequences read, total length is %d.\n",(int) sequences.size(), total_len);
	return true;
}

inline bool fasta::enumerate(const string filename, int read_length) {
	FILE *file=fopen(filename.c_str(), "wt");
	if (!file){
		cout << "error opening file for writing:" << filename << endl;
		return false;
	}
	if (verbose) printf("writing fasta file %s.\n", filename.c_str());
	int total_sequences = 0, total_len = 0;
	for (int i = 0; i < (int)tags.size(); i++) {
		int l = (int)sequences[i].length();
		for (int j = 0; j <= l - read_length; j++) {
			if (!write_to_file(file, tags[i] + ":" + int2str(j+1), sequences[i].substr(j, read_length))) return false;
			total_len += read_length;
			total_sequences++;
		}
	}
	if (verbose) printf("%d sequences wrote, total length is %d.\n", total_sequences, total_len);
	fclose(file);
	return true;
}

#endif //FASTA_H
