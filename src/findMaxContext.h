#ifndef FINDMAXCONTEXT_H
#define FINDMAXCONTEXT_H
#include "stl.h"

string SP_reverseSeq(const string &seq)
{
        string new_seq(seq);

        for (unsigned int i=0;i<seq.size();i++){
                new_seq[i] = seq[seq.size()-1-i];
        }
        return (new_seq);
}


map<string, vector<unsigned int> > buildGenomeIndex(string & genomeSeq, int seed_len=8)
{
	map<string, vector<unsigned int> > genomeIndex;
	for (unsigned int i=0;i<genomeSeq.size()-seed_len+1;i++){
		string cur_context = genomeSeq.substr(i,seed_len);
		genomeIndex[cur_context].push_back(i);
	}
	return genomeIndex;
}
vector<int> findMaxContextIndex(string cur_context, vector<string> &ex_ref, vector<vector<double> > &ipd_ref)
{
	int cur_pos = 30;
	int max_left = 0;
	int max_right = 0;
	int max_idx = 0;
	int max_cvg = 0;
	for (int i=0;i<(int)ex_ref.size();i++){
		int j=cur_pos;
		while(cur_context[j]==ex_ref[i][j]){j--;if(j<0) break;}
		if (cur_pos-j>max_left){
			max_left = cur_pos-j;
			j=cur_pos+1;
			while(cur_context[j]==ex_ref[i][j]){j++;if(j>=(int)cur_context.size()) break;}
			max_idx = i;
			max_right = j-cur_pos-1;	
			max_cvg = ipd_ref[i].size();
		}else{
			if (cur_pos-j==max_left){
	                        j=cur_pos+1;
        	                while(cur_context[j]==ex_ref[i][j]){j++;if(j>=(int)cur_context.size()) break;}
                        	if (j-cur_pos-1 > max_right){
                                	max_idx = i;
					max_right = j-cur_pos-1;
					max_cvg = ipd_ref[i].size();
				}else{
					if (j-cur_pos-1 == max_right && ipd_ref[i].size()>max_cvg ){
						max_idx = i;
	                                        max_right = j-cur_pos-1;
        	                                max_cvg = ipd_ref[i].size();
					}			
				}
			}
		}
	}
	//Rprintf("max_left = %d\nmax_right = %d\nmax_idx = %d\n", max_left, max_right, max_idx);
	vector<int> result;
	result.push_back(max_idx);
	result.push_back(max_left);
	result.push_back(max_right);
	result.push_back(max_cvg);
	return result;	
}

void findMaxContext(string & genomeSeq, string & genomeSeqRef,  map<string, list<unsigned int> > & genomeIndex,
			SEXP Ripd, SEXP Rstart, int strand, int strand_ref, int min_cvg, int left_len, int right_len, 
			vector<int> & control_len, vector<vector<double> > & control_data,
			vector<vector<vector<double> > > &ref_data)
{
	int start = INTEGER(Rstart)[0] - 1;
	int seed_len = left_len + right_len + 1;
	list<unsigned int> index_all;
	for (int i=0;i<(int)genomeSeq.size();i++){
		string cur_context;
		if (strand==0){ 
			if (i-left_len<0 || i+right_len>genomeSeq.size()-1) continue;
			if (strand_ref==0)
				cur_context = genomeSeq.substr(i-left_len, left_len+right_len+1);	
			else
				cur_context = SP_reverseSeq(genomeSeq.substr(i-left_len, left_len+right_len+1));
		}else{	
			if (i-right_len<0 || i+left_len>genomeSeq.size()-1) continue;
			if (strand_ref==0)
				cur_context = SP_reverseSeq(genomeSeq.substr(i-right_len, left_len+right_len+1));
			else
				cur_context = genomeSeq.substr(i-right_len, left_len+right_len+1);
		}
		
		// get positions with the same context and find the maximal one
		map<string, list<unsigned int> >::iterator it_idx;
		it_idx = genomeIndex.find(cur_context);
		if (it_idx==genomeIndex.end()) continue;
		
		index_all = it_idx->second;
		list<unsigned int>::iterator it;
		vector <int> cur_len(index_all.size(),0);
		int max_len = 0;
		int j=0;
		for (it=index_all.begin();it!=index_all.end();it++){
			int cur_idx = (int)*it;
			int k=0;
			if (strand==0){
				if (strand_ref==0){
					while(genomeSeq[i-left_len-k]==genomeSeqRef[cur_idx-k]){k++;if (i-left_len-k<0 || cur_idx-k<0 || k>20) break;}
				}else{
					while(genomeSeq[i-left_len-k]==genomeSeqRef[cur_idx+seed_len-1+k]){
						k++;
						if (i-left_len-k<0 || cur_idx+seed_len-1+k>genomeSeqRef.size()-1 || k>20) break;
					}	
				}
			}else{
				if (strand_ref==0){
                                        while(genomeSeq[i+left_len+k]==genomeSeqRef[cur_idx-k]){k++;if (i+left_len+k>genomeSeq.size()-1 || cur_idx-k<0 || k>20) break;}
                                }else{
                                        while(genomeSeq[i+left_len+k]==genomeSeqRef[cur_idx+seed_len-1+k]){
						k++;
                                                if (i+left_len+k>genomeSeq.size()-1 || cur_idx+seed_len-1+k>genomeSeqRef.size()-1 || k>20) break;
					}
                                }
			}
			k--;
			cur_len[j] = seed_len + k;
			if (cur_len[j]>max_len) max_len = cur_len[j];
			j++;				
		}
		// get historical data
		j=0;
                for (it=index_all.begin();it!=index_all.end();it++){
			int cur_idx = (int) *it;
			double * ipd;
			int cvg;
			if (strand_ref==0){
				ipd = REAL(VECTOR_ELT(Ripd, cur_idx+left_len-start));
                                cvg = Rf_length(VECTOR_ELT(Ripd, cur_idx+left_len-start));
			}else{
				ipd = REAL(VECTOR_ELT(Ripd, cur_idx+right_len-start));
                                cvg = Rf_length(VECTOR_ELT(Ripd, cur_idx+right_len-start));	
			}

			if (cvg<min_cvg) continue;
			if (cur_len[j]==max_len){
                        	if (max_len > control_len[i]) control_data[i].clear();
                               	if (max_len >= control_len[i]){
					control_len[i] = max_len;
					for (int t=0;t<cvg;t++) control_data[i].push_back(ipd[t]);
				}else{
					vector<double> ipd_vec;
                                	for (int t=0;t<cvg;t++) ipd_vec.push_back(ipd[t]);
                                	ref_data[i].push_back(ipd_vec);
				}
				// control_len[i] = max_len;
                                //for (int t=0;t<cvg;t++) control_data[i].push_back(ipd[t]);
                       	}else{
				vector<double> ipd_vec;
				for (int t=0;t<cvg;t++) ipd_vec.push_back(ipd[t]);
                                ref_data[i].push_back(ipd_vec);
                        }
			j++;
		}
	}
					
}

void getGenomeIndex(SEXP RgenomeIndex, map<string, list<unsigned int> > &genomeIndex)
{
	SEXP Rcontext = Rf_getAttrib(RgenomeIndex, R_NamesSymbol);	

        for (int i=0;i<Rf_length(RgenomeIndex);i++){
                string cur_context = CHAR(STRING_ELT(Rcontext,i));
		list<unsigned int> index;
        	int *index_tmp = INTEGER(VECTOR_ELT(RgenomeIndex,i));
		int index_len = Rf_length(VECTOR_ELT(RgenomeIndex,i));
		for (int j=0;j<index_len;j++) index.push_back(index_tmp[j]);
		genomeIndex[cur_context] = index;
	}
}

#endif //FINDMAXCONTEXT_H



