#include <Rcpp.h>
#include <R.h> 
#include <Rinternals.h> 
#include "headers.h"
#include "mathopt.h"
#include "fasta.h"
#include "sampling.h"
#include "hieModelEB_alt.h"
#include "pooling.h"
#include "findMaxContext.h"
#include <math.h>
#include "statistics.h"


using namespace std;
using namespace Rcpp;
string SP_reverseSeq(const string &seq);

struct DataPos
{
	vector<double> IPD;
	vector<double> PW;	
	vector<double> moleculeID;
};

void constructDataByGenome(int *tStart, int *tEnd, SEXP RrefSeq, int * strand, int nRead, 
			 map<string, vector<DataPos> > & ipdByGenome_pos,
			 map<string, vector<DataPos> > & ipdByGenome_neg,
			 map<string, int> & genome_start_pos, map<string, int> & genome_start_neg)
{
	map<string, int> genome_end_pos;
	int genome_size;

        map<string, int> genome_end_neg;
	
	map<string, int>::iterator it_start;
	map<string, int>::iterator it_end;
	
	// get genome size
	for (int i=0;i<nRead;i++){
		const char *refSeq = CHAR(STRING_ELT(RrefSeq,i));
		// for positive reads
		if (strand[i] == 0){
			it_start = genome_start_pos.find(refSeq);
			it_end = genome_end_pos.find(refSeq);
			if (it_start == genome_start_pos.end()){
				genome_start_pos[refSeq] = tStart[i];
			}else{ 	
				if(tStart[i] < genome_start_pos[refSeq])
				genome_start_pos[refSeq] = tStart[i];
			}
			
			if (it_end == genome_end_pos.end()){
				genome_end_pos[refSeq] = tEnd[i];
			}else{
				if(tEnd[i] > genome_end_pos[refSeq])
                        		genome_end_pos[refSeq] = tEnd[i];
			}
		}else{ 
		// for negative reads 
			it_start = genome_start_neg.find(refSeq);
                        it_end = genome_end_neg.find(refSeq);
                        if (it_start == genome_start_neg.end()){
                                genome_start_neg[refSeq] = tStart[i];
                        }else{
                                if(tStart[i] < genome_start_neg[refSeq])
                                genome_start_neg[refSeq] = tStart[i];
                        }
                        
                        if (it_end == genome_end_neg.end()){
                                genome_end_neg[refSeq] = tEnd[i];
                        }else{
                                if(tEnd[i] > genome_end_neg[refSeq])
                                        genome_end_neg[refSeq] = tEnd[i];
                        }	
		}
	}

	// construct data (positive reads)	
	it_start = genome_start_pos.begin();
	it_end = genome_end_pos.begin();
	while(it_start!=genome_start_pos.end()){
		genome_size = it_end->second - it_start->second + 1;
		vector<DataPos> new_data(genome_size);
		ipdByGenome_pos[it_start->first] = new_data;
		it_start++;
		it_end++;
	}
	
	// construct data (negative reads)      
        it_start = genome_start_neg.begin();
        it_end = genome_end_neg.begin();
        while(it_start!=genome_start_neg.end()){
                genome_size = it_end->second - it_start->second + 1;
                vector<DataPos> new_data(genome_size);
                ipdByGenome_neg[it_start->first] = new_data;
                it_start++;
                it_end++;
        }
	
}

vector<int> getAlnsFColIdx(SEXP alnsFNames)
{
	vector<int> colIdx(4,0);
	int len = Rf_length(alnsFNames);
	for (int i=0;i<len;i++){
		string cur_name = CHAR(STRING_ELT(alnsFNames,i));
		if (cur_name=="IPD")
			colIdx[0] = i;
		if (cur_name=="PulseWidth")
                        colIdx[1] = i;
		if (cur_name=="read")
                        colIdx[2] = i;
		if (cur_name=="reference")
                        colIdx[3] = i;
	}
	return colIdx;
}

vector<int> getAlnsIdxColIdx(SEXP RalnsIdxNames)
{
	vector<int> colIdx(4,0);
        int len = Rf_length(RalnsIdxNames);
        for (int i=0;i<len;i++){
                string cur_name = CHAR(STRING_ELT(RalnsIdxNames,i));
                if (cur_name=="fullRefName")
                        colIdx[0] = i;
                if (cur_name=="tStart")
                        colIdx[1] = i;
                if (cur_name=="tEnd")
                        colIdx[2] = i;
                if (cur_name=="alignedStrand")
                        colIdx[3] = i;
		if (cur_name=="moleculeID")
                        colIdx[4] = i;
        }
        return colIdx;

}

int getIdxByName (SEXP data, string query)
{
	SEXP names = Rf_getAttrib(data,R_NamesSymbol);
	int len = Rf_length(names);
	for (int i=0;i<len;i++){
		string cur_name = CHAR(STRING_ELT(names,i));
		if (cur_name == query) return i;
	}
	return -1;
}
RcppExport SEXP getFeaturesAlongGenome(SEXP RalnsF, SEXP RalnsIdx, SEXP RalnsFNames, SEXP RalnsIdxNames, SEXP RuseCCS)
{

	vector<int> alnsF_colIdx = getAlnsFColIdx(RalnsFNames);
	vector<int> alnsIdx_colIdx = getAlnsIdxColIdx(RalnsIdxNames);
	SEXP RrefSeq = VECTOR_ELT(RalnsIdx,alnsIdx_colIdx[0]);
	int *tStart = INTEGER(VECTOR_ELT(RalnsIdx,alnsIdx_colIdx[1]));
	int *tEnd = INTEGER(VECTOR_ELT(RalnsIdx,alnsIdx_colIdx[2]));
	int *strand = INTEGER(VECTOR_ELT(RalnsIdx,alnsIdx_colIdx[3]));
	int *moleculeID = INTEGER(VECTOR_ELT(RalnsIdx,alnsIdx_colIdx[4]));

	int useCCS = INTEGER(RuseCCS)[0];
	int nRead = Rf_length(RalnsF);
	map<string, vector<DataPos> > ipdByGenome_pos;
	map<string, vector<DataPos> > ipdByGenome_neg;
	map<string, int> genome_start_pos;
	map<string, int> genome_start_neg;
	constructDataByGenome(tStart, tEnd, RrefSeq, strand, nRead, ipdByGenome_pos, ipdByGenome_neg, genome_start_pos, genome_start_neg);

	// map IPD in alignment to genome 	
	int cur_pos = 0;
	int cur_idx = 0;
	for (int i=0;i<nRead;i++){
		if ((i+1)%1000==0)
			Rprintf("processed %d reads\r",i+1);
		if (strand[i]==0){
			cur_pos = tStart[i];	
		}else{
			cur_pos = tEnd[i];
		}
		string refSeq = CHAR(STRING_ELT(RrefSeq,i));
		SEXP read_data = VECTOR_ELT(RalnsF,i);
		double *IPD = REAL(VECTOR_ELT(read_data,alnsF_colIdx[0]));
		//double *PW = REAL(VECTOR_ELT(read_data,alnsF_colIdx[1]));
		SEXP R_read = VECTOR_ELT(read_data,alnsF_colIdx[2]);
		SEXP R_ref = VECTOR_ELT(read_data,alnsF_colIdx[3]);
		int read_len = Rf_length(VECTOR_ELT(read_data,0));
		int temp_len = 0;
		int signal_len = 0;
		for (int j=0;j<read_len;j++){
			string read = CHAR(STRING_ELT(R_read,j));
			string ref = CHAR(STRING_ELT(R_ref,j));
			if (ref!="-")
				temp_len++;
			// remove the data besides indels and mismatachs 
			if (j>=1){
				string read_pre = CHAR(STRING_ELT(R_read,j-1));
                                string ref_pre = CHAR(STRING_ELT(R_ref,j-1));
                                if (read=="-"||ref=="-"||read_pre=="-"||ref_pre=="-"||read!=ref)
                                	continue;
			}
			/*if (j==0){
				string read_next = CHAR(STRING_ELT(R_read,j+1));
				string ref_next = CHAR(STRING_ELT(R_ref,j+1));
				if (read=="-"||ref=="-"||read_next=="-"||ref_next=="-")
					continue;
			}else{
				if (j==read_len-1){
					string read_pre = CHAR(STRING_ELT(R_read,j-1));
                               		string ref_pre = CHAR(STRING_ELT(R_ref,j-1));
					if (read=="-"||ref=="-"||read_pre=="-"||ref_pre=="-")
                                       		continue;
				}else{
					string read_next = CHAR(STRING_ELT(R_read,j+1));
                               		string ref_next = CHAR(STRING_ELT(R_ref,j+1));
					string read_pre = CHAR(STRING_ELT(R_read,j-1));
                                       	string ref_pre = CHAR(STRING_ELT(R_ref,j-1));
					if (read=="-"||ref=="-"||read_next=="-"||ref_next=="-"||read_pre=="-"||ref_pre=="-")
                                       		continue;
				}
			}*/
			// remove IPD or PW equals NAN ( at start of a read)
			//if (IPD[j]!=IPD[j] || PW[j]!=PW[j])
			if (IPD[j]!=IPD[j]) continue;
			// record data
			if (strand[i]==0){
				cur_idx = cur_pos + temp_len - 1 - genome_start_pos[refSeq];
				ipdByGenome_pos[refSeq][cur_idx].IPD.push_back(IPD[j]);
        	                //ipdByGenome_pos[refSeq][cur_idx].PW.push_back(PW[j]);
				if (useCCS!=0) ipdByGenome_pos[refSeq][cur_idx].moleculeID.push_back(moleculeID[i]);
	
			}else{
				cur_idx = cur_pos - temp_len + 1 - genome_start_neg[refSeq];
				ipdByGenome_neg[refSeq][cur_idx].IPD.push_back(IPD[j]);
                                //ipdByGenome_neg[refSeq][cur_idx].PW.push_back(PW[j]);
				if (useCCS!=0) ipdByGenome_neg[refSeq][cur_idx].moleculeID.push_back(moleculeID[i]);	
			}
			signal_len++;
		}
	}	
	Rprintf("processed %d reads\n",nRead);
	/*------------return results--------------*/	
	Rprintf("output result\n");
        map<string, map<string, vector<vector<double> > > > result_final;

	map<string, vector<vector<double> > > result_pos_moleculeID;
	map<string, vector<vector<double> > > result_neg_moleculeID;
        // postive ipd
        map<string, vector<DataPos> >::iterator it;
        map<string, vector<vector<double> > > result_pos_ipd;
        for (it=ipdByGenome_pos.begin();it!=ipdByGenome_pos.end();it++){
                vector<vector<double> >tmp;
                vector<vector<double> >tmp_moleculeID;

		string chr = it->first;
                for (unsigned int i=0;i<it->second.size();i++){
                        tmp.push_back(it->second[i].IPD);
			if (useCCS!=0) tmp_moleculeID.push_back(it->second[i].moleculeID);
                }
                result_pos_ipd[chr] = tmp ;
        	result_pos_moleculeID[chr] = tmp_moleculeID;
	}

        // negative ipd
        map<string, vector<vector<double> > > result_neg_ipd;
        for (it=ipdByGenome_neg.begin();it!=ipdByGenome_neg.end();it++){
                vector<vector<double> >tmp;
		vector<vector<double> >tmp_moleculeID;

                string chr = it->first;
                for (unsigned int i=0;i<it->second.size();i++){
                        tmp.push_back(it->second[i].IPD);
			if (useCCS!=0) tmp_moleculeID.push_back(it->second[i].moleculeID); 
               }
                result_neg_ipd[chr] = tmp ;
		result_neg_moleculeID[chr] = tmp_moleculeID;
        }

        // positive pw
        /*map<string, vector<vector<double> > > result_pos_pw;
        for (it=ipdByGenome_pos.begin();it!=ipdByGenome_pos.end();it++){
                vector<vector<double> >tmp;
                string chr = it->first;
                for (unsigned int i=0;i<it->second.size();i++){
                        tmp.push_back(it->second[i].PW);
                }
                result_pos_pw[chr] = tmp ;
        }*/

        // negative pw
        /*map<string, vector<vector<double> > > result_neg_pw;
        for (it=ipdByGenome_neg.begin();it!=ipdByGenome_neg.end();it++){
                vector<vector<double> >tmp;
                string chr = it->first;
                for (unsigned int i=0;i<it->second.size();i++){
                        tmp.push_back(it->second[i].PW);
                }
                result_neg_pw[chr] = tmp ;
        }*/

	if (useCCS==TRUE){
		result_final["moleculeID_pos"] = result_pos_moleculeID;
		result_final["moleculeID_neg"] = result_neg_moleculeID;
	}
        result_final["ipd_pos"] = result_pos_ipd;
        result_final["ipd_neg"] = result_neg_ipd;
        //result_final["pw_pos"] = result_pos_pw;
        //result_final["pw_neg"] = result_neg_pw;

        return Rcpp::wrap(result_final);

}

RcppExport SEXP getGenomeSeq(SEXP Rfastafilename)
{

	string fastafilename = CHAR(STRING_ELT(Rfastafilename,0));
        //Rprintf("load sequences from \"%s\"\n",fastafilename.c_str());
        fasta fastaObj;
	fastaObj.read_from_file(fastafilename);
	
	map<string,string> data;
        Rprintf("%u\n", fastaObj.tags.size() );
	for (unsigned int i=0;i<fastaObj.tags.size();i++){
		Rprintf("%s\n",fastaObj.tags[i].c_str());
		data[fastaObj.tags[i]] = fastaObj.sequences[i];
	}
        map<string,string> data_comp;
	map<string,map<string,string> > result;
	map<string,string>::iterator it;
	for (it=data.begin();it!=data.end();it++){
		string seq_comp(it->second.size(),'N');
		for (unsigned int i=0;i<it->second.size();i++){
			//Rprintf("%d\n",i);
			switch(it->second[i]){
				case 'A':
					seq_comp[i]='T';
					break;
				case 'C':
					seq_comp[i]='G';
					break;
				case 'G':
					seq_comp[i]='C';
					break;
				case 'T':
					seq_comp[i]='A';
					break;
				case 'N':
					seq_comp[i]='N';
					break;
				default:
					Rprintf("error: shouldn't contain letters other than A C G T and N");
					result.clear();
					return Rcpp::wrap(result);			
			};
		}
		data_comp[it->first] = seq_comp;
	}	
	result["pos"] = data;
	result["neg"] = data_comp;
	SEXP rl = Rcpp::wrap(result);
        return rl;
}


RcppExport SEXP getFeaturesByContext(SEXP Rsignal, SEXP RstartPos, SEXP RgenomeSeq, SEXP RleftLen, SEXP RrightLen, SEXP Rstrand)
{
	double *signal = REAL(Rsignal);
        int sig_size = Rf_length(Rsignal);
	int start_pos = INTEGER(RstartPos)[0] - 1;
	string genomeSeq = CHAR(STRING_ELT(RgenomeSeq,0));
	int left_len = INTEGER(RleftLen)[0];
	int right_len = INTEGER(RrightLen)[0];
	int strand = INTEGER(Rstrand)[0];
	map<string, vector<double> > signal_bycontext;
	
	int num=0;	
	for (int i=0;i<sig_size;i++){
		if (signal[i]!=signal[i]){
			continue;
		}
		if (strand==0){
			if (start_pos + i - left_len < 0  || start_pos + i + right_len >= sig_size)
				continue;
			string context = genomeSeq.substr(start_pos + i - left_len, left_len + right_len + 1);
			signal_bycontext[context].push_back(signal[i]);	
		}else{
			if (start_pos + i - right_len < 0  || start_pos + i + left_len >= sig_size)
                                continue;
			string context = genomeSeq.substr(start_pos + i - right_len, left_len + right_len + 1);
			// reverse referece sequence
			string context_rev = reverse_dna_seq(context);
			signal_bycontext[context_rev].push_back(signal[i]);	
		}
		num++;
		if ((i+1)%10000==0)
			Rprintf("processed %d positions\r",i+1);	
	}
	Rprintf("processed %d positions\n",sig_size);
	Rprintf("non NA is %d\n",num);	
	
	SEXP rl = Rcpp::wrap(signal_bycontext);
	return rl;
}

RcppExport SEXP correctContextEffect(SEXP RalnsF, SEXP RcontextEffect, SEXP Rcontexts, SEXP Rleft_len, SEXP Rright_len)
{
	// load context effect table
	double *contextEffet_IPD = REAL(RcontextEffect);
	int context_size = Rf_length(RcontextEffect);
	map<string, double> contextEffect;
	for (int i=0;i<context_size;i++){
		string tmp_context = CHAR(STRING_ELT(Rcontexts,i));
		contextEffect.insert( pair<string,double> (tmp_context, contextEffet_IPD[i]) );
	}
	int left_len = INTEGER(Rleft_len)[0];
	int right_len = INTEGER(Rright_len)[0];
	int context_len = left_len + right_len + 1;	
	// load alignments
	int read_num = Rf_length(RalnsF);	
	for (int i=0;i<read_num;i++){
		if ((i+1)%1000==0)
			Rprintf("processed %d reads\r",i+1);
		double *IPD = REAL(VECTOR_ELT( VECTOR_ELT(RalnsF,i),0));
		int len = Rf_length(VECTOR_ELT( VECTOR_ELT(RalnsF,i),0));
		SEXP Rcur_context = VECTOR_ELT( VECTOR_ELT(RalnsF,i),3);
		string cur_context(Rf_length(Rcur_context),'N');
		for (unsigned int j=0;j<cur_context.size();j++){
			cur_context[j] = CHAR(STRING_ELT(Rcur_context,j))[0];
		}
		
		// correct IPD by context
		for (int j=0;j<len;j++){
			if (j<left_len||j>len-1-right_len){
				IPD[j] = sqrt(-1);
				continue;
			}
			if (IPD[j]!=IPD[j] || cur_context[j]=='-' || cur_context[j]=='N'){
				IPD[j] = sqrt(-1);
				continue;
			}
			// find context 
			string the_context(context_len,'N');
			int signal_pos = left_len;
			the_context[signal_pos] = cur_context[j];
			int k=1;
			int tmp_count=1;
			bool is_N = false;
			// search upstream
			while (tmp_count<=left_len){
				if (cur_context[j-k]=='-'){
					k++;
					continue;
				}
				if (cur_context[j-k]=='N'){
					is_N = true;
					break;
				}
				if (j-k<0){
                                        is_N = true;
                                        break;
                                }
				the_context[signal_pos-tmp_count] = cur_context[j-k];
				tmp_count++;
				k++;
			}
			// search downstream
			k=1;
			tmp_count=1;
			while (tmp_count<=right_len){
                                if (cur_context[j+k]=='-'){
                                        k++;
                                        continue;
                                }
                                if (cur_context[j+k]=='N'){
                                        is_N = true;
                                        break;
                                }
				if (j+k>=len){
                                        is_N = true;
                                        break;
                                }
                                the_context[signal_pos+tmp_count] = cur_context[j+k];
                                tmp_count++;
                                k++;
                        }
	
			//Rprintf("%s\n",the_context.c_str());
			if (is_N==true){
				IPD[j] = sqrt(-1);
				continue;
			}
			// correct context effect
			double IPD_pred = contextEffect[the_context];
			if (IPD_pred==0)
				IPD[j] = sqrt(-1);
			else
				IPD[j] = IPD[j]/IPD_pred;		
		}
	}
	Rprintf("processed %d reads\n",read_num);

	//return Rcpp::wrap(contextEffect); 
	SEXP rl=R_NilValue;
	return rl;
	 
}
string findContextForRead(const string & cur_context, int idx, int left_len, int right_len)
{
	int j = idx;
	string the_context(left_len + right_len + 1, 'N');
	int signal_pos = left_len;
        the_context[signal_pos] = cur_context[j];
        int k=1;
        int tmp_count=1;
        bool is_N = false;
	
	// search upstream
	while (tmp_count<=left_len){
        	if (j-k<0){
                        is_N = true;
                        break;
                }

		if (cur_context[j-k]=='-'){
                	k++;
                        continue;
                }
                if (cur_context[j-k]=='N'){
                        is_N = true;
                        break;
                }
                the_context[signal_pos-tmp_count] = cur_context[j-k];
                tmp_count++;
                k++;
        }

	// search downstream
	k=1;
        tmp_count=1;
        while (tmp_count<=right_len){
        	if (j+k>=int(cur_context.size())){
                        is_N = true;
                        break;
                }
		if (cur_context[j+k]=='-'){
                	k++;
                        continue;
                }
                if (cur_context[j+k]=='N'){
                	is_N = true;
                        break;
                }
                the_context[signal_pos+tmp_count] = cur_context[j+k];
                tmp_count++;
                k++;
	}
	if (is_N==true){
        	return "";
	}else{	
		return the_context;	
	}
}

RcppExport SEXP getContextEffectForRead(SEXP RalnsF, SEXP RcontextEffect, SEXP Rcontexts, SEXP Rleft_len, SEXP Rright_len)
{
	int left_len = INTEGER(Rleft_len)[0];
        int right_len = INTEGER(Rright_len)[0];
	int read_num = Rf_length(RalnsF);
	double *contextEffect_IPD = REAL(RcontextEffect);
	int context_size = Rf_length(RcontextEffect);
	map<string, double> contextEffect;
	for (int i=0;i<context_size;i++){
		string tmp_context = CHAR(STRING_ELT(Rcontexts,i));
		contextEffect.insert( pair<string,double> (tmp_context, contextEffect_IPD[i]) );
	}
	
	// scan each subread and get context effect for each position, [-1,1] around indels are omitted 	
	vector<vector<double> > contextEffectForRead_all;
	for (int i=0;i<read_num;i++){
                if ((i+1)%1000==0)
                        Rprintf("processed %d reads\r",i+1);
                double *IPD = REAL(VECTOR_ELT( VECTOR_ELT(RalnsF,i),0));
                int len = Rf_length(VECTOR_ELT( VECTOR_ELT(RalnsF,i),0));
                SEXP Rcur_read = VECTOR_ELT( VECTOR_ELT(RalnsF,i),2);
                string cur_read(Rf_length(Rcur_read),'N');
                for (unsigned int j=0;j<cur_read.size();j++){
                        cur_read[j] = CHAR(STRING_ELT(Rcur_read,j))[0];
                }
                SEXP Rcur_context = VECTOR_ELT( VECTOR_ELT(RalnsF,i),3);
                string cur_context(Rf_length(Rcur_context),'N');
                for (unsigned int j=0;j<cur_context.size();j++){
                        cur_context[j] = CHAR(STRING_ELT(Rcur_context,j))[0];
                }
		
		vector<double> contextEffectForRead(len,-1);		
		for (int j=0;j<len;j++){
        		if (j<left_len||j>len-1-right_len){
				contextEffectForRead[j] = sqrt(-1);
     	        		continue;
                	}
                	if (IPD[j]!=IPD[j] || cur_context[j]=='-' || cur_context[j]=='N' ||
 					cur_context[j-1]=='-' || cur_context[j+1]=='-' || 
					cur_read[j]=='-' || cur_read[j-1]=='-' || cur_read[j+1]=='-'){
				contextEffectForRead[j] = sqrt(-1);
                        	continue;
                	}
			string the_context = findContextForRead(cur_context,j,left_len, right_len);
			if (the_context==""){
				contextEffectForRead[j] = sqrt(-1);
                        	continue;
			}
			//Rprintf("%d : %s\n", j+1, the_context.c_str());
			double IPD_pred = contextEffect[the_context];
                	if (IPD_pred==0)
                		contextEffectForRead[j] = sqrt(-1);
                	else
                        	contextEffectForRead[j] = IPD_pred;	
		}
		contextEffectForRead_all.push_back(contextEffectForRead);
	}
	Rprintf("processed %d reads\n", read_num);
	SEXP rl = Rcpp::wrap(contextEffectForRead_all);
				
	//SEXP rl=R_NilValue;
        return rl;		
}


RcppExport SEXP getContextEffect(SEXP RalnsF, SEXP Rleft_len, SEXP Rright_len)
{
	int left_len = INTEGER(Rleft_len)[0];
	int right_len = INTEGER(Rright_len)[0];
	map<string, vector<double> > contextEffect;
	int read_num = Rf_length(RalnsF);	
	int context_len = left_len + right_len + 1;	
	for (int i=0;i<read_num;i++){
                if ((i+1)%1000==0)
                        Rprintf("processed %d reads\r",i+1);
                double *IPD = REAL(VECTOR_ELT( VECTOR_ELT(RalnsF,i),0));
                int len = Rf_length(VECTOR_ELT( VECTOR_ELT(RalnsF,i),0));
                SEXP Rcur_read = VECTOR_ELT( VECTOR_ELT(RalnsF,i),2);
        	string cur_read(Rf_length(Rcur_read),'N');
        	for (unsigned int j=0;j<cur_read.size();j++){
                	cur_read[j] = CHAR(STRING_ELT(Rcur_read,j))[0];
        	}
		SEXP Rcur_context = VECTOR_ELT( VECTOR_ELT(RalnsF,i),3);
                string cur_context(Rf_length(Rcur_context),'N');
                for (unsigned int j=0;j<cur_context.size();j++){
                        cur_context[j] = CHAR(STRING_ELT(Rcur_context,j))[0];
                }
			
		// get context effect
		for (int j=0;j<len;j++){
                        if (j<left_len||j>len-1-right_len){
                                continue;
                        }
                        if (IPD[j]!=IPD[j] || cur_context[j]=='-' || cur_context[j]=='N' ||
                                	cur_context[j-1]=='-' || cur_context[j+1]=='-' ||
                                	cur_read[j]=='-' || cur_read[j-1]=='-' || cur_read[j+1]=='-'){
                                
				continue;
                        }
			
			// find context
			string the_context = findContextForRead(cur_context,j,left_len, right_len);
			if (the_context==""){
				continue;
			}
			// get context effect
			contextEffect[the_context].push_back(IPD[j]);
		}	
	}
	Rprintf("processed %d reads\n",read_num);
	SEXP rl=Rcpp::wrap(contextEffect);
	//SEXP rl=R_NilValue;
        return rl;	
}

RcppExport SEXP getContextEffectByPos(SEXP Ripd_pos, SEXP Ripd_neg, SEXP RgenomeSeq_pos, SEXP RgenomeSeq_neg, SEXP Rstart_pos, SEXP Rstart_neg, SEXP Rleft_len, SEXP Rright_len)
{
	// parse arguments
	int start_pos = INTEGER(Rstart_pos)[0];
	int start_neg = INTEGER(Rstart_neg)[0];
	start_pos--;
	start_neg--;
	int left_len = INTEGER(Rleft_len)[0];
        int right_len = INTEGER(Rright_len)[0];
	const string genomeSeq_pos = CHAR(STRING_ELT(RgenomeSeq_pos,0));
	const string genomeSeq_neg = CHAR(STRING_ELT(RgenomeSeq_neg,0));	

	map<const string, vector<vector<double> > > context_effect;
	// get data from forward strand	
	Rprintf("get ipd from forward strand.\n");
	int Ripd_pos_len = Rf_length(Ripd_pos);
	for (int i=0;i<Ripd_pos_len;i++){
		if (i%10000==0){
			Rprintf("processed %d positions\r",i);
		}
		if (start_pos + i - left_len<0 || start_pos + i + right_len >= Ripd_pos_len 
			|| start_pos + i - 30 <0 || start_pos + i + 10 >=Ripd_pos_len)
			continue;
		double * cur_ipd_T =  REAL( VECTOR_ELT(Ripd_pos,i) );
		int cur_ipd_T_len = Rf_length(VECTOR_ELT(Ripd_pos,i));
		vector<double> cur_ipd(cur_ipd_T_len,sqrt(-1));
		for (int j=0;j<cur_ipd_T_len;j++)
			cur_ipd[j] = cur_ipd_T[j];
		const string context = genomeSeq_pos.substr(start_pos+i-left_len, left_len + right_len + 1);	
		context_effect[context].push_back(cur_ipd); 
	}				
	Rprintf("processed %d positions\n",Ripd_pos_len);

	// get data from backward strand
	Rprintf("get ipd from backward strand.\n");
        int Ripd_neg_len = Rf_length(Ripd_neg);
        for (int i=0;i<Ripd_neg_len;i++){
                if (i%10000==0){
                        Rprintf("processed %d positions\r",i);
                }
                if (start_neg + i - right_len<0 || start_neg + i + left_len >= Ripd_neg_len || 
			start_neg + i - 10 <0 || start_neg + i + 30 >= Ripd_neg_len)
                        continue;
                double * cur_ipd_T =  REAL( VECTOR_ELT(Ripd_neg,i) );
                int cur_ipd_T_len = Rf_length(VECTOR_ELT(Ripd_neg,i));
                vector<double> cur_ipd(cur_ipd_T_len,sqrt(-1));
                for (int j=0;j<cur_ipd_T_len;j++)
                        cur_ipd[j] = cur_ipd_T[j];
                const string context = SP_reverseSeq(genomeSeq_neg.substr(start_neg + i-right_len, left_len + right_len + 1));
                context_effect[context].push_back(cur_ipd);
        }
        Rprintf("processed %d positions\n",Ripd_neg_len);
	
	SEXP RcontextEffect = Rcpp::wrap(context_effect);
	return RcontextEffect;	
}

RcppExport SEXP getContextExByPos(SEXP Ripd_pos, SEXP Ripd_neg, SEXP RgenomeSeq_pos, SEXP RgenomeSeq_neg, SEXP Rstart_pos, SEXP Rstart_neg, SEXP Rleft_len, SEXP Rright_len)
{
        // parse arguments
        int start_pos = INTEGER(Rstart_pos)[0];
        int start_neg = INTEGER(Rstart_neg)[0];
        start_pos--;
        start_neg--;
        int left_len = INTEGER(Rleft_len)[0];
        int right_len = INTEGER(Rright_len)[0];
        const string genomeSeq_pos = CHAR(STRING_ELT(RgenomeSeq_pos,0));
        const string genomeSeq_neg = CHAR(STRING_ELT(RgenomeSeq_neg,0));

        map<const string, vector<string> > context_extended;
        // get data from forward strand
        Rprintf("get ipd from forward strand.\n");
        int Ripd_pos_len = Rf_length(Ripd_pos);
        for (int i=0;i<Ripd_pos_len;i++){
                if (i%10000==0){
                        Rprintf("processed %d positions\r",i);
                }
                if (start_pos + i - left_len<0 || start_pos + i + right_len >= Ripd_pos_len
                        || start_pos + i - 30 <0 || start_pos + i + 10 >=Ripd_pos_len)
                        continue;
                double * cur_ipd_T =  REAL( VECTOR_ELT(Ripd_pos,i) );
                int cur_ipd_T_len = Rf_length(VECTOR_ELT(Ripd_pos,i));
                vector<double> cur_ipd(cur_ipd_T_len,sqrt(-1));
                for (int j=0;j<cur_ipd_T_len;j++)
                        cur_ipd[j] = cur_ipd_T[j];
                const string context = genomeSeq_pos.substr(start_pos+i-left_len, left_len + right_len + 1);
                context_extended[context].push_back(genomeSeq_pos.substr(start_pos+i-30, 41));
        }
        Rprintf("processed %d positions\n",Ripd_pos_len);

        // get data from backward strand
        Rprintf("get ipd from backward strand.\n");
        int Ripd_neg_len = Rf_length(Ripd_neg);
        for (int i=0;i<Ripd_neg_len;i++){
                if (i%10000==0){
                        Rprintf("processed %d positions\r",i);
                }
                if (start_neg + i - right_len<0 || start_neg + i + left_len >= Ripd_neg_len ||
                        start_neg + i - 10 <0 || start_neg + i + 30 >= Ripd_neg_len)
                        continue;
                double * cur_ipd_T =  REAL( VECTOR_ELT(Ripd_neg,i) );
                int cur_ipd_T_len = Rf_length(VECTOR_ELT(Ripd_neg,i));
                vector<double> cur_ipd(cur_ipd_T_len,sqrt(-1));
                for (int j=0;j<cur_ipd_T_len;j++)
                        cur_ipd[j] = cur_ipd_T[j];
                const string context = SP_reverseSeq(genomeSeq_neg.substr(start_neg + i-right_len, left_len + right_len + 1));
                context_extended[context].push_back(SP_reverseSeq(genomeSeq_neg.substr(start_neg+i-10, 41)));
        }
        Rprintf("processed %d positions\n",Ripd_neg_len);

        SEXP RcontextExtended = Rcpp::wrap(context_extended);
        return RcontextExtended;

}
RcppExport SEXP getContextEffectByPosMerge(SEXP RcontextEffect, SEXP RcontextEffectRefSeq, SEXP RcontextEffect_names, SEXP RcontextEffectRefSeq_names)
{
	map<string, vector<vector<double> > > contextEffect;
	map<string, vector<vector<double> > > contextEffectRefSeq;
	int contextEffect_len = Rf_length(RcontextEffect);
	int contextEffectRefSeq_len = Rf_length(RcontextEffectRefSeq);
	
	// parse arguments	
	for (int i=0;i<contextEffect_len;i++){
	const string motif = CHAR(STRING_ELT(RcontextEffect_names,i));
	int cur_motif_len = Rf_length(VECTOR_ELT(RcontextEffect,i));
		for (int j=0;j<cur_motif_len;j++){
			double * cur_ipd_T = REAL(VECTOR_ELT(VECTOR_ELT(RcontextEffect,i),j));
			int cur_ipd_T_len = Rf_length(VECTOR_ELT(VECTOR_ELT(RcontextEffect,i),j));
			vector<double> cur_ipd(cur_ipd_T_len, sqrt(-1));
			for (int k=0;k<cur_ipd_T_len;k++){
				cur_ipd[k] = cur_ipd_T[k];
			}
			contextEffect[motif].push_back(cur_ipd);
		}
	}
	for (int i=0;i<contextEffectRefSeq_len;i++){
                const string motif = CHAR(STRING_ELT(RcontextEffectRefSeq_names,i));
                int cur_motif_len = Rf_length(VECTOR_ELT(RcontextEffectRefSeq,i));
                for (int j=0;j<cur_motif_len;j++){
                        double * cur_ipd_T = REAL(VECTOR_ELT(VECTOR_ELT(RcontextEffectRefSeq,i),j));
                        int cur_ipd_T_len = Rf_length(VECTOR_ELT(VECTOR_ELT(RcontextEffectRefSeq,i),j));
                        vector<double> cur_ipd(cur_ipd_T_len, sqrt(-1));
                        for (int k=0;k<cur_ipd_T_len;k++){
                                cur_ipd[k] = cur_ipd_T[k];
                        }
                        contextEffectRefSeq[motif].push_back(cur_ipd);
                }
        }

	// merge context effect
	for (int i=0;i<contextEffectRefSeq_len;i++){
                const string motif = CHAR(STRING_ELT(RcontextEffectRefSeq_names,i));
                vector<vector<double> > cur_contextEffectRefSeq = contextEffectRefSeq[motif];        
		for (int j=0; j<cur_contextEffectRefSeq.size(); j++){
			contextEffect[motif].push_back(cur_contextEffectRefSeq[j]);
		}
        }
	
	SEXP rl = Rcpp::wrap(contextEffect);
        return rl;
	
}
RcppExport SEXP getContextExByPosMerge(SEXP RcontextEx, SEXP RcontextExRefSeq, SEXP RcontextEx_names, SEXP RcontextExRefSeq_names)
{
        map<string, vector<string> > contextEx;
        map<string, vector<string> > contextExRefSeq;

        // parse arguments
        for (int i=0;i<Rf_length(RcontextEx);i++){
                string cur_context = CHAR(STRING_ELT(RcontextEx_names,i));
                SEXP Rcur_Ex = VECTOR_ELT(RcontextEx,i);
                for (int j=0;j<Rf_length(Rcur_Ex);j++){
                        contextEx[cur_context].push_back(CHAR(STRING_ELT(Rcur_Ex,j)));
                }
        }
	
	for (int i=0;i<Rf_length(RcontextExRefSeq);i++){
                string cur_context = CHAR(STRING_ELT(RcontextExRefSeq_names,i));
                SEXP Rcur_Ex = VECTOR_ELT(RcontextExRefSeq,i);
                for (int j=0;j<Rf_length(Rcur_Ex);j++){
                        contextExRefSeq[cur_context].push_back(CHAR(STRING_ELT(Rcur_Ex,j)));
                }
        }

        // merge context effect
        for (int i=0;i<Rf_length(RcontextExRefSeq);i++){
                const string motif = CHAR(STRING_ELT(RcontextExRefSeq_names,i));
                vector<string > cur_contextExRefSeq = contextExRefSeq[motif];
                for (int j=0; j<cur_contextExRefSeq.size(); j++){
                        contextEx[motif].push_back(cur_contextExRefSeq[j]);
                }
        }

        SEXP rl = Rcpp::wrap(contextEx);
        return rl;

}

RcppExport SEXP mergeContextEffectSummary(SEXP RcontextEffect_1, SEXP RcontextEffect_2, SEXP Rcontext_1, SEXP Rcontext_2)
{
	// load context effect 1
	map<string, map<string, vector<double> > > contextEffect_1;
        if (!Rf_isNull(RcontextEffect_1)){
                for (int i=0;i<Rf_length(RcontextEffect_1);i++){
                        string cur_context = CHAR(STRING_ELT(Rcontext_1,i));
                        map<string, vector<double > > cur_data;
                        vector<double> tmp;
                        cur_data["mean"] = tmp; cur_data["var"] = tmp; cur_data["len"] = tmp;
                        SEXP cur_Rdata = VECTOR_ELT(RcontextEffect_1,i);
                        double *dat;
                        // get mean
                        dat = REAL(VECTOR_ELT(cur_Rdata, getIdxByName(cur_Rdata, "mean") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,0)); j++){cur_data["mean"].push_back(dat[j]);}
                        // get var
                        dat = REAL(VECTOR_ELT(cur_Rdata, getIdxByName(cur_Rdata, "var") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,1)); j++){cur_data["var"].push_back(dat[j]);}
                        // get len
                        dat = REAL(VECTOR_ELT(cur_Rdata, getIdxByName(cur_Rdata, "len") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,2)); j++){cur_data["len"].push_back(dat[j]);}
                        contextEffect_1[cur_context] = cur_data;
                }
        }
	
	// load context effect 2
        map<string, map<string, vector<double> > > contextEffect_2;
        if (!Rf_isNull(RcontextEffect_2)){
                for (int i=0;i<Rf_length(RcontextEffect_2);i++){
                        string cur_context = CHAR(STRING_ELT(Rcontext_2,i));
                        map<string, vector<double > > cur_data;
                        vector<double> tmp;
                        cur_data["mean"] = tmp; cur_data["var"] = tmp; cur_data["len"] = tmp;
                        SEXP cur_Rdata = VECTOR_ELT(RcontextEffect_2,i);
                        double *dat;
                        // get mean
                        dat = REAL(VECTOR_ELT(cur_Rdata, getIdxByName(cur_Rdata, "mean") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,0)); j++){cur_data["mean"].push_back(dat[j]);}
                        // get var
                        dat = REAL(VECTOR_ELT(cur_Rdata, getIdxByName(cur_Rdata, "var") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,1)); j++){cur_data["var"].push_back(dat[j]);}
                        // get len
                        dat = REAL(VECTOR_ELT(cur_Rdata, getIdxByName(cur_Rdata, "len") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,2)); j++){cur_data["len"].push_back(dat[j]);}
                        contextEffect_2[cur_context] = cur_data;
                }
        }
	
	// merge contextEffect_1 and contextEffect_2	
	int Rcontext_2_len = Rf_length(Rcontext_2);
	
	for (int i=0;i<Rcontext_2_len;i++){
                if ((i+1)%1000==0) Rprintf("merged %d contexts\r", i+1);
		const string motif = CHAR(STRING_ELT(Rcontext_2,i));
                map<string, vector<double> > &cur_contextEffect = contextEffect_2[motif];
                // merge mean
		for (int j=0; j<cur_contextEffect["mean"].size(); j++){
			contextEffect_1[motif]["mean"].push_back(cur_contextEffect["mean"][j]);
                }
		// merge var
		for (int j=0; j<cur_contextEffect["var"].size(); j++){
                        contextEffect_1[motif]["var"].push_back(cur_contextEffect["var"][j]);
                }
		// merge len
		for (int j=0; j<cur_contextEffect["len"].size(); j++){
                        contextEffect_1[motif]["len"].push_back(cur_contextEffect["len"][j]);
                }
        }
	Rprintf("merged %d contexts\n", Rcontext_2_len);
	SEXP rl = Rcpp::wrap(contextEffect_1);
	return rl;			
}

RcppExport SEXP correctTrendEffect(SEXP RalnsF, SEXP RtrendEffect)
{
	if (Rf_length(RalnsF) != Rf_length(RtrendEffect)){
		Rprintf("Error: length(alnsF and trendEffect should be the same.\n)");
		return R_NilValue;
	}
	int read_num = Rf_length(RalnsF);
	for (int i=0;i<read_num;i++){
		double *IPD = REAL(VECTOR_ELT( VECTOR_ELT(RalnsF,i),0));
                double *trendEffect = REAL(VECTOR_ELT(RtrendEffect,i));
		int len = Rf_length(VECTOR_ELT( VECTOR_ELT(RalnsF,i),0));
		if (len!= Rf_length(VECTOR_ELT(RtrendEffect,i))){
			Rprintf("Error: length of signal and trend effect should be the same.\n");
			return R_NilValue;
		}
		for (int j=0;j<len;j++){
			IPD[j] = IPD[j]/trendEffect[j];
		}	
	}
	SEXP rl=R_NilValue;
	return rl;
}

RcppExport SEXP getTimeBias(SEXP RalnsF, SEXP RstartTime, SEXP RstepLen, SEXP RtimeBinLen)
{
	double stepLen = REAL(RstepLen)[0];
	int timeBinLen = INTEGER(RtimeBinLen)[0]; 
	int numRead = Rf_length(RalnsF); 
	vector <vector<double> > IPDbyBin(timeBinLen, vector<double>(0,0));
	//Rprintf("%d\n",timeBinLen);	

	for (int i=0;i<numRead;i++){
		if ((i+1)%1000==0)
			Rprintf("processed %d reads\r", i+1);
		int len = Rf_length(VECTOR_ELT(VECTOR_ELT(RalnsF,i),0));
		if (len != Rf_length(VECTOR_ELT(RstartTime,i))){
			Rprintf("error: length of start time and IPD should be the same.\n");
			return R_NilValue;
		}
		//Rprintf("%d\n",len);
		double *IPD = REAL(VECTOR_ELT(VECTOR_ELT(RalnsF,i),0));
		double *sTime = REAL(VECTOR_ELT(RstartTime,i));
		for (int j=0;j<len;j++){
			if (sTime[j]!=sTime[j] || IPD[j]!=IPD[j])
				continue;
			int idx_bin = int(sTime[j]/stepLen);
			if (idx_bin >= timeBinLen)
				idx_bin = timeBinLen-1;
			//Rprintf("%d %f %f\n",idx_bin,sTime[j],stepLen);
			IPDbyBin[idx_bin].push_back(IPD[j]);
		} 
	}
	Rprintf("processed %d reads\n", numRead);
	SEXP rl = Rcpp::wrap(IPDbyBin);	
	//SEXP rl=R_NilValue;
        return rl;

}

RcppExport SEXP buildGenomeIndexR(SEXP RgenomeSeq, SEXP Rseed_len)
{
	string genomeSeq = CHAR(STRING_ELT(RgenomeSeq,0));	
	int seed_len = INTEGER(Rseed_len)[0];
        return Rcpp::wrap(buildGenomeIndex(genomeSeq, seed_len));
}

RcppExport SEXP buildControl_NC(SEXP RgenomeSeq, SEXP RgenomeSeqRef_pos, SEXP RgenomeSeqRef_neg, SEXP Ripd_pos, SEXP Ripd_neg, SEXP Rstart_pos, SEXP Rstart_neg, SEXP RgenomeSeqRefIndex_pos, SEXP RgenomeSeqRefIndex_neg, SEXP Rstrand, SEXP Rmin_cvg, SEXP Rleft_len, SEXP Rright_len)
{	
	// pasre arguments
	string genomeSeq = CHAR(STRING_ELT(RgenomeSeq,0));
	int strand = INTEGER(Rstrand)[0];
	int min_cvg = INTEGER(Rmin_cvg)[0];
	int left_len = INTEGER(Rleft_len)[0];
	int right_len = INTEGER(Rright_len)[0];	
		
	/*---------- get psuedo control and historical data ---------*/
	// init
	vector<int> control_len(genomeSeq.size(),0);
	vector<vector<double> > control_data;
	vector<vector<vector<double> > > ref_data;
	for (unsigned int i=0; i<genomeSeq.size(); i++) {vector<double> tmp1; vector<vector<double> > tmp2; control_data.push_back(tmp1); ref_data.push_back(tmp2); }
	
	// get control and ref data from historical data	
	map<string, list<unsigned int> > genomeIndex;
	string genomeSeqRef;
	int num_chr = Rf_length(RgenomeSeqRef_pos);
	for (int i=0;i<num_chr;i++){
		// forward strand
		genomeIndex.clear();
		getGenomeIndex( VECTOR_ELT(RgenomeSeqRefIndex_pos,i), genomeIndex );
		genomeSeqRef = CHAR( STRING_ELT(VECTOR_ELT(RgenomeSeqRef_pos,i),0) );
		findMaxContext(genomeSeq, genomeSeqRef, genomeIndex, VECTOR_ELT(Ripd_pos,i), 
				VECTOR_ELT(Rstart_pos,i),  strand, 0, min_cvg, left_len, right_len,
				 control_len, control_data, ref_data);
		// backward strand
		genomeIndex.clear();
                getGenomeIndex( VECTOR_ELT(RgenomeSeqRefIndex_neg,i), genomeIndex );
                genomeSeqRef = CHAR( STRING_ELT(VECTOR_ELT(RgenomeSeqRef_neg,i),0) );
                findMaxContext(genomeSeq, genomeSeqRef, genomeIndex, VECTOR_ELT(Ripd_neg,i),
                                VECTOR_ELT(Rstart_neg,i),  strand, 1, min_cvg, left_len, right_len,
                                 control_len, control_data, ref_data);
	}

	return Rcpp::List::create(Rcpp::Named("control_len")=Rcpp::wrap(control_len), Rcpp::Named("control_data")=Rcpp::wrap(control_data),
				 Rcpp::Named("ref_data")=Rcpp::wrap(ref_data));	
}

RcppExport SEXP hieModel_NC (SEXP Ripd_native, SEXP Rgenome_start_native, SEXP RgenomeSeq,  SEXP RcontextEffect,
			 SEXP Rcontext, SEXP Rstrand, SEXP Rleft_len, SEXP Rright_len, SEXP Rmethod, SEXP Ris_cvg_all)
{
	// load context effect
	map<string, map<string, vector<double> > > contextEffect;
        if (!Rf_isNull(RcontextEffect)){
                for (int i=0;i<Rf_length(RcontextEffect);i++){
                        string cur_context = CHAR(STRING_ELT(Rcontext,i));
                        map<string, vector<double > > cur_data;
                        vector<double> tmp;
                        cur_data["mean"] = tmp; cur_data["var"] = tmp; cur_data["len"] = tmp;
                        SEXP cur_Rdata = VECTOR_ELT(RcontextEffect,i);
                        double *dat;
                        // get mean
                        dat = REAL(VECTOR_ELT(cur_Rdata,getIdxByName(cur_Rdata, "mean") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,0)); j++){cur_data["mean"].push_back(dat[j]);}
                        // get var
                        dat = REAL(VECTOR_ELT(cur_Rdata,getIdxByName(cur_Rdata, "var") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,1)); j++){cur_data["var"].push_back(dat[j]);}
                        // get len
                        dat = REAL(VECTOR_ELT(cur_Rdata,getIdxByName(cur_Rdata, "len") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,2)); j++){cur_data["len"].push_back(dat[j]);}
                        contextEffect[cur_context] = cur_data;
                }
        }

	/*map<string, vector<vector<double> > > contextEffect;
        if (!Rf_isNull(RcontextEffect)){
                for (int i=0;i<Rf_length(RcontextEffect);i++){
                        string cur_context = CHAR(STRING_ELT(Rcontext,i));
                        vector<vector<double> > ipd;
                        for (int j=0;j<Rf_length(VECTOR_ELT(RcontextEffect,i)); j++) {vector<double> tmp; ipd.push_back(tmp);}
                        for (int j=0;j<Rf_length(VECTOR_ELT(RcontextEffect,i)); j++){
                                double * cur_ipd = REAL(VECTOR_ELT(VECTOR_ELT(RcontextEffect,i),j));
                                int cur_cvg = Rf_length(VECTOR_ELT(VECTOR_ELT(RcontextEffect,i),j));
                                for (int k=0; k<cur_cvg; k++ ){
                                        ipd[j].push_back(cur_ipd[k]);
                                }
                        }
                        contextEffect[cur_context] = ipd;
                }
        }*/

	// load data and parameters
	int genome_start_native = INTEGER(Rgenome_start_native)[0];
	string genomeSeq = CHAR(STRING_ELT(RgenomeSeq,0));
	int strand = INTEGER(Rstrand)[0];
	int left_len = 	INTEGER(Rleft_len)[0];
	int right_len = INTEGER(Rright_len)[0];
	string method = CHAR(STRING_ELT(Rmethod,0));
	int is_cvg_all = INTEGER(Ris_cvg_all)[0];
	Rprintf("method : %s\n",method.c_str());

	/*-----------------start to detect modifications----------------------*/	
	vector<double> detection(Rf_length(Ripd_native), sqrt(-1));	
	vector<double> ipd_ratio(Rf_length(Ripd_native), sqrt(-1));
	vector<double> is_findContext(Rf_length(Ripd_native), sqrt(-1));	
	
	vector<double> mu0(Rf_length(Ripd_native), sqrt(-1));
	vector<double> kappa0(Rf_length(Ripd_native), sqrt(-1));
	vector<double> upsilon0(Rf_length(Ripd_native), sqrt(-1));
	vector<double> sigma0(Rf_length(Ripd_native), sqrt(-1));

	vector<double> mu_native_est(Rf_length(Ripd_native), sqrt(-1));
	vector<double> sigma2_native_est(Rf_length(Ripd_native), sqrt(-1));
	vector<double> sampleSize_native(Rf_length(Ripd_native), sqrt(-1));

	vector<double> mu_ctrl_est(Rf_length(Ripd_native), sqrt(-1));
        vector<double> sigma2_ctrl_est(Rf_length(Ripd_native), sqrt(-1));
        vector<double> sampleSize_ctrl(Rf_length(Ripd_native), sqrt(-1));

	vector<double> effect_size(Rf_length(Ripd_native), sqrt(-1));

	vector<double> L_log_null(Rf_length(Ripd_native), sqrt(-1));
        vector<double> L_log_alt(Rf_length(Ripd_native), sqrt(-1));
        vector<double> LR_log(Rf_length(Ripd_native), sqrt(-1));

	vector<double> sampleSize_history(Rf_length(Ripd_native), sqrt(-1));
	
	for (int i=0;i<Rf_length(Ripd_native);i++){
		if ((i+1)%10000==0) Rprintf("detected %d positions\r", i+1);
		// get ipd for test
		double *ipd_native_tmp = REAL(VECTOR_ELT(Ripd_native,i));
		vector<double> ipd_native(ipd_native_tmp, ipd_native_tmp + Rf_length(VECTOR_ELT(Ripd_native,i)) );
		int cur_pos = i + genome_start_native - 1;
		string cur_context;
		if (strand == 0){
			if (cur_pos - left_len<0 || cur_pos + right_len>Rf_length(Ripd_native)-1 || cur_pos + right_len > (int) genomeSeq.size()-1)
				continue;
			cur_context = genomeSeq.substr(cur_pos - left_len , right_len + left_len + 1);		
		}else{
			if (cur_pos - right_len<0 || cur_pos + left_len>Rf_length(Ripd_native)-1 || cur_pos + left_len > (int) genomeSeq.size()-1)
                                continue;
			cur_context = SP_reverseSeq (genomeSeq.substr(cur_pos - right_len, right_len + left_len + 1) ) ;
		}
		
		if (ipd_native.size()<3){
                        sampleSize_native[i] = ipd_native.size();
                        continue;
                }

		map<string, map<string, vector<double> > >::iterator it = contextEffect.find(cur_context);
               	if (it!=contextEffect.end())
			is_findContext[i] = 1;
		else
			is_findContext[i] = 0;
		
		double mu_native = sqrt(-1);
		double sigma_native = sqrt(-1);
		double n_native = sqrt(-1);

		double mu_ctrl = sqrt(-1);
                double sigma_ctrl = sqrt(-1);
                double n_ctrl = sqrt(-1);

		if (method=="hieModel"){
			if (it!=contextEffect.end()){
				//vector<vector<double> > ipd_ref = it->second;	
				vector<double> ipd_ref_mean = it->second["mean"];
                                vector<double> ipd_ref_var = it->second["var"];
                                vector<double> ipd_ref_len = it->second["len"];

				// fit hierarchical model
				map<string, vector<double> > hie_fit = hieModelEB_alt(ipd_native, ipd_ref_mean, ipd_ref_var, ipd_ref_len);
				mu_native = mean(ipd_native);
				sigma_native = var(ipd_native);
				n_native = ipd_native.size();

				double theta = hie_fit["hyperPara"][0];
				double kappa = hie_fit["hyperPara"][1];
				double upsilon = hie_fit["hyperPara"][2];
				double tau2 = hie_fit["hyperPara"][3];	
				
				mu_ctrl = theta;
				sigma_ctrl = tau2*upsilon/(upsilon-2);
				n_ctrl = kappa;
				detection[i] = get_t_stat(mu_native, sigma_native, n_native, mu_ctrl, sigma_ctrl, n_ctrl);
				
				mu0[i] = hie_fit["hyperPara"][0];
				kappa0[i] = hie_fit["hyperPara"][1];
				upsilon0[i] = hie_fit["hyperPara"][2];
				sigma0[i] = hie_fit["hyperPara"][3];
				
				// liklihood ratio test
				int k = ipd_ref_mean.size();
                                double cur_L_log_alt = getLogLikelihood_marginal(ipd_ref_mean, ipd_ref_var, ipd_ref_len,
							hie_fit["hyperPara"][0],hie_fit["hyperPara"][1]
                                                        ,hie_fit["hyperPara"][2], hie_fit["hyperPara"][3] );
				cur_L_log_alt += lgamma((upsilon0[i] + n_native)/2) + upsilon0[i]*log(upsilon0[i]*sigma0[i])/2 -
                                        lgamma(upsilon0[i]/2) - (upsilon0[i]+n_native)*log( (n_native-1)*sigma_native + upsilon0[i]*sigma0[i] )/2
                                        -n_native*log(my_pi)/2;

                                ipd_ref_mean.push_back(mean(ipd_native));
				ipd_ref_var.push_back(var(ipd_native));
				ipd_ref_len.push_back(ipd_native.size());
			
                                map<string, vector<double> > hie_fit_null = hieModelEB(ipd_ref_mean, ipd_ref_var, ipd_ref_len);
                                double cur_L_log_null = getLogLikelihood_marginal(ipd_ref_mean, ipd_ref_var, ipd_ref_len,
							hie_fit_null["hyperPara"][0],hie_fit_null["hyperPara"][1]
                                                        ,hie_fit_null["hyperPara"][2], hie_fit_null["hyperPara"][3] );
                                /*cur_L_log_alt += sum(ldnorm(ipd_native, mean(ipd_native), hie_fit_null["sigma"][k]));
                              	double sigma_native_null = (upsilon0[i]*sigma0[i] + (n_native-1)*sigma_native +
                                                        kappa0[i]*n_native*(mu_native-mu0[i])*(mu_native-mu0[i])/(kappa0[i]+n_native))
                                                         / (upsilon0[i] + n_native - 2);*/
				/*double sigma_native_null = (upsilon0[i]*sigma0[i] + (n_native-1)*sigma_native)/ (upsilon0[i] + n_native - 2);
				cur_L_log_alt += sum(ldnorm(ipd_native, mu_native, sigma_native_null));*/
				
				L_log_null[i] = cur_L_log_null;
                                L_log_alt[i] = cur_L_log_alt;
                                LR_log[i] = cur_L_log_alt - cur_L_log_null;
				
				/*double mu_native_null = (kappa0[i]*mu0[i] + n_native*mu_native)/(kappa0[i] + n_native);
	                        double sigma_native_null = (upsilon0[i]*sigma0[i] + (n_native-1)*sigma_native +
                        	                        kappa0[i]*n_native*(mu_native-mu0[i])*(mu_native-mu0[i])/(kappa0[i]+n_native))
                	                                 / (upsilon0[i] + n_native - 2);
        	                double cur_L_log_alt = getLogLikelihood_marginal(ipd_ref,hie_fit["hyperPara"][0],hie_fit["hyperPara"][1]
                                	                ,hie_fit["hyperPara"][2], hie_fit["hyperPara"][3] )
                                        	        + sum(ldnorm(ipd_native, mu_native, sigma_native_null));
                        	ipd_ref.push_back(ipd_native);
				double cur_L_log_null = getLogLikelihood_marginal(ipd_ref,hie_fit["hyperPara"][0],hie_fit["hyperPara"][1]
                                	                ,hie_fit["hyperPara"][2], hie_fit["hyperPara"][3] );
                        	L_log_null[i] = cur_L_log_null;
                        	L_log_alt[i] = cur_L_log_alt;
                        	LR_log[i] = cur_L_log_alt - cur_L_log_null;*/
			}
		}
		
		mu_native_est[i] = mu_native;
		sigma2_native_est[i] = sigma_native;
		sampleSize_native[i] = n_native;
		
		mu_ctrl_est[i] = mu_ctrl;
                sigma2_ctrl_est[i] = sigma_ctrl;
                sampleSize_ctrl[i] = n_ctrl;
	
		effect_size[i] = (mu_native-mu_ctrl)/sqrt(sigma_native+sigma_ctrl);
	}
	Rprintf("detected %d positions\n", Rf_length(Ripd_native));
	// load results
	map<string, vector<double> > result;	
	result["t_stat"] = detection;
	result["ipd_ratio"] = ipd_ratio; 
	result["is_findContext"] = is_findContext;
	
	result["theta"] = mu0;
	result["kappa"] = kappa0;
	result["upsilon"] = upsilon0;
	result["tau2"] = sigma0;

	result["mu_native"] = mu_native_est;
	result["sigma2_native"] = sigma2_native_est;
	result["sampleSize_native"] = sampleSize_native;

	result["mu_ctrl"] = mu_ctrl_est;
        result["sigma2_ctrl"] = sigma2_ctrl_est;
        result["sampleSize_ctrl"] = sampleSize_ctrl;

	result["effect.size"] = effect_size;

	result["L_log_null"] = L_log_null;
	result["L_log_alt"] = L_log_alt;
	result["LR_log"] = LR_log;	
	SEXP rl = Rcpp::wrap(result);
	return rl;

}


RcppExport SEXP hieModel_CC (SEXP Ripd_native, SEXP Ripd_ctrl, SEXP Rgenome_start_native, SEXP Rgenome_start_ctrl,
			 SEXP RgenomeSeq, SEXP RcontextEffect, SEXP Rcontext, SEXP Rstrand, SEXP Rleft_len, SEXP Rright_len, SEXP Rmethod)
{

	// load context effect
	map<string, map<string, vector<double> > > contextEffect;
        if (!Rf_isNull(RcontextEffect)){
                for (int i=0;i<Rf_length(RcontextEffect);i++){
                        string cur_context = CHAR(STRING_ELT(Rcontext,i));
                        map<string, vector<double > > cur_data;
                        vector<double> tmp;
                        cur_data["mean"] = tmp; cur_data["var"] = tmp; cur_data["len"] = tmp;
                        SEXP cur_Rdata = VECTOR_ELT(RcontextEffect,i);
                        double *dat;
                        // get mean
                        dat = REAL(VECTOR_ELT(cur_Rdata, getIdxByName(cur_Rdata, "mean") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,0)); j++){cur_data["mean"].push_back(dat[j]);}
                        // get var
                        dat = REAL(VECTOR_ELT(cur_Rdata, getIdxByName(cur_Rdata, "var") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,1)); j++){cur_data["var"].push_back(dat[j]);}
                        // get len
                        dat = REAL(VECTOR_ELT(cur_Rdata, getIdxByName(cur_Rdata, "len") ));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,2)); j++){cur_data["len"].push_back(dat[j]);}
                        contextEffect[cur_context] = cur_data;
                }
        }

	/*map<string, vector<vector<double> > > contextEffect;
	if (!Rf_isNull(RcontextEffect)){
		for (int i=0;i<Rf_length(RcontextEffect);i++){
			string cur_context = CHAR(STRING_ELT(Rcontext,i));
			vector<vector<double> > ipd;
			for (int j=0;j<Rf_length(VECTOR_ELT(RcontextEffect,i)); j++) {vector<double> tmp; ipd.push_back(tmp);}
			for (int j=0;j<Rf_length(VECTOR_ELT(RcontextEffect,i)); j++){
				double * cur_ipd = REAL(VECTOR_ELT(VECTOR_ELT(RcontextEffect,i),j));
				int cur_cvg = Rf_length(VECTOR_ELT(VECTOR_ELT(RcontextEffect,i),j));
				for (int k=0; k<cur_cvg; k++ ){
					ipd[j].push_back(cur_ipd[k]);
				}
			}
			contextEffect[cur_context] = ipd;
		}
	}*/

	// load data and parameters
	int genome_start_native = INTEGER(Rgenome_start_native)[0];
	int genome_start_ctrl = INTEGER(Rgenome_start_ctrl)[0];
	string genomeSeq = CHAR(STRING_ELT(RgenomeSeq,0));
	int strand = INTEGER(Rstrand)[0];
	int left_len = 	INTEGER(Rleft_len)[0];
	int right_len = INTEGER(Rright_len)[0];
	string method = CHAR(STRING_ELT(Rmethod,0));
	Rprintf("method : %s\n",method.c_str());
		
	/*-----------------start to detect modifications----------------------*/	
	vector<double> detection(Rf_length(Ripd_native), sqrt(-1));	
	vector<double> ipd_ratio(Rf_length(Ripd_native), sqrt(-1));
	vector<double> is_findContext(Rf_length(Ripd_native), sqrt(-1));	
	
	vector<double> mu0(Rf_length(Ripd_native), sqrt(-1));
	vector<double> kappa0(Rf_length(Ripd_native), sqrt(-1));
	vector<double> upsilon0(Rf_length(Ripd_native), sqrt(-1));
	vector<double> sigma0(Rf_length(Ripd_native), sqrt(-1));

	vector<double> mu_native_est(Rf_length(Ripd_native), sqrt(-1));
	vector<double> sigma2_native_est(Rf_length(Ripd_native), sqrt(-1));
	vector<double> sampleSize_native(Rf_length(Ripd_native), sqrt(-1));

	vector<double> mu_ctrl_est(Rf_length(Ripd_native), sqrt(-1));
        vector<double> sigma2_ctrl_est(Rf_length(Ripd_native), sqrt(-1));
        vector<double> sampleSize_ctrl(Rf_length(Ripd_native), sqrt(-1));

	vector<double> effect_size(Rf_length(Ripd_native), sqrt(-1));

        vector<double> L_log_alt(Rf_length(Ripd_native), sqrt(-1));
	vector<double> L_log_null(Rf_length(Ripd_native), sqrt(-1));
	vector<double> LR_log (Rf_length(Ripd_native), sqrt(-1));

	vector<double> sampleSize_history(Rf_length(Ripd_native), sqrt(-1));

	for (int i=0;i<Rf_length(Ripd_native);i++){
		if ((i+1)%10000==0) Rprintf("detected %d positions\r", i+1);
		// get ipd for test
		if (i + genome_start_native - genome_start_ctrl<0 || i + genome_start_native - genome_start_ctrl > Rf_length(Ripd_ctrl)-1)
			continue;
		double *ipd_native_tmp = REAL(VECTOR_ELT(Ripd_native,i));
		double *ipd_ctrl_tmp = REAL(VECTOR_ELT(Ripd_ctrl, i + genome_start_native - genome_start_ctrl ));
		vector<double> ipd_native(ipd_native_tmp, ipd_native_tmp + Rf_length(VECTOR_ELT(Ripd_native,i)) );
		vector<double> ipd_ctrl(ipd_ctrl_tmp, ipd_ctrl_tmp + Rf_length(VECTOR_ELT(Ripd_ctrl, i + genome_start_native - genome_start_ctrl)) );
		
		int cur_pos = i + genome_start_native - 1;
		string cur_context;
		if (strand == 0){
			if (cur_pos - left_len<0 || cur_pos + right_len>Rf_length(Ripd_native)-1 || cur_pos + right_len > (int) genomeSeq.size()-1)
				continue;
			cur_context = genomeSeq.substr(cur_pos - left_len , right_len + left_len + 1);		
		}else{
			if (cur_pos - right_len<0 || cur_pos + left_len>Rf_length(Ripd_native)-1 || cur_pos + left_len > (int) genomeSeq.size()-1)
                                continue;
			cur_context = SP_reverseSeq (genomeSeq.substr(cur_pos - right_len, right_len + left_len + 1) ) ;
		}
		
		if (ipd_native.size()<3 || ipd_ctrl.size() <3){
                        sampleSize_native[i] = ipd_native.size();
                        sampleSize_ctrl[i] = ipd_ctrl.size();
                        continue;
                }
		map<string, map<string, vector<double> > >::iterator it = contextEffect.find(cur_context);
               	if (it!=contextEffect.end())
			is_findContext[i] = 1;
		else
			is_findContext[i] = 0;
		
		double mu_native = sqrt(-1);
		double sigma_native = sqrt(-1);
		double n_native = sqrt(-1);

		double mu_ctrl = sqrt(-1);
                double sigma_ctrl = sqrt(-1);
                double n_ctrl = sqrt(-1);

		if (method=="hieModel"){
			if (it!=contextEffect.end()){
				vector<double> ipd_ref_mean = it->second["mean"];
				vector<double> ipd_ref_var = it->second["var"];
				vector<double> ipd_ref_len = it->second["len"];
				ipd_ref_mean.push_back(mean(ipd_ctrl));
				ipd_ref_var.push_back(var(ipd_ctrl));
				ipd_ref_len.push_back(ipd_ctrl.size());
				//ipd_ref.push_back(ipd_ctrl);
					
				// fit hierarchical model
				//map<string, vector<double> > hie_fit = hieModelEB_alt(ipd_ref,ipd_native);
				map<string, vector<double> > hie_fit = hieModelEB_alt(ipd_native, ipd_ref_mean, ipd_ref_var, ipd_ref_len);

				// get t-statistic and ipd ratio
				mu_native = mean(ipd_native);
				sigma_native = var(ipd_native);
				n_native = ipd_native.size();
				int k = (int) hie_fit["mu"].size()-1;
				mu_ctrl = hie_fit["mu"][k];
				sigma_ctrl = hie_fit["sigma"][k] * hie_fit["upsilon"][k]/(hie_fit["upsilon"][k]-2);
				n_ctrl = ipd_ctrl.size();
				detection[i] = get_t_stat(mu_native, sigma_native, n_native, mu_ctrl, sigma_ctrl, n_ctrl+hie_fit["hyperPara"][1]);
				ipd_ratio[i] = exp(mu_native - mu_ctrl);
				// get hyper parameters
				mu0[i] = hie_fit["hyperPara"][0];
				kappa0[i] = hie_fit["hyperPara"][1];
				upsilon0[i] = hie_fit["hyperPara"][2];
				sigma0[i] = hie_fit["hyperPara"][3];
				
				// liklihood ratio test
				double cur_L_log_alt = getLogLikelihood_marginal(ipd_ref_mean, ipd_ref_var, ipd_ref_len,
							hie_fit["hyperPara"][0],hie_fit["hyperPara"][1]
                                                        ,hie_fit["hyperPara"][2], hie_fit["hyperPara"][3] );
				cur_L_log_alt += lgamma((upsilon0[i] + n_native)/2) + upsilon0[i]*log(upsilon0[i]*sigma0[i])/2 -
                                        lgamma(upsilon0[i]/2) - (upsilon0[i]+n_native)*log( (n_native-1)*sigma_native + upsilon0[i]*sigma0[i] )/2
                                        -n_native*log(my_pi)/2;

				//ipd_ref.pop_back();
				ipd_ref_mean.pop_back();
				ipd_ref_var.pop_back();
				ipd_ref_len.pop_back();

				vector<double> ipd_merge (ipd_ctrl);
				for (int j=0;j<(int)ipd_native.size();j++) ipd_merge.push_back(ipd_native[j]);
				
				ipd_ref_mean.push_back(mean(ipd_merge));
				ipd_ref_var.push_back(var(ipd_merge));
				ipd_ref_len.push_back(ipd_merge.size());
				
				map<string, vector<double> > hie_fit_null = hieModelEB(ipd_ref_mean, ipd_ref_var, ipd_ref_len);
				double cur_L_log_null = getLogLikelihood_marginal(ipd_ref_mean, ipd_ref_var, ipd_ref_len,
							hie_fit_null["hyperPara"][0],hie_fit_null["hyperPara"][1]
							,hie_fit_null["hyperPara"][2], hie_fit_null["hyperPara"][3] );
				/*cur_L_log_alt += sum(ldnorm(ipd_native, mean(ipd_native), hie_fit_null["sigma"][k]));
				double sigma_native_null = (upsilon0[i]*sigma0[i] + (n_native-1)*sigma_native +
                                                        kappa0[i]*n_native*(mu_native-mu0[i])*(mu_native-mu0[i])/(kappa0[i]+n_native))
                                                         / (upsilon0[i] + n_native - 2);*/
				/*double sigma_native_null = (upsilon0[i]*sigma0[i] + (n_native-1)*sigma_native)/ (upsilon0[i] + n_native - 2);
				cur_L_log_alt += sum(ldnorm(ipd_native, mu_native, sigma_native_null));*/
				L_log_null[i] = cur_L_log_null;
				L_log_alt[i] = cur_L_log_alt;
				LR_log[i] = cur_L_log_alt - cur_L_log_null;

				/*double mu_native_null = (kappa0[i]*mu0[i] + n_native*mu_native)/(kappa0[i] + n_native);
	                        double sigma_native_null = (upsilon0[i]*sigma0[i] + (n_native-1)*sigma_native +
        	                                        kappa0[i]*n_native*(mu_native-mu0[i])*(mu_native-mu0[i])/(kappa0[i]+n_native))
                	                                 / (upsilon0[i] + n_native - 2);
                        	double cur_L_log_alt = getLogLikelihood_marginal(ipd_ref,hie_fit["hyperPara"][0],hie_fit["hyperPara"][1]
                                	                ,hie_fit["hyperPara"][2], hie_fit["hyperPara"][3] )
                                        	        + sum(ldnorm(ipd_native, mu_native, sigma_native_null));
                        	for (int j=0;j<(int)ipd_native.size();j++) ipd_ref[k].push_back(ipd_native[j]);
                        	double cur_L_log_null = getLogLikelihood_marginal(ipd_ref,hie_fit["hyperPara"][0],hie_fit["hyperPara"][1]
                                	                ,hie_fit["hyperPara"][2], hie_fit["hyperPara"][3] );
                        	L_log_null[i] = cur_L_log_null;
                        	L_log_alt[i] = cur_L_log_alt;
                        	LR_log[i] = cur_L_log_alt - cur_L_log_null;*/
			
				/*double sigma_native_null = (upsilon0[i]*sigma0[i] + (n_native-1)*sigma_native +
                                                        kappa0[i]*n_native*(mu_native-mu0[i])*(mu_native-mu0[i])/(kappa0[i]+n_native))
                                                         / (upsilon0[i] + n_native - 2);
				LR_log[i] = get_t_stat(mu_native, sigma_native_null, n_native, mu_ctrl, sigma_ctrl, n_ctrl+hie_fit["hyperPara"][1]);*/

			}else{
				mu_native = mean(ipd_native);
                	        sigma_native = var(ipd_native);
                	        n_native = ipd_native.size();
			
				mu_ctrl = mean(ipd_ctrl);
                	        sigma_ctrl = var(ipd_ctrl);
                	        n_ctrl = ipd_ctrl.size();
				
				detection[i] = get_t_stat(mu_native, sigma_native, n_native, mu_ctrl, sigma_ctrl, n_ctrl);
				ipd_ratio[i] = exp(mu_native - mu_ctrl);
			
				// get likelyhood 
				vector<double> ipd_merge (ipd_ctrl);
                                for (int j=0;j<(int)ipd_native.size();j++) ipd_merge.push_back(ipd_native[j]);
				double sigma_merge = var(ipd_merge);
				L_log_null[i] = sum(ldnorm(ipd_merge, mean(ipd_merge), sigma_merge));	
				L_log_alt[i] = sum(ldnorm(ipd_native, mu_native, sigma_merge)) + sum(ldnorm(ipd_ctrl, mu_ctrl, sigma_merge));
				LR_log[i] = L_log_alt[i] - L_log_null[i];
			}
		}else{
			if (method=="CC"){
				mu_native = mean(ipd_native);
                                sigma_native = var(ipd_native);
                                n_native = ipd_native.size();

                                mu_ctrl = mean(ipd_ctrl);
                                sigma_ctrl = var(ipd_ctrl);
                                n_ctrl = ipd_ctrl.size();
				detection[i] = get_t_stat(mu_native, sigma_native, n_native, mu_ctrl, sigma_ctrl, n_ctrl);
                                //detection[i] = - fabs(get_t_stat(mu_native, sigma_native, n_native, mu_ctrl, sigma_ctrl, n_ctrl));
                                //detection[i] = - fabs(get_t_stat(mu_native, sigma_native, 1, mu_ctrl, sigma_ctrl, 1));
				ipd_ratio[i] = exp(mu_native - mu_ctrl);
			}
		}
		mu_native_est[i] = mu_native;
		sigma2_native_est[i] = sigma_native;
		sampleSize_native[i] = n_native;
		
		mu_ctrl_est[i] = mu_ctrl;
                sigma2_ctrl_est[i] = sigma_ctrl;
                sampleSize_ctrl[i] = n_ctrl;

		effect_size[i] = (mu_native-mu_ctrl)/sqrt(sigma_native+sigma_ctrl);	
	}
	Rprintf("detected %d positions\n", Rf_length(Ripd_native));
	// load results
	map<string, vector<double> > result;	
	result["t_stat"] = detection;
	result["ipd_ratio"] = ipd_ratio; 
	result["is_findContext"] = is_findContext;
	
	result["theta"] = mu0;
	result["kappa"] = kappa0;
	result["upsilon"] = upsilon0;
	result["tau2"] = sigma0;

	result["mu_native"] = mu_native_est;
	result["sigma2_native"] = sigma2_native_est;
	result["sampleSize_native"] = sampleSize_native;

	result["mu_ctrl"] = mu_ctrl_est;
        result["sigma2_ctrl"] = sigma2_ctrl_est;
        result["sampleSize_ctrl"] = sampleSize_ctrl;

	result["effect.size"] = effect_size;
		
	result["L_log_null"] = L_log_null;
	result["L_log_alt"] = L_log_alt;
	result["LR_log"] = LR_log;

        SEXP rl = Rcpp::wrap(result);
	return rl;
}



RcppExport SEXP collapseContextEffectByPos(SEXP RcontextEffectByPos, SEXP RcvgCutoff)
{
	SEXP context = Rf_getAttrib(RcontextEffectByPos, R_NamesSymbol);
	int cvgCutoff = INTEGER(RcvgCutoff)[0];
	int nContext = Rf_length(RcontextEffectByPos);
	vector <vector<double> > result;
	for (int i=0; i<nContext; i++){
		if ((i+1)%100000==0) Rprintf("collapsed %d contexts\r",i+1);
		int nSite = Rf_length(VECTOR_ELT(RcontextEffectByPos,i));
		string cur_context = CHAR(STRING_ELT(context,i));
		for (int j=0;j<nSite;j++){
			vector<double> cur_result;
			int cur_cvg = Rf_length(VECTOR_ELT(VECTOR_ELT(RcontextEffectByPos,i),j));
			if (cur_cvg < cvgCutoff) continue;
			// get log mean ipd for each position that pass filtering
			double *ipd = REAL(VECTOR_ELT(VECTOR_ELT(RcontextEffectByPos,i),j));
			vector<double> ipd_vec;  for (int k=0;k<cur_cvg;k++) ipd_vec.push_back(log(ipd[k]+0.01));	
			ipd_vec = filter_outlier(ipd_vec);
			double ipd_mean_log = mean(ipd_vec); 
			//double ipd_mean_log = 0;
			//for (int k=0;k<cur_cvg;k++) ipd_mean_log += log(ipd[k]+0.01);
			//ipd_mean_log = ipd_mean_log/cur_cvg;

			// record ipd and context
			cur_result.push_back(ipd_mean_log);	
			for (int k=0;k<cur_context.size();k++){ 
				switch(cur_context[k]){
					case 'A': cur_result.push_back(1); break;
					case 'C': cur_result.push_back(2); break;
					case 'G': cur_result.push_back(3); break;	
					case 'T': cur_result.push_back(4); break;
					default : cur_result.push_back(0);
				}	
			}
			result.push_back(cur_result);
		}
			
	}
	Rprintf("collapsed %d contexts\n",nContext);

	// output result as matrix : coding is A:1 C:2 G:3 T:4
	int nrows = (int) result.size();
	int ncols = (int) result[0].size();
	SEXP Rrl;
        double * rl;
        PROTECT(Rrl = Rf_allocMatrix(REALSXP,nrows,ncols));
        rl = REAL(Rrl);
        for (int i=0;i<nrows;i++)
                for (int j=0;j<ncols;j++)
                        rl[i+nrows*j] = result[i][j];
        UNPROTECT(1);
	return Rcpp::wrap(Rrl);
	//return Rcpp::wrap(result);
}
RcppExport SEXP permutation_genomeF(SEXP Ripd_native, SEXP Ripd_ctrl, SEXP Rstart_native, SEXP Rstart_ctrl)
{
        int n_native = Rf_length(Ripd_native);
        int n_ctrl = Rf_length(Ripd_ctrl);
        int start_native = INTEGER(Rstart_native)[0];
        int start_ctrl = INTEGER(Rstart_ctrl)[0];

        for (int i=0;i<n_native;i++){
                //Rprintf("permuated %d positions \r",i+1);
                if ((i+1)%10000==0) Rprintf("permuated %d positions \r",i+1);
                int idx_ctrl = i + start_native - start_ctrl;
                if (idx_ctrl < 0 || idx_ctrl>=n_ctrl) continue;
                double *cur_ipd_native = REAL(VECTOR_ELT(Ripd_native,i));
                double *cur_ipd_ctrl = REAL(VECTOR_ELT(Ripd_ctrl,idx_ctrl));
                int cur_cvg_native = Rf_length(VECTOR_ELT(Ripd_native,i));
                int cur_cvg_ctrl = Rf_length(VECTOR_ELT(Ripd_ctrl,idx_ctrl));
                vector<double> cur_ipd_pool;
                for(int j=0;j<cur_cvg_native;j++) cur_ipd_pool.push_back(cur_ipd_native[j]);
                for(int j=0;j<cur_cvg_ctrl;j++) cur_ipd_pool.push_back(cur_ipd_ctrl[j]);
                vector<double> cur_ipd_pool_B = sample_SELECT(cur_ipd_pool);

                for(int j=0;j<cur_cvg_native;j++) cur_ipd_native[j] = cur_ipd_pool_B[j];
                for(int j=cur_cvg_native;j<(int)cur_ipd_pool.size();j++) cur_ipd_ctrl[j-cur_cvg_native]= cur_ipd_pool_B[j];
        }
        Rprintf("permuated %d positions \n",n_native);
	
	SEXP rl=R_NilValue;
        return rl;
}
RcppExport SEXP mixContextEffect(SEXP RcontextEffect_1, SEXP RcontextEffect_2)
{
	map<string, vector<vector<double> > > contextEffect_1;
	SEXP RcontextEffect =RcontextEffect_1 ;
	SEXP Rcontext = Rf_getAttrib(RcontextEffect_1, R_NamesSymbol);	
	Rprintf("load context effect 1\n");
	for (int i=0;i<Rf_length(RcontextEffect);i++){
		string cur_context = CHAR(STRING_ELT(Rcontext,i));
		vector<vector<double> > ipd;
		for (int j=0;j<Rf_length(VECTOR_ELT(RcontextEffect,i)); j++) {vector<double> tmp; ipd.push_back(tmp);}
		for (int j=0;j<Rf_length(VECTOR_ELT(RcontextEffect,i)); j++){
			double * cur_ipd = REAL(VECTOR_ELT(VECTOR_ELT(RcontextEffect,i),j));
			int cur_cvg = Rf_length(VECTOR_ELT(VECTOR_ELT(RcontextEffect,i),j));
			for (int k=0; k<cur_cvg; k++ ){
				ipd[j].push_back(cur_ipd[k]);
			}
		}
		contextEffect_1[cur_context] = ipd;
	}
	Rprintf("load context effect 2\n");
	map<string, vector<vector<double> > > contextEffect_2;
        RcontextEffect =RcontextEffect_2 ;
        Rcontext = Rf_getAttrib(RcontextEffect_2, R_NamesSymbol);
        for (int i=0;i<Rf_length(RcontextEffect);i++){
                string cur_context = CHAR(STRING_ELT(Rcontext,i));
                vector<vector<double> > ipd;
                for (int j=0;j<Rf_length(VECTOR_ELT(RcontextEffect,i)); j++) {vector<double> tmp; ipd.push_back(tmp);}
                for (int j=0;j<Rf_length(VECTOR_ELT(RcontextEffect,i)); j++){
                        double * cur_ipd = REAL(VECTOR_ELT(VECTOR_ELT(RcontextEffect,i),j));
                        int cur_cvg = Rf_length(VECTOR_ELT(VECTOR_ELT(RcontextEffect,i),j));
                        for (int k=0; k<cur_cvg; k++ ){
                                ipd[j].push_back(cur_ipd[k]);
                        }
                }
                contextEffect_2[cur_context] = ipd;
        }
	// merge contextEffect_2
	Rprintf("merge context effect\n");
	map<string, vector<vector<double> > >::iterator it;
	int n_context=0;
	for (it=contextEffect_2.begin();it!=contextEffect_2.end();it++){
		n_context++;
		if (n_context%1000==0)Rprintf("combined %d contexts\r", n_context);
		string cur_context = it->first;
		for (int i=0;i<(int)it->second.size();i++){
			contextEffect_1[cur_context].push_back(it->second[i]);
		}	
	}
	Rprintf("combined %d contexts\n", n_context);
	return Rcpp::wrap(contextEffect_1);					
}

// bootstrap to get motif frequency in background
RcppExport SEXP bootstrap_bg_motif(SEXP Rmotif, SEXP RgenomeSeq, SEXP Ridx_sel, SEXP Rstrand, SEXP Rtol, SEXP Rmotif_len, SEXP RB, SEXP Rn_context)
{
	/*--------parse arguments----------*/
	const string genomeSeq = CHAR(STRING_ELT(RgenomeSeq,0)); 
	const vector<int> idx_sel_all( INTEGER(Ridx_sel), INTEGER(Ridx_sel)+ Rf_length(Ridx_sel) );
	int strand = INTEGER(Rstrand)[0];	// 0: forward , 1: backward
	int tol = INTEGER(Rtol)[0];
	int motif_len = INTEGER(Rmotif_len)[0];
	int B = INTEGER(RB)[0];
	int n_context = INTEGER(Rn_context)[0];		

	// remove marginal positions 
	vector<int> idx_sel;
	for (int i=0;i<(int)idx_sel_all.size();i++){ if(idx_sel_all[i]>=tol && idx_sel_all[i]<genomeSeq.size()-tol) idx_sel.push_back(idx_sel_all[i]); } 
			
	// bootstrap	
	//vector<map<string, int> > motif_freq;
	map<string, vector<double> > motif_freq;
	map<string, vector<double> > motif_count;
	vector<int> cur_idx;
	vector<string> context_B;
	for (int i=0;i<B;i++){
		map<string, int> cur_motif_freq;	
		map<string, int>::iterator it;
		cur_idx = RNG_sample(idx_sel, n_context);
		for (int j=0;j<(int)cur_idx.size();j++){
			string cur_context = genomeSeq.substr(cur_idx[j]-tol, 2*tol+1);
			if(i==0) context_B.push_back(cur_context);
			map<string, bool> is_find;
			for (int k=0;k<(int)cur_context.size()-motif_len+1;k++){
				string the_motif = cur_context.substr(k, motif_len);
				if (is_find.find(the_motif)==is_find.end()) is_find[the_motif] = true;
				else continue;	
				it = cur_motif_freq.find(the_motif);
				if (it==cur_motif_freq.end()) cur_motif_freq[the_motif] = 1;	
				else cur_motif_freq[the_motif]++;	
			}				
		}
		
		for (it = cur_motif_freq.begin(); it!=cur_motif_freq.end(); it++){
			motif_freq[it->first].push_back(1.0*it->second/n_context);
			motif_count[it->first].push_back(it->second);
		}
		//motif_freq.push_back(cur_motif_freq);		
		if (i%100==0) Rprintf("%d\r",i);	
	}	
	Rprintf("%d\n",B);
	return Rcpp::List::create(Rcpp::Named("motif_freq")=Rcpp::wrap(motif_freq), Rcpp::Named("motif_count")=Rcpp::wrap(motif_count),
				Rcpp::Named("n_context")=Rcpp::wrap(n_context) );
	//return Rcpp::wrap(motif_freq);
	//return Rcpp::List::create(Rcpp::Named("motif_freq")=Rcpp::wrap(motif_freq), Rcpp::Named("context_B")=Rcpp::wrap(context_B) );
}










/*RcppExport SEXP hieModelEB_core(SEXP data, SEXP Rmax_iter)
{
	vector<vector<double> > ipd;
	for (int i=0;i<Rf_length(data);i++){
		SEXP Rcur_ipd = VECTOR_ELT(data,i);
		vector<double> cur_ipd(REAL(Rcur_ipd), REAL(Rcur_ipd) + Rf_length(Rcur_ipd));
		ipd.push_back(cur_ipd);
	}

	int max_iter = INTEGER(Rmax_iter)[0];
	return Rcpp::wrap(hieModelEB(ipd, max_iter));
}
RcppExport SEXP hieModelEB_core_alt(SEXP data, SEXP Ripd_native, SEXP Rmax_iter)
{
        vector<vector<double> > ipd;
        for (int i=0;i<Rf_length(data);i++){
                SEXP Rcur_ipd = VECTOR_ELT(data,i);
                vector<double> cur_ipd(REAL(Rcur_ipd), REAL(Rcur_ipd) + Rf_length(Rcur_ipd));
                ipd.push_back(cur_ipd);
        }
	
	vector<double> ipd_native(REAL(Ripd_native), REAL(Ripd_native) + Rf_length(Ripd_native));
        int max_iter = INTEGER(Rmax_iter)[0];
	//return Rcpp::wrap(ipd_native);        
	return Rcpp::wrap(hieModelEB_alt(ipd, ipd_native, max_iter));
}*/


RcppExport SEXP filter_outlier_wrap(SEXP data, SEXP Rk, SEXP Rtail)
{
	double * cur_ipd_ay = REAL(data);
	int cur_cvg = Rf_length(data);
	vector<double> cur_ipd(cur_ipd_ay, cur_ipd_ay + cur_cvg);

	double k = REAL(Rk)[0];
	
	string tail = CHAR(STRING_ELT(Rtail,0));
	if (cur_cvg<=3) return Rcpp::wrap(cur_ipd);
        else return Rcpp::wrap(filter_outlier(cur_ipd, k, tail));
	
}
RcppExport SEXP filter_outlier_by_genomeF(SEXP data)
{
	vector<vector<double> > ipd_ft;
	for (int i=0;i<Rf_length(data);i++){
		if ((i+1)%10000==0) Rprintf("filtered %d positions\r",i+1);
		double * cur_ipd_ay = REAL(VECTOR_ELT(data,i));
		int cur_cvg = Rf_length(VECTOR_ELT(data,i));
		vector<double> cur_ipd(cur_ipd_ay, cur_ipd_ay + cur_cvg);
		if (cur_cvg<=3) ipd_ft.push_back(cur_ipd);
		else ipd_ft.push_back(filter_outlier(cur_ipd));
	}
	Rprintf("filtered %d positions\n",Rf_length(data));
	return Rcpp::wrap(ipd_ft);
}



RcppExport SEXP reverseSeq(SEXP Rseq)
{
	const string seq = CHAR(STRING_ELT(Rseq,0));
	string new_seq(seq);
	int lock = 0;	
	for (unsigned int i=0;i<seq.size();i++){
		new_seq[i] = seq[seq.size()-1-i];
		if (new_seq[i]=='[' || new_seq[i]==']'){
			if (lock == 0){
				new_seq[i] = '[';
				lock = 1;
			}else{
				new_seq[i] = ']';
				lock = 0;
			}
		}
	}
	SEXP rl = Rcpp::wrap(new_seq);
	//SEXP rl=R_NilValue;
	return rl;  

}

RcppExport SEXP mergeGenomeF(SEXP Rdata_list, SEXP Rgenome_start_list, SEXP Rgenome_start, SEXP Rgenome_end)
{
	int genome_start = INTEGER(Rgenome_start)[0];
	int genome_end = INTEGER(Rgenome_end)[0];
	int data_size = genome_end - genome_start + 1;

	int *genome_start_list = INTEGER(Rgenome_start_list);
	int len = Rf_length(Rdata_list);

	// initialization
	vector<vector<double> > data_merged;	
	
	for (int i=0;i<data_size;i++){
		vector<double> data_elt;
		data_merged.push_back(data_elt);	
	}
	//return Rcpp::wrap(data_merged);	
	
	// write Rdata_list into data_merged
	for (int i=0;i<len;i++){
		int idx_shift = genome_start_list[i] - genome_start;
		SEXP cur_Rdata_list = VECTOR_ELT(Rdata_list,i);	
		int cur_len = Rf_length(cur_Rdata_list);	
		for (int j=0;j<cur_len;j++){
			double *cur_vec = REAL(VECTOR_ELT(cur_Rdata_list,j));
			int cur_size = Rf_length(VECTOR_ELT(cur_Rdata_list,j));
			for (int k=0;k<cur_size;k++){
				data_merged[j+idx_shift].push_back(cur_vec[k]);
			}
		}
	}		
	return Rcpp::wrap(data_merged);
	return R_NilValue;
}



RcppExport SEXP testFindMaxContextIndex(SEXP Rcur_context, SEXP Rex_ref, SEXP Ripd_ref)
{
	string cur_context = CHAR(STRING_ELT(Rcur_context,0));
	vector<string> ex_ref;	
	for (int i=0;i<Rf_length(Rex_ref);i++){
		ex_ref.push_back(CHAR(STRING_ELT(Rex_ref,i)));
	}
	vector<vector<double> > ipd_ref;
	return Rcpp::wrap(findMaxContextIndex(cur_context,ex_ref, ipd_ref));	
	return Rcpp::wrap(ex_ref);
}

RcppExport SEXP test_filter(SEXP data)
{
	double *rl = REAL(data);
	int size = Rf_length(data);
	vector<double> rl_vec; for (int i=0;i<size;i++) rl_vec.push_back(rl[i]);	
	rl_vec = filter_outlier(rl_vec, 1.5);
	
        return Rcpp::wrap(rl_vec);
}
RcppExport SEXP testNA(SEXP number)
{
	/*double x = REAL(number)[0];
	if (x!=x)
		Rprintf("na\n");*/
	
	double *tmp = REAL(number);
	
	tmp[0] = sqrt(-1);
	return Rcpp::wrap(tmp[0]); 
	//SEXP rl=R_NilValue;
        //return rl;
}

RcppExport SEXP testRNG(SEXP Rdata, SEXP Rsize)
{
	/*GetRNGstate();
	int N = 1000;
	vector<double> x;
	for (int i=0;i<N;i++){
		//x.push_back ( exp_rand() );
		//x.push_back ( norm_rand() );
		x.push_back ( unif_rand() );
	}
	PutRNGstate();*/
	//vector<double> x = RNG_gamma(10000,10,5);
	//vector<double> data = RNG_normal(20,1.5,0.9);	
	//int n = Rf_length(Rdata);
	int size = INTEGER(Rsize)[0];
	vector<double> data(REAL(Rdata), REAL(Rdata)+Rf_length(Rdata));
	vector<double> x = RNG_sample(data, size);
	return Rcpp::wrap(x);	
	/*vector<vector<double> > data;
	for (int i=0;i<Rf_length(Rdata);i++){
		vector<double> cur_data(REAL(VECTOR_ELT(Rdata,i)), REAL(VECTOR_ELT(Rdata,i)) + Rf_length(VECTOR_ELT(Rdata,i)));
		data.push_back(cur_data);
	}
	data = sample_group(data);	
	return Rcpp::wrap(data);*/
}
RcppExport SEXP test_rbinom(SEXP Rn)
{
	int n = INTEGER(Rn)[0];
	return Rcpp::wrap(RNG_rbinom(n, 100, 0.5));
}

RcppExport SEXP testList()
{
	map<string, vector<double> > result_1;
	vector <double> result_2(4,100);
	return Rcpp::List::create(Rcpp::Named("context")=Rcpp::wrap(result_1), Rcpp::Named("genomeF")=Rcpp::wrap(result_2)); 	
}
RcppExport SEXP testGetNames(SEXP data)
{
	SEXP names = Rf_getAttrib(data,R_NamesSymbol);
	string first_names = CHAR(STRING_ELT(names,0));
	SEXP rl= Rcpp::wrap(first_names);
        return rl;
}

RcppExport SEXP testMatrix(SEXP data)
{
	SEXP Rrl;
	int * rl;
	PROTECT(Rrl = Rf_allocMatrix(INTSXP,2000000,20));
	rl = INTEGER(Rrl);	
	for (int i=0;i<2000000;i++)
		for (int j=0;j<20;j++)
			rl[i+2000000*j] = (i+1)+j+1;	
	UNPROTECT(1);
        return Rrl;
}

double f (double x, void * params) {
       double alpha = *(double *) params;
       double f = log(alpha*x) / sqrt(x);
       return f;
}

/*RcppExport SEXP testhieModel(SEXP Rdata)
{
	vector<vector<double> > data;
	int n_ref = Rf_length(Rdata);
	for (int i=0;i<n_ref;i++){
		vector<double> cur_data(REAL(VECTOR_ELT(Rdata,i)), REAL(VECTOR_ELT(Rdata,i)) + Rf_length(VECTOR_ELT(Rdata,i)));
		data.push_back(cur_data);	
	}	
	
	map<string, vector<double> > result = hieModelEB(data);		
	SEXP rl=Rcpp::wrap(result);
        return rl;
}*/

RcppExport SEXP testString(SEXP RcontextEffect, SEXP RcontextEx)
{
	// load contextEx
	map<string, vector<string> >contextEx;
	SEXP Rcontext = Rf_getAttrib(RcontextEx,R_NamesSymbol);
	for (int i=0;i<Rf_length(RcontextEx);i++){
		string cur_context = CHAR(STRING_ELT(Rcontext,i)); 
		SEXP Rcur_Ex = VECTOR_ELT(RcontextEx,i);
		for (int j=0;j<Rf_length(Rcur_Ex);j++){
			contextEx[cur_context].push_back(CHAR(STRING_ELT(Rcur_Ex,j)));	
		}
	}
	// set names
	Rcontext = Rf_getAttrib(RcontextEffect, R_NamesSymbol);
	for (int i=0;i<Rf_length(RcontextEffect);i++){
		Rprintf("%d\n",i);
		string cur_context = CHAR(STRING_ELT(Rcontext,i));
		SEXP cur_effect = VECTOR_ELT(RcontextEffect,i);
		Rf_setAttrib(cur_effect, R_NamesSymbol, Rcpp::wrap(contextEx[cur_context]) );	
	}
	
	return RcontextEffect;
}

RcppExport SEXP testlgamma(SEXP Rx)
{
	vector<double> x(REAL(Rx), REAL(Rx) + Rf_length(Rx));
	vector<double> y;
	for (int i=0;i<(int)x.size();i++)
		y.push_back(lgamma(x[i]));
	return Rcpp::wrap(y);
}

/*RcppExport SEXP testGetLikelihood_marginal(SEXP data, SEXP Rtheta, SEXP Rkappa, SEXP Rupsilon, SEXP Rtau2)
{
	double theta = REAL(Rtheta)[0];
	double kappa = REAL(Rkappa)[0];
	double upsilon = REAL(Rupsilon)[0];
	double tau2 = REAL(Rtau2)[0];
	vector<vector<double> > y;
	for (int i=0;i<Rf_length(data);i++){
		SEXP cur_data = VECTOR_ELT(data,i);
		vector<double> cur_y(REAL(cur_data), REAL(cur_data) + Rf_length(cur_data));
		y.push_back(cur_y);
	}
	return Rcpp::wrap(getLogLikelihood_marginal(y, theta, kappa, upsilon, tau2));	
}*/

RcppExport SEXP testLogNormalDensity(SEXP Rx, SEXP Rmu, SEXP Rsigma2)
{
	vector<double> x(REAL(Rx), REAL(Rx)+Rf_length(Rx));
	double mu = REAL(Rmu)[0];
	double sigma2 = REAL(Rsigma2)[0];
	return Rcpp::wrap(ldnorm(x, mu, sigma2));
}


RcppExport SEXP testloadContextEffect(SEXP RcontextEffect, SEXP Rcontext)
{
	map<string, map<string, vector<double> > > contextEffect;
        if (!Rf_isNull(RcontextEffect)){
                for (int i=0;i<Rf_length(RcontextEffect);i++){
                        string cur_context = CHAR(STRING_ELT(Rcontext,i));
                        map<string, vector<double > > cur_data;
                        vector<double> tmp;
                        cur_data["mean"] = tmp; cur_data["var"] = tmp; cur_data["len"] = tmp;
                        SEXP cur_Rdata = VECTOR_ELT(RcontextEffect,i);
                        double *dat;
                        // get mean
                        dat = REAL(VECTOR_ELT(cur_Rdata,0));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,0)); j++){cur_data["mean"].push_back(dat[j]);}
                        // get var
                        dat = REAL(VECTOR_ELT(cur_Rdata,1));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,1)); j++){cur_data["var"].push_back(dat[j]);}
                        // get len
                        dat = REAL(VECTOR_ELT(cur_Rdata,2));
                        for (int j=0;j<Rf_length(VECTOR_ELT(cur_Rdata,2)); j++){cur_data["len"].push_back(dat[j]);}
                        contextEffect[cur_context] = cur_data;
                }
        }
	return Rcpp::wrap(contextEffect);
}

RcppExport SEXP testGetIdxByName(SEXP data, SEXP Rquery)
{
	string query = CHAR(STRING_ELT(Rquery,0) );
	return Rcpp::wrap(getIdxByName(data,query));
}


