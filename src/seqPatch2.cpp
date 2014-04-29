#include "seqPatch2_headers.h"


/*-------------------internal functions---------------------*/
bool loadContextHyperPara(SEXP R_context_hyperPara, SEXP R_context, map<string, map<string, double> > &context_hyperPara)
{
	if (Rf_length(R_context_hyperPara)!=Rf_length(R_context))
                printf("loadContextHyperPara : length of R_context_hyperPara and R_context should be the same.\n");
	
	for (int i=0;i<Rf_length(R_context_hyperPara);i++){
                        string cur_context = CHAR(STRING_ELT(R_context,i));
                        map<string, double> cur_data;
			SEXP R_cur_data = VECTOR_ELT(R_context_hyperPara,i);
                        
			// get data
                        cur_data["theta"] = REAL(VECTOR_ELT(R_cur_data,getIdxByName(R_cur_data, "theta")))[0];
               		cur_data["kappa"] = REAL(VECTOR_ELT(R_cur_data,getIdxByName(R_cur_data, "kappa")))[0];
			cur_data["upsilon"] = REAL(VECTOR_ELT(R_cur_data,getIdxByName(R_cur_data, "upsilon")))[0];
			cur_data["tau2"] = REAL(VECTOR_ELT(R_cur_data,getIdxByName(R_cur_data, "tau2")))[0];
			 
			// load to context_hyperPara 
			context_hyperPara[cur_context] = cur_data;
                }
	
	return true;
}

double getTstat_var_equal(double * x, int n_x, double * y, int n_y)
{
	double t = sqrt(-1);
	if (n_x + n_y <= 2) return t;
	 
	// s2 is summation of variance 
	double m_x = 0;
	double s2_x = 0;
	double m_y = 0;
	double s2_y = 0;

	for (int i=0;i<n_x;i++){
		m_x += x[i];
		s2_x += x[i]*x[i];
	} 		
	m_x = m_x / n_x;
	s2_x = s2_x - n_x*m_x*m_x;
 
	for (int i=0;i<n_y;i++){
                m_y += y[i];
                s2_y += y[i]*y[i];
        }
        m_y = m_y / n_y;
        s2_y = s2_y - n_y*m_y*m_y;
	
	double S = sqrt( (s2_x + s2_y)*(n_x+n_y)/( (n_x+n_y-2)*n_x*n_y ) );
	t = (m_x - m_y) / S;
	return t;
}
double getTstat_var_equal_ref (double x_avg, double x_var, double n_x, double y_avg, double y_var, double n_y)
{
	double t = sqrt(-1);
        if (n_x + n_y <= 2) return t;
        
	double S = sqrt( (n_x*x_var + n_y*y_var)*(n_x+n_y) /( (n_x+n_y-2)*n_x*n_y ) );
        t = (x_avg - y_avg) / S;
        return t;
}


vector<int> bin_search(double query, double *temp, int temp_len)
{
	vector<int> idx_pair(2,-1);
	if (temp_len<2) return idx_pair;
	temp_len--;	
	idx_pair[0] = 0;
	idx_pair[1] = temp_len;
	
	if (query<temp[0]){
		idx_pair[1] = 0;
		return idx_pair;
	}	
	if (query>temp[temp_len]){
		idx_pair[0] = temp_len;
                return idx_pair;
	}
	// assume temp is sorted (small to large)
	while (1){
		//Rprintf("<%d,%d>\n",idx_pair[0], idx_pair[1]);
		int cur_idx = idx_pair[0] + int((idx_pair[1] - idx_pair[0])/2.0);
		//Rprintf("cur_idx %d\n",cur_idx);	
		if (query <= temp[cur_idx]) 
			idx_pair[1] = cur_idx;
		else 
			idx_pair[0] = cur_idx;		

		if (idx_pair[1]-idx_pair[0]==1 || idx_pair[0]==temp_len || idx_pair[1]==0 )
			break;	
	
	}
	//Rprintf("finish\n");	
	return idx_pair;
}


/*-----------------------R API-----------------------*/
//RcppExport SEXP R_API_getNullDist_ID(SEXP R_genomeSeq, SEXP R_seed_len)
RcppExport SEXP R_API_getNullDist_ID(SEXP R_genomeSeq, SEXP R_left, SEXP R_right)
{
	string genomeSeq = CHAR(STRING_ELT(R_genomeSeq,0));
	//int seed_len = INTEGER(R_seed_len)[0];
	int left = INTEGER(R_left)[0];
	int right = INTEGER(R_right)[0];
	int seed_len = left + right + 1;
	EBmixture EBmixtureObj;
	map<string, vector<int> > genome_index = EBmixtureObj.buildGenomeIndex(genomeSeq, seed_len);	
	
	vector<vector<int> > max_index_all;
	for (int i=0; i<(int)genomeSeq.size(); i++){
		vector<int> max_index = EBmixtureObj.findMaxContext(i, genomeSeq, genome_index, left, right);
		max_index_all.push_back(max_index);
		if ((i+1)%10000 == 0) Rprintf("%d\r",i+1);
	}
	Rprintf("%d\n",genomeSeq.size());
	
	return Rcpp::wrap(max_index_all);
}

RcppExport SEXP R_API_getNullDist(SEXP R_IPD, SEXP R_mol_id, SEXP R_start, SEXP R_max_index)
{	
	if (Rf_length(R_IPD) != Rf_length(R_mol_id) ){
		Rprintf("length of IPD and molecule ID are not equal.\n");
		return R_NilValue;
	}
	map<string, vector<vector<double> > > rl;
			
	int start = INTEGER(R_start)[0];	
	int len = Rf_length(R_IPD);	
	for (int i=0; i<len; i++){
		int loci = start - 1 + i;
		SEXP R_cur_max_index = VECTOR_ELT(R_max_index, loci);
		vector<int> cur_max_index(INTEGER(R_cur_max_index), INTEGER(R_cur_max_index) + Rf_length(R_cur_max_index));
		
		if (cur_max_index.size()>0){
			int loci_homo = RNG_sample(cur_max_index, 1)[0];
			int i_homo = loci_homo - start + 1;
			SEXP R_IPD_homo = VECTOR_ELT(R_IPD, i_homo);
			SEXP R_mol_id_homo = VECTOR_ELT(R_mol_id, i_homo);
			vector<double> IPD_homo(REAL(R_IPD_homo), REAL(R_IPD_homo) + Rf_length(R_IPD_homo));
	                vector<double> mol_id_homo(REAL(R_mol_id_homo), REAL(R_mol_id_homo) + Rf_length(R_mol_id_homo));
			vector<double> rl_i_homo;
			rl_i_homo.push_back(i_homo);

			rl["IPD"].push_back(IPD_homo);
                        rl["mol_id"].push_back(mol_id_homo);
			rl["i_homo"].push_back(rl_i_homo);

		}else{
			vector<double> IPD_homo;
	                vector<double> mol_id_homo;
			vector<double> rl_i_homo;
			rl["IPD"].push_back(IPD_homo);	
			rl["mol_id"].push_back(mol_id_homo);
			rl["i_homo"].push_back(rl_i_homo);
		}
		if ((i+1)%10000 == 0)
			Rprintf("%d\r",i+1);
	}
	Rprintf("%d\n",len);
	return Rcpp::wrap(rl);
	//return R_NilValue;
}




RcppExport SEXP R_API_getIPDandMolID(SEXP R_prop, SEXP R_n_ccs_all, SEXP R_n_mol_all, SEXP R_mu_0, SEXP R_sd_0, SEXP R_d)
{
	map<string, vector<vector<double> > > rl;
	vector<double> prop(REAL(R_prop), REAL(R_prop) + Rf_length(R_prop));	
	vector<double> n_ccs_all(REAL(R_n_ccs_all), REAL(R_n_ccs_all) + Rf_length(R_n_ccs_all));
	vector<double> n_mol_all(REAL(R_n_mol_all), REAL(R_n_mol_all) + Rf_length(R_n_mol_all));		
	vector<double> mu_0(REAL(R_mu_0), REAL(R_mu_0) + Rf_length(R_mu_0));
	vector<double> sd_0(REAL(R_sd_0), REAL(R_sd_0) + Rf_length(R_sd_0));
	double d = REAL(R_d)[0];
		
	int len = Rf_length(R_prop);	
	vector<double> n_mol = RNG_sample(n_mol_all, len);

	for (int i=0; i<len; i++){
		int cur_n_mol = int(n_mol[i]);
		vector<double> n_ccs = RNG_sample(n_ccs_all, cur_n_mol);
		vector<double> mol_id;
		for (int j=0; j<cur_n_mol; j++)
			for (int k=0; k<n_ccs[j]; k++)
				mol_id.push_back(j+1);
		
		int cur_n_mol_mod =  int(1.0*cur_n_mol * prop[i] + 0.5);
		int cur_n_mol_unmod = cur_n_mol - cur_n_mol_mod;
	
		int cur_n_ipd_mod = 0;
		for (int j=0; j<cur_n_mol_mod; j++){
			cur_n_ipd_mod += n_ccs[j];
		}	
		int cur_n_ipd_unmod = mol_id.size() - cur_n_ipd_mod;
		vector<double> IPD = RNG_normal(cur_n_ipd_mod, mu_0[i] + d, sd_0[i]);
		vector<double> IPD_unmod = RNG_normal(cur_n_ipd_unmod, mu_0[i], sd_0[i]);			
		for (int j=0; j<cur_n_ipd_unmod; j++)
			IPD.push_back(IPD_unmod[j]);
		rl["mol_id"].push_back(mol_id);
		rl["IPD"].push_back(IPD);
		if ((i+1) % 10000 == 0) 
			Rprintf("%d\r",i+1);	
	}	
	Rprintf("%d\n",len);
	//return Rcpp::wrap(n_mol);
	return Rcpp::wrap(rl);
}

RcppExport SEXP R_API_getNumMol(SEXP R_ipd, SEXP R_mol_id)
{
	vector<double> num_mol;
	for (int i=0;i<Rf_length(R_mol_id);i++){
                double * ipd = REAL(VECTOR_ELT(R_ipd,i));
                double * mol_id = REAL(VECTOR_ELT(R_mol_id,i));
                int ipd_len = Rf_length(VECTOR_ELT(R_ipd,i));
                int mol_id_len =  Rf_length(VECTOR_ELT(R_mol_id,i));
                EBmixture EBmixtureObj;
                EBmixtureObj.getMoleculeMeanIPD(ipd, mol_id, ipd_len, mol_id_len);
                num_mol.push_back(EBmixtureObj.get_n_mol());
        	if ((i+1)%10000==0) Rprintf("%d\r",i+1);
        }
        Rprintf("%d\n", Rf_length(R_mol_id));
	return Rcpp::wrap(num_mol);

}


RcppExport SEXP R_API_getNumCCS(SEXP R_ipd, SEXP R_mol_id, SEXP R_is_collapse)
{
	int is_collapse = INTEGER(R_is_collapse)[0];
	if (is_collapse == 0){
		vector<vector<double> > num_ccs;	
		for (int i=0;i<Rf_length(R_mol_id);i++){
			double * ipd = REAL(VECTOR_ELT(R_ipd,i));
			double * mol_id = REAL(VECTOR_ELT(R_mol_id,i));
			int ipd_len = Rf_length(VECTOR_ELT(R_ipd,i));
			int mol_id_len =  Rf_length(VECTOR_ELT(R_mol_id,i));
			EBmixture EBmixtureObj;
			EBmixtureObj.getMoleculeMeanIPD(ipd, mol_id, ipd_len, mol_id_len);
			num_ccs.push_back(EBmixtureObj.get_ipd_n());
			if ((i+1)%10000==0) Rprintf("%d\r",i+1);
		}
		Rprintf("%d\n", Rf_length(R_mol_id));
		return Rcpp::wrap(num_ccs);
	}else{
		vector<double> num_ccs;
		for (int i=0;i<Rf_length(R_mol_id);i++){
                        double * ipd = REAL(VECTOR_ELT(R_ipd,i));
                        double * mol_id = REAL(VECTOR_ELT(R_mol_id,i));
                        int ipd_len = Rf_length(VECTOR_ELT(R_ipd,i));
                        int mol_id_len =  Rf_length(VECTOR_ELT(R_mol_id,i));
                        EBmixture EBmixtureObj;
                        EBmixtureObj.getMoleculeMeanIPD(ipd, mol_id, ipd_len, mol_id_len);
                        vector<double> cur_num_ccs = EBmixtureObj.get_ipd_n();
			for (int j=0; j<(int)cur_num_ccs.size(); j++){
				num_ccs.push_back(cur_num_ccs[j]);
			}
			if ((i+1)%10000==0) Rprintf("%d\r",i+1);
                }
                Rprintf("%d\n", Rf_length(R_mol_id));
                return Rcpp::wrap(num_ccs);
	}
}

RcppExport SEXP R_API_sample_subreads(SEXP R_IPD, SEXP R_idx, SEXP R_rate)
{
	double *IPD = REAL(R_IPD);	
	int IPD_len = Rf_length(R_IPD);
	double *idx = REAL(R_idx);
	int idx_len = Rf_length(R_idx); 
	double rate = REAL(R_rate)[0]; 
	
	EBmixture EBmixtureObj;
	return Rcpp::wrap(EBmixtureObj.subsample(IPD, idx, IPD_len, idx_len, rate));
}

RcppExport SEXP R_API_detectModProp_EB(SEXP R_IPD_native, SEXP R_idx_native, SEXP R_IPD_wga, SEXP R_start_native, 
			SEXP R_start_wga, SEXP R_f1_x, SEXP R_f1_y, SEXP R_is_f1_variable, SEXP R_max_iter, SEXP R_is_EM)
{
	int start_native = INTEGER(R_start_native)[0];
	int start_wga = INTEGER(R_start_wga)[0];
	int shift = start_native - start_wga;
	vector<double> f1_x(REAL(R_f1_x), REAL(R_f1_x) + Rf_length(R_f1_x));
        vector<double> f1_y(REAL(R_f1_y), REAL(R_f1_y) + Rf_length(R_f1_y));

	int is_f1_variable = INTEGER(R_is_f1_variable)[0];
	int max_iter = INTEGER(R_max_iter)[0];
	int is_EM = INTEGER(R_is_EM)[0];
	if (is_EM==0)
		Rprintf("fit mixture by variational inference\n");
	else 
		Rprintf("fit mixture by EM algorithm\n");

	int len = Rf_length(R_IPD_native);
	int len_wga = Rf_length(R_IPD_wga);

	map<string, vector<double> > rl;
        vector<double> prop(len, sqrt(-1));
        vector<double> n_mol(len, sqrt(-1));
	vector<double> N_0(len, sqrt(-1));
	vector<double> N_1(len, sqrt(-1));
	vector<double> avg_n(len, sqrt(-1));
	vector<double> cvg_wga(len, sqrt(-1));
	vector<double> n_iter(len, sqrt(-1));

	vector<double> mu_d(len, sqrt(-1));
	vector<double> mu_1(len, sqrt(-1));

	vector<double> log_likely_0(len, sqrt(-1));
	vector<double> log_likely_1(len, sqrt(-1));
	
	vector<double> BIC_0(len, sqrt(-1));	
	vector<double> BIC_1(len, sqrt(-1));

	vector<double> AIC_0(len, sqrt(-1));
        vector<double> AIC_1(len, sqrt(-1));	

	for (int i=0; i<len; i++){
		if ((i+1)%1000==0) Rprintf("processed %d sites\r",i+1);
		double * IPD_native = REAL(VECTOR_ELT(R_IPD_native,i));
		int n_IPD_native = Rf_length(VECTOR_ELT(R_IPD_native,i));			
		double * idx_native = REAL(VECTOR_ELT(R_idx_native,i));
		int n_idx_native = Rf_length(VECTOR_ELT(R_idx_native,i));

		int i_wga = i + shift;
		if (i_wga<0 || i_wga>=len_wga)
			continue;
		
		double * ipd_wga = REAL(VECTOR_ELT(R_IPD_wga, i_wga));
                int n_ipd_wga = Rf_length(VECTOR_ELT(R_IPD_wga,i_wga));

		double mu_0;
		double sigma_0;
		if (n_ipd_wga <= 2){
			mu_0 = sqrt(-1);
			sigma_0 = sqrt(-1);	
			cvg_wga[i] = n_ipd_wga;
			continue;
		}else{
			mu_0 = mean_c(ipd_wga, n_ipd_wga);
			sigma_0 = sqrt((n_ipd_wga - 1)*var_c(ipd_wga, n_ipd_wga)/n_ipd_wga);
		}
		
		if (is_EM==0){
			EBmixture EBmixtureObj;
        		EBmixtureObj.setParameters(mu_0, sigma_0, f1_x, f1_y, max_iter);
        		EBmixtureObj.getMoleculeMeanIPD(IPD_native, idx_native, n_IPD_native, n_idx_native);
			if (is_f1_variable!=0)
        			EBmixtureObj.run(true, true, 0.2);	
			else 
				EBmixtureObj.run(true, false, 0.2);
			prop[i] = EBmixtureObj.get_prop();
	                n_mol[i] = EBmixtureObj.get_n_mol();
        	        N_0[i] = EBmixtureObj.get_N_0();
                	N_1[i] = EBmixtureObj.get_N_1();
                	avg_n[i] = EBmixtureObj.get_avg_n();
                	cvg_wga[i] = n_ipd_wga;
                	n_iter[i] = EBmixtureObj.get_n_iter();
                	mu_d[i] = EBmixtureObj.get_mu_d();
                	mu_1[i] = EBmixtureObj.get_mu_1();

		}else{
			EBmixture_EM_naive EBmixtureObj;
                        EBmixtureObj.setParameters(mu_0, sigma_0, max_iter);
                        EBmixtureObj.getMoleculeMeanIPD(IPD_native, idx_native, n_IPD_native, n_idx_native);
                        //EBmixtureObj.run();
			EBmixtureObj.run_reinit();
	
			prop[i] = EBmixtureObj.get_prop();
                	n_mol[i] = EBmixtureObj.get_n_mol();
                	N_0[i] = EBmixtureObj.get_N_0();
                	N_1[i] = EBmixtureObj.get_N_1();
                	avg_n[i] = EBmixtureObj.get_avg_n();
                	cvg_wga[i] = n_ipd_wga;
                	n_iter[i] = EBmixtureObj.get_n_iter();
                	mu_1[i] = EBmixtureObj.get_mu_1();

			log_likely_0[i] = EBmixtureObj.get_log_likely_0();
			log_likely_1[i] = EBmixtureObj.get_log_likely_1();
			BIC_0[i] = EBmixtureObj.get_BIC_0();
			BIC_1[i] = EBmixtureObj.get_BIC_1();
			AIC_0[i] = EBmixtureObj.get_AIC_0();
                        AIC_1[i] = EBmixtureObj.get_AIC_1();

		}
	}	
	Rprintf("processed %d sites\n", len);

	rl["prop"] = prop;
	rl["n_mol"] = n_mol;	
	rl["N_0"] = N_0;
	rl["N_1"] = N_1;
	rl["avg_n"] = avg_n;
	rl["cvg_wga"] = cvg_wga;
	rl["n_iter"] = n_iter;
	rl["mu_d"] = mu_d;
	rl["mu_1"] = mu_1;
	
	rl["log_likely_0"] = log_likely_0;
	rl["log_likely_1"] = log_likely_1;
	rl["BIC_0"] = BIC_0;
	rl["BIC_1"] = BIC_1;
	rl["AIC_0"] = AIC_0;
	rl["AIC_1"] = AIC_1;
	return Rcpp::wrap(rl);
}
RcppExport SEXP R_API_EBmixture_EM_naive(SEXP R_IPD, SEXP R_idx, SEXP R_mu_0, SEXP R_sigma_0, SEXP R_max_iter)
{
        double * IPD = REAL(R_IPD);
        int len_IPD = Rf_length(R_IPD);
        double * idx = REAL(R_idx);
        int len_idx = Rf_length(R_idx);

        double mu_0 = REAL(R_mu_0)[0];
        double sigma_0 = REAL(R_sigma_0)[0];
        int max_iter = INTEGER(R_max_iter)[0];

        EBmixture_EM_naive EBmixtureObj;
        EBmixtureObj.setParameters(mu_0, sigma_0, max_iter);
        EBmixtureObj.getMoleculeMeanIPD(IPD, idx, len_IPD, len_idx);
        //EBmixtureObj.run();
	EBmixtureObj.run_reinit();
        return Rcpp::List::create(Rcpp::Named("ipd_avg")=Rcpp::wrap(EBmixtureObj.get_ipd_avg()),
                                Rcpp::Named("ipd_n")=Rcpp::wrap(EBmixtureObj.get_ipd_n()),
                                Rcpp::Named("ipd_var")=Rcpp::wrap(EBmixtureObj.get_ipd_var()),
                                Rcpp::Named("n_mol")=Rcpp::wrap(EBmixtureObj.get_n_mol()),
                                Rcpp::Named("prop")=Rcpp::wrap(EBmixtureObj.get_prop_track()),
                                Rcpp::Named("N_0")=Rcpp::wrap(EBmixtureObj.get_N_0_track()),
                                Rcpp::Named("N_1")=Rcpp::wrap(EBmixtureObj.get_N_1_track()),
                                Rcpp::Named("gamma_0")=Rcpp::wrap(EBmixtureObj.get_gamma_0()),
                                Rcpp::Named("gamma_1")=Rcpp::wrap(EBmixtureObj.get_gamma_1()) );

}

RcppExport SEXP R_API_EBmixture_EM(SEXP R_IPD, SEXP R_idx, SEXP R_mu_0, SEXP R_sigma_0,
                SEXP R_f1_x, SEXP R_f1_y, SEXP R_max_iter)
{
        double * IPD = REAL(R_IPD);
        int len_IPD = Rf_length(R_IPD);
        double * idx = REAL(R_idx);
        int len_idx = Rf_length(R_idx);

        double mu_0 = REAL(R_mu_0)[0];
        double sigma_0 = REAL(R_sigma_0)[0];
        vector<double> f1_x(REAL(R_f1_x), REAL(R_f1_x) + Rf_length(R_f1_x));
        vector<double> f1_y(REAL(R_f1_y), REAL(R_f1_y) + Rf_length(R_f1_y));
        int max_iter = INTEGER(R_max_iter)[0];

        EBmixture_EM EBmixtureObj;
        EBmixtureObj.setParameters(mu_0, sigma_0, f1_x, f1_y, max_iter);
        EBmixtureObj.getMoleculeMeanIPD(IPD, idx, len_IPD, len_idx);
        EBmixtureObj.run(true, 0.2);

        return Rcpp::List::create(Rcpp::Named("ipd_avg")=Rcpp::wrap(EBmixtureObj.get_ipd_avg()),
                                Rcpp::Named("ipd_n")=Rcpp::wrap(EBmixtureObj.get_ipd_n()),
                                Rcpp::Named("ipd_var")=Rcpp::wrap(EBmixtureObj.get_ipd_var()),
                                Rcpp::Named("n_mol")=Rcpp::wrap(EBmixtureObj.get_n_mol()),
                                Rcpp::Named("prop")=Rcpp::wrap(EBmixtureObj.get_prop_track()),
                                Rcpp::Named("N_0")=Rcpp::wrap(EBmixtureObj.get_N_0_track()),
                                Rcpp::Named("N_1")=Rcpp::wrap(EBmixtureObj.get_N_1_track()),
                                Rcpp::Named("f1_x")=Rcpp::wrap(EBmixtureObj.get_f1_x()),
                                Rcpp::Named("gamma_0")=Rcpp::wrap(EBmixtureObj.get_gamma_0()),
                                Rcpp::Named("gamma_1")=Rcpp::wrap(EBmixtureObj.get_gamma_1()) );

}


RcppExport SEXP R_API_EBmixture(SEXP R_IPD, SEXP R_idx, SEXP R_mu_0, SEXP R_sigma_0, 
		SEXP R_f1_x, SEXP R_f1_y, SEXP R_max_iter)
{
	double * IPD = REAL(R_IPD);
	int len_IPD = Rf_length(R_IPD);
	double * idx = REAL(R_idx);
	int len_idx = Rf_length(R_idx);
	
	double mu_0 = REAL(R_mu_0)[0];
        double sigma_0 = REAL(R_sigma_0)[0];
        vector<double> f1_x(REAL(R_f1_x), REAL(R_f1_x) + Rf_length(R_f1_x));
        vector<double> f1_y(REAL(R_f1_y), REAL(R_f1_y) + Rf_length(R_f1_y));
	int max_iter = INTEGER(R_max_iter)[0];

	EBmixture EBmixtureObj;
	EBmixtureObj.setParameters(mu_0, sigma_0, f1_x, f1_y, max_iter);
	EBmixtureObj.getMoleculeMeanIPD(IPD, idx, len_IPD, len_idx);		
	EBmixtureObj.run();

	return Rcpp::List::create(Rcpp::Named("ipd_avg")=Rcpp::wrap(EBmixtureObj.get_ipd_avg()),
				Rcpp::Named("ipd_n")=Rcpp::wrap(EBmixtureObj.get_ipd_n()),
				Rcpp::Named("ipd_var")=Rcpp::wrap(EBmixtureObj.get_ipd_var()),
				Rcpp::Named("n_mol")=Rcpp::wrap(EBmixtureObj.get_n_mol()),
                                Rcpp::Named("prop")=Rcpp::wrap(EBmixtureObj.get_prop_track()),
				Rcpp::Named("N_0")=Rcpp::wrap(EBmixtureObj.get_N_0_track()),
				Rcpp::Named("N_1")=Rcpp::wrap(EBmixtureObj.get_N_1_track()),
				Rcpp::Named("f_mu_1")=Rcpp::wrap(EBmixtureObj.get_f_mu_1()),
				Rcpp::Named("f1_x")=Rcpp::wrap(EBmixtureObj.get_f1_x()),
				Rcpp::Named("gamma_0")=Rcpp::wrap(EBmixtureObj.get_gamma_0()),
				Rcpp::Named("gamma_1")=Rcpp::wrap(EBmixtureObj.get_gamma_1()) );

}

RcppExport SEXP R_API_detectModProp_EB_pooling(SEXP R_z_score, SEXP R_mu_0, SEXP R_sigma_0, SEXP R_f1_x, SEXP R_f1_y, SEXP R_max_iter)
{
	double mu_0 = REAL(R_mu_0)[0];
        double sigma_0 = REAL(R_sigma_0)[0];
        vector<double> f1_x(REAL(R_f1_x), REAL(R_f1_x) + Rf_length(R_f1_x));
        vector<double> f1_y(REAL(R_f1_y), REAL(R_f1_y) + Rf_length(R_f1_y));
        int max_iter = INTEGER(R_max_iter)[0];
	
	map<string, vector<double> > rl;
	int genome_size = Rf_length(R_z_score);
	for (int i=0;i<genome_size;i++){
		if ((i+1)%1000==0) Rprintf("%d\r", i+1); 
		SEXP cur_R_z_score = VECTOR_ELT(R_z_score, i);
		vector<double> z_score(REAL(cur_R_z_score), REAL(cur_R_z_score) + Rf_length(cur_R_z_score));
		
		EBmixture_pooling EBmixtureObj;
        	EBmixtureObj.setParameters(mu_0 , sigma_0, f1_x, f1_y);
        	EBmixtureObj.setData(z_score);
        	EBmixtureObj.setMaxIter(max_iter);
        	EBmixtureObj.run();
	
		rl["N_0"].push_back(EBmixtureObj.get_N_0());
		rl["N_1"].push_back(EBmixtureObj.get_N_1());
		rl["prop"].push_back(EBmixtureObj.get_prop());
	}	
	Rprintf("%d\n", genome_size);
	return Rcpp::wrap(rl);
}

RcppExport SEXP R_API_EBmixture_pooling(SEXP R_z_score, SEXP R_mu_0, SEXP R_sigma_0, SEXP R_f1_x, SEXP R_f1_y, SEXP R_max_iter)
{
	vector<double> z_score(REAL(R_z_score), REAL(R_z_score) + Rf_length(R_z_score));
	double mu_0 = REAL(R_mu_0)[0];
	double sigma_0 = REAL(R_sigma_0)[0];
	vector<double> f1_x(REAL(R_f1_x), REAL(R_f1_x) + Rf_length(R_f1_x));
        vector<double> f1_y(REAL(R_f1_y), REAL(R_f1_y) + Rf_length(R_f1_y));	
	int max_iter = INTEGER(R_max_iter)[0];	

	EBmixture_pooling EBmixtureObj;
	EBmixtureObj.setParameters(mu_0 , sigma_0, f1_x, f1_y);
	EBmixtureObj.setData(z_score);
	EBmixtureObj.setMaxIter(max_iter);
	EBmixtureObj.run();
	
	
	return Rcpp::List::create(Rcpp::Named("rho_0")=Rcpp::wrap(EBmixtureObj.get_rho_0()),
				Rcpp::Named("rho_1")=Rcpp::wrap(EBmixtureObj.get_rho_1()),
				Rcpp::Named("gamma_0")=Rcpp::wrap(EBmixtureObj.get_gamma_0()),
				Rcpp::Named("gamma_1")=Rcpp::wrap(EBmixtureObj.get_gamma_1()),
				Rcpp::Named("N_0")=Rcpp::wrap(EBmixtureObj.get_N_0_track()),
				Rcpp::Named("N_1")=Rcpp::wrap(EBmixtureObj.get_N_1_track()),
				Rcpp::Named("prop")=Rcpp::wrap(EBmixtureObj.get_prop_track()));
}


RcppExport SEXP R_API_bin_search(SEXP R_query, SEXP R_temp, SEXP R_temp_len)
{
	double query = REAL(R_query)[0];
	double *temp = REAL(R_temp);
	int temp_len = INTEGER(R_temp_len)[0];			
	//vector<int> rl = bin_search(query, temp, temp_len);
	EBmixture_pooling EBmixtureObj;
	vector<int> rl = EBmixtureObj.bin_search(query, temp, temp_len);
	Rprintf("<%d,%d>\n",rl[0],rl[1]);
	//return Rcpp::wrap(rl);
	return R_NilValue;	
}

RcppExport SEXP R_API_unlist(SEXP x)
{
	vector<double> z;
	vector<double> idx;
	map<string, vector<double> > rl;
	rl["z"] = z;
	rl["idx"] = idx;
	int len = Rf_length(x);
	for (int i=0;i<len;i++){
		double * cur_x = REAL(VECTOR_ELT(x,i));
		int cur_x_len = Rf_length(VECTOR_ELT(x,i));
		for (int j=0;j<cur_x_len;j++){
			rl["z"].push_back(cur_x[j]);
			rl["idx"].push_back(i+1);
		}	
	}	
	return Rcpp::wrap(rl);	
}

RcppExport SEXP R_API_getTstat_var_equal(SEXP R_x, SEXP R_y)
{
	double * x = REAL(R_x);
	int x_len = Rf_length(R_x);
	double * y = REAL(R_y);
	int y_len = Rf_length(R_y);

	return Rcpp::wrap(getTstat_var_equal(x, x_len, y, y_len));
}

RcppExport SEXP R_API_getLocfdrSingleMolecule(SEXP R_z, SEXP R_x, SEXP R_fdr)
{	
	double *x = REAL(R_x);
	int x_len = Rf_length(R_x);
	double *fdr = REAL(R_fdr);
	int fdr_len = Rf_length(R_fdr);
	if (x_len != fdr_len){
		Rprintf("length of x and fdr should be the same\n");
		return R_NilValue;
	}
	
	vector<vector<double> > z_fdr;
	int z_len = Rf_length(R_z);	
	for (int i=0;i<z_len;i++){
		if ((i+1)%100000==0) Rprintf("processed %d sites\r", i+1);
		vector<double> cur_z_fdr;
		double *cur_z = REAL(VECTOR_ELT(R_z,i));
		int cur_z_len = Rf_length(VECTOR_ELT(R_z,i));
		for (int j=0;j<cur_z_len;j++){
			// remove NaN
			if (cur_z[j]!=cur_z[j]){
				cur_z_fdr.push_back(sqrt(-1));
				continue;
			}
			// remove Inf
			if (cur_z[j]>1e12){
				cur_z_fdr.push_back(sqrt(-1));
				continue;
			}
			// get fdr
			if (cur_z[j]<0){
				cur_z_fdr.push_back(1);
				continue;
			}
			vector<int> idx_pair = bin_search(cur_z[j], x, x_len);
			if (idx_pair[1]==0){
				cur_z_fdr.push_back(fdr[0]);
				continue;
			}
			if (idx_pair[0]==fdr_len-1){
                                cur_z_fdr.push_back(fdr[fdr_len-1]);
                                continue;
                        }
			double fdr_pred = fdr[idx_pair[0]]+(cur_z[j] - x[idx_pair[0]])* 
					(fdr[idx_pair[1]] - fdr[idx_pair[0]]) / (x[idx_pair[1]] - x[idx_pair[0]]);
			cur_z_fdr.push_back(fdr_pred);
		}
		z_fdr.push_back(cur_z_fdr);	
	}
	Rprintf("processed %d sites\n", z_len);
	return Rcpp::wrap(z_fdr);
	/*double *z = REAL(R_z);
	Rprintf("%lf\n",z[0]);
	if (z[0]>=1e12) Rprintf("yes");*/
	//return R_NilValue;
}

RcppExport SEXP R_API_getZscoreSingleMolecule_CC(SEXP R_IPD_native, SEXP R_IPD_ctrl, SEXP R_mol_id_native, SEXP R_mol_id_ctrl,
					SEXP R_start_native, SEXP R_start_ctrl, SEXP R_is_z, SEXP R_is_delta)
{
	int start_native = INTEGER(R_start_native)[0];
	int start_ctrl = INTEGER(R_start_ctrl)[0];
	int shift = start_native - start_ctrl;
	int len = Rf_length(R_IPD_native);	
	int is_z = INTEGER(R_is_z)[0];	
	int is_delta = INTEGER(R_is_delta)[0];

	vector<vector<double> > t_score;
	for (int i=0;i<len;i++){
		if ((i+1)%10000==0) Rprintf("processed %d sites\r",i+1);
		double * ipd_native = REAL(VECTOR_ELT(R_IPD_native,i));
		int n_ipd_native = Rf_length(VECTOR_ELT(R_IPD_native,i));			
		double * mol_idx_native = REAL(VECTOR_ELT(R_mol_id_native,i));
		int n_mol_idx_native = Rf_length(VECTOR_ELT(R_mol_id_native,i));

		int i_ctrl = i + shift;
		if (i_ctrl<0){
			vector<double> cur_t_score;
			cur_t_score.push_back(sqrt(-1));
			t_score.push_back(cur_t_score);
			continue;
		}
		double * ipd_ctrl = REAL(VECTOR_ELT(R_IPD_ctrl,i_ctrl));
                int n_ipd_ctrl = Rf_length(VECTOR_ELT(R_IPD_ctrl,i_ctrl));
                //double * mol_idx_ctrl = REAL(VECTOR_ELT(R_mol_id_ctrl,i_ctrl));
                //int n_mol_idx_ctrl = Rf_length(VECTOR_ELT(R_mol_id_ctrl,i_ctrl));
		
		double ipd_avg_ctrl;
		double ipd_var_ctrl;
		if (n_ipd_ctrl < 2){
			ipd_avg_ctrl = sqrt(-1);
			ipd_var_ctrl = sqrt(-1);	
		}else{
			ipd_avg_ctrl = mean_c(ipd_ctrl, n_ipd_ctrl);
			ipd_var_ctrl = (n_ipd_ctrl - 1)*var_c(ipd_ctrl, n_ipd_ctrl)/n_ipd_ctrl;
		}
			
	

		BayesianMixtureModel BayesianMixtureModel_obj;
		BayesianMixtureModel_obj.getMoleculeMeanIPD(ipd_native, mol_idx_native, n_ipd_native, n_mol_idx_native);
		vector <double> ipd_avg_native = BayesianMixtureModel_obj.get_ipd_avg();
		vector <double> ipd_var_native = BayesianMixtureModel_obj.get_ipd_var();
		vector <double> ipd_n_native = BayesianMixtureModel_obj.get_ipd_n();
	
		int cur_n_mol = ipd_avg_native.size();
		vector<double> cur_t_score;
		for (int j=0;j<cur_n_mol;j++){
			double t;
			if (is_delta == 1){
			 	t = ipd_avg_native[j] - ipd_avg_ctrl;
			}else{
				t = getTstat_var_equal_ref(ipd_avg_native[j], ipd_var_native[j], ipd_n_native[j],
						ipd_avg_ctrl, ipd_var_ctrl, n_ipd_ctrl);
				if (is_z!=0){
					double df = n_ipd_ctrl + ipd_n_native[j] - 2;
					if (df < 0) t = sqrt(-1);
					else t = Rf_qnorm5(Rf_pt(t, df, 1, 1), 0, 1, 1, 1);
				}
			}	
			cur_t_score.push_back(t);
		}
		t_score.push_back(cur_t_score);
	}
	Rprintf("processed %d sites\n",len);
	return Rcpp::wrap(t_score);
	//return R_NilValue;
}

RcppExport SEXP R_API_BayesianMixtureModel (SEXP R_IPD, SEXP R_idx, SEXP R_theta_0, SEXP R_kappa_0, SEXP R_upsilon_0, SEXP R_tau2_0, 
				SEXP R_theta_1, SEXP R_kappa_1, SEXP R_upsilon_1, SEXP R_tau2_1, SEXP R_max_iter)
{
	BayesianMixtureModel_NC BayesianMixtureModelObj;
	//BayesianMixtureModelObj.setHyperParametersNull(REAL(R_theta_0)[0], REAL(R_kappa_0)[0], REAL(R_upsilon_0)[0], REAL(R_tau2_0)[0]);
	BayesianMixtureModelObj.setHyperParameters(REAL(R_theta_0)[0], REAL(R_kappa_0)[0], REAL(R_upsilon_0)[0], REAL(R_tau2_0)[0],
						REAL(R_theta_1)[0], REAL(R_kappa_1)[0], REAL(R_upsilon_1)[0], REAL(R_tau2_1)[0]);
	BayesianMixtureModelObj.getMoleculeMeanIPD(REAL(R_IPD), REAL(R_idx), Rf_length(R_IPD), Rf_length(R_idx));
	BayesianMixtureModelObj.run(INTEGER(R_max_iter)[0]);		
	

	return Rcpp::List::create(Rcpp::Named("theta_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_theta_0_t_track()),
				Rcpp::Named("kappa_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_kappa_0_t_track()),
	 			Rcpp::Named("upsilon_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_upsilon_0_t_track()),
				Rcpp::Named("tau2_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_tau2_0_t_track()),
				Rcpp::Named("N_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_N_0_t_track()),
				Rcpp::Named("N_gamma_0_t")=Rcpp::wrap(BayesianMixtureModelObj.get_N_gamma_0_t_track()),

				Rcpp::Named("theta_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_theta_1_t_track()),
                                Rcpp::Named("kappa_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_kappa_1_t_track()),
                                Rcpp::Named("upsilon_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_upsilon_1_t_track()),
                                Rcpp::Named("tau2_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_tau2_1_t_track()),
                                Rcpp::Named("N_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_N_1_t_track()),
                                Rcpp::Named("N_gamma_1_t")=Rcpp::wrap(BayesianMixtureModelObj.get_N_gamma_1_t_track()),
		
				Rcpp::Named("gamma_0")=Rcpp::wrap(BayesianMixtureModelObj.get_gamma_0()),
				Rcpp::Named("gamma_1")=Rcpp::wrap(BayesianMixtureModelObj.get_gamma_1()),
				
				Rcpp::Named("ipd_avg")=Rcpp::wrap(BayesianMixtureModelObj.get_ipd_avg()),
				Rcpp::Named("ipd_n")=Rcpp::wrap(BayesianMixtureModelObj.get_ipd_n()),
				Rcpp::Named("ipd_var")=Rcpp::wrap(BayesianMixtureModelObj.get_ipd_var()));		

      	return R_NilValue;
}

RcppExport SEXP R_API_DetectModProp_CC(SEXP R_IPD, SEXP R_idx, SEXP R_IPD_ctrl, SEXP R_genome_start, SEXP R_genome_start_ctrl, SEXP R_max_iter)
{
        /*---------------------check inputs---------------------*/

	if (Rf_length(R_IPD)!=Rf_length(R_idx)){
                Rprintf("inconsistent IPD and moleculeID.\n");
                return R_NilValue;
        }
	/*---------------- load data and parameters -----------------------*/

	int genome_start = INTEGER(R_genome_start)[0];
	int genome_start_ctrl = INTEGER(R_genome_start_ctrl)[0];
	int max_iter = INTEGER(R_max_iter)[0];
	
	int shift = genome_start - genome_start_ctrl;

	/*------------ initialized results that need to be recorded ----------*/

	// prior distribution parameters ( " *_0 " for null distribution and " *_1 " for alternative(modified) distribution )
	vector<double> theta_0(Rf_length(R_IPD), sqrt(-1)); 
	vector<double> kappa_0(Rf_length(R_IPD), sqrt(-1)); 
	vector<double> upsilon_0(Rf_length(R_IPD), sqrt(-1)); 
	vector<double> tau2_0(Rf_length(R_IPD), sqrt(-1));
	
	// posterior distribution parameters
	vector<double> theta_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> theta_1_t(Rf_length(R_IPD), sqrt(-1));
        vector<double> kappa_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> kappa_1_t(Rf_length(R_IPD), sqrt(-1));
        vector<double> upsilon_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> upsilon_1_t(Rf_length(R_IPD), sqrt(-1));
        vector<double> tau2_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> tau2_1_t(Rf_length(R_IPD), sqrt(-1));

	vector<double> N_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> N_1_t(Rf_length(R_IPD), sqrt(-1));
	vector<double> N_gamma_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> N_gamma_1_t(Rf_length(R_IPD), sqrt(-1));

	vector<double> prop_mod(Rf_length(R_IPD), sqrt(-1));
	
	// number of iteration, coverge 
	vector<double> n_steps(Rf_length(R_IPD), sqrt(-1));
	vector<double> cvg(Rf_length(R_IPD), sqrt(-1));	
	
	/*----------------- run ------------------*/
	Rprintf("run.\n");
	BayesianMixtureModel_NC BayesianMixtureModelObj;
	int genome_size = Rf_length(R_IPD);
	for (int i=0;i<genome_size;i++){
		if ((i+1)%10000==0) Rprintf("detected %d positions\r", i+1);
		if (i + shift<0 || i+shift>=Rf_length(R_IPD_ctrl)) continue;
		
		// check data and load IPD and moleculeID
		if (Rf_length(VECTOR_ELT(R_IPD,i)) != Rf_length(VECTOR_ELT(R_idx,i))){
			Rprintf("in the %dth position: inconsistent IPD and moleculeID length.\n",i+1);
			return R_NilValue;
		}
		vector<double> cur_IPD(REAL(VECTOR_ELT(R_IPD,i)), REAL(VECTOR_ELT(R_IPD,i)) + Rf_length(VECTOR_ELT(R_IPD,i)));
		vector<double> cur_idx(REAL(VECTOR_ELT(R_idx,i)), REAL(VECTOR_ELT(R_idx,i)) + Rf_length(VECTOR_ELT(R_idx,i)));		
		vector<double> cur_IPD_ctrl(REAL(VECTOR_ELT(R_IPD_ctrl,i+shift)), REAL(VECTOR_ELT(R_IPD_ctrl,i+shift)) + Rf_length(VECTOR_ELT(R_IPD_ctrl,i+shift)));
		
		// get prior from control data 
		if (cur_IPD.size() < 3 || cur_IPD_ctrl.size() < 3 ) continue;
		theta_0[i] = mean(cur_IPD_ctrl);
		kappa_0[i] = 1000000;
		upsilon_0[i] = 1000000;
		tau2_0[i] = var(cur_IPD_ctrl);
		
		// fit Bayesian Mixture Model
		BayesianMixtureModelObj.setHyperParameters(theta_0[i], kappa_0[i], upsilon_0[i], tau2_0[i],
							theta_0[i],1/10000,upsilon_0[i],tau2_0[i]);
		BayesianMixtureModelObj.getMoleculeMeanIPD(&cur_IPD[0], &cur_idx[0], cur_IPD.size(), cur_idx.size());
		BayesianMixtureModelObj.run(max_iter);

		// get results	
		theta_0_t[i] = BayesianMixtureModelObj.get_theta_0_t();
		kappa_0_t[i] = BayesianMixtureModelObj.get_kappa_0_t();							
		upsilon_0_t[i] = BayesianMixtureModelObj.get_upsilon_0_t();
		tau2_0_t[i] = BayesianMixtureModelObj.get_tau2_0_t();

		theta_1_t[i] = BayesianMixtureModelObj.get_theta_1_t();
                kappa_1_t[i] = BayesianMixtureModelObj.get_kappa_1_t();
                upsilon_1_t[i] = BayesianMixtureModelObj.get_upsilon_1_t();
                tau2_1_t[i] = BayesianMixtureModelObj.get_tau2_1_t();

		N_0_t[i] = BayesianMixtureModelObj.get_N_0_t();	
		N_1_t[i] = BayesianMixtureModelObj.get_N_1_t();

		N_gamma_0_t[i] = BayesianMixtureModelObj.get_N_gamma_0_t();
                N_gamma_1_t[i] = BayesianMixtureModelObj.get_N_gamma_1_t();

		prop_mod[i] = N_1_t[i] / (N_0_t[i] + N_1_t[i]);

		n_steps[i] = BayesianMixtureModelObj.get_n_steps();

		cvg[i] = cur_IPD.size();	
	}	
	Rprintf("detected %d positions\n", genome_size);

	/*--------------------- output --------------------------*/
	map<string,vector<double> > rl;
	rl["theta_0"] = theta_0;
	rl["kappa_0"] = kappa_0;
	rl["upsilon_0"] = upsilon_0;
	rl["tau2_0"] = tau2_0;

	rl["theta_0_t"] = theta_0_t;
        rl["kappa_0_t"] = kappa_0_t;
        rl["upsilon_0_t"] = upsilon_0_t;
        rl["tau2_0_t"] = tau2_0_t;

	rl["theta_1_t"] = theta_1_t;
        rl["kappa_1_t"] = kappa_1_t;
        rl["upsilon_1_t"] = upsilon_1_t;
        rl["tau2_1_t"] = tau2_1_t;

	rl["N_0_t"] = N_0_t;
        rl["N_1_t"] = N_1_t;
        rl["N_gamma_0_t"] = N_gamma_0_t;
        rl["N_gamma_1_t"] = N_gamma_1_t;
        
	rl["prop_mod"] = prop_mod;
       	rl["n_steps"] = n_steps;
        rl["cvg"] = cvg;

	return Rcpp::wrap(rl);
	//return R_NilValue;
}

RcppExport SEXP R_API_DetectModProp_NC(SEXP R_IPD, SEXP R_idx, SEXP R_genome_start, SEXP R_genomeSeq, SEXP R_context_hyperPara,
				 SEXP R_context, SEXP R_strand, SEXP R_left_len, SEXP R_right_len, SEXP R_max_iter)
{
	/*---------------------check inputs---------------------*/
	if (Rf_length(R_IPD)!=Rf_length(R_idx)){
		Rprintf("inconsistent IPD and moleculeID.\n");
		return R_NilValue;
	}	

	/*---------------- load data and parameters -----------------------*/
	
	// load context_hyperPara
	Rprintf("load context effect hyperparameters.\n");
	map<string, map<string, double> > context_hyperPara;
	loadContextHyperPara(R_context_hyperPara, R_context, context_hyperPara);
	
	// load parameters
	int genome_start = INTEGER(R_genome_start)[0];
	string genomeSeq = CHAR(STRING_ELT(R_genomeSeq,0));
	int strand = INTEGER(R_strand)[0];
	int left_len = 	INTEGER(R_left_len)[0];
	int right_len = INTEGER(R_right_len)[0];
	int max_iter = INTEGER(R_max_iter)[0];

	/*------------ initialized results that need to be recorded ----------*/
	
	vector<double> is_findContext(Rf_length(R_IPD), sqrt(-1));	
	
	// prior distribution parameters ( " *_0 " for null distribution and " *_1 " for alternative(modified) distribution )		
	vector<double> theta_0(Rf_length(R_IPD), sqrt(-1)); 
	vector<double> kappa_0(Rf_length(R_IPD), sqrt(-1)); 
	vector<double> upsilon_0(Rf_length(R_IPD), sqrt(-1)); 
	vector<double> tau2_0(Rf_length(R_IPD), sqrt(-1));

	// posterior distribution parameters 
	vector<double> theta_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> theta_1_t(Rf_length(R_IPD), sqrt(-1));
        vector<double> kappa_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> kappa_1_t(Rf_length(R_IPD), sqrt(-1));
        vector<double> upsilon_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> upsilon_1_t(Rf_length(R_IPD), sqrt(-1));
        vector<double> tau2_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> tau2_1_t(Rf_length(R_IPD), sqrt(-1));

	vector<double> N_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> N_1_t(Rf_length(R_IPD), sqrt(-1));
	vector<double> N_gamma_0_t(Rf_length(R_IPD), sqrt(-1)); vector<double> N_gamma_1_t(Rf_length(R_IPD), sqrt(-1));

	vector<double> prop_mod(Rf_length(R_IPD), sqrt(-1));
	
	// number of iteration, coverge 
	vector<double> n_steps(Rf_length(R_IPD), sqrt(-1));
	vector<double> cvg(Rf_length(R_IPD), sqrt(-1));	

	/*----------------- run ------------------*/
	Rprintf("run.\n");
	BayesianMixtureModel_NC BayesianMixtureModelObj;
	int genome_size = Rf_length(R_IPD);
	for (int i=0;i<genome_size;i++){
		if ((i+1)%10000==0) Rprintf("detected %d positions\r", i+1);

		// check data and load IPD and moleculeID
		if (Rf_length(VECTOR_ELT(R_IPD,i)) != Rf_length(VECTOR_ELT(R_idx,i))){
			Rprintf("in the %dth position: inconsistent IPD and moleculeID length.\n",i+1);
			return R_NilValue;
		}
		vector<double> cur_IPD(REAL(VECTOR_ELT(R_IPD,i)), REAL(VECTOR_ELT(R_IPD,i)) + Rf_length(VECTOR_ELT(R_IPD,i)));
		vector<double> cur_idx(REAL(VECTOR_ELT(R_idx,i)), REAL(VECTOR_ELT(R_idx,i)) + Rf_length(VECTOR_ELT(R_idx,i)));
		
		// find context of current position
		int cur_pos = i + genome_start - 1;
		string cur_context;
		if (strand == 0){
			if (cur_pos - left_len<0 || cur_pos + right_len>Rf_length(R_IPD)-1 || cur_pos + right_len > (int) genomeSeq.size()-1)
				continue;
			cur_context = genomeSeq.substr(cur_pos - left_len , right_len + left_len + 1);		
		}else{
			if (cur_pos - right_len<0 || cur_pos + left_len>Rf_length(R_IPD)-1 || cur_pos + left_len > (int) genomeSeq.size()-1)
                                continue;
			cur_context = SP_reverseSeq (genomeSeq.substr(cur_pos - right_len, right_len + left_len + 1) ) ;
		}
		
		if (cur_IPD.size()<3){
                	cvg[i] = cur_IPD.size();
                	continue;
        	}
		
		// find context hyperparameters of current position
		map<string, map<string, double> >::iterator it = context_hyperPara.find(cur_context);
               	if (it != context_hyperPara.end()){
			is_findContext[i] = 1;
		}else{
			is_findContext[i] = 0;
			continue;	
		}

		// load hyperparameters of current position
		theta_0[i] = it->second["theta"];
		kappa_0[i] = it->second["kappa"];
		upsilon_0[i] = it->second["upsilon"];
		tau2_0[i] = it->second["tau2"];
		
		// fit Bayesian Mixture Model	
		//BayesianMixtureModelObj.setHyperParametersNull(theta_0[i], kappa_0[i], upsilon_0[i], tau2_0[i]);
		BayesianMixtureModelObj.setHyperParameters(theta_0[i], kappa_0[i]+10000, upsilon_0[i]+1000, tau2_0[i],
							theta_0[i],0.0001,0.0001,tau2_0[i]);
		BayesianMixtureModelObj.getMoleculeMeanIPD(&cur_IPD[0], &cur_idx[0], cur_IPD.size(), cur_idx.size());
		BayesianMixtureModelObj.run(max_iter);
		
		// get results
		theta_0_t[i] = BayesianMixtureModelObj.get_theta_0_t();
		kappa_0_t[i] = BayesianMixtureModelObj.get_kappa_0_t();							
		upsilon_0_t[i] = BayesianMixtureModelObj.get_upsilon_0_t();
		tau2_0_t[i] = BayesianMixtureModelObj.get_tau2_0_t();

		theta_1_t[i] = BayesianMixtureModelObj.get_theta_1_t();
                kappa_1_t[i] = BayesianMixtureModelObj.get_kappa_1_t();
                upsilon_1_t[i] = BayesianMixtureModelObj.get_upsilon_1_t();
                tau2_1_t[i] = BayesianMixtureModelObj.get_tau2_1_t();

		N_0_t[i] = BayesianMixtureModelObj.get_N_0_t();	
		N_1_t[i] = BayesianMixtureModelObj.get_N_1_t();
		
		N_gamma_0_t[i] = BayesianMixtureModelObj.get_N_gamma_0_t();
                N_gamma_1_t[i] = BayesianMixtureModelObj.get_N_gamma_1_t();
		
		prop_mod[i] = N_1_t[i] / (N_0_t[i] + N_1_t[i]);
	
		n_steps[i] = BayesianMixtureModelObj.get_n_steps();
	
		cvg[i] = cur_IPD.size();

	}
	Rprintf("detected %d positions\n", genome_size);

	/*--------------------- output --------------------------*/
	map<string,vector<double> > rl;
	rl["theta_0"] = theta_0;
	rl["kappa_0"] = kappa_0;
	rl["upsilon_0"] = upsilon_0;
	rl["tau2_0"] = tau2_0;

	rl["theta_0_t"] = theta_0_t;
        rl["kappa_0_t"] = kappa_0_t;
        rl["upsilon_0_t"] = upsilon_0_t;
        rl["tau2_0_t"] = tau2_0_t;

	rl["theta_1_t"] = theta_1_t;
        rl["kappa_1_t"] = kappa_1_t;
        rl["upsilon_1_t"] = upsilon_1_t;
        rl["tau2_1_t"] = tau2_1_t;

	rl["N_0_t"] = N_0_t;
        rl["N_1_t"] = N_1_t;
        rl["N_gamma_0_t"] = N_gamma_0_t;
        rl["N_gamma_1_t"] = N_gamma_1_t;
        
	rl["prop_mod"] = prop_mod;
       	rl["n_steps"] = n_steps;
        rl["cvg"] = cvg;
        rl["is_findContext"] = is_findContext;

	return Rcpp::wrap(rl);

	/*return Rcpp::List::create(Rcpp::Named("theta_0")=Rcpp::wrap(theta_0),
				Rcpp::Named("kappa_0")=Rcpp::wrap(kappa_0),
				Rcpp::Named("upsilon_0")=Rcpp::wrap(upsilon_0),	
				Rcpp::Named("tau2_0")=Rcpp::wrap(tau2_0),
				Rcpp::Named("theta_0_t")=Rcpp::wrap(theta_0_t),
				Rcpp::Named("kappa_0_t")=Rcpp::wrap(kappa_0_t),
				Rcpp::Named("upsilon_0_t")=Rcpp::wrap(upsilon_0_t),
				Rcpp::Named("tau2_0_t")=Rcpp::wrap(tau2_0_t),
				Rcpp::Named("theta_1_t")=Rcpp::wrap(theta_1_t),
                                Rcpp::Named("kappa_1_t")=Rcpp::wrap(kappa_1_t),
                                Rcpp::Named("upsilon_1_t")=Rcpp::wrap(upsilon_1_t),
                                Rcpp::Named("tau2_1_t")=Rcpp::wrap(tau2_1_t),
				Rcpp::Named("N_0_t")=Rcpp::wrap(N_0_t),
				Rcpp::Named("N_1_t")=Rcpp::wrap(N_1_t),
				Rcpp::Named("N_gamma_0_t")=Rcpp::wrap(N_gamma_0_t),
				Rcpp::Named("N_gamma_1_t")=Rcpp::wrap(N_gamma_1_t),
				Rcpp::Named("prop_mod")=Rcpp::wrap(prop_mod),
				Rcpp::Named("n_steps")=Rcpp::wrap(n_steps),
				Rcpp::Named("cvg")=Rcpp::wrap(cvg),
				Rcpp::Named("is_findContext")= Rcpp::wrap(is_findContext));	
	*/
	//return Rcpp::wrap(context_hyperPara);	
        //return R_NilValue;
}


RcppExport SEXP R_API_hieModel_core(SEXP R_IPD_mean, SEXP R_IPD_var, SEXP R_IPD_n, SEXP R_max_iter)
{
	vector<double> IPD_mean(REAL(R_IPD_mean), REAL(R_IPD_mean) + Rf_length(R_IPD_mean) );
        vector<double> IPD_var(REAL(R_IPD_var), REAL(R_IPD_var) + Rf_length(R_IPD_var) );
        vector<double> IPD_n(REAL(R_IPD_n), REAL(R_IPD_n) + Rf_length(R_IPD_n) );
	int max_iter = INTEGER(R_max_iter)[0];
	
	return Rcpp::wrap(hieModelEB(IPD_mean, IPD_var, IPD_n, max_iter) );	
}

 

RcppExport SEXP R_API_writeGenomeFtoCSV(SEXP R_IPD, SEXP R_mol_ID, SEXP R_ref_pos, SEXP R_file_name)
{
	string file_name = CHAR(STRING_ELT(R_file_name,0));
	double * ref_pos = REAL(R_ref_pos);
	int genome_size = Rf_length(R_ref_pos);
	
	FILE *outfile = fopen(file_name.c_str(), "w");
	for (int i=0;i<genome_size;i++){
		double * IPD = REAL(VECTOR_ELT(R_IPD,i));
		double * mol_ID = REAL(VECTOR_ELT(R_mol_ID,i));
		int cvg = Rf_length(VECTOR_ELT(R_IPD,i));
		if (cvg!=Rf_length(VECTOR_ELT(R_mol_ID,i))){
			Rprintf("inconsistent length of IPD and moleculeID in %d.\n", ref_pos[i]);
			return R_NilValue;
		}
		for (int j=0;j<cvg;j++) fprintf(outfile, "%.0lf\t%.0lf\t%lf\n", mol_ID[j], ref_pos[i], IPD[j]);	
		if ((i+1)%10000==0) Rprintf("converted %d positions.\r", i+1);
	}
	Rprintf("converted %d positions.\n", genome_size);
	fclose(outfile);
	return R_NilValue;

}

RcppExport SEXP R_API_write_z_score_to_CSV(SEXP R_z, SEXP R_ref_pos, SEXP R_file_name)
{
	string file_name = CHAR(STRING_ELT(R_file_name,0));
	double * ref_pos = REAL(R_ref_pos);
	int genome_size = Rf_length(R_ref_pos);
	
	FILE *outfile = fopen(file_name.c_str(), "w");
        for (int i=0;i<genome_size;i++){
                double * z_score = REAL(VECTOR_ELT(R_z,i));
                int cvg = Rf_length(VECTOR_ELT(R_z,i));
                for (int j=0;j<cvg;j++) fprintf(outfile, "%.0lf\t%lf\n", ref_pos[i], z_score[j]);
                if ((i+1)%10000==0) Rprintf("converted %d positions.\r", i+1);
        }
        Rprintf("converted %d positions.\n", genome_size);
        fclose(outfile);
		
	return R_NilValue;
}


/*------------------- test -----------------------*/

RcppExport SEXP test_f0(SEXP R_x_avg, SEXP R_x_var, SEXP R_x_n)
{
	double x_avg = REAL(R_x_avg)[0];
	double x_var = REAL(R_x_var)[0];	
	double x_n = REAL(R_x_n)[0];
	
	EBmixture EBmixtureObj;
	return Rcpp::wrap(EBmixtureObj.f0(x_avg, x_var, x_n));
}

RcppExport SEXP test_f0_log(SEXP R_x_avg, SEXP R_x_var, SEXP R_x_n, SEXP R_mu_0, SEXP R_sigma_0, SEXP R_f1_x, SEXP R_f1_y)
{
        double x_avg = REAL(R_x_avg)[0];
        double x_var = REAL(R_x_var)[0];
        double x_n = REAL(R_x_n)[0];

	double mu_0 = REAL(R_mu_0)[0];
        double sigma_0 = REAL(R_sigma_0)[0];
        vector<double> f1_x(REAL(R_f1_x), REAL(R_f1_x) + Rf_length(R_f1_x));
        vector<double> f1_y(REAL(R_f1_y), REAL(R_f1_y) + Rf_length(R_f1_y));
        EBmixture EBmixtureObj;
        EBmixtureObj.setParameters(mu_0, sigma_0, f1_x, f1_y, 1);
	

        return Rcpp::wrap(EBmixtureObj.f0_log(x_avg, x_var, x_n));
}

RcppExport SEXP test_f1_log_int(SEXP R_x_avg, SEXP R_x_var, SEXP R_x_n, SEXP R_mu_0, SEXP R_sigma_0, SEXP R_f1_x, SEXP R_f1_y)
{
        double x_avg = REAL(R_x_avg)[0];
        double x_var = REAL(R_x_var)[0];
        double x_n = REAL(R_x_n)[0];

	double mu_0 = REAL(R_mu_0)[0];
	double sigma_0 = REAL(R_sigma_0)[0];
	vector<double> f1_x(REAL(R_f1_x), REAL(R_f1_x) + Rf_length(R_f1_x));
	vector<double> f1_y(REAL(R_f1_y), REAL(R_f1_y) + Rf_length(R_f1_y));
        EBmixture EBmixtureObj;
	EBmixtureObj.setParameters(mu_0, sigma_0, f1_x, f1_y, 0);	
        EBmixtureObj.run();
	return Rcpp::wrap(EBmixtureObj.f1_log_int(x_avg, x_var, x_n));
}


RcppExport SEXP test_getMoleculeMeanIPD_subsample(SEXP R_IPD, SEXP R_idx, SEXP R_rate)
{
	EBmixture EBmixtureObj;
	double *ipd = REAL(R_IPD);
	double *idx = REAL(R_idx);
	int ipd_len = Rf_length(R_IPD);
	int idx_len = Rf_length(R_idx);
	double rate = REAL(R_rate)[0];

	map<string, vector<double> > tmp = EBmixtureObj.subsample(ipd, idx, ipd_len, idx_len, rate);
				
	map<int, vector<double> > ipd_map = EBmixtureObj.get_ipd_map();
	map<string, vector<double> > rl;
	map<int, vector<double> >::iterator it = ipd_map.begin();
	while(it!=ipd_map.end()){
		char idx_c[256];
		sprintf(idx_c,"%d",it->first);
		string idx_s = idx_c;
		rl[idx_s] = it->second;
		it++;
	}
	return Rcpp::wrap(tmp);
	//return Rcpp::wrap(rl);
	//return R_NilValue;
}


