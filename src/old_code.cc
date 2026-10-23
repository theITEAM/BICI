/*
/// Adding node to new branch (accounting for genetic observation)
void State::change_add_node_with_gobs(double t, const NodeRef &nr_add, TransTreeChange &gc, const GenObsFunc &gof, unsigned int mutBC_num, RegraftStore rs) const 
{
	const auto &inf_node = genetic_value.inf_node;
	auto mut_rate = genetic_value.mut_rate*model.details.dt;
	
	string warn;
	
	vector <unsigned int> val_samp;
	
	RegraftConnect ty;
	unsigned int num;
	double t1, t2, t3, frac, lam; 
	
	if(nr_add.node != OUTSIDE_INF && nr_add.node != ENTER_INF){
		const auto &in_add = inf_node[nr_add.node];
		const auto &inf_ev_add = in_add.inf_ev;
		auto j_add = nr_add.index;
		
		if(j_add > 0) t1 = inf_ev_add[j_add-1].tdiv;
		else t1 = in_add.tdiv_start;
		t2 = t;
		
		if(j_add == inf_ev_add.size()){                // Adds at the end
			if(rs == REGRAFT_RETRIEVE){
				for(auto k = 0u; k < N_REGRAFT_TRY; k++) val_samp.push_back(gc.rev_store[1+k]);
			}
			else{
				for(auto k = 0u; k < N_REGRAFT_TRY; k++) val_samp.push_back(poisson_sample(mut_rate*(t2-t1),warn));
			}
			
			ty = RG_ADD_END;
		}
		else{
			t3 = inf_ev_add[j_add].tdiv;
			frac = (t2-t1)/(t3-t1);
			num = inf_ev_add[j_add].mut_num;
			
			if(rs == REGRAFT_RETRIEVE){
				for(auto k = 0u; k < N_REGRAFT_TRY; k++) val_samp.push_back(gc.rev_store[1+k]);
			}
			else{
				for(auto k = 0u; k < N_REGRAFT_TRY; k++) val_samp.push_back(binomial_sample(frac,num,warn));
			}
			
			ty = RG_ADD_MID;
		}
	}
	else{
		auto seq_var = genetic_value.seq_var;
		
		lam = seq_var+mut_rate*t;
		if(rs == REGRAFT_RETRIEVE){
			for(auto k = 0u; k < N_REGRAFT_TRY; k++) val_samp.push_back(gc.rev_store[1+k]);
		}
		else{
			for(auto k = 0u; k < N_REGRAFT_TRY; k++) val_samp.push_back(poisson_sample(lam,warn));
		}
		
		ty = RG_ORIGIN;
	}
	
	auto mut_sum = gof.mut_sum;
	int mutA_num = gc.mutA_num;
			
	double pr;
	
	vector <unsigned int> mutA_st;
	
	auto sum = 0.0;
	vector <double> pr_st, sum_st;
	for(auto k = 0u; k < N_REGRAFT_TRY; k++){
		auto val = val_samp[k];
					
		auto ms = mut_sum;
		if(mutBC_num == UNSET) ms += val;
		else ms += mutBC_num-val;

		int sh = int(ran()*(2*ms+1))-(int)ms;
				cout << sh << " " << -sh << " " << mutA_num << " shh\n";
		if(-sh > (int)mutA_num){
			pr = 0;
			mutA_st.push_back(UNSET);
		}
		else{
			auto mutA = mutA_num+sh;
			mutA_st.push_back(mutA);
			
			cout << calc_genetic_obs_fast(val,mutA,gof) << " va\n";
			pr = exp(calc_genetic_obs_fast(val,mutA,gof));
		}
		
		sum += pr;
		
		pr_st.push_back(pr);
		sum_st.push_back(sum);
		cout << pr << " prob\n";
	}
	
	if(sum == 0){
		gc.type = TRANS_TREE_FAIL;
		return;
	}
	
	auto z = ran()*sum;
	auto k = 0u; while(k < N_REGRAFT_TRY && z > sum_st[k]) k++;
	if(k == N_REGRAFT_TRY) emsg("end wrong");
	
	gc.probif = log(pr_st[k])/sum;
	
	auto val = val_samp[k];
	gc.mut_num = val;
	gc.mutA_num = mutA_st[k];
		
	switch(ty){
	case RG_ADD_END:
		{
			auto prob = poisson_probability(val,mut_rate*(t2-t1));
			gc.dlike_genetic_process += prob;
			gc.probif += prob;
		}
		break;
			
	case RG_ADD_MID:
		gc.dlike_genetic_process += poisson_probability(val,mut_rate*(t2-t1)) 
																+ poisson_probability(num-val,mut_rate*(t3-t2))
																- poisson_probability(num,mut_rate*(t3-t1));                   
		gc.probif += binomial_probability(gc.mut_num,frac,num);
		break;
	
	case RG_ORIGIN:
		{
			auto prob = poisson_probability(gc.mut_num,lam);
			gc.dlike_genetic_process += prob;
			gc.probif += prob;
		}
		break;
		
	default: emsg("pro op"); break;
	}
}

			
/// Removing node from branch
void State::change_remove_node_with_gobs(double t, const NodeRef &nr_from, TransTreeChange &gc, const GenObsFunc &gof, RegraftStore rs) const 		
{	
	const auto &inf_node = genetic_value.inf_node;
	auto mut_rate = genetic_value.mut_rate*model.details.dt;

	//RegraftConnect ty;
	vector <unsigned int> val_samp;
	
	if(nr_from.node != OUTSIDE_INF && nr_from.node != ENTER_INF){ // Removing current node position
		const auto &in_from = inf_node[nr_from.node];
		const auto &inf_ev_from = in_from.inf_ev;
		auto j_from = nr_from.index;
		
		double t1; 
		if(j_from > 0) t1 = inf_ev_from[j_from-1].tdiv; 
		else t1 = in_from.tdiv_start;
		auto t2 = inf_ev_from[j_from].tdiv;
		
		if(j_from+1 == inf_ev_from.size()){            // End event
			auto num = inf_ev_from[j_from].mut_num;
			val_samp.push_back(num);
			
			//if(rs == REGRAFT_RETRIEVE){
				//for(auto k = 0u; k < N_REGRAFT_TRY; k++) val_samp.push_back(gc.rev_store[1+k]);
			//}
			//else{
			//for(auto k = 1u; k < N_REGRAFT_TRY; k++) val_samp.push_back(poisson_sample(mut_rate*(t2-t1),warn));

			if(rs == REGRAFT_STORE) gc.rev_store.push_back(num);
			
			auto prob = poisson_probability(num,mut_rate*(t2-t1));
			gc.dlike_genetic_process -= prob; 
			gc.probfi += prob;
		}
		else{
			auto t3 = inf_ev_from[j_from+1].tdiv;
			auto numA = inf_ev_from[j_from].mut_num;
			auto numB = inf_ev_from[j_from+1].mut_num;

			if(rs == REGRAFT_STORE) gc.rev_store.push_back(numA);
			if(rs == REGRAFT_STORE) gc.rev_store.push_back(gc.mutA_num);
	
			gc.dlike_genetic_process += poisson_probability(numA+numB,mut_rate*(t3-t1)) 
																- poisson_probability(numA,mut_rate*(t2-t1)) 
																- poisson_probability(numB,mut_rate*(t3-t2)); 
															
			auto frac = (t2-t1)/(t3-t1);
			gc.probfi += binomial_probability(numA,frac,numA+numB);
		}
	}
	else{
		auto seq_var = genetic_value.seq_var;
	
		auto lam = seq_var+mut_rate*t;
		auto num = genetic_value.inf_origin[nr_from.index].mut_num;
		if(rs == REGRAFT_STORE) gc.rev_store.push_back(num);
	
		auto prob = poisson_probability(num,lam);
		gc.dlike_genetic_process -= prob;
		gc.probfi += prob;
	}
	
	}
*/





	for(auto spli : spline_up_list){
			const auto &spl = model.spline[spli];
			unsigned int ind_start;
			if(ti == 0) ind_start = spl.div[0].index;
			else ind_start = spl.div[ti-1].index+1;
			
			for(auto ind = ind_start; ind <= spl.div[ti_next-1].index; ind++){
				const auto &pr = spl.param_ref[ind];
		
		for(auto spli : spline_up_list){
			const auto &spl = model.spline[spli];
			unsigned int ind_start;
			if(ti == 0) ind_start = spl.div[0].index;
			else ind_start = spl.div[ti-1].index+1;
			
			for(auto ind = ind_start; ind <= spl.div[ti_next-1].index; ind++){
				const auto &pr = spl.param_ref[ind];
				
				if(!pr.cons){
					auto th = pr.index;
					const auto &pv = model.param_vec[th];
					auto ti = pv.reparam_spl_ti;
					
					if(ti != UNSET){
						const auto &par = model.param[pv.th];
						auto eq_ref = par.get_eq_ref(pv.index);

						model.precalc_eqn.calculate(pv.spec_precalc_before,param_val,true);
									
						param_val.value_change(th);
						cout << "cal" << endl;
						value[th] = model.eqn[eq_ref].calculate_all_time(ti,popnum_t,precalc);
					
						model.precalc_eqn.calculate(pv.set_param_spec_precalc,param_val,true);
					
						const auto &spre = pv.spec_precalc_after;
						model.precalc_eqn.calculate(spre,param_val,true);
					
						if(spre.list_time.size() != 1) emsg("cannot find list time");
						
						const auto &pct = spre.list_time[0];
						if(pct.size() == 0) emsg("zero list time");
						
						auto ti_start = pct[0];
						auto ti_end = pct[pct.size()-1]+1;
						
						for(const auto &mer : spl.markov_eqn_ref){
							species[mer.p].likelihood_ib_spline_section(mer.e,ti_start,ti_end,popnum_t,like_ch);
						}
						
						for(const auto &trar : spl.trans_ref){	
							species[trar.p].likelihood_pop_spline_section(trar.tr,ti_start,ti_end,popnum_t,like_ch);
						}
					}
				}
			}
		
			spline_up_map[spli] = false;
		}
		
		/*
	// After parameters have been defined
	auto dif = false;
	if(param_list.size() > 1){
		const auto &pv = model.param_vec[param_list[0]];
		for(auto i = 1u; i < param_list.size(); i++){
			const auto &pv2 = model.param_vec[param_list[i]];
			if(equal_vec(pv.list_precalc_time,pv2.list_precalc_time) == false){ dif = true; break;}
		}
	}
	
	//cout << name << " na";
	//zz
	
	if(dif == true){  // If there is a difference in the times of the proposals then do seperately
		for(auto j : param_list){
			const auto &pv = model.param_vec[j];
			auto list_new = remove_precalc_done(pv.list_precalc,mapl);
			if(list_new.size() > 0){				
				SpecPrecalc up_pre;
				up_pre.list_precalc = list_new;
				up_pre.list_precalc_time = pv.list_precalc_time;			
				spec_precalc.push_back(up_pre);
			}
		}
		
		emsg("dif time");
	}
	else{
		const auto &pv = model.param_vec[param_list[0]];
		
		SpecPrecalc up_pre;
		up_pre.list_precalc_time = pv.list_precalc_time;			
		
		vector <unsigned int> list_precalc;
		
		auto C = param_list.size();
		vector <unsigned int> index(C,0);

		do{
			auto imin = LARGE;
			for(auto j = 0u; j < C; j++){
				const auto &pc = model.param_vec[param_list[j]].list_precalc;
				if(index[j] < pc.size()){
					auto i = pc[index[j]];
					if(i < imin) imin = i;
				}
			}
			
			if(imin == LARGE) break;
			
			if(mapl[imin] == false) list_precalc.push_back(imin);
			
			for(auto j = 0u; j < C; j++){
				const auto &pc = model.param_vec[param_list[j]].list_precalc;
				if(index[j] < pc.size() && pc[index[j]] == imin) index[j]++;
			}
		}while(true);
		
		if(false && C > 1){
			for(auto j = 0u; j < C; j++){
				const auto &pc = model.param_vec[param_list[j]].list_precalc;
				for(auto i : pc) cout << i << ",";
				cout << endl;			
			}
			
			cout << "after" << endl;
			for(auto i : list_precalc) cout << i << ",";
			cout << endl;
			emsg("comb");
		}
		
		if(list_precalc.size() > 0){
			up_pre.list_precalc = list_precalc;
			spec_precalc.push_back(up_pre);
		}
	}
	*/
	
	/*
	// For dependent parameters
	auto M = model.precalc_eqn.calcu.size();
	vector <bool> mapl(M,false);
	{
		dependent_spec_precalc.clear();
		
		for(auto k = 0u; k < dependent.size(); k++){
			auto j = dependent[k];
			const auto &pv = model.param_vec[j];
			
			SpecPrecalc up_pre;
			up_pre.list_precalc_time = pv.list_precalc_time;
			for(auto i : pv.list_precalc_before){
				if(mapl[i] == false){
					up_pre.list_precalc.push_back(i);
					mapl[i] = true;
				}
			}
			
			dependent_spec_precalc.push_back(up_pre);
		}
	}
	
	// After parameters have been defined
	auto dif = false;
	if(param_list.size() > 1){
		const auto &pv = model.param_vec[param_list[0]];
		for(auto i = 1u; i < param_list.size(); i++){
			const auto &pv2 = model.param_vec[param_list[i]];
			if(equal_vec(pv.list_precalc_time,pv2.list_precalc_time) == false){ dif = true; break;}
		}
	}
	
	//cout << name << " na";
	//zz
	
	if(dif == true){  // If there is a difference in the times of the proposals then do seperately
		for(auto j : param_list){
			const auto &pv = model.param_vec[j];
			auto list_new = remove_precalc_done(pv.list_precalc,mapl);
			if(list_new.size() > 0){				
				SpecPrecalc up_pre;
				up_pre.list_precalc = list_new;
				up_pre.list_precalc_time = pv.list_precalc_time;			
				spec_precalc.push_back(up_pre);
			}
		}
		
		emsg("dif time");
	}
	else{
		const auto &pv = model.param_vec[param_list[0]];
		
		SpecPrecalc up_pre;
		up_pre.list_precalc_time = pv.list_precalc_time;			
		
		vector <unsigned int> list_precalc;
		
		auto C = param_list.size();
		vector <unsigned int> index(C,0);

		do{
			auto imin = LARGE;
			for(auto j = 0u; j < C; j++){
				const auto &pc = model.param_vec[param_list[j]].list_precalc;
				if(index[j] < pc.size()){
					auto i = pc[index[j]];
					if(i < imin) imin = i;
				}
			}
			
			if(imin == LARGE) break;
			
			if(mapl[imin] == false) list_precalc.push_back(imin);
			
			for(auto j = 0u; j < C; j++){
				const auto &pc = model.param_vec[param_list[j]].list_precalc;
				if(index[j] < pc.size() && pc[index[j]] == imin) index[j]++;
			}
		}while(true);
		
		if(false && C > 1){
			for(auto j = 0u; j < C; j++){
				const auto &pc = model.param_vec[param_list[j]].list_precalc;
				for(auto i : pc) cout << i << ",";
				cout << endl;			
			}
			
			cout << "after" << endl;
			for(auto i : list_precalc) cout << i << ",";
			cout << endl;
			emsg("comb");
		}
		
		if(list_precalc.size() > 0){
			up_pre.list_precalc = list_precalc;
			spec_precalc.push_back(up_pre);
		}
	}
	*/

	const auto &calcu = model.precalc_eqn.calcu;
	auto M = calcu.size();
	//vector <bool> mapl_first(M,false);
	vector <bool> mapl(M,false);
		
		
			model.precalc_eqn.sv_add_map(dependent_spec_precalc,pv.
			
			SpecPrecalc up_pre;
			up_pre.list_precalc_time = pv.list_precalc_time;
			for(auto i : pv.list_precalc_before){
				const auto &ca = calcu[i];
				if(ca.time_dep){
					auto fl = false;
					for(auto ti : pv.list_precalc_time){
						if(mapl[i+ti] == false){
							mapl[i+ti] = true;
							fl = true;
						}
					}
					if(fl == true) up_pre.list_precalc.push_back(i);
				}
				else{
					if(mapl[i] == false){
						up_pre.list_precalc.push_back(i);
						mapl[i] = true;
					}
				}
			}
			
			dependent_spec_precalc.push_back(up_pre);
			
			/*
											else{
												const auto &el = par.element[ele.index];
												auto prr = el.prior_ref;
												if(prr != UNSET){
													const auto &pri = prior[prr];	
													if(pri.type == FIX_PR){
														const auto &dp = pri.dist_param[0];
														cout << dp.eq_ref << " "<< dp.te_raw << " " << par.name << " " << dp.value << " dp";
														
														
														if(dp.value != UNSET){
															item.type = NUMERIC; 
															item.num = constant.add(dp.value);
															op.push_back(item); 
															done = true;
														}
													}
												}
											}
											*/
											
											This is parameterised through its diagonals and correlation coefficients . For reasons of numerical stability, it is inadvisable for the diagonals to exceed around four , or for the correlation coefficients to go outside the range -0.9 to 0.9. If simulating, ensure the values are within this range. If performing inference, appropriate priors could be “uniform(0,4)” for diagonals and “uniform(-0.9,0.9)” for correlations.

/*		
		auto tab = load_table(dist_split);
		if(tab.error == true) return;	
		
		auto col_name = pp.dep;
		col_name.push_back("Dist");
		
		auto subtab = get_subtable(tab,col_name); if(subtab.error == true) return;
		
		auto ncol = subtab.ncol;
		
		for(auto r = 0u; r < subtab.nrow; r++){
			vector <unsigned int> ind(ncol-1);
			
			auto fl = false;
			for(auto i = 0u; i < ncol-1; i++){
				auto ele = subtab.ele[r][i];
				
				ind[i] = par.dep[i].hash_list.find(ele);
	
				if(ind[i] == UNSET){ 
					fl = true;
					
					if(par.dep[i].hash_list_out.find(ele) == UNSET){	
						alert_import("The table element '"+ele+"' is not valid (column '"+subtab.heading[i]+"', row "+tstr(r+2)+")");
						return;
					}
				}
			}
			
			if(fl == false){
				auto pri = convert_text_to_prior(subtab.ele[r][ncol-1],line_num,par.full_name,true);
			
				if(pri.error != ""){
					alert_import("The table element '"+subtab.ele[r][ncol-1]+"' is not a valid distribution specification: "+pri.error+" (col '"+subtab.heading[ncol-1]+"', row "+tstr(r+2)+").");
					return;
				}
				
				set_prior_element(par,ind,pri);
			}
		}
		*/
		
		for(auto ti = 0u; ti < T; ti++){
						for(auto c = 0u; c < ssp.cpop_st[ti].size(); c++){
							if(cpop_st_st[ti][c] != ssp.cpop_st[ti][c]){
								cout << ti << " " << c << " dif";
							}
						}
					}
					
					
	const auto &par = param[pv.th];
			switch(par.variety){
			case CONST_PARAM:
			case REPARAM_PARAM: 
			



		/*
		if(al.type == PRECALC_AFFECT){  // Combines together precalculation
			vector <unsigned int> list_precalc_new;
			auto k = 0u;
			const auto &pc = al.list_precalc;
			for(auto i : vec[i].list_precalc){
				while(k < pc.size() && pc[k] < i){
					list_precalc_new.push_back(pc[k]);
					k++;
				}
				if(k < pc.size() && pc[k] == i) k++;
				list_precalc_new.push_back(i);
			}
			
			while(k < pc.size()){
				list_precalc_new.push_back(pc[k]);
				k++;
			}
			
			if(false){
				for(auto i : vec[i].list_precalc) cout << i << ","; 
				cout << " before" << endl;
				for(auto i : pc) cout << i << ","; 
				cout << " add" << endl;
				for(auto i : list_precalc_new) cout << i << ","; 
				cout << " new" <<  endl;
				emsg("f");
			}
	
			vec[i].list_precalc = list_precalc_new;
		}
		*/
		
		/*
	// Stores spline time maps for each of the model variables
	vector < vector <bool> > spline_map;
	spline_map.resize(model.param_vec.size());
	for(auto k = 0u; k < model.param_vec.size(); k++){
		const auto &pv = model.param_vec[k];
		
		auto th = pv.th, ind = pv.index;
	
		const auto &par = model.param[th];
		
		if(par.spline_info.on == true){	
			auto nknot = par.spline_info.knot_tdiv.size();
		
			auto s = 0u;
			while(s < model.spline.size() && !(th == model.spline[s].th && ind/nknot == model.spline[s].index)) s++;
			if(s == model.spline.size()) emsg_input("Cannot find spline");
		
			const auto &spl = model.spline[s];
			
			vector <bool> map(T,false);
			for(auto ti = 0u; ti < T; ti++){ 
				if(spl.div[ti].th1 == k || spl.div[ti].th2 == k) map[ti] = true;
			}
	
			spline_map[k] = map;
			AffectLike al; al.type = SPLINE_AFFECT; al.num = s; al.num2 = UNSET; al.map = map;
			param_vec_add_affect(model.param_vec[k].affect_like,al);		
		}		
	}
	*/
	
		/*
	/// The effect of parameters on Markov equations
	for(auto p = 0u; p < model.species.size(); p++){
		auto &sp =  model.species[p];
		switch(sp.type){
		case INDIVIDUAL:
			{
				for(auto i = 0u; i < sp.markov_eqn.size(); i++){
					const auto &me = sp.markov_eqn[i];
				
					const auto &eqn = model.eqn[me.eqn_ref];
				
					for(auto &pref : eqn.param_ref){
						auto th = pref.th, ind = pref.index;
						auto &par = model.param[th];
					
						if(par.variety != CONST_PARAM){
							auto k = par.get_param_vec(ind);	
							if(k != UNSET){
								if(par.spline_info.on == true){	
									auto nknot = par.spline_info.knot_tdiv.size();
									
									for(auto t = 0u; t < nknot; t++){ // Goes down the spline
										k = par.get_param_vec(ind+t);
										AffectLike al; al.map = spline_map[k];
										al.type = DIV_VALUE_AFFECT; al.num = p; al.num2 = i;
										param_vec_add_affect(model.param_vec[k].affect_like,al);		
								
										al.type = MARKOV_LIKE_AFFECT;
										param_vec_add_affect(model.param_vec[k].affect_like,al);
									}
								}
								else{
									AffectLike al; 
									if(me.time_vari == true) al.map.resize(T,true);
									else al.map.resize(1,true);
						
									al.type = DIV_VALUE_AFFECT; al.num = p; al.num2 = i;
									param_vec_add_affect(model.param_vec[k].affect_like,al);
								
									al.type = MARKOV_LIKE_AFFECT;
									param_vec_add_affect(model.param_vec[k].affect_like,al);
								}
							}
						}
					}
				}
			}
			break;
			
		case POPULATION:
			{
				for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
					const auto &tra = sp.tra_gl[tr];
					const auto &eqn = model.eqn[tra.dist_param[0].eq_ref];
					for(auto &pref : eqn.param_ref){
						auto th = pref.th, ind = pref.index;
					
						auto &par = model.param[th];
						if(par.variety != CONST_PARAM){
							auto k = par.get_param_vec(ind);	
							if(k != UNSET){
								if(par.spline_info.on == true){
									auto nknot = par.spline_info.knot_tdiv.size();
								
									for(auto t = 0u; t < nknot; t++){ // Goes down the spline
										k = par.get_param_vec(ind+t);
										AffectLike al; al.map = spline_map[k];
										al.type = MARKOV_POP_AFFECT; al.num = p; al.num2 = tr; 
										param_vec_add_affect(model.param_vec[k].affect_like,al);
									}
								}
								else{
									AffectLike al; al.map.resize(T,true);
									al.type = MARKOV_POP_AFFECT; al.num = p; al.num2 = tr;
									param_vec_add_affect(model.param_vec[k].affect_like,al);
								}
							}
						}
					}
				}
			}
			break;
		}
	}
	*/
	
	
	/*
	cout << model.npopcomb << " "<< model.popcomb.size() << " jjj"; 
	
	for(auto &der : model.derive){
		for(auto eqi : der.eq){
			const auto &eq = model.eqn[eqi.eq_ref];
 			cout << eq.te_raw << " DERIVEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE\n";
			eq.print_calculation();
			
			auto imax = eq.integral.size(); if(imax > 2) imax = 2;
			for(auto i = 0u; i < imax; i++){
				auto it = eq.integral[i];
				eq.print_calc("integral",it.calc);
			}
			cout << "\n\n";
		}
		
	}
	emsg("derive print");
	*/
	
	/*
			for(auto ti = 0u; ti < ti_end; ti++){
				vector <double> tnum(sp.tra_gl.size());
				for(auto tr = 0u; tr < sp.tra_gl.size(); tr++){
					tnum[tr] = trans_num[tr][ti];
				}
			
				cpop_st.push_back(cpop);
				update_cpop(ti,cpop,tnum);
			}
			cpop_st.push_back(cpop);
			*/