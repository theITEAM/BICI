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
	