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
	