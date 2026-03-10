// Checks that inputs are correctly specified

#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <algorithm>

using namespace std;

#include "input.hh"
#include "utils.hh"

/// Works out where memory is being used
void Input::profile_memory() const
{
	cout << endl << "MEMORY USAGE" << endl;
	auto sum_tot = 0.0;
	sum_tot += model_profile();
	sum_tot += species_profile();
	sum_tot += equation_profile();
	sum_tot += precalc_profile();
	
	auto f = 100000.0/(1024*total_memory());
	cout << "TOTAL MEMORY: " << (unsigned int)(sum_tot*f) << endl;
}


/// Removes anything that is no longer required
void Input::reduce_memory()
{
	switch(model.mode){
	case DATA_SIM: case DATA_SHOW: case DATA_DEL: case DATA_CLEAR: 
	case TORNADO_SETUP: case TORNADO_RESULT:
	case SCAN_SETUP: case SCAN_RESULT:
		return;
	
	default: break;
	}
	//return;
	// Removes hash tables
	
	if(false){
		model.spline.clear();
	}
	else{
		for(auto &spl : model.spline){
			spl.name.clear();
			//spl.param_ref.clear();
			spl.div.clear();
			spl.const_val.clear();
			spl.markov_eqn_ref.clear();
			spl.hash_trans_ref.off();
		}
	}

	model.constant.hash.off();
	model.hash_pop.off();
	model.hash_spline.off();

	for(auto &sp : model.species){
		sp.sim_linear_speedup.lin_form.hash_po.off();
		
		for(auto &claa : sp.cla){
			claa.hash_comp.off();
			for(auto &isl : claa.island) isl.hash_comp.off();
			//claa.hash_swap.off(); // This is used for local proposals
		}
		
		for(auto &val : sp.ind_effect){
			val.hash_markov_eqn_ref.off();
			val.hash_nm_trans_incomp_ref.off();
			val.hash_pop_ref.off();
			val.hash_nm_trans_ref.off();
		};
	
		//for(auto &val : sp.ind_eff_group){
			//val.A_matrix.hash_ind_list.off();
		//}
		
		for(auto &val : sp.fix_effect){
			val.X_vector.hash_ind_list.off();
			val.hash_markov_eqn_ref.off();
			val.hash_nm_trans_ref.off();
			val.hash_nm_trans_incomp_ref.off();
			val.hash_pop_ref.off();		
		}
		
		sp.hash_obs_eqn.off();
		sp.hash_obs_trans_eqn.off();
		sp.hash_ind.off();
		sp.hash_enter.off();
	}

	for(auto &par : model.param){
		if(!par.trace_output && !par.state_output){
			par.element.clear();
			par.element_ref.clear();
			par.weight.clear();
			for(auto &dep : par.dep){
				if(dep.index != "t") clear_dep(dep);
			}
		}
	}

	for(auto &pv : model.param_vec){
		for(auto &al : pv.affect_like){
			al.lin_form.hash_po.off();
			//al.map.clear();
		}
	}
	
	
	hash_off(model.spec_precalc);
	hash_off(model.spec_precalc_derive);
	hash_off(model.spec_precalc_sample);
	hash_off(model.spec_precalc_all);
		
	for(auto &spl : model.spec_precalc_list) hash_off(spl.spec_precalc);
	for(auto &po : model.pop) po.hash_spline_update.off();

	model.precalc_eqn.hash_off();
}


/// Turns off hash table in SpecPrecalc
void Input::hash_off(SpecPrecalc &spre) const
{
	spre.hash.off();
	spre.hash_time.off();
}


/// Clears a dependency
void Input::clear_dep(Dependency &dep)
{
	dep.list.clear();
	dep.list_out.clear();
	dep.hash_list.off();
	dep.hash_list_out.off();
}


/// Profiles memory usage for equations
double Input::equation_profile() const
{
	cout << "EQUATION" << endl;
	auto f = 100000.0/(1024*total_memory());
	auto sum_min = 5.0/f;
	
	auto sum_tot = 0.0;
	
	{	
		auto sum = 0.0;
		for(const auto &eqn : model.eqn){
			const auto &lin = eqn.linearise;
			sum += mem(lin.no_pop_calc_store);
			for(const auto &v : lin.pop_grad_calc_store) sum += mem(v);
			sum += mem(lin.factor_calc);
			sum += lin.pop_grad_precalc.size()*8;
			
			sum += sizeof(vector < vector <PopRefFromPo> >);
			for(const auto &val : lin.pop_ref_from_po){
				sum += sizeof(vector <PopRefFromPo>);
				sum += val.size()*sizeof(PopRefFromPo);
			}
		}		
		sum_tot += sum;
		if(sum > sum_min) cout << "lin: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(const auto &eqn : model.eqn){
			sum += sizeof(Equation);
			
			sum += eqn.timer.size()*sizeof(double);
			sum += mem(eqn.calcu);
			sum += sizeof(vector <Integral>);
			for(auto &in : eqn.integral){
				sum += mem(in.calc)+2*sizeof(unsigned int);
			}
			sum += eqn.param_ref.size()*sizeof(ParamRef);
			sum += eqn.derive_ref.size()*sizeof(DeriveRef);
			sum += eqn.pop_time_ref.size()*sizeof(PopTimeRef);
			sum += eqn.pop_ref.size()*sizeof(unsigned int);
			sum += eqn.source_tr_gl.size()*sizeof(unsigned int);
			sum += eqn.ind_eff_mult.size()*sizeof(unsigned int);
			sum += eqn.ind_eff_exist.size()*sizeof(bool);
			sum += eqn.fix_eff_mult.size()*sizeof(unsigned int);
			sum += sizeof(eqn.te);
			sum += sizeof(eqn.te_raw);
			sum += mem(eqn.comp_pref_convert);
		
			sum += sizeof(eqn.warn);
			sum += mem(eqn.time);
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "eqn: " << (unsigned int)(sum*f) << endl;
	}
	
	return sum_tot;
}


double Input::model_profile() const
{
	cout << "MODEL" << endl;
	auto f = 100000.0/(1024*total_memory());
	auto sum_min = 5.0/f;
	
	auto sum_tot = 0.0;
	
	{
		auto sum = 0.0; 
		for(const auto &pr : model.prior) sum += mem(pr);
		sum_tot += sum;
		if(sum > sum_min) cout << "prior: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0; 
		for(const auto &par : model.param){
			sum += sizeof(Param);
			sum += mem(par.full_name);
			sum += mem(par.name);
			for(const auto &de : par.dep) sum += mem(de);	
			sum += sizeof(SplineInfo)+mem(par.spline_info.knot_tdiv);
			sum += sizeof(SplineOut)+mem(par.spline_out.list);
			sum += par.element_ref.size()*sizeof(ElementRef);
			for(const auto &va : par.element){
				sum += sizeof(ParamElement);
				sum += mem(va.value);
				sum += va.parent.size()*sizeof(ParamRef);
				sum += va.child.size()*sizeof(ParamRef);
			}
			sum += mem(par.weight) + mem(par.reparam_eqn);
			sum += par.ieg_ref.size()*sizeof(IEGref);
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "param: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(const auto &pv : model.param_vec){
			sum += sizeof(ParamVecEle);
			sum += mem(pv.name)+mem(pv.affect_like);
			sum += mem(pv.spec_precalc_before);
			sum += mem(pv.set_param_spec_precalc);
			sum += mem(pv.spec_precalc_after);
		}
		sum_tot += sum;
		
		auto m1 = 0.0, m2 = 0.0, m3 = 0.0, m4 = 0.0, m5 = 0.0, m6 = 0.0;
		for(const auto &pv : model.param_vec){
			m1 += sizeof(ParamVecEle);
			m2 += mem(pv.name);
			m6 += mem(pv.affect_like);
			m3 += mem(pv.spec_precalc_before);
			m4 += mem(pv.set_param_spec_precalc);
			m5 += mem(pv.spec_precalc_after);
		}
		cout << m1 << " " << m2 << " " << m3 << " "<< m4 << " "<<m5 << " "<<m6 << " m" << endl;
		if(sum > sum_min) cout << "param_vec: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.param_vec_prop);
		sum_tot += sum;
		if(sum > sum_min) cout << "param_vec_prop: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.param_vec_tvreparam);
		sum_tot += sum;
		if(sum > sum_min) cout << "param_vec_tvreparam: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.pop_reparam_th);
		sum_tot += sum;
		if(sum > sum_min) cout << "pop_reparam_th: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(const auto &po : model.pop){
			sum += sizeof(Population);
			sum += mem(po.name);
			sum += mem(po.ind_eff_mult);
			sum += mem(po.fix_eff_mult);
			sum += po.term.size()*sizeof(PopulationTerm);
			sum += po.markov_eqn_ref.size()*sizeof(PopMarkovEqnRef);
			sum += po.trans_ref.size()*sizeof(PopTransRef);
			sum += mem(po.hash_spline_update);
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "pop: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = model.hash_pop.mem();
		sum_tot += sum;
		if(sum > sum_min) cout << "hash_pop: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		const auto &gen = model.genetic_data;
		auto sum = 0.0;
		sum += mem(gen.mut_rate)+mem(gen.seq_var);
		for(const auto &ob : gen.obs){
			sum += sizeof(ObsGeneticData)+mem(ob.name)+ob.snp.size()*sizeof(GenChar);
		}
	
		for(const auto &va : gen.ind_gen_obs){
			for(const auto &va2 : va){
				sum += va2.size()*sizeof(ObsGenRef);
			}
		}
		sum += mem(gen.gen_dif);
		
		sum_tot += sum;
		if(sum > sum_min) cout << "genetic_data: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(const auto &va : model.inf_cause){
			for(const auto &va2 : va){
				for(const auto &va3 : va2){
					sum += va3.size()*sizeof(InfCause);
				}
			}
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "inf_cause: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		sum += model.ieg_ref.size()*sizeof(IEGref);
		sum_tot += sum;
		if(sum > sum_min) cout << "ieg_ref: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = model.hash_spline.mem();
		for(const auto &spl : model.spline) sum += mem(spl);
		
		auto sum1 = 0.0, sum2 = 0.0, sum3 = 0.0, sum4 = 0.0, sum5 = 0.0, sum6 = 0.0, sum7 = 0.0, sum8 = 0.0; 
	
		for(const auto &spl : model.spline){ 
			sum1 += sizeof(Spline);
			sum2 += sizeof(spl.name);
			sum3 += spl.param_ref.size()*sizeof(ElementRef);
			sum4 += spl.div.size()*sizeof(SplineDiv);
			sum5 += spl.cubic_div.size()*sizeof(CubicDiv);
			sum6 += spl.const_val.size()*sizeof(double);
			sum7 += spl.markov_eqn_ref.size()*sizeof(PopMarkovEqnRef);
			sum8 += mem(spl.hash_trans_ref);
		}
		
		sum_tot += sum;
		if(sum > sum_min) cout << "spline: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.timepoint);
		sum_tot += sum;
		if(sum > sum_min) cout << "timepoint: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(const auto &der : model.derive){
			sum += sizeof(Derive);
			sum += mem(der.name);
			sum += mem(der.full_name);
			sum += mem(der.eq_raw);
			for(const auto &de : der.dep) sum += mem(de);
			for(const auto &eq : der.eq) sum += mem(eq);
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "derive: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(const auto &def : model.define){
			sum += mem(def.name);
			sum += mem(def.value);
			sum += mem(def.full_name);
			for(const auto &de : def.dep) sum += mem(de);
			sum += mem(def.swap_temp);
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "define: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(const auto &den : model.density){
			sum += sizeof(Density)+mem(den.value);
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "density: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(auto &tab : model.param_samp_store){
			sum += mem(tab);
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "param_samp_store: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		sum += mem(model.sample);
		sum_tot += sum;
		if(sum > sum_min) cout << "sample: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(const auto &ti : model.terminal_info){
			sum += sizeof(TerminalInfo);
			for(const auto &va : ti.prop_info_store) sum += mem(va);
			sum += mem(ti.av);
			sum += mem(ti.av2);
		}
		sum_tot += sum;
		cout << "terminal_info: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.constant.value)+model.constant.hash.mem();
		sum_tot += sum;
		if(sum > sum_min) cout << "constant: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.param_vec_ref);
		sum_tot += sum;
		if(sum > sum_min) cout << "param_vec_ref: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.spline_ref);
		sum_tot += sum;
		if(sum > sum_min) cout << "spline_ref: " << (unsigned int)(sum*f) << endl;
	}
	
	return sum_tot;
}
	

double Input::precalc_profile() const
{
	cout << "PRECALC" << endl;
	auto f = 100000.0/(1024*total_memory());
	auto sum_min = 5.0/f;
	
	auto sum_tot = 0.0;
	
	const auto &pe = model.precalc_eqn;
	
	{
		auto sum = 0.0;
		sum += sizeof(vector <PreCalc>);
		for(const auto &va : pe.pcalcu){
			sum += sizeof(PreCalc);
			sum += va.item.size()*sizeof(EqItem);
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "pcalcu: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.precalc_eqn.pcalcu_ref);
		sum_tot += sum;
		if(sum > sum_min) cout << "pcalcu_ref: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = model.precalc_eqn.hash_ca_mem();
		sum_tot += sum;
		if(sum > sum_min) cout << "hash_ca: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.spec_precalc);
		sum_tot += sum;
		if(sum > sum_min) cout << "spec_precalc: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.spec_precalc_derive);
		sum_tot += sum;
		if(sum > sum_min) cout << "spec_precalc_derive: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.spec_precalc_sample);
		sum_tot += sum;
		if(sum > sum_min) cout << "spec_precalc_sample: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.spec_precalc_all);
		sum_tot += sum;
		if(sum > sum_min) cout << "spec_precalc_all: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		sum += mem(model.spec_precalc_time_ref);
		sum_tot += sum;
		if(sum > sum_min) cout << "spec_precalc_time_ref: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = 0.0;
		for(const auto &sp : model.spec_precalc_list){
			sum += mem(sp.pv);
			sum += mem(sp.spec_precalc);
		}
		sum_tot += sum;
		if(sum > sum_min) cout << "spec_precalc_time: " << (unsigned int)(sum*f) << endl;
	}
	
	{
		auto sum = mem(model.precalc_init);
		sum_tot += sum;
		if(sum > sum_min) cout << "precalc_init: " << (unsigned int)(sum*f) << endl;
	}
	
	return sum_tot;
}


double Input::species_profile() const
{
	cout << "SPECIES" << endl;
	auto f = 100000.0/(1024*total_memory());
	auto sum_min = 5.0/f;
	
	auto sum_tot = 0.0;
	
	for(const auto &sp : model.species){
		cout << sp.name << endl;
		
		{
			
			for(const auto &claa : sp.cla){
				cout << "CLASSIFICATION " << claa.name << endl;
				
				{
					auto sum = 0.0;
					for(const auto &co : claa.comp){
						sum += sizeof(Compartment);
						sum += mem(co.name);
						sum += mem(co.erlang_source);
						sum += mem(co.tr_leave);
						sum += mem(co.tr_enter);
					}
					sum_tot += sum;
					if(sum > sum_min) cout << "comp: " << (unsigned int)(sum*f) << endl;
				}
				
				{
					auto sum = claa.hash_comp.mem();
					sum_tot += sum;
					if(sum > sum_min) cout << "hash_comp: " << (unsigned int)(sum*f) << endl;
				}
				
				{
					auto sum = 0.0;
					for(const auto &tr : claa.tra){
						sum += sizeof(Transition);
						sum += mem(tr.name);
						sum += mem(tr.bp);
						sum += mem(tr.dist_param);
					}
					sum_tot += sum;
					if(sum > sum_min) cout << "tra: " << (unsigned int)(sum*f) << endl;
				}
				
				{
					auto sum = 0.0;
					for(const auto &isl : claa.island){
						for(const auto &co : isl.comp){
							sum += sizeof(IslandComp);
							for(const auto &le : co.leave){
								sum += sizeof(IslandTrans);
								for(const auto &val : le.markov_eqn_ref){
									sum += val.size()*sizeof(MEref);
								}
								sum += mem(le.nm_trans_ref);
							}
						}
						sum += mem(isl.hash_comp);
					}
					sum_tot += sum;
					if(sum > sum_min) cout << "island: " << (unsigned int)(sum*f) << endl;
				}
				
				{
					auto sum = 0.0;
					for(const auto &val : claa.tr_swap_possible){
						sum += val.size()*sizeof(bool);
					}
					sum_tot += sum;
					if(sum > sum_min) cout << "tr_swap_possible: " << (unsigned int)(sum*f) << endl;
				}
				
				{
					auto sum = claa.hash_swap.mem();
					sum_tot += sum;
					if(sum > sum_min) cout << "hash_swap: " << (unsigned int)(sum*f) << endl;
				}
				
				{
					auto sum = 0.0;
					for(const auto &val : claa.swap){
						sum += sizeof(Swap);
						sum += mem(val.name);
						sum += val.start.size()*sizeof(TrSwap);
						for(const auto &val2 : val.swap){
							sum += val2.size()*sizeof(TrSwap);
						}
						sum += mem(val.swap_rep_ref);
					}
					sum_tot += sum;
					if(sum > sum_min) cout << "swap: " << (unsigned int)(sum*f) << endl;
				}
				
				{
					auto sum = 0.0;
					for(const auto &val : claa.swap_rep){
						sum += sizeof(SwapRep);
						sum += mem(val.name);
					}
					sum_tot += sum;
					if(sum > sum_min) cout << "swap_rep: " << (unsigned int)(sum*f) << endl;
				}
				cout << "END" << endl;
			}
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.ind_effect){
				sum += sizeof(IndEffect);
				sum += mem(val.name);
				sum += mem(val.markov_eqn_ref);
				sum += mem(val.hash_markov_eqn_ref);
				sum += mem(val.nm_trans_ref);
				sum += mem(val.hash_nm_trans_ref);
				sum += mem(val.nm_trans_incomp_ref);
				sum += mem(val.hash_nm_trans_incomp_ref);
				sum += mem(val.pop_ref);
				sum += mem(val.hash_pop_ref);
				sum_tot += sum;
			}
			if(sum > sum_min) cout << "ind_effect: " << (unsigned int)(sum*f) << endl;
		}
	
		{
			auto sum = 0.0;
			for(const auto &val : sp.ind_eff_group){
				sum += sizeof(IEgroup);
				sum += mem(val.name);
				for(const auto &val2 : val.list){
					sum += sizeof(IEname)+ mem(val2.name);
				}
				
				sum += mem(val.A_matrix.value);
				sum += mem(val.A_matrix.ind_list);
				sum += val.A_matrix.hash_ind_list.mem();
			
				sum += mem(val.omega_pv);
				sum += mem(val.markov_eqn_ref);
				sum += mem(val.nm_trans_ref);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "ind_effect_group: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.fix_effect){
				sum += sizeof(FixedEffect);
				sum += mem(val.name);
				sum += mem(val.X_vector.ind_list);
				sum += val.X_vector.hash_ind_list.mem();
				sum += mem(val.X_vector.value);
				sum += mem(val.markov_eqn_ref);
				sum += mem(val.hash_markov_eqn_ref);
				sum += mem(val.nm_trans_ref);
				sum += mem(val.hash_nm_trans_ref);
				sum += mem(val.nm_trans_incomp_ref);
				sum += mem(val.hash_nm_trans_incomp_ref);
				sum += mem(val.pop_ref);
				sum += mem(val.hash_pop_ref);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "fix_effect: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.comp_mult);
			sum_tot += sum;
			if(sum > sum_min) cout << "comp_mult: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.comp_gl){
				sum += sizeof(CompGlobal);
				sum += mem(val.name);
				sum += mem(val.cla_comp);
				sum += val.pop_ref.size()*sizeof(PopRef);
				sum += mem(val.pop_ref_simp);
				sum += mem(val.tr_enter);
				sum += mem(val.tr_leave);
				sum += mem(val.tr_leave_markov);
				for(const auto &val2 : val.tra_leave_group){
					sum += sizeof(CompGlTransGroup);
					sum += mem(val2.tr_list);
				}
				sum += mem(val.me_ref);
				sum += mem(val.nmtransincomp_ref);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "comp_gl: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.tra_gl){
				sum += sizeof(TransGlobal);
				sum += mem(val.name);
				sum += sizeof(TransInfection)+mem(val.infection.pop_ref);
				sum += mem(val.bp);
				sum += mem(val.bp_other_eq);
				sum += mem(val.dist_param);
				sum += mem(val.tform);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "tra_gl: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.tr_shift);
			sum_tot += sum;
			if(sum > sum_min) cout << "tr_shift: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.trg_from_tr);
			sum_tot += sum;
			if(sum > sum_min) cout << "trg_from_tr: " << (unsigned int)(sum*f) << endl;
		}
	
		{
			auto sum = mem(sp.comp_global_convert);
			sum_tot += sum;
			if(sum > sum_min) cout << "comp_global_convert: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.tra_ie_affect);
			sum_tot += sum;
			if(sum > sum_min) cout << "tra_ie_affect: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.tra_ie_affect);
			sum_tot += sum;
			if(sum > sum_min) cout << "tra_fe_affect: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.markov_eqn){
				sum += sizeof(MarkovEqn);
				sum += mem(val.source_tr_gl);
				sum += mem(val.ind_eff_mult);
				sum += mem(val.fix_eff_mult);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "markov_eqn: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.markov_tree.node){
				sum += sizeof(MarkovNode) + mem(val.child);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "markov_tree: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.cgl_begin_nm);
			sum_tot += sum;
			if(sum > sum_min) cout << "cgl_begin_nm: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.cgl_incomp_nmtrans_cl);
			sum_tot += sum;
			if(sum > sum_min) cout << "cgl_incomp_nmtrans_cl: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.cgl_tr_source);
			sum_tot += sum;
			if(sum > sum_min) cout << "cgl_tr_source: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.cgl_tr_sink);
			sum_tot += sum;
			if(sum > sum_min) cout << "cgl_tr_sink: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.individual){
				sum += sizeof(IndData);
				sum += mem(val.name);
				sum += val.ev.size()*sizeof(EventData);
				for(const auto &val2 : val.obs){
					sum += sizeof(ObsData);
					sum += mem(val2.c_obs_prob_eqn);
					sum += mem(val2.eqn_zero);
					sum += mem(val2.obs_eqn_ref);
					//sum += mem(val2.Se_eqn);
					//sum += mem(val2.Sp_eqn);
				}
				
				for(const auto &val2 : val.cl_trig_ev_ref){
					sum += val2.size()*sizeof(TrigEventRef);
				}
				for(const auto &val2 : val.fixed_trans_ev){
					sum += val2.size()*sizeof(TrigEventRef);
				}
				sum += val.trig_ev_ref.size()*sizeof(TrigEventRef);
				sum += mem(val.sample_needed);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "individual: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.enter){
				sum += sizeof(Enter);
				sum += mem(val.name);
				for(const auto &val2 : val.cla){
					sum += sizeof(EnterCla);
					sum += mem(val2.eqn);
					sum += mem(val2.obs_eqn_ref);	
				}
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "enter: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.source){
				sum += sizeof(DataSource);
				sum += mem(val.name);
				for(const auto &val2 : val.tags){
					sum += sizeof(Tag);
					sum += mem(val2.name);
					sum += mem(val2.value);
				}
				sum += mem(val.command_name);
				sum += mem(val.table);
				sum += mem(val.pop_prior);
				sum += mem(val.filter_str);
				sum += mem(val.filter_trans_str);
				sum += mem(val.obs_model);
				sum += mem(val.SNP_root);
				sum += mem(val.mut_rate_str);
				sum += mem(val.seq_var_str);
				for(const auto &val2 : val.load_col){
					sum += sizeof(LoadCol);
					sum += mem(val2.heading);
					sum += mem(val2.desc);
				}
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "source: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			const auto &ic = sp.init_cond;
			sum += sizeof(InitCond);
			sum += mem(ic.cnum);
			for(const auto &val : ic.comp_prior) sum += mem(val);
			sum += mem(ic.alpha_focal);
			sum += mem(ic.comp_reduce);
			sum += ic.comp_reduce_ref.size()*sizeof(CompRedRef);
			sum += mem(ic.pop_prior);
			sum += mem(ic.alpha);
			sum_tot += sum;
			if(sum > sum_min) cout << "init_cond: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.nm_trans){
				sum += sizeof(NMTrans);
				sum += mem(val.name);
				sum += mem(val.bp_other_eq);
				sum += mem(val.bp_all_eq);
				sum += mem(val.dist_param_eq_ref);
				sum += val.ind_fac_rate.size()*sizeof(IndFacRate);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "nm_trans: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.nm_trans_incomp){
				sum += sizeof(NMTransIncomp);
				sum += mem(val.name);
				sum += mem(val.nmtrans_ref);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "nm_trans_incomp: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.obs_eqn);
			sum_tot += sum;
			if(sum > sum_min) cout << "obs_eqn: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.hash_obs_eqn);
			sum_tot += sum;
			if(sum > sum_min) cout << "hash_obs_eqn: " << (unsigned int)(sum*f) << endl;
		}
	
		{
			auto sum = mem(sp.obs_trans_eqn_ref);
			sum_tot += sum;
			if(sum > sum_min) cout << "obs_trans_eqn_ref: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.obs_trans_eqn);
			sum_tot += sum;
			if(sum > sum_min) cout << "obs_trans_eqn: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.hash_obs_trans_eqn);
			sum_tot += sum;
			if(sum > sum_min) cout << "hash_obs_trans_eqn: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.obs_trans){
				sum += mem(val.name);
				sum += mem(val.tra_prob_eqn);
				sum += mem(val.obs_eqn_ref);
				sum += mem(val.is_one);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "obs_trans: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.comp_period);
			sum_tot += sum;
			if(sum > sum_min) cout << "comp_period: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.comp_terminal_cl);
			sum_tot += sum;
			if(sum > sum_min) cout << "comp_terminal_cl: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.sink_exist);
			sum_tot += sum;
			if(sum > sum_min) cout << "sink_exist: " << (unsigned int)(sum*f) << endl;
		}
	
		{
			auto sum = mem(sp.comp_terminal);
			sum_tot += sum;
			if(sum > sum_min) cout << "comp_terminal: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.add_rem_pop);
			sum_tot += sum;
			if(sum > sum_min) cout << "add_rem_pop: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.add_rem_pop_change);
			sum_tot += sum;
			if(sum > sum_min) cout << "add_rem_pop_change: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.pop_filter){
				sum += sizeof(PopFilter);
				sum += mem(val.comp_prob_eqn);
				sum += mem(val.eqn_zero);
				sum += mem(val.comp_obs_mod_ref);
				sum += mem(val.c_nonzero);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "pop_filter: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.pop_data){
				sum += sizeof(PopData);
				sum += mem(val.comp_obs_mod_ref);     
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "pop_data: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.pop_data_ref);
			sum_tot += sum;
			if(sum > sum_min) cout << "pop_data_ref: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.pop_trans_filter){
				sum += sizeof(PopTransFilter);
				sum += mem(val.name);
				sum += mem(val.trans_prob_eqn);
				sum += mem(val.eqn_zero);
				sum += mem(val.trans_obs_mod_ref);
				sum += mem(val.tr_nonzero);     				
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "pop_trans_filter: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.pop_trans_data){
				sum += sizeof(PopTransData);
				sum += mem(val.trans_obs_mod_ref);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "pop_trans_data: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.pop_trans_ref);
			sum_tot += sum;
			if(sum > sum_min) cout << "pop_trans_ref: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.par_event_joint){
				sum += sizeof(ParEventJointProp);
				sum += mem(val.tr_list);
			}
			sum_tot += sum;
			if(sum > sum_min) cout << "par_event_joint: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			for(const auto &val : sp.warn){
				sum += sizeof(WarnData);
				sum += mem(val.te);
			}
			if(sum > sum_min) cout << "warn: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = 0.0;
			sum += mem(sp.sim_linear_speedup.calc);
			sum += mem(sp.sim_linear_speedup.lin_form);
			if(sum > sum_min) cout << "sim_linear_speedup: " << (unsigned int)(sum*f) << endl;
		}
	
		{
			auto sum = mem(sp.tr_after);
			sum_tot += sum;
			if(sum > sum_min) cout << "tr_after: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.tr_before);
			sum_tot += sum;
			if(sum > sum_min) cout << "tr_before: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = sp.hash_ind.mem();
			sum_tot += sum;
			if(sum > sum_min) cout << "hash_ind: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = sp.hash_enter.mem();
			sum_tot += sum;
			if(sum > sum_min) cout << "hash_enter: " << (unsigned int)(sum*f) << endl;
		}
		
		{
			auto sum = mem(sp.data_warning);
			sum_tot += sum;
			if(sum > sum_min) cout << "data_warning: " << (unsigned int)(sum*f) << endl;
		}
	}
	
	return sum_tot;
}
