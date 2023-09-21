"use strict";

// Contains functions which update the data in response to changes in the model


///////// CHANGES TO CLASSIFICATION /////////////

/// Updates data based on a classification changing name
function data_update_rename_classification(p,cl,old_name,new_name)
{
	let sp = model.species[p];
	
	for(let loop = 0; loop < 2; loop++){
		let source;
		switch(loop){
		case 0: source = sp.sim_source; break;
		case 1: source = sp.inf_source; break;
		}

		for(let i = 0; i < source.length; i++){
			let so = source[i];
			if(so.error != true){
				let tab = so.table;
				
				for(let c = 0; c < so.load_col.length; c++){
					let lc = so.load_col[c];
					let ty = lc.type;
					if((ty == "compartment" || ty == "multi_comp" || ty == "compartment_prob") &&  lc.heading == old_name){
						lc.heading = new_name;
						if(tab.ncol > c) tab.heading[c] = new_name;
					}
				}			
							
				switch(so.type){
				case "Init. Pop.":
					{
						switch(so.spec.radio.value){
						case "Graphical":
							switch(so.spec.radio2.value){
							case "Focal": 
								//if(so.cla[cl].name_store == old_name) so.cla[cl].name_store = new_name;
								break;

							case "All":	break;
							}
							break;

						case "File":
							switch(so.spec.radio2.value){
							case "Focal":	
								if(so.spec.focal.te == old_name) so.spec.focal.te = new_name;
								break;

							case "All":
								if(tab.heading[cl] == old_name) tab.heading[cl] = new_name;
								break;
							}
							break;
						}
					}
					break;

				case "Add Ind.":
					if(tab.heading[cl+2] == old_name) tab.heading[cl+2] = new_name;
					break;

				case "Move Ind.":
					if(so.spec.cl_drop.te == old_name) so.spec.cl_drop.te = new_name;
					break;

				case "Init. Pop. Prior":
					switch(so.spec.radio.value){
					case "Flat":
						break;

					case "Dirichlet":
						if(tab.heading[cl] == old_name) tab.heading[cl] = new_name;
						break;
					}
					break;

				case "Compartment":
					if(so.spec.cl_drop.te == old_name){
						so.spec.cl_drop.te = new_name;
						if(tab.heading[2] == old_name) tab.heading[2] = new_name;
					}
					break;

				case "Diag. Test":
					{
						let cb = so.spec.check_box;
						if(cb.name == old_name) cb.name = new_name;
						if(so.spec.cl_drop.te == old_name) so.spec.cl_drop.te = new_name;
					}
					break;

				case "Population":
					{
						//let clz = so.spec.filter.cla[cl];
						//if(clz.name_store == old_name) clz.name_store = new_name;
					}
					break;

				case "Pop. Trans.":
					{
						if(so.spec.cl_drop.te == old_name) so.spec.cl_drop.te = new_name;
						if(so.spec.filter.te == old_name) so.spec.filter.te = new_name;
					}
					break;

				case "Set Traps":
					{

					}
					break;
				}	
			}
			add_source_description(so);
		}
	}
}


/// Updates data based on a classification being added 
function data_update_add_classification(p)
{
	let sp = model.species[p];
	let cl = sp.ncla-1;
	let claa = sp.cla[cl];

	for(let loop = 0; loop < 2; loop++){
		let source;
		switch(loop){
		case 0: source = sp.sim_source; break;
		case 1: source = sp.inf_source; break;
		}

		for(let i = 0; i < source.length; i++){
			let so = source[i];
			if(so.error != true){
				let tab = so.table;
						
				switch(so.type){
				case "Init. Pop.":
					{
						switch(so.spec.radio.value){
						case "Graphical":
							switch(so.spec.radio2.value){
							case "Focal": break;
							case "All": break;
							}
							break;

						case "File":
							switch(so.spec.radio2.value){
							case "Focal": break;
							case "All": break;
							}
							break;
						}
					}
					break;
							
				case "Add Ind.":break;

				case "Population": case "Pop. Trans.": case "Set Traps":
					let ob = {cl_name_store:claa.name, comp:[], radio: {value: "All"}};
					for(let c = 0; c < claa.ncomp; c++){
						ob.comp.push({check:false, comp_name_store:claa.comp[c].name});
					}
					so.spec.filter.cla.push(ob);
					break;
				}
			}
			add_source_description(so);
		}
	}
}


/// Updates data based on a classification changing name
function data_update_delete_classification(p,cl,name,comp_list)
{
	let sp = model.species[p];
	
	for(let loop = 0; loop < 2; loop++){
		let source;
		switch(loop){
		case 0: source = sp.sim_source; break;
		case 1: source = sp.inf_source; break;
		}

		for(let i = 0; i < source.length; i++){
			let so = source[i];
			
			if(so.error != true){
				let tab = so.table;
						
				switch(so.type){
				case "Init. Pop.":
					{
						switch(so.spec.radio.value){
						case "Graphical":
							switch(so.spec.radio2.value){
							case "Focal": 
								let focal = so.spec.focal.te;
								let cl2 = find(sp.cla,"name",focal);
								if(cl2 != undefined) so.cla.splice(cl,1);
								break;

							case "All":
								let glob_comp = model.get_glob_comp(p);
								for(let j = 0; j < glob_comp.length; j++){
									let gc = glob_comp[j].cla;

									let sum = 0;
									for(let jj = 0; jj < so.glob_comp.length; jj++){
										let gc2 = so.glob_comp[jj].cla;

										let cl2;
										for(cl2 = 0; cl2 < sp.ncla; cl2++){
											let cl3 = cl2; if(cl3 >= cl) cl3++;
											if(gc[cl2] != gc2[cl3]) break;
										}
										if(cl2 == sp.ncla) sum += so.glob_comp[jj].pop;
									}
									glob_comp[j].pop = sum;
								}
								so.glob_comp = glob_comp;
								break;
							}
							break;

						case "File":
							switch(so.spec.radio2.value){
							case "Focal":
								let focal = so.spec.focal.te;
								let cl2 = find(sp.cla,"name",focal);
								if(cl2 != cl){
									let r = 0;
									while(r < tab.nrow){
										if(find_in(comp_list,tab.ele[r][0]) != undefined){
											tab.ele.splice(r,1);
											tab.nrow--;
										}
										else r++;
									}
								}
								break;

							case "All":
								let ncla = sp.ncla+1;
								for(let r = 0; r < tab.nrow; r++){
									let rr = r+1;
									while(rr < tab.nrow){
										let cl2 = 0; 
										while(cl2 < ncla && 
											(tab.ele[r][cl2] == tab.ele[rr][cl2] || cl2 == cl)) cl2++;

										if(cl2 == ncla){
											tab.ele[r][ncla] = Number(tab.ele[r][ncla])+Number(tab.ele[rr][ncla]);
											tab.ele.splice(rr,1);
											tab.nrow--;
										}
										else rr++;
									}
									tab.ele[r].splice(cl,1);
								}
								if(tab.col_used != undefined) tab.col_used.splice(cl,1);
								tab.heading.splice(cl,1);
								tab.ncol--;
								break;
							}
							break;
						}
					}
					break;
							
				case "Add Ind.":	
					{
						for(let r = 0; r < tab.nrow; r++){
							tab.ele[r].splice(cl+2,1);
						}
						if(tab.col_used != undefined) tab.col_used.splice(cl+2,1);
						tab.heading.splice(cl+2,1);
						tab.ncol--;
					}
					break;

				case "Population": case "Pop. Trans.": case "Set Traps":
					so.spec.filter.cla.splice(cl,1);
					break;
				}
			}
			add_source_description(so);
		}
	}
}


/////////////// CHANGES TO TRANSITIONS ////////////////

/// Updates data sources when a transition is deleted
function data_update_add_transition(p,cl)
{
	let sp = model.species[p];
	let claa = sp.cla[cl];

	for(let loop = 0; loop < 2; loop++){
		let source;
		switch(loop){
		case 0: source = sp.sim_source; break;
		case 1: source = sp.inf_source; break;
		}

		for(let i = 0; i < source.length; i++){
			let so = source[i];
			if(so.error != true){
				let tab = so.table;
							
				switch(so.type){
				case "Pop. Trans.":
					{
						let cl_sel = find(sp.cla,"name",so.spec.cl_drop.te);
						if(cl == cl_sel){
							let ob = {tra_name_store:claa.tra[claa.ntra-1].name, check:false}
							so.spec.filter.tra.push(ob);
						}
					}
					break;
				}
			}
			add_source_description(so);
		}
	}
}


/// Updates data sources when a transition is deleted
function data_update_delete_transition(p,cl,ind)
{
	let sp = model.species[p];
	let claa = sp.cla[cl];

	for(let loop = 0; loop < 2; loop++){
		let source;
		switch(loop){
		case 0: source = sp.sim_source; break;
		case 1: source = sp.inf_source; break;
		}

		for(let i = 0; i < source.length; i++){
			let so = source[i];
			if(so.error != true){
				let tab = so.table;
							
				switch(so.type){
				case "Pop. Trans.":
					{
						let cl_sel = find(sp.cla,"name",so.spec.cl_drop.te);
						if(cl == cl_sel){
							so.spec.filter.tra.splice(ind,1);
						}
					}
					break;
				}
			}
			add_source_description(so);
		}
	}
}


/////////////// CHANGES TO COMPARTMENTS ///////////////

// Updates data sources corresponding to change in compartment name
function data_update_rename_compartment(p,cl,old_name,new_name)
{
	let sp = model.species[p];
	let claa = sp.cla[cl];

	for(let loop = 0; loop < 2; loop++){
		let source;
		switch(loop){
		case 0: source = sp.sim_source; break;
		case 1: source = sp.inf_source; break;
		}

		for(let i = 0; i < source.length; i++){
			let so = source[i];
			if(so.error != true){
				let tab = so.table;
							
				for(let c = 0; c < tab.ncol; c++){
					let ty = so.load_col[c].type;
					if((ty == "compartment" || ty == "multi_comp" || ty == "compartment_prob") && tab.heading[c] == claa.name){
						for(let r = 0; r < tab.nrow; r++){
							if(tab.ele[r][c] == old_name) tab.ele[r][c] = new_name;
						}
					}
				}
						
				switch(so.type){
				case "Init. Pop.":
					{
						switch(so.spec.radio.value){
						case "Graphical":
							switch(so.spec.radio2.value){
							case "Focal":
							/*
								let clz = so.cla[cl];
								for(let c = 0; c < claa.ncomp; c++){
									if(clz.comp_init_pop[c].name_store == old_name){
										clz.comp_init_pop[c].name_store = new_name;
									}	
								}
								*/
								break;

							case "All":	
								for(let j = 0; j < so.glob_comp.length; j++){
									if(so.glob_comp[j].cla[cl] == old_name){
										so.glob_comp[j].cla[cl] = new_name;
									}
								}
								break;
							}
							break;

						case "File":
							switch(so.spec.radio2.value){
							case "Focal":	
								for(let r = 0; r < tab.nrow; r++){
									if(tab.ele[r][0] == old_name) tab.ele[r][0] = new_name;
								}
								break;

							case "All":	
								for(let r = 0; r < tab.nrow; r++){
									if(tab.ele[r][cl] == old_name) tab.ele[r][cl] = new_name;
								}
								break;
							}
							break;
						}
					}
					break;

				case "Add Ind.":
					{
						for(let r = 0; r < tab.nrow; r++){
							if(tab.ele[r][cl+2] == old_name) tab.ele[r][cl+2] = new_name;
						}
					}
					break;

				case "Move Ind.":
					{
						let cl_sel = find(model.species[p].cla,"name",so.spec.cl_drop.te);
						if(cl_sel == cl){
							for(let r = 0; r < tab.nrow; r++){
								if(tab.ele[r][2] == old_name) tab.ele[r][2] = new_name;
								if(tab.ele[r][3] == old_name) tab.ele[r][3] = new_name;
							}
						}
					}
					break;

				case "Init. Pop. Prior":
					{
						switch(so.spec.radio.value){
						case "Flat":
							break;
		
						case "Dirichlet":
							for(let r = 0; r < tab.nrow; r++){
								if(tab.ele[r][cl] == old_name) tab.ele[r][cl] = new_name;
							}
							break;
						}
					}
					break;

				case "Transition": case "Source": case "Sink":
					{
						let tr_name = so.spec.tr_drop.te;
						so.spec.tr_drop.te = convert_transition_name(tr_name,old_name,new_name);
					}
					break;

				case "Diag. Test":
					{

						/*
						let cb = so.spec.check_box;
						let cl_sel = find(sp.cla,"name",cb.name);
						if(cl_sel == cl){
							for(let c = 0; c < claa.ncomp; c++){
								if(cb.value[c].name_store == old_name) cb.value[c].name_store = new_name;
							}
						}
						*/
					}
					break;

				case "Population":
					{
					}
					break;

				case "Pop. Trans.":
					{
						let cl_sel = find(sp.cla,"name",so.spec.cl_drop.te);
						let filt = so.spec.filter;
						if(cl_sel == cl){
							for(let j = 0; j < claa.ntra; j++){
								let tr_name = filt.tra[j].tra_name_store;
								tr_name = convert_transition_name(tr_name,old_name,new_name);
								filt.tra[j].tra_name_store = tr_name;
							}
						}
						//else{
							//let clz = filt.cla[cl];
							//for(let c = 0; c < claa.ncomp; c++){
								//if(clz.comp[c].name_store == old_name) clz.comp[c].name_store = new_name;
							//}
						//}		
					}
					break;
		
				case "Set Traps":
					{
						/*
						let clz = so.spec.filter.cla[cl];
						for(let c = 0; c < claa.ncomp; c++){
							if(clz.comp[c].name_store == old_name) clz.comp[c].name_store = new_name;
						}
						*/
					}
					break;
				}
			}
			add_source_description(so);
		}
	}
}


/// Updates data sources when a compartment is added
function data_update_add_compartment(p,cl)
{
	let sp = model.species[p];
	let claa = sp.cla[cl];
	let name = claa.comp[claa.ncomp-1].name;

	for(let loop = 0; loop < 2; loop++){
		let source;
		switch(loop){
		case 0: source = sp.sim_source; break;
		case 1: source = sp.inf_source; break;
		}

		for(let i = 0; i < source.length; i++){
			let so = source[i];
			if(so.error != true){
				let tab = so.table;
							
				switch(so.type){
				case "Init. Pop.":
					switch(so.spec.radio.value){
					case "Graphical":
						switch(so.spec.radio2.value){
						case "Focal":
							let clz = so.cla[cl];
							clz.comp_init_pop.push({comp_name_store:name, pop:0, pop_per:0});
							break;

						case "All":
							{
								let glob_comp = model.get_glob_comp(p);
								for(let k = 0; k < glob_comp.length; k++){
									let gc = glob_comp[k].cla;

									let j;
									for(j = 0; j < so.glob_comp.length; j++){
										let gc2 = so.glob_comp[j].cla;

										let cl2 = 0; while(cl2 < sp.ncla && gc2[cl2] == gc[cl2]) cl2++;
										if(cl2 == sp.ncla) break;
									}
									if(j < so.glob_comp.length){
										glob_comp[k].pop = so.glob_comp[j].pop;
									}
									else{
										glob_comp[k].pop = 0;
									}
								}	
								
								so.glob_comp = glob_comp;
							}
							break;
						}
						break;

					case "File":
						switch(so.spec.radio2.value){
						case "Focal":
							{
								let cl_sel = find(sp.cla,"name",so.spec.focal.te);
		
								let r = tab.nrow;
								tab.ele[r] = [];
								tab.ele[r][0] = name;

								let val = "0"; if(cl_sel != cl) val += "%";
								tab.ele[r][1] = val;
								tab.nrow++;
							}
							break;

						case "All":
							{
								let glob_comp = model.get_glob_comp(p);
								for(let r = 0; r < glob_comp.length; r++){
									if(glob_comp[r].cla[cl] == name){
										let rr = tab.nrow;
										tab.ele[rr] = [];
										for(let cl2 = 0; cl2 < sp.ncla; cl2++){
											tab.ele[rr][cl2] = glob_comp[r].cla[cl2];
										}
										tab.ele[rr][sp.ncla] = 0;
										tab.nrow++;
									}
								}
							}
							break;
						}
						break;
					}
					break;

				case "Init. Pop. Prior":
					switch(so.spec.radio.value){
					case "Flat":
						break;

					case "Dirichlet":
						break;
					}
					break;

				case "Diag. Test":
					{
						let cb = so.spec.check_box;
						let cl_sel = find(sp.cla,"name",cb.name);
						if(cl_sel == cl){
							cb.value.push({comp_name_store:name, check:false});
						}
					}
					break;

				case "Population":
					so.spec.filter.cla[cl].comp.push({check:false, comp_name_store:name});
					break;

				case "Set Traps":
					so.spec.filter.cla[cl].comp.push({check:false, comp_name_store:name});
					break;
		
				case "Pop. Trans.":
					let cl_sel = find(sp.cla,"name",so.spec.cl_drop.te);
					if(cl != cl_sel){
						so.spec.filter.cla[cl].comp.push({check:false, comp_name_store:name});
					}
					break;
				}
			}
			add_source_description(so);
		}
	}
}


/// Updates data sources when a compartment is deleted
function data_update_delete_compartment(p,cl,ind,name)
{
	let sp = model.species[p];
	let claa = sp.cla[cl];

	for(let loop = 0; loop < 2; loop++){
		let source;
		switch(loop){
		case 0: source = sp.sim_source; break;
		case 1: source = sp.inf_source; break;
		}

		for(let i = 0; i < source.length; i++){
			let so = source[i];
			if(so.error != true){
				let tab = so.table;
							
				switch(so.type){
				case "Init. Pop.":
					switch(so.spec.radio.value){
					case "Graphical":
						switch(so.spec.radio2.value){
						case "Focal":
							let clz = so.cla[cl];
							clz.comp_init_pop.splice(ind,1);
							break;

						case "All":
							{
								let j = 0;
								while(j < so.glob_comp.length){
									if(so.glob_comp[j].cla[cl] == name) so.glob_comp.splice(j,1);
									else j++;
								}	
							}
							break;
						}
						break;

					case "File":
						switch(so.spec.radio2.value){
						case "Focal":
							{
								let r = 0;
								while(r < tab.nrow){
									if(tab.ele[r][0] == name){ 
										tab.ele.splice(r,1);
										tab.nrow--;
									}
									else r++;
								}
							}
							break;

						case "All":
							{
								let r = 0;
								while(r < tab.nrow){
									if(tab.ele[r][cl] == name){ 
										tab.ele.splice(r,1);
										tab.nrow--;
									}
									else r++;
								}
							}
							break;
						}
						break;
					}
					break;

				case "Init. Pop. Prior":
					switch(so.spec.radio.value){
					case "Flat":
						break;

					case "Dirichlet":
						{
							let r = 0;
							while(r < tab.nrow){
								if(tab.ele[r][cl] == name){ 
									tab.ele.splice(r,1);
									tab.nrow--;
								}
								else r++;
							}
						}
						break;
					}
					break;

				case "Diag. Test":
					{
						let cb = so.spec.check_box;
						let cl_sel = find(sp.cla,"name",cb.name);
						if(cl_sel == cl){
							cb.value.splice(ind,1);
						}
					}
					break;

				case "Population":
					so.spec.filter.cla[cl].comp.splice(ind,1);
					break;

				case "Set Traps":
					so.spec.filter.cla[cl].comp.splice(ind,1);
					break;
		
				case "Pop. Trans.":
					let cl_sel = find(sp.cla,"name",so.spec.cl_drop.te);
					if(cl != cl_sel){
						so.spec.filter.cla[cl].comp.splice(ind,1);
					}
					break;
				}
			}

			add_source_description(so);
		}
	}
}


/// Based on a compartment changing name
function convert_transition_name(tr_name,old_name,new_name)
{
	let spl = tr_name.split("→");
	switch(spl.length){
	case 1:
		if(tr_name == old_name) return new_name;
		break;
	
	case 2:
		if(spl[0] == old_name) return new_name+"→"+spl[1];
		if(spl[1] == old_name) return spl[0]+"→"+new_name;
		break;

	default: error("Error option not recognised"); break;
	} 

	return tr_name;
}


/// Ensures that the data sources have the correct value for p
function reset_info_p()
{
	for(let p = 0; p < model.species.length; p++){
		let sp = model.species[p];
		for(let i = 0; i < sp.sim_source.length; i++){
			sp.sim_source[i].info.p = p;
		}
		for(let i = 0; i < sp.inf_source.length; i++){
			sp.inf_source[i].info.p = p;
		}
	}
}

			