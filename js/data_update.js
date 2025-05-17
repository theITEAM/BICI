"use strict";
// Functions which update the data in response to changes in the model

///////// CHANGES TO CLASSIFICATION /////////////

/// Updates data based on a classification changing name
function data_update_rename_classification(p,cl,old_name,new_name)
{
	for(let loop = 0; loop < 3; loop++){
		let source;
		switch(loop){
		case 0: source = model.species[p].sim_source; break;
		case 1: source = model.species[p].inf_source; break;
		case 2:
			{		
				let pf = model.inf_res.plot_filter;
				if(pf){
					let sp = pf.species[p];
					if(sp) source = sp.ppc_source;
				}
			}
			break;
		}

		if(source){
			for(let i = 0; i < source.length; i++){
				let so = source[i];
				if(so.error != true){
					let tab = so.table;
					for(let c = 0; c < so.load_col.length; c++){
						let lc = so.load_col[c];
						let ty = lc.type;
						if((ty == "compartment" || ty == "compartment_prob") && lc.heading == old_name){
							lc.heading = new_name;
							if(tab && tab.ncol > c) tab.heading[c] = new_name;
						}
					}			
					
					switch(so.type){
					case "Init. Pop.":
						{
							if(so.spec.focal.te == old_name) so.spec.focal.te = new_name;
													
							switch(so.spec.radio.value){
							case "Graphical":		
								switch(so.spec.radio2.value){
								case "Focal": break;
								case "All":	break;
								}
								break;

							case "File":
								switch(so.spec.radio2.value){
								case "Focal":	
									break;

								case "All":
									if(tab.heading[cl] == old_name) tab.heading[cl] = new_name;
									break;
								}
								break;
							}
						}
						break;

					case "Add Pop.": case "Remove Pop.":
						if(tab.heading[cl+1] == old_name) tab.heading[cl+1] = new_name;
						break;
					
					case "Add Ind.":
						if(tab.heading[cl+2] == old_name) tab.heading[cl+2] = new_name;
						break;

					case "Move Ind.":
						if(so.spec.cl_drop.te == old_name) so.spec.cl_drop.te = new_name;
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

					case "Transition":
						if(so.spec.cl_drop.te == old_name) so.spec.cl_drop.te = new_name;
						if(so.spec.filter.te ==  old_name) so.spec.filter.te = new_name;
						break;
						
					case "Population":
						{
						}
						break;

					case "Pop. Trans.":
						{
							if(so.spec.cl_drop.te == old_name) so.spec.cl_drop.te = new_name;
							if(so.spec.filter.te == old_name) so.spec.filter.te = new_name;
						}
						break;
					}	
				}
			}
		}
	}
}


/// Updates data based on a classification being added 
function data_update_add_classification(p)
{
	for(let loop = 0; loop < 2; loop++){
		let sp = model.species[p];
	
		let cl = sp.ncla-1;
		let claa = sp.cla[cl];
		
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
				
				case "Add Pop.": case "Remove Pop.": break;
					
				case "Add Ind.": break;

				case "Population": case "Pop. Trans.": 
					let ob = {cl_name_store:claa.name, comp:[], radio: {value: "All"}};
					for(let c = 0; c < claa.ncomp; c++){
						ob.comp.push({check:false, comp_name_store:claa.comp[c].name});
					}
					so.spec.filter.cla.push(ob);
					break;
				}
			}
		}
	}
}


/// Updates data based on a classification changing name
function data_update_delete_classification(p,cl,name,comp_list)
{
	for(let loop = 0; loop < 2; loop++){
		let sp = model.species[p];
		
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
								let glob_comp = get_glob_comp(p);
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
										if(cl2 == sp.ncla) sum += Number(so.glob_comp[jj].pop);
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
				
				default:
					{
						let c = find_in(tab.heading,name);
						if(c != undefined){
							for(let r = 0; r < tab.nrow; r++){
								tab.ele[r].splice(c,1);
							}
							tab.heading.splice(c,1);
							tab.ncol--;
							so.load_col.splice(c,1);		
						}
					}
					break;
				}
				
				if(so.load_col){
					for(let k = 0; k < so.load_col.length; k++){
						let co = so.load_col[k];
						if(co.cl){
							if(co.cl > cl) co.cl--;
						}
					}
				}
			}
		}
	}
}


/////////////// CHANGES TO COMPARTMENTS ///////////////

// Updates data sources corresponding to change in compartment name
function data_update_rename_compartment(p,cl,old_name,new_name)
{
	let sp = model.species[p];
	
	for(let loop = 0; loop < 2; loop++){
		let sp = model.species[p];
	
		let claa = sp.cla[cl];

		let source;
		switch(loop){
		case 0: source = sp.sim_source; break;
		case 1: source = sp.inf_source; break;
		}

		for(let i = 0; i < source.length; i++){
			let so = source[i];
			
			let tab = so.table;
			if(so.error != true){
				if(tab != undefined){
					for(let c = 0; c < tab.ncol; c++){
						let ty = so.load_col[c].type;
						if((ty == "compartment" || ty == "compartment_prob") && tab.heading[c] == claa.name){
							for(let r = 0; r < tab.nrow; r++){
								if(tab.ele[r][c] == old_name) tab.ele[r][c] = new_name;
							}
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
					
				case "Add Pop.": case "Remove Pop.":
					{
						for(let r = 0; r < tab.nrow; r++){
							if(tab.ele[r][cl+1] == old_name) tab.ele[r][cl+1] = new_name;
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

				case "Diag. Test":
					{
					}
					break;

				case "Population":
					{
					}
					break;

				case "Pop. Trans.": case "Transition": 
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
					}
					break;
				}
			}
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
	
	let pf = model.inf_res.plot_filter;
	if(pf){
		for(let p = 0; p < pf.species.length; p++){
			let sp = pf.species[p];
			for(let i = 0; i < sp.ppc_source.length; i++){
				sp.ppc_source[i].info.p = p;
			}
		}
	}
}
			