"use strict";
// Functions for creating a SIRE file (used only for diagnostics)

/// Creates a .sire file from BICI 
function create_sire()
{
	let dt = "";
	let col = [];
	let groups = [];
		
	let par_list= [];
		
	par_list.push({name:"beta"});
	
	let name_list = model.get_all_data_individual(0);
	let N = name_list.length;
	
	let claa = model.species[0].cla[0];
	
	let tr_inf = claa.tra[0].value.rate_eqn;

	{
		let hash = new Hash();
			
		let list = [];
		
		for(let j = 0; j < N; j++){
			list.push(name_list[j].name);
			hash.add(name_list[j].name,j);
		}
		
		col.push({title:"id", data:list});
		
		let source = model.species[0].inf_source;
	
		let temp = []; for(let j = 0; j < N; j++) temp[j] = "NA";
		
		for(let j = 0; j < source.length; j++){
			let so = source[j];
			let tab = so.table;
			switch(so.type){
			case "Add Ind.":
				{
					let list = copy(temp);
					let list2 = copy(temp);
					for(let k = 0; k < tab.nrow; k++){
						let j = hash.find(tab.ele[k][0]);
						if(j == undefined) pr("Problem getting individual");
						list[j] = tab.ele[k][2];
						list2[j] = tab.ele[k][3];
					}
					
					col.push({title:"initial_comp", data:list});
					col.push({title:"group", data:list2});
				}
				break;
			}
		}
			
		for(let j = 0; j < source.length; j++){
			let so = source[j];
			let tab = so.table;
			
			switch(so.type){
			case "Add Ind.":
				break;
				
			case "Transition":
				{		
					let claa = model.species[0].cla[0];
					
					let title;
					for(let k = 0; k < claa.ntra; k++){
						if(so.spec.filter.tra[k].check == true){
							claa.tra[k].obs = true;
							title = claa.tra[k].name;
							break;
						}
					}
					
					let list = copy(temp);
					let c = 0; while(c < col.length && col[c].title != "group") c++;
					if(c == col.length) pr("Cannto find group:");
		
					for(let k = 0; k < N; k++){
						if(col[c].data[k] != ".") list[k] = "no"; 
					}
					
					for(let k = 0; k < tab.nrow; k++){
						let j = hash.find(tab.ele[k][0]);
						if(j == undefined) pr("Problem getting individual");
						list[j] = tab.ele[k][1];
					}
					
					col.push({title:title, data:list});
				}
				break;
				
			case "Ind. Eff.":
				{
					let na = so.spec.drop.te;
					let av = 0.0;
					for(let k = 0; k < tab.nrow; k++){
						av += Math.log(Number(tab.ele[k][1]));
					}
					av /= tab.nrow;
					
					let list = copy(temp);
					for(let k = 0; k < tab.nrow; k++){
						let j = hash.find(tab.ele[k][0]);
						if(j == undefined) pr("Problem getting individual");
						list[j] = Math.log(Number(tab.ele[k][1]))-av;
					}
					col.push({title:na, data:list});
				}
				break;
				
			case "Ind. Group":
				{
					let st = "";
					for(let k = 0; k < tab.nrow; k++){
						if(k != 0) st += ",";
						st += tab.ele[k][0];
					}
					groups.push({name:so.spec.gname, te:st});
				}
				break; 
			}
		}
		
		pr("fix");
		pr(model.species[0].fix_eff);
		
		// Data from fixed effects
		for(let i = 0; i < model.species[0].fix_eff.length; i++){
			let fe = model.species[0].fix_eff[i];
			
			par_list.push({name:fe.name});
			
			let list = [];
			
			for(let j = 0; j < fe.X_vector.ind_list.length; j++){
				let ii = hash.find(fe.X_vector.ind_list[j]);
				if(ii != undefined) list[ii] = Number(fe.X_vector.X_value[j]);
			}
			for(let j = 0; j < N; j++){
				if(list[j] == undefined) pr("X problem");
			}
			
			col.push({title:fe.name, data:list});
		}	
				
		for(let c = 0; c < col.length; c++){
			let co = col[c];
			if(c != 0) dt += ",";
			dt += co.title;
		}
		dt += "\n";
		
		for(let j = 0; j < N; j++){
			for(let c = 0; c < col.length; c++){
				let co = col[c];
				if(c != 0) dt += ",";
				dt += co.data[j];
			}
			dt += "\n";
		}
	}
	
	let te = '<?xml version="1.0" encoding="UTF-8"?>\n';
	te += '<SIRE version="2.2">\n';

	te += '<!--MCMC options-->\n'

	//	te += '<pas output_dir="out" phi_final="1" nsample_per_gen="50" nsample="2000" burnin="400" thin="2" sample_states="0" ie_output="true"/>\n';

	let nsamp = model.inf_details.sample;
	te += '<mcmc output_dir="out" nsample="'+nsamp+'" burnin="'+Math.floor(nsamp/5)+'" thin="2" ie_output="true"/>\n';

	te += '<!--Model compartments-->\n';
	
	for(let c = 0; c < claa.comp.length; c++){
		let na = claa.comp[c].name;
		te += '<comp name="'+na+'"';
		if(na != "S" && na != "E" && na != "R") te += ' relative_infectivity="1"';
		te += '/>\n';
	}
	
	te += '<!--Model transitions-->\n';
	
	for(let j = 0; j < claa.tra.length; j++){
		let tr = claa.tra[j];
		pr("tr");
		pr(tr);
		
		te += '<trans ';
		te += 'from="'+claa.comp[tr.i].name+'" ';
		te += 'to="'+claa.comp[tr.f].name+'" ';
		//if(tr.type == "exp(rate)"){
		if(j == 0){	
			te += 'type="infection" beta="beta" inf_model="density dependent" ';
		}
		else{
			let pa = tr.value.mean_eqn.param;
			pr("TR");
			pr(tr.value);
			
			if(pa.length != 1) pr("ERROR shoud be one");
			let na = pa[0].name;
			
			te += 'type="exp" mean="'+na+'" ';
			par_list.push({name:na});
		}
		
		if(j == 0){
			let ind_eff = tr.value.rate_eqn.ind_eff;
			if(ind_eff.length > 0){
				te += 'individual_effect="';
				for(let k = 0; k < ind_eff.length; k++){
					let na = ind_eff[k].name;
					if(na.substr(0,1) != "f"){
						if(k != 0) te += ',';
						te += na;
					}
				}
				te += '" ';
			}
			
			let fix_eff = tr.value.rate_eqn.fix_eff;
			if(fix_eff.length > 0){
				te += 'fixed_effect="';
				for(let k = 0; k < fix_eff.length; k++){
					let na = fix_eff[k].name;
					if(na.substr(0,1) != "f"){
						if(k != 0) te += ',';
						te += na;
					}
				}
				te += '" ';
			}
		}
		else{
			let ind_eff = tr.value.mean_eqn.ind_eff;
			if(ind_eff.length > 0){
				te += 'individual_effect="';
				for(let k = 0; k < ind_eff.length; k++){
					let na = ind_eff[k].name;
					if(k != 0) te += ',';
					te += na;
				}
				te += '" ';
			}
			
			let fix_eff = tr.value.rate_eqn.fix_eff;
			if(fix_eff.length > 0){
				te += 'fixed_effect="';
				for(let k = 0; k < fix_eff.length; k++){
					let na = fix_eff[k].name;
					if(k != 0) te += ',';
					te += na;
				}
				te += '" ';
			}
		}
		
		if(tr.obs == true)  te += ' data_column="'+tr.name+'"';
		te += '/>\n';
	}

	te +=  '<!--Add FEs or IEs to infectivity here-->\n';
	
	{
		for(let i = 0; i < tr_inf.ind_eff.length; i++){
			let na = tr_inf.ind_eff[i].name;
		}
		
		te += '<infectivity individual_effect="';

		let te2 = "";
		for(let i = 0; i < tr_inf.ind_eff.length; i++){
			let na = tr_inf.ind_eff[i].name;
			if(na.substr(0,1) == "f"){
				if(te2 != "") te2 += ",";
				te2 += na;
			}	
		}
		te += te2;
		te += '" ';
	
		te += 'fixed_effect="';
		for(let i = 0; i < tr_inf.fix_eff.length; i++){
			let na = tr_inf.fix_eff[i].name;
			if(na.substr(0,1) == "f"){
				te += na;
			}	
		}
		te += '"/>\n';
	}
	
	{
		let k = 0; while(k < tr_inf.param.length &&  tr_inf.param[k].name != "G") k++;
		if(k < tr_inf.param.length){
			te += '<!--Group effect-->\n';
			te += '<group_effect sigma="sigma"/>\n';
			par_list.push({name:"sigma", min:0,max:1});
		}
	}
		
	let ind_eff_group = model.species[0].ind_eff_group;
	
	te += '<!--Genetic covariance between different individual effects-->\n';
	
	for(let k = 0; k < ind_eff_group.length; k++){
		let ieg = ind_eff_group[k];
		
		te += '<covariance individual_effect="';
		let N = ieg.ie_list.length;
		for(let j = 0; j < N; j++){
			if(j != 0) te += ',';
			te += ieg.ie_list[j].name;
		}
		te += '" relationship_matrix="';
		if(ieg.A_matrix.check) te += "A";
		else te += "I";
		te += '">\n';
		te += '<variance>\n';
		for(let j = 0; j < N; j++){
			let na = ieg.ie_list[j].name;
			let na2 = 'Ω^'+na+","+na;
			te += na2+'\n';
			par_list.push({name:na2});
		}
		te += '</variance>\n';
		
		if(N > 1){
			te += '<correlation>\n';
			for(let j = 0; j < N; j++){
				for(let i = 0; i < N; i++){
					if(i != 0) te += "\t";
					if(j == i) te += '1';
					else{
						let na2;
						if(j < i) na2 = "ω_"+ieg.ie_list[j].name+","+ieg.ie_list[i].name;
						else na2 = "ω_"+ieg.ie_list[i].name+","+ieg.ie_list[j].name;
						te += na2;
						if(i < j) par_list.push({name:na2, min:-0.9,max:0.9});
					}
				}
				te += endl;
			}
			te += '</correlation>\n';
		}
		te += '</covariance>\n';
	}
	
	let Ate = "";
	for(let i = 0; i < ind_eff_group.length; i++){
		let ieg = ind_eff_group[i];
		if(ieg.A_matrix.check == true){
			let ha = new Hash();
			let ind_list = ieg.A_matrix.ind_list;
			for(let k = 0; k < ind_list.length; k++){
				ha.add(ind_list[k],k);
			}
			
			for(let k = 0; k < ind_list.length; k++){
				if(k != 0) Ate += ",";
				Ate += name_list[k].name;
			}
			Ate += '\n';
			
			let val = ieg.A_matrix.A_value;
			for(j = 0; j < N; j++){
				let jj = ha.find(name_list[j].name);
				if(jj == undefined) pr("prob");
				for(k = 0; k < N; k++){
					let kk = ha.find(name_list[k].name);
					if(kk == undefined) pr("prob");
					
					if(k != 0) Ate += ",";
					//pr(jj+" "+
					Ate += val[jj][kk];
				}
				Ate += "\n";
			}
			
			te += '<A name="A" file="A.csv"/>\n';
	
			break;
		}
	}
	
	let claz = model.species[0].cla[1];
 
	te += '<!--Inference and observation periods-->\n';
	te += '<inference group="';
	for(let c = 0; c < claz.comp.length; c++){
		if(c != 0) te += ',';
		te += claz.comp[c].name;
	}
	te += '" tmin="'+model.sim_details.t_start+'" tmax="'+model.sim_details.t_end+'"/>\n';
	
	te += '<observation group="';
	for(let c = 0; c < claz.comp.length; c++){
		if(c != 0) te += ',';
		te += claz.comp[c].name;
	}
	te += '" tmin="'+model.sim_details.t_start+'" tmax="'+model.sim_details.t_end+'"/>\n';
	
	te += '<!--Model priors-->\n';
	for(let j = 0; j < par_list.length; j++){
		let par = par_list[j];
		if(!par.min){
			let th = find(model.param,"name",par.name);
			if(th == undefined) pr("Cannot find prior:"+par.name);
			let pa = model.param[th];
			
			switch(pa.prior.type.te){
			case "fix":
				par.val = pa.prior.value.mean_eqn.te;
				//par.min = pa.prior.value.mean_eqn.te;
				//par.max = pa.prior.value.mean_eqn.te;
				break;
				
			case "uniform":
				par.min = pa.prior.value.min_eqn.te;
				par.max = pa.prior.value.max_eqn.te;
				break;
				
			default: pr("ERROR prior type not recognised"); break;
			}	
		}
		
		
		te += '<prior parameter="'+par.name+'" ';
		
		if(par.val != undefined) te += 'type="Fixed" val="'+par.val+'"/>\n';
		else te += 'type="Flat" val1="'+par.min+'" val2="'+par.max+'"/>\n';
	}
	
	te += '<!--Information about individuals-->\n';


	te += '<datatable file="table.csv" ';
	for(let c= 0; c < col.length; c++){
		te += col[c].title+'="'+col[c].title+'" ';
	}
	te += '/>\n';
	
	for(let k = 0; k < groups.length; k++){
		te += '<prediction_accuracy name="'+groups[k].name+'" ind="'+groups[k].te+'"/>\n';
	}
	
	te += '</SIRE>\n'

	post({te:te, Ate:Ate, dt:dt});
}
