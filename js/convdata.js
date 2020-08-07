function converttoobs(ty)                                 // Converts data sources into individual measurements
{
	var nindtotal, ind=[];
	var dat, probexp=[], tdatamin, tdatamax, d, i, j, te="", cl2, cl;

	if(ty == "sim") da = simdata;
	else da = data;
	
	tdatamin = large; tdatamax = -large;
	for(d = 0; d < da.length; d++){
		dat = da[d];
		if(!isNaN(dat.tmin)){
			if(dat.tmin < tdatamin && dat.tmin != -large) tdatamin = dat.tmin;
			if(dat.tmax > tdatamax && dat.tmax != large) tdatamax = dat.tmax;
		}
	}

	nindtotal = 0;
	ind=[];
	
	for(d = 0; d < da.length; d++){
		dat = da[d];

		if(dat.variety != "cap" && dat.variety != "capid" && dat.variety != "cappd" &&
 		   dat.variety != "pop" && dat.variety != "der" && dat.variety != "presence"){
			for(i = 0; i < dat.id.length; i++){
				id = dat.id[i]; 
				
				if(dat.t[i] == "All") t = tdatamin;
				else{ dat.t[i] = parseFloat(dat.t[i]); t = dat.t[i];}
		
				j = 0; while(j < nindtotal && ind[j].id != id) j++;	
				if(j == nindtotal){
					ind[nindtotal] = {id:id, cl:[]};
					for(cl = 0; cl < ncla; cl++) ind[nindtotal].cl[cl] = {ev:[]};
					nindtotal++;
				}
		
				cl = dat.cl;
				
				switch(dat.variety){
				case "state":
					val = dat.val[i];
					p = 0; pmax = dat.pos.length; while(p < pmax && dat.pos[p] != val) p++;
					if(p == pmax){ alertp("Error code EC26");}
			
					for(k = 0; k < cla[cl].ncomp; k++){
						switch(dat.type){
						case "simple":
							if(k != dat.posref[p]) probexp[k] = "0";
							else probexp[k] = "1";	
							break;
						
						case "binary":
							if(dat.posbinary[p][k] == 0) probexp[k] = "0";
							else probexp[k] = "1";	
							break;
							
						case "test":
							na = dat.testname;
							switch(dat.postestres[p]){
							case 1:
								if(dat.sensitive[k] == 1)  probexp[k] = "[Se("+na+")]";
								else probexp[k] = "1-[Sp("+na+")]";
								break;
						
							case 0:
								if(dat.sensitive[k] == 1)  probexp[k] = "1-[Se("+na+")]";
								else probexp[k] = "[Sp("+na+")]";
								break;
						
							case 2:
								probexp[k] = "1";
								break;
							}	
							break;
						
						case "expression":
							probexp[k] = dat.posexpression[p][k];
							break;
						}
					}
					addstatedata(d,ind,i,j,t,cl,probexp);
					break;
				
				case "trans":
					for(cl2 = 0; cl2 < 2; cl2++){
						if(dat.filt[cl2] != "All"){
							for(c = 0; c < cla[cl2].ncomp; c++){
								if(cla[cl2].comp[c].name == dat.filt[cl2]) probexp[c] = "1"; else probexp[c] = "0";
							}
							addstatedata(d,ind,i,j,t,cl2,probexp);
						}
					}
					if(dat.transty == "+" || dat.transty == "-") cl = ncla-2;
					
					k = 0; kmax = ind[j].cl[cl].ev.length; while(k < kmax && t > ind[j].cl[cl].ev[k].t) k++;
					if(k < kmax && t == ind[j].cl[cl].ev[k].t){
						te = "The transtions on individual "+ind[j].id+" in data sources '"+dat.name+"' and '";
						if(ind[j].cl[cl].ev[k].variety == "state") te += data[ind[j].cl[cl].ev[k].obsdata[0]].name;
						else te += data[ind[j].cl[cl].ev[k].obsdata].name;
						te += "' are at exactly the same time t="+t+".";
					}
					else{
						newev = {t:t, variety:"trans", obsdata:d, transi:dat.transi, transf:dat.transf, transty:dat.transty, col:[]};
						ind[j].cl[cl].ev.splice(k,0,newev);
					}
					
					t = dat.tmin;
					k = 0; kmax = ind[j].cl[cl].ev.length; while(k < kmax && t > ind[j].cl[cl].ev[k].t) k++;
					var evtmin = {t:t, variety:"transtmin", name:dat.name};
					ind[j].cl[cl].ev.splice(k,0,evtmin);	
					
					t = dat.tmax;
					k = 0; kmax = ind[j].cl[cl].ev.length; while(k < kmax && t > ind[j].cl[cl].ev[k].t) k++;
					var evtmax = {t:t, variety:"transtmax", name:dat.name};
					ind[j].cl[cl].ev.splice(k,0,evtmax);	
					break;
					
				case "move":
					if(dat.transty == "+" || dat.transty == "-") cl = ncla-2;
					k = 0; kmax = ind[j].cl[cl].ev.length; while(k < kmax && t > ind[j].cl[cl].ev[k].t) k++;
					if(k < kmax && t == ind[j].cl[cl].ev[k].t){
						te = "The transtions on individual "+ind[j].id+" in data sources '"+dat.name+"' and '";
						if(ind[j].cl[cl].ev[k].variety == "state") te += data[ind[j].cl[cl].ev[k].obsdata[0]].name;
						else te += data[ind[j].cl[cl].ev[k].obsdata].name;
						te += "' are at exactly the same time t="+t+".";
					}
					else{
						newev = {t:t, variety:"move", obsdata:d, transi:dat.transi, transf:dat.transf, transty:dat.transty, col:[]};
						ind[j].cl[cl].ev.splice(k,0,newev);
					}
					break;
				}
			}
		}
	}
	
	for(d = 0; d < da.length; d++){
		dat = da[d];
		if(dat.variety == "presence"){
			for(i = 0; i < dat.id.length; i++){
				id = dat.id[i];
				
				if(dat.t[i] == "All") t =  tdatamin;
				else{ dat.t[i] = parseFloat(dat.t[i]); t = dat.t[i];}
					
				j = 0; while(j < nindtotal && ind[j].id != id) j++;	
				if(j == nindtotal){
					ind[nindtotal] = {id:id, cl:[]};
					for(cl = 0; cl < ncla; cl++) ind[nindtotal].cl[cl] = {ev:[]};
					nindtotal++;
				}
				
				cl = 0;
				k = 0; kmax = ind[j].cl[cl].ev.length; 
				while(k < kmax && (t > ind[j].cl[cl].ev[k].t || 
					    (t == ind[j].cl[cl].ev[k].t && ind[j].cl[cl].ev[k].variety != "state"))) k++;
				if(k == kmax || t > ind[j].cl[cl].ev[k].t){ 
					newev = {t:t, variety:"presence"};
					ind[j].cl[cl].ev.splice(k,0,newev);
				}		
			}				
		}
	}
	
	ind.sort( function(a, b){ return orderstring(a.id,b.id);});
	
	for(i = 0; i < nindtotal; i++){  	// Works out colours
		for(cl = 0; cl < ncla; cl++){
			for(e = 0; e < ind[i].cl[cl].ev.length; e++){
				ev = ind[i].cl[cl].ev[e];
				switch(ev.variety){
				case "state":
					var testres = -1;
				
					for(ee = 0; ee < ev.obsdata.length; ee++){  // works out if it could be a diagnostic test result
						d = ev.obsdata[ee];
						if(da[d].type != "test") break;
					}
					
					if(ee == ev.obsdata.length){  // Test results
						for(ee = 0; ee < ev.obsdata.length; ee++){
							d = ev.obsdata[ee];
							val = da[d].val[ev.obsdatai[ee]];
							
							for(j = 0; j < da[d].pos.length; j++) if(da[d].pos[j] == val) break;
							if(j == da[d].pos.length) alert("not pos");
							switch(da[d].postestres[j]){
							case 0: ev.col.push(-2); break;
							case 1: ev.col.push(-1); break;
							case 2: ev.col.push(-3); break;
							default: alert("Prob"); break;
							}
						}
					}
					else{   // Normal results
						for(j = 0; j < cla[cl].ncomp; j++){
							if(ev.probexp[j] != "0") ev.col.push(j);
						}
					}
					break;
					
				case "trans": case "move":
					switch(ev.transty){
					case "+": ev.col.push(-1); break;
					case "-": ev.col.push(-2); break;
					case "trans":
						for(j = 0; j < cla[cl].ncomp; j++) if(cla[cl].comp[j].name == ev.transi){ ev.col.push(j); break;}
						if(j == cla[cl].ncomp) alertp("Error code EC28");
						for(j = 0; j < cla[cl].ncomp; j++) if(cla[cl].comp[j].name == ev.transf){ ev.col.push(j); break;}
						if(j == cla[cl].ncomp) alertp("Error code EC29");
						break;
					}
					break;
				}
			}			
		}
	}
	
	if(ty == "sim"){
		indsim = {nindtotal:nindtotal, ind:ind};
	}
	else{
		var indlist=[]; for(i = 0; i < nindtotal; i++) indlist.push(ind[i].id);
		datares={indlist:indlist, tmin:tdatamin, tmax:tdatamax, ncla:ncla, cla:copy(cla), agecl:ncla-2, settimecl:ncla-1};
		
		inddata = {nindtotal:nindtotal, ind:ind};
		dt = tdatamax-tdatamin;
		if(tdatamin < tpostmin) tpostmin = tdatamin;
		if(tdatamax > tpostmax) tpostmax = tdatamax;
	}
	
	if(te != ""){ changepage(INFERENCEPAGE,0,0); alertp(te); return 0;}
	return 1;
}

function addstatedata(d,ind,i,j,t,cl,probexp)    // Adds information about state
{
	var k, kmax, jj, st, expold, expnew;
	
	k = 0; kmax = ind[j].cl[cl].ev.length; 
	while(k < kmax && (t > ind[j].cl[cl].ev[k].t || 
		(t == ind[j].cl[cl].ev[k].t && ind[j].cl[cl].ev[k].variety != "state"))) k++;
		
	if(k < kmax && t == ind[j].cl[cl].ev[k].t){  // Multiple times
		ind[j].cl[cl].ev[k].obsdata.push(d);
		ind[j].cl[cl].ev[k].obsdatai.push(i);
	
		for(jj = 0; jj < cla[cl].ncomp; jj++){
			expold = ind[j].cl[cl].ev[k].probexp[jj];
			expnew = probexp[jj];
			if(expold == "0" || expnew == "0") st = "0";
			else{
				if(expold == "1") st = expnew;
				else{
					if(expnew == "1") st = expold;
					else{
						st = "("+expold+")*("+expnew+")";
					}
				}
			}
			ind[j].cl[cl].ev[k].probexp[jj] = st;
		}				
	}
	else{
		newev = {t:t, variety:"state", probexp:[], obsdata:[], obsdatai:[], col:[]};
		ind[j].cl[cl].ev.splice(k,0,newev);
		ind[j].cl[cl].ev[k].obsdata[0] = d; 
		ind[j].cl[cl].ev[k].obsdatai[0] = i;
		for(jj = 0; jj < cla[cl].ncomp; jj++) ind[j].cl[cl].ev[k].probexp[jj] = probexp[jj];
	}
}
