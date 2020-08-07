function priorbuts()                                      // Prior page                      
{
	x = menux+tab; y = 100; 

	cornx = x+20; corny = y;
	addbutton("",cornx,corny,priorwidth,priorheight,CANVASBUT,CANVASBUT,-1,-1);

	if(drawpriors() == 1){
		addbutton("This page defines the prior distributions for all the parameters in the model.",x+15,70,150,0,-1,NTEXTBUT,-1,-1);
	}
	else{
		y = 70;
		if(page == MODELPAGE) addbutton("There are currently no distributions.",x+18,y,900,0,-1,PARAGRAPHBUT,1,-1);
		else addbutton("There are no model parameters.",x+18,y,900,0,-1,PARAGRAPHBUT,1,-1);
	}

	if(tableyfrac < 1) addbutton("",menux+tab+20+priorwidth+10,corny,13,priorheight,SLIDEAC,YSLIDEBUT,-1,-1);
}

function drawpriors()                                     // Draws the prior buttons
{
	var x, y, fl = 0, list=[], j, dy1 = 19, dy2 = 25, dy1 = priordy, dy2 = priordy;
	
	x = 0; nob = 0;
	
	switch(page){
	case INFERENCEPAGE:
		y = -8;
		list=[];
		for(p = 0; p < param.length; p++){ 
			if(param[p].classname != "Obs. Mod." && param[p].classname != "Initial"){
				if(param[p].dist == 1){ list.push(p); fl = 1;}
			}
		}
		if(list.length > 0){
			y += dy1; addob(x,y,OBPRIORTITLE,"Defined",217); y += dy2;
			for(j = 0; j < list.length; j++){ addob(x,y,OBPRIOR,list[j]); y += priordy;}
		}
		
		list=[];
		for(p = 0; p < param.length; p++){ 
			if(param[p].classname != "Obs. Mod." && param[p].classname != "Initial"){
				if(param[p].dist == 0){ list.push(p); fl = 1;}
			}
		}
		if(list.length > 0){
			y += dy1; addob(x,y,OBPRIORTITLE,"Compartmental model",218); y += dy2;
			for(j = 0; j < list.length; j++){ addob(x,y,OBPRIOR,list[j]); y += priordy;}
		}

		list=[];
		for(p = 0; p < param.length; p++){ 
			if(param[p].classname == "Obs. Mod."){ list.push(p); fl = 1;}
		}
		if(list.length > 0){
			y += dy1; addob(x,y,OBPRIORTITLE,"Observation model",219); y += dy2;
			for(j = 0; j < list.length; j++){ addob(x,y,OBPRIOR,list[j]); y += priordy;}
		}
		
		list=[];
		for(p = 0; p < param.length; p++){ 
			if(param[p].classname == "Initial"){ list.push(p); fl = 1;}
		}
		if(list.length > 0){
			y += dy1; addob(x,y,OBPRIORTITLE,"Initial state",220); y += dy2;
			for(j = 0; j < list.length; j++){ addob(x,y,OBPRIOR,list[j]); y += priordy;}
		}
		
		if(paramagesmooth.length > 0){
			y += dy1; addob(x,y,OBPRIORTITLE,"Age smoothing",221); y += dy2;
			for(j = 0; j < paramagesmooth.length; j++){ addob(x,y,OBAGESMOOTH,j); y += priordy;}
		}
		
		if(paramtimesmooth.length > 0){
			y += dy1; addob(x,y,OBPRIORTITLE,"Time smoothing",222); y += dy2;
			for(j = 0; j < paramtimesmooth.length; j++){ addob(x,y,OBTIMESMOOTH,j); y += priordy;}
		}
		break;
		
	case MODELPAGE:
		y = 15;
		for(p = 0; p < param.length; p++){ 
			if(param[p].classname != "Obs. Mod." && param[p].classname != "Initial"){
				if(adddepfl == 1){
					if(param[p].dist == 0){ addob(x,y,OBPRIOR2,p); fl = 1; y += priordy;}
				}
				else{
					if(param[p].dist == 1){ addob(x,y,OBPRIOR,p); fl = 1; y += priordy;}
				}
			}
		}
		break;
	}
	y += 70;
	ytot = y;
	
	tableyfrac = priorheight/ytot;
	placeob();	
	return fl;
}

function priorsubbut(param,p,k)                           // Expands a subscript on a prior variable
{
	q = param[p].dep[k];
	cl = 0; while(cl < ncla && cla[cl].name != q) cl++;
	if(cl < ncla){ // A class name
		for(j = 0; j < cla[cl].ncomp-1; j++){
			var v = copy(param[p]);
			param.splice(p,0,v);
		}
		for(j = 0; j < cla[cl].ncomp; j++) param[p+j].dep[k] = cla[cl].comp[j].name;
	}
	else{
		for(cl = 0; cl < ncla; cl++){
			for(j = 0; j < cla[cl].ncomp; j++) if(cla[cl].comp[j].name == q) break;
			if(j < cla[cl].ncomp) break;
		}
		if(cl < ncla){
			dif = 0;
			oth=[];
			for(pp = 0; pp < param.length; pp++){
				if(param[pp].name == param[p].name){
					kkmax = param[pp].dep.length; 
					if(kkmax == param[p].dep.length){
						for(kk = 0; kk < kkmax; kk++){
							if(kk != k){
								if(param[pp].dep[kk] != param[p].dep[kk]) break;
							}
						}
						if(kk == kkmax){
							if(param[pp].prior){
								if(param[pp].prior != param[p].prior) dif = 1;
								for(j = 0; j < param[p].val.length; j++) if(param[p].val[j] != param[pp].val[j]) dif = 1;
							}
							else{
								if(param[pp].sim != param[p].sim) dif = 1;
							}
							
							oth.push(pp);
						}
					}
				}
			}
			
			if(dif == 0 || (dif == 1 && confirm("Are you sure?"))){
				param[p].dep[k] = cla[cl].name;
				for(j = oth.length-1; j >= 0; j--) if(oth[j] != p) param.splice(oth[j],1);
			}
		}
		else alertp("Error code EC36");
	}
	return param;
}

function findparam(na)                                   // Finds a paramter from a name
{
	var vpri, dep=[], p, k, kk, cl, clname;

	sp = na.split("_");
	if(sp.length > 1){
		if(sp.length == 2) dep = sp[1].split(",");
		else alert("There is a problem with this expression:"+na);
	}
	
	for(p = 0; p < param.length; p++){
		if(param[p].name == sp[0] && param[p].dep.length == dep.length){
			for(k = 0; k < dep.length; k++){
				st = dep[k]; 
				clname = "";
				for(cl = 0; cl < ncla; cl++){
					for(kk = 0; kk < cla[cl].ncomp; kk++) if(cla[cl].comp[kk].name == st) clname = cla[cl].name;
				}
				
				for(kk = 0; kk < param[p].dep.length; kk++){ if(param[p].dep[kk] == st || param[p].dep[kk] == clname) break;}
				if(kk == param[p].dep.length) break; 
			}
			if(k == dep.length) return p;
		}
	}
	
	return param.length;
}

function collectvariables()                               // Collects together all the variables in the model
{
	var cl, i, n, d;

	paramnew = []; paramsimnew = [];
	for(cl = 0; cl < ncla; cl++){   // makes a list of parameters based on transitions
		for(i = 0; i < cla[cl].ntra; i++){
			for(n = 0; n < cla[cl].tra[i].nratevar; n++){
				addvar(cla[cl].name,cla[cl].tra[i].ratevar[n],cla[cl].tra[i].ratevardep[n]);
			}
		}
	}

	for(d = 0; d < data.length; d++){
		dat = data[d];
		switch(dat.variety){
		case "cap":
			if(dat.obspd == "set" && dat.pdsame == "same"){
				for(n = 0; n < dat.npdvar; n++){
					if("["+dat.pdvar[n]+"]" == dat.pd) addvar("Obs. Mod.",dat.pdvar[n],dat.pdvardep[n],"flat");
					else addvar("Obs. Mod.",dat.pdvar[n],dat.pdvardep[n]);
				}
			}
			break;
			
		case "trans":
			if(dat.obspd == "set"){
				for(n = 0; n < dat.npdvar; n++) addvar("Obs. Mod.",dat.pdvar[n],dat.pdvardep[n]);
			}
			break;
		
		case "cappd":
			for(n = 0; n < dat.npdvar; n++) addvar("Obs. Mod.",dat.pdvar[n],dat.pdvardep[n]);	
			break;
		
		case "state":
			switch(dat.type){
			case "test":
				na = dat.testname;
				//addvar("Obs. Mod.","Se("+na+")","","flat");
				//addvar("Obs. Mod.","Sp("+na+")","","flat2");
				addvar("Obs. Mod.","Se("+na+")","");
				addvar("Obs. Mod.","Sp("+na+")","");
				break;
			}			
			break;
		}
	}
			
	for(d = 0; d < derive.length; d++){
		for(n = 0; n < derive[d].ndervar; n++) addvar("Derived",derive[d].dervar[n],[]);
	}	
		
	st = ""; for(cl = 0; cl < ncla-1; cl++){ if(cla[cl].ncomp > 1){ if(st != "") st += ","; st += cla[cl].name;}}
	if(st != "") addvar("Initial","ξ",st);

	for(p = 0; p < param.length; p++){
		if(param[p].dist == 1 && onparamnew(p,param,paramnew) >= 0){
			switch(param[p].prior){			
			case "Gamma": case "Normal": case "Log-Normal": case "Beta": case "Weibull":
				if(param[p].val[0] != ""){
					nvarlist = 0; varlist=[]; vardeplist=[];   
					checkeqn(param[p].val[0]);
					for(j = 0; j < nvarlist; j++) addvar("Hyper",varlist[j],"");
				}
				
				if(param[p].val[1] != ""){
					nvarlist = 0; varlist=[]; vardeplist=[];   
					checkeqn(param[p].val[1]);
					for(j = 0; j < nvarlist; j++) addvar("Hyper",varlist[j],"");
				}
				
				break;
		
			case "Exponential":
				if(param[p].val[0] != ""){
					nvarlist = 0; varlist=[]; vardeplist=[];   
					checkeqn(param[p].val[0]);
					for(j = 0; j < nvarlist; j++) addvar("Hyper",varlist[j],"");
				}
				break;
			}
		}
	}
	paramnew = paramexpand(paramnew);
	param = paramexpand(param);	
	
	param = paramcombine(param,paramnew);
	
	for(p = 0; p < param.length; p++){
		if(param[p].dist == 0 && param[p].classname != "Initial" && param[p].classname != "Obs. Mod."){
			paramsimnew.push({name:param[p].name, dep:copy(param[p].dep), classname:param[p].classname, sim:"Unspecified"});
		}
	}
	param = paramcontract(param);

	paramsim = paramexpand(paramsim);
	paramsim = paramcombine(paramsim,paramsimnew);
	paramsim = paramcontract(paramsim);	

	smoothing();
}

function addvar(clname,st,dep,prior)                           // Adds a variables to the model
{
	if(dep == "") depa=[]; else depa = dep.split(",");
	
	v = 0; while(v < paramnew.length && paramnew[v].name != st) v++;
	if(v == paramnew.length){
		switch(clname){
		case  "Initial":
			paramnew.push({name:st, dep:depa, classname:clname, prior:"Dirichlet", val:[1], dist:0});
			break;
			
		case "Obs. Mod.":
			if(prior == "flat"){
				paramnew.push({name:st, dep:depa, classname:clname, prior:"Flat", val:[0,1], dist:0});
			}
			else{
				if(prior == "flat2"){
					paramnew.push({name:st, dep:depa, classname:clname, prior:"Flat", val:[0.5,1], dist:0});
				}
				else{
					paramnew.push({name:st, dep:depa, classname:clname, prior:"Unspecified", val:[], dist:0});	
				}
			}
			break;
			
		default:
			paramnew.push({name:st, dep:depa, classname:clname, prior:"Unspecified", val:[], dist:0});
			break;
		}
	}
	else{
		fl = 0;
		if(paramnew[v].dep.length != depa.length) fl = 1;
		else{ for(j = 0; j < depa.length; j++) if(paramnew[v].dep[j] != depa[j]) fl = 1;}
		if(fl == 1){ alertp("Error: Parameter "+st+" has different dependancies"); return;}
	}
}

function getparamname(param,p)                           // Gives the name of a parameter
{
	st = param[p].name;
	dmax = param[p].dep.length;
	if(dmax > 0){
		st += "_";
		for(d = 0; d < dmax; d++){
			st += param[p].dep[d]; if(d < dmax-1) st += ",";
		}
	}
	return st;
}

function paramchname(param,fr,to)                        // Changes the name of a parameter
{
	var p, j;
	
	for(p = 0; p < param.length; p++){
		if(param[p].classname == fr) param[p].classname = to;
		
		for(j = 0; j < param[p].dep.length; j++){
			if(param[p].dep[j] == fr) param[p].dep[j] = to;
		}
	}
}
					
function paramcombine(param,paramnew)                    // Combines equal parameters together
{
	var p, pp;
	
	for(p = 0; p < param.length; p++){
		pp = onparamnew(p,param,paramnew);
		if(pp != -1) paramnew[pp] = copy(param[p]);
	}
	param = copy(paramnew);

	return param;
}

function onparamnew(p,param,paramnew)                    // Works out if parameter already exists
{
	var pp, kmax, k;
	for(pp = 0; pp < paramnew.length; pp++){
		if(param[p].name == paramnew[pp].name){
			kmax = param[p].dep.length;
			if(kmax == paramnew[pp].dep.length){
				for(k = 0; k < kmax; k++) if(param[p].dep[k] != paramnew[pp].dep[k]) break;
				if(k == kmax) break;
			}
		}	
	}
	if(pp < paramnew.length) return pp; else return -1;
}

function paramexpand(paramin)                            // Exapands subscripts in parameters
{
	var k, kmax, param=[], nparam;
	
	kmax = 1; for(p = 0; p < paramin.length; p++){ nd = paramin[p].dep.length; if(nd > kmax) kmax = nd;}
	
	for(k = 0; k < kmax; k++){
		param=[];
		nparam = 0;
		for(p = 0; p < paramin.length; p++){
			nd = paramin[p].dep.length;
			don = 0;
			if(k < nd){
				na = paramin[p].dep[k];
				cl = 0; while(cl < ncla && na != cla[cl].name) cl++;
				if(cl < ncla){
					don = 1;
					for(j = 0; j < cla[cl].ncomp; j++){
						param[nparam] = copy(paramin[p]);
						param[nparam].dep[k] = cla[cl].comp[j].name;
						nparam++;
					}
				}
			}
			
			if(don == 0){ param[nparam] = copy(paramin[p]); nparam++;}
		}
		paramin = param;
	}
	
	return param;
}

function copy(inp){ return JSON.parse(JSON.stringify(inp));} // Copies an object

function paramcontract(param)                                // Contracts a parameter list
{
	var oth=[];
	kmax = 0; for(p = 0; p < param.length; p++){ nd = param[p].dep.length; if(nd > kmax) kmax = nd;}

	for(k = kmax-1; k >= 0; k--){
		for(p = 0; p < param.length; p++){
			if(k < param[p].dep.length){ 
				na = param[p].dep[k];
				for(cl = 0; cl < ncla; cl++){
					for(j = 0; j < cla[cl].ncomp; j++) if(cla[cl].comp[j].name == na) break;
					if(j < cla[cl].ncomp) break;
				}

				oth=[];
				for(pp = p+1; pp < param.length; pp++){
					same = 0;
					if(param[pp].name == param[p].name){
						if(param[p].prior){
							if(param[pp].prior == param[p].prior){
								for(j = 0; j < param[p].val.length; j++) if(param[p].val[j] != param[pp].val[j]) break;
								if(j == param[p].val.length) same = 1;
							}
						}
						else{
							if(param[pp].sim == param[p].sim) same = 1;
						}
					}
					
					if(same == 1){
						kkmax = param[pp].dep.length; 
						if(kkmax == param[p].dep.length){
							for(kk = 0; kk < kkmax; kk++){
								if(kk != k){
									if(param[pp].dep[kk] != param[p].dep[kk]) break;
								}
							}
							if(kk == kkmax) oth.push(pp);
						}
					}
				}
				
				if(cl < ncla){
					if(oth.length > cla[cl].ncomp-1) alertp("Error code EC35");
					if(oth.length == cla[cl].ncomp-1){
						param[p].dep[k] = cla[cl].name;
						for(j = oth.length-1; j >= 0; j--) param.splice(oth[j],1);
					}
				}					
			}
		}
	}
	return param;
}

function checkprior(dist)                               // Checks the values in the priors
{
	var p, val1, val2, pa, su, ty;

	if(dist == 0){ pa = INFERENCEPAGE; su = 1; ty = "Prior";} else{ pa = MODELPAGE; su = -3; ty = "Distribution";}

	for(p = 0; p < param.length; p++){
		if(param[p].dist == dist){
			if(param[p].prior == "Unspecified"){ changepage(pa,su,0); alertp2(ty+" unspecified for "+param[p].name+"."); return 1;}
		
			if(param[p].val.length == 0){ changepage(pa,su,0); alertp2(ty+" not set for "+param[p].name+"."); return 1;}
		
			if(dist == 0){
				val1 = parseFloat(param[p].val[0]);
				if(isNaN(val1)){ changepage(pa,su,0); alertp2(ty+" for "+param[p].name+" must be set."); return 1;}

				if(typeof param[p].val[1] != "undefined"){
					val2 = parseFloat(param[p].val[1]);
					if(isNaN(val2)){ changepage(pa,su,0); alertp2(ty+" for "+param[p].name+" must be set."); return 1;}
				}
			}
			else{
				param[p].val[0] = param[p].val[0].trim();
				if(param[p].val[0].length == 0){ alertp2(ty+" for "+param[p].name+" must be set."); return 1;}
				
				if(typeof param[p].val[1] != "undefined"){
					param[p].val[1] = param[p].val[1].trim();
					if(param[p].val[1].length == 0){ alertp2(ty+" for "+param[p].name+" must be set."); return 1;}
				}				
			}
			
			switch(param[p].prior){
			case "Unbounded":
				break;
				
			case "Flat":
				if(dist == 0 && val1 >= val2){
					changepage(pa,su,0);
					alertp2("For "+param[p].name+" the maximum must be larger than the minimum.");					
					return 1;
				}
				break;
				
			case "Gamma":
				if(dist == 0 && val1 <= 0){
					changepage(INFERNECEPAGE,1,0);
					alertp2("For "+param[p].name+" the mean must be greater than zero.");					
					return 1;
				}

				if(dist == 0 && val2 <= 0){
					changepage(pa,su,0);
					alertp2("For "+param[p].name+" the standard deviation must be greater than zero.");					
					return 1;
				}
				break;
		
			case "Normal": case "Log-Normal":
				if(dist == 0 && val2 <= 0){
					changepage(pa,su,0);
					alertp2("For "+param[p].name+" the standard deviation must be greater than zero.");					
					return 1;
				}
				break;
				
			case "Fix":
				break;
			}
		}
	}

	if(dist == 0){ 
		for(p = 0; p < paramagesmooth.length; p++){
			if(paramagesmooth[p].type != "None"){
				val = paramagesmooth[p].val;
				if(val == ""){
					changepage(pa,su,0);
					alertp2("The value for the smoothing prior for "+paramagesmooth[p].name+" must be set.");					
					return 1;
				}
				else{
					if(isNaN(val)){
						changepage(pa,su,0);
						alertp2("The value for the smoothing prior for "+paramagesmooth[p].name+" must be a number.");					
						return 1;
					}
				}
			}
		}
		
		for(p = 0; p < paramtimesmooth.length; p++){
			if(paramtimesmooth[p].type != "None"){
				val = paramtimesmooth[p].val;
				if(val == ""){
					changepage(pa,su,0);
					alertp2("The value for the smoothing prior for "+paramtimesmooth[p].name+" must be set.");					
					return 1;
				}
				else{
					if(isNaN(val)){
						changepage(pa,su,0);
						alertp2("The value for the smoothing prior for "+paramtimesmooth[p].name+" must be a number.");					
						return 1;
					}
				}
			}
		}
	}
		
	return 0;
}

function setsimparam(na,val)                              // Sets the values used for simulation
{
	var dep=[], sp = na.split("_"), n;
	nam = sp[0];
	if(sp.length > 1){
		if(sp.length == 2) dep = sp[1].split(",");
		else{ alertimp("There is a problem with this expression"); return 1;}
	}

	for(n = 0; n < paramsim.length; n++) if(paramsim[n].name == nam) break;
	if(n == paramsim.length) collectvariables();
	for(n = 0; n < paramsim.length; n++) if(paramsim[n].name == nam) break;
	if(n == paramsim.length){ alertimp("Cannot find parameter"); return 1;}
		
	paramsim = paramexpand(paramsim);
	for(n = 0; n < paramsim.length; n++){
		if(paramsim[n].name == nam){
			if(dep.length != paramsim[n].dep.length){ alertimp("The dependencies do not match"); return 1;}
			for(nn = 0; nn < dep.length; nn++){
				if(dep[nn] != paramsim[n].dep[nn]){
					cl2 = getclfromcomp(paramsim[n].dep[nn]);   
					if(dep[nn] != cla[cl2].name) break;
				}
			}
			if(nn == dep.length) paramsim[n].sim = val;	
		}
	}
	return 0;
}

function setinfprior(na,prior,val,dist)                   // Sets the prior used for inference       
{
	var dep=[], sp = na.split("_"), n;
	
	nam = sp[0];
	if(sp.length > 1){
		if(sp.length == 2) dep = sp[1].split(",");
		else{ alertimp("There is a problem with this expression"); return 1;}
	}

	for(n = 0; n < param.length; n++) if(param[n].name == nam) break;
	if(n == param.length) collectvariables();
	for(n = 0; n < param.length; n++) if(param[n].name == nam) break;
	if(n == param.length){ alertimp("Cannot find parameter"); return 1;}
			
	param = paramexpand(param);

	for(n = 0; n < param.length; n++){
		if(param[n].name == nam){
			if(dep.length != param[n].dep.length){ alertimp("The dependencies do not match"); return 1;}
			for(nn = 0; nn < dep.length; nn++){
				if(dep[nn] != param[n].dep[nn]){
					cl2 = getclfromcomp(param[n].dep[nn]);   
					if(dep[nn] != cla[cl2].name) break;
				}
			}
			if(nn == dep.length){ param[n].prior = prior; param[n].val = val; param[n].dist = dist;}	
		}
	}
	return 0;
}

function smoothing()
{
	var paramagesmoothnew=[], paramtimesmoothnew=[];
	
	for(p = 0; p < param.length; p++){
		if(param[p].classname != "Initial"){
			for(j = 0; j < param[p].dep.length; j++){
				if(param[p].dep[j] == "Age"){
					for(pp = 0; pp < paramagesmooth.length; pp++){ if(paramagesmooth[pp].name == param[p].name) break;}
					if(pp == paramagesmooth.length){
						paramagesmoothnew.push({name:param[p].name, type:"None", val:""});	
					}
					else{
						paramagesmoothnew.push(paramagesmooth[pp]);
					}
				}
			}
		}
	}
	
	for(p = 0; p < param.length; p++){
		if(param[p].classname != "Initial"){
			for(j = 0; j < param[p].dep.length; j++){
				if(param[p].dep[j] == "Time"){
					for(pp = 0; pp < paramtimesmooth.length; pp++){ if(paramtimesmooth[pp].name == param[p].name) break;}
					if(pp == paramtimesmooth.length){
						paramtimesmoothnew.push({name:param[p].name, type:"None", val:""});	
					}
					else{
						paramtimesmoothnew.push(paramtimesmooth[pp]);
					}
				}
			}
		}
	}
	
	paramagesmooth = paramagesmoothnew; paramtimesmooth = paramtimesmoothnew;
}

