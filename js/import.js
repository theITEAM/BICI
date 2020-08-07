
function importfile(text)                                 // Import a script to define (or partially define) model and data 
{
	var i, j, ist, type, cl, cl2, tr, trr, type, te, x, y, quote, quotetype;
	var lines = text.split('\n'), linesref=[];
	var content=[], frag=[], clcompch=[], cltrach=[];
		
	textFromFile = save(0);
	
	for(j = 0; j < lines.length; j++) linesref[j] = j;
	
	j = 0;     // Gets rid of comments
	while(j < lines.length){ 
		lines[j] = lines[j].trim();
		if(lines[j].length >= 2 && lines[j].substr(0,2) == "//"){ lines.splice(j,1); linesref.splice(j,1);}
		else{
			if(lines[j].length >= 1 && lines[j].substr(0,1) == "#"){ lines.splice(j,1); linesref.splice(j,1);}
			else j++;
		}
	}
	
	for(j = 0; j < lines.length; j++){
		trr = lines[j];
		if(trr != ""){
			frag=[]; jsto = linesref[j];
			i = 0; quote = 0;
			do{
				ist = i; i++; while(i == trr.length && quote==1 && j < lines.length-1){ j++; trr += lines[j].trim()+"\n";}
				while(i < trr.length && trr.substr(i,1) != "=" && 
						!(quote == 0 && trr.substr(i,1) == '"')&& 
						!(quote == 0 && trr.substr(i,1) == "'") && 
						!(quote == 1 && trr.substr(i,1) == '"' && quotetype == 0) && 
						!(quote == 1 && trr.substr(i,1) == "'" && quotetype == 1) && 
				      !(quote == 0 && trr.substr(i,1) == " ")){
					i++; while(i == trr.length && quote==1 && j < lines.length-1){ j++; trr += lines[j].trim()+"\n";}
				}
				//if(i < trr.length && trr.substr(i,1) == '"'
				
				te = trr.substr(ist,i-ist);
				while(te.substr(0,1) == "\n") te = te.substring(1);
			    while(te.substr(te.length-1,1) == "\n") te = te.substring(0,te.length-1);
				te = te.trim(); if(quote == 0) te = te.toLowerCase();
				if(te != "") frag.push({text:trr.substr(ist,i-ist).trim(), quote:quote});
				while(i < trr.length && quote == 0 && trr.substr(i,1) == " ") i++;
				if(i < trr.length){
					if(trr.substr(i,1) == '"' || trr.substr(i,1) == "'"){ 
						if(trr.substr(i,1) == '"') quotetype = 0; else quotetype = 1;
						quote = 1-quote; i++;
					}
				
/*				
					if(quote == 0){
						if(trr.substr(i,1) == '"' || trr.substr(i,1) == "'"){ 
							if(trr.substr(i,1) == '"') quotetype = 0; else quotetype = 1;
							quote = 1; i++;
						}
					}
					else{
						if((trr.substr(i,1) == '"' && quotetype == 0) || (trr.substr(i,1) == "'" && quotetype == 1)){
							quote = 0; i++;
						}						
					}
					*/
				}
				while(i == trr.length && quote==1 && j < lines.length-1){ j++; trr += lines[j].trim()+"\n";}
			}while(i < trr.length);
			
			if(frag.length == 0){ alertimp("Does not contain any content"); return;}
			
			if(frag[0].quote == 1){ alertimp("Should not start with quotes"); return;}
			type = frag[0].text; typest = {type:type, j:j};
			num = (frag.length-1)/3;

			if(num != Math.floor(num)){ alertimp("Syntax error"); return;}
			tags=[];
			for(n = 0; n < num; n++){
				if(frag[1+3*n+0].quote != 0){ alertimp("Syntax error1"); return;}
				if(frag[1+3*n+1].text != "="){ alertimp("Syntax error2"); return;}
				if(frag[1+3*n+2].quote != 1){ alertimp("Syntax error3"); return;}
				if(frag[1+3*n+0].text == ""){ alertimp("Syntax error4"); return;}
				if(frag[1+3*n+2].text == ""){ alertimp("Property "+frag[1+3*n+0]+" must have content"); return;}
				tags.push({name:frag[1+3*n+0].text,value:frag[1+3*n+2].text,done:0});
			}
			
			for(n = 0; n < tags.length-1; n++){
				for(nn = n+1; nn < tags.length; nn++){
					if(tags[n].name == tags[nn].name){ alertimp("The tag '"+tags[n].name+"' is set more than once"); return;}
				}
			}
			
			switch(type){
			case "compartment":
				claname = gettag("classification"); if(claname == ""){ cannotfind(); return;}
				if(claname == "Age"){ alertimp("Classification cannot have the name 'Age'"); return;}
				if(claname == "Time"){ alertimp("Classification cannot have the name 'Time'"); return;}
			
				name = gettag("name"); if(name == ""){ cannotfind(); return;}
				x = gettag("x"); if(x == ""){ cannotfind(); return;}
				y = gettag("y"); if(y == ""){ cannotfind(); return;}
				x = Number(x); if(isNaN(x)){ alertimp("x must be a number"); return;}
				y = Number(y); if(isNaN(y)){ alertimp("y must be a number"); return;}
				
				color = gettag("color"); if(color == ""){ cannotfind(); return;}
				
				cl = 0; while(cl < ncla-2 && cla[cl].name != claname) cl++;
				if(cl == ncla-2) addclass(claname);
				else{
					if(clcompch[cl] != 1){ while(cla[cl].ncomp > 0) deletecomp(cl,0);}
				}
				clcompch[cl] = 1;
			
				for(cl2 = 0; cl2 < ncla; cl2++){
					c = 0; while(c < cla[cl2].ncomp && cla[cl2].comp[c].name != name) c++;
					if(c < cla[cl2].ncomp){ alertimp("The compartment "+name+" already exists"); return;}
				}
				
				addcompartment(cl,name,color,x,y,0,0,0);
				break;
			
			case "transition":
				type = gettag("type"); if(type == ""){ cannotfind(); return;}
				
				if(type == "source"){
					fr = "Source"; 
					te = gettag("from"); if(te != ""){ alertimp("A source transition should not have a 'from' tag"); return;}
				}
				else{ fr = gettag("from"); if(fr == ""){ cannotfind(); return;}}
				
				if(type == "sink"){
					to = "Sink"; 
					te = gettag("to"); if(te != ""){ alertimp("A sink transition should not have a 'to' tag"); return;}
				}
				else{ to = gettag("to"); if(to == ""){ cannotfind(); return;}}
				
				if(fr == "Source"){
					ci = -2;
					for(cl = 0; cl < ncla-2; cl++){
						for(cf = 0; cf < cla[cl].ncomp; cf++){ if(cla[cl].comp[cf].name == to) break;}
						if(cf < cla[cl].ncomp) break;
					}
					if(cl == ncla-2){ alertimp("Cannot find compartment '"+to+"'"); return;}
				}
				else{
					if(to == "Sink"){
						cf = -2;
						for(cl = 0; cl < ncla-2; cl++){
							for(ci = 0; ci < cla[cl].ncomp; ci++){ if(cla[cl].comp[ci].name == fr) break;}
							if(ci < cla[cl].ncomp) break;
						}
						if(cl == ncla-2){ alertimp("Cannot find compartment '"+fr+"'"); return;}
					}
					else{
						for(cl = 0; cl < ncla-2; cl++){
							for(ci = 0; ci < cla[cl].ncomp; ci++){ if(cla[cl].comp[ci].name == fr) break;}
							if(ci < cla[cl].ncomp) break;
						}
						if(cl == ncla-2){ alertimp("Cannot find compartment '"+fr+"'"); return;}
						
						for(cl2 = 0; cl2 < ncla-2; cl2++){
							for(cf = 0; cf < cla[cl2].ncomp; cf++){ if(cla[cl2].comp[cf].name == to) break;}
							if(cf < cla[cl2].ncomp) break;
						}
						if(cl2 == ncla-2){ alertimp("Cannot find compartment '"+to+"'"); return;}
						
						if(cl != cl2){ alertimp("Transition must be within a classification"); return;}
					}
				}
				
				if(cltrach[cl] != 1){
					while(cla[cl].ntra > 0) deletetrans(cl,0);
					cltrach[cl] = 1;
				}
				
				for(tr = 0; tr < cla[cl].ntra; tr++){
					if(cla[cl].tra[tr].i == ci && cla[cl].tra[tr].f == cf){
						alertimp("Transition from '"+fr+"' to '"+to+"' already exists"); return;
					}
				}
				
				tr = cla[cl].ntra;
				switch(type){
				case "markovian": addtrans(cl,"Exponential"); break;
				case "source": addtrans(cl,"Source"); break;
				case "sink": addtrans(cl,"Sink"); break;
				case "gamma": addtrans(cl,"Gamma"); break;
				case "weibull": addtrans(cl,"Weibull"); break;
				default: alertimp("No transition type '"+type+"'"); return;
				}
				addtransmode = 0;
				
				cla[cl].tra[tr].i = ci;
				cla[cl].tra[tr].f = cf;
				
				x = gettag("x");
				if(x != ""){
					var xlist = x.split(",");
					y = gettag("y");
					if(y == ""){ alertimp("Must set x and y points"); return;}
					var ylist = y.split(",");
					if(xlist.length != ylist.length){ alertimp("Must have the same number of x and y points"); return;}
			
					for(n = 0; n < xlist.length; n++){
						xlist[n] = Number(xlist[n]); if(isNaN(xlist[n])){ alertimp("x must be a number"); return;}
						ylist[n] = Number(ylist[n]); if(isNaN(ylist[n])){ alertimp("y must be a number"); return;}
					}
					
					for(n = 0; n < xlist.length; n++) cla[cl].tra[tr].p.push({x:xlist[n],y:ylist[n]});
				}
				else{
					if(type == "source"){ alertimp("The position of the source must be set"); return;}
					if(type == "sink"){ alertimp("The position of the sink must be set"); return;}
				}
				
				switch(type){
				case "markovian": case "source": case "sink":
					rate = gettag("rate"); 
					if(rate != "") cla[cl].tra[tr].rate = rate;
					else{
						time = gettag("time"); if(time == ""){ tagst = "rate/time"; cannotfind(); return;}
						cla[cl].tra[tr].rate = time;
						cla[cl].tra[tr].ratetime = 1;
					}
					break;
					
				case "gamma":
					mean = gettag("mean"); if(mean == ""){ cannotfind(); return;}
					shape = gettag("shape"); if(shape == ""){ cannotfind(); return;}
					cla[cl].tra[tr].mean = mean;
					cla[cl].tra[tr].shape = shape;
					break;
					
				case "weibull":
					lam = gettag("lambda"); if(lam == ""){ cannotfind(); return;}
					k = gettag("k"); if(k == ""){ cannotfind(); return;}
					cla[cl].tra[tr].lam = lam;
					cla[cl].tra[tr].k = k;
					break;
				}
				checkrate(cl,tr);
				if(errmsg != ""){ alertimp(errmsg); return;}
				break;
				
			case "agetransitions":
				age=[];
				val = gettag("values");
				if(val != ""){
					var vallist = val.split(",");
					for(n = 0; n < vallist.length; n++){
						vallist[n] = Number(vallist[n]); if(isNaN(vallist[n])){ alertimp("Age transitions must be numbers"); return;}
						if(vallist[n] <= 0){ alertimp("Age transitions must be greater than zero"); return;}
					}
					for(n = 0; n < vallist.length-1; n++){
						if(vallist[n] > vallist[n+1]){ alertimp("Age transitions must be in order"); return;}
					}
					for(n = 0; n < vallist.length; n++) age[n] = vallist[n];
				}
				setageon();
				break;
				
			case "timetransitions":
				time=[];
				val = gettag("values");
				if(val != ""){
					var vallist = val.split(",");
					for(n = 0; n < vallist.length; n++){
						vallist[n] = Number(vallist[n]); if(isNaN(vallist[n])){ alertimp("Time transitions must be numbers"); return;}
					}
					for(n = 0; n < vallist.length-1; n++){
						if(vallist[n] > vallist[n+1]){ alertimp("Time transitions must be in order"); return;}
					}
					for(n = 0; n < vallist.length; n++) time[n] = vallist[n];
				}
				settimeon();
				break;
				
			case "setparam":
				par = gettag("param"); if(par == ""){ cannotfind(); return;}
				val = gettag("value"); if(val == ""){ cannotfind(); return;}
				if(setsimparam(par,val) == 1) return;
				break;
				
			case "setprior": case "setdistribution":
				if(type == "setprior") dist = 0; else dist = 1;
				
				par = gettag("param"); if(par == ""){ cannotfind(); return;}
				pri = gettag("prior"); if(pri == ""){ cannotfind(); return;}
				pri = pri.substr(0,1).toUpperCase()+pri.substring(1).toLowerCase();
				if(pri == "Log-normal") pri = "Log-Normal";
				
				switch(pri){
				case "Flat":
					min = gettag("min"); if(min == ""){ cannotfind(); return;}
					max = gettag("max"); if(max == ""){ cannotfind(); return;}
					val=[min,max];
					break;
					
				case "Gamma": case "Normal": case "Log-Normal":
					mean = gettag("mean"); if(mean == ""){ cannotfind(); return;}
					sd = gettag("sd"); if(sd == ""){ cannotfind(); return;}
					val=[mean,sd];
					break;
					
				case "Exponential":
					rate = gettag("rate"); if(rate == ""){ cannotfind(); return;}
					val=[rate];
					break;
					
				case "Beta":
					alpha = gettag("alpha"); if(alpha == ""){ cannotfind(); return;}
					beta = gettag("beta"); if(beta == ""){ cannotfind(); return;}
					val=[alpha,beta];
					break;
				
				case "Weibull":
					lambda = gettag("lambda"); if(lambda == ""){ cannotfind(); return;}
					k = gettag("k"); if(k == ""){ cannotfind(); return;}
					val=[lambda,k];
					break;
					
				case "Fix":
					value = gettag("value"); if(value == ""){ cannotfind(); return;}
					val=[value];
					break;
				}
				
				for(n = 0; n < val.length; n++){
					if(dist == 0){
						te = val[n];
						val[n] = Number(val[n]); if(isNaN(val[n])){ alertimp("'"+te+"' must be a number"); return;}
					}
					else{
						te = Number(val[n]); if(!isNaN(te)){ alertimp("'"+te+"' must not be a number"); return;}
					}
				}
				if(setinfprior(par,pri,val,dist) == 1) return;
				break;
			
			case "setinitprior":
				state = gettag("state"); if(state == ""){ cannotfind(); return;}
				value = gettag("value"); if(value == ""){ cannotfind(); return;}
				setinfprior("ξ_"+state,"Dirichlet",[value],dist);
				break;
			
			case "addderived":
				par = gettag("param"); if(par == ""){ cannotfind(); return;}
				exp = gettag("expression"); if(exp == ""){ cannotfind(); return;}
				if(addderive(par) == 1){ alertimp("The derived name cannot contain any of the following characters: space, underscore, comma, '[', ']', '{', '}'"); return;}
				if(setderive(derive.length-1,exp) == 1){ alertimp(errmsg); return;}
				break;
					
			case "siminitpop":
				clfl = -1;
				for(cl = 0; cl < ncla-1; cl++){
					for(c = 0; c < cla[cl].ncomp; c++){
						value = gettag(cla[cl].comp[c].name); 
						if(value != ""){
							if(clfl == -1) clfl = cl;
							else{
								if(clfl != cl){ alertimp("Can only specify initial populations in a single classification"); return;}
							}
							value = Number(value); if(isNaN(value)){ alertimp("Value must be a number"); return;}
							if(value != Math.floor(value)){ alertimp("Value must be an integer"); return;}
							if(value < 0){alertimp("Value must be non-negative"); return;}
							cla[cl].comp[c].simpopinit = value;
						}
					}
					if(cl == clfl){
						siminit = "Manualpop"; simpopinitset = cl; viewpopinitset = cl;
					}
				}			
				break;
				
			case "siminitpercent":
				clfl = -1;
				for(cl = 0; cl < ncla-1; cl++){
					num = 0; sum = 1;
					for(c = 0; c < cla[cl].ncomp; c++){
						value = gettag(cla[cl].comp[c].name); 
						if(value != ""){
							if(clfl == -1){ clfl = cl; siminit = "Manualpop";}
							else{
								if(clfl != cl){ alertimp("Can only specify percentages in a single classification"); return;}
							}
							num++;
							
							value = Number(value); if(isNaN(value)){ alertimp("Value must be a number"); return;}
							if(value < 0){alertimp("Value must be non-negative"); return;}
							cla[cl].comp[c].simfracinit = value/100;
							sum -= value/100;
						}
						else cst = c;
					}
					if(cl == clfl){
						siminit = "Manualpop";
						if(num != cla[cl].ncomp-1){ alertimp("All but one compartment in '"+cla[cl].name+"' must be set"); return;} 
						
						if(sum < 0){ alertimp("The percentages must add up to less than 100"); return;}
						cla[cl].comp[cst].simfracinit = sum;
					}
				}			
				break;
				
			case "simtimerange": case "inftimerange":
				min = gettag("min"); if(min == ""){ cannotfind(); return;}
				max = gettag("max"); if(max == ""){ cannotfind(); return;}
				min = Number(min); if(isNaN(min)){ alertimp("Minimum must be a number"); return;}
				max = Number(max); if(isNaN(max)){ alertimp("Maximum must be a number"); return;}
				if(min >= max){ alertimp("Minimum must be smaller than maximum"); return;}
				if(type == "simtimerange"){ tsimmin = min; tsimmax = max;}
				else{ tpostmin = min; tpostmax = max;}
				break;
			
			case "runnumber": 
				val = gettag("value"); if(val == ""){ cannotfind(); return;}
				val = Number(val); if(isNaN(val)){ alertimp("Value must be a number"); return;}
				if(val != Math.floor(val)){ alertimp("Value must be an integer"); return;}
				if(val < 1 || val > 20){ alertimp("Value out of range"); return;}
				nchain = val;
				break;
				
			case "simnumber": 
				val = gettag("value"); if(val == ""){ cannotfind(); return;}
				val = Number(val); if(isNaN(val)){ alertimp("Value must be a number"); return;}
				if(val != Math.floor(val)){ alertimp("Value must be an integer"); return;}
				if(val < 1){ alertimp("Value must be positive"); return;}
				if(val == 1) simty = 0;
				else{ simty = 1; simnumber = val;}
				break;
				
			case "runnumber": 
				val = gettag("value"); if(val == ""){ cannotfind(); return;}
				val = Number(val); if(isNaN(val)){ alertimp("Value must be a number"); return;}
				if(val != Math.floor(val)){ alertimp("Value must be an integer"); return;}
				if(val < 1 || val > 20){ alertimp("Value out of range"); return;}
				nchain = val;
				break;
				
			case "paramsampmax": case "eventsampmax":
				val = gettag("value"); if(val == ""){ cannotfind(); return;}
				val = Number(val); if(isNaN(val)){ alertimp("Value must be a number"); return;}
				if(val != Math.floor(val)){ alertimp("Value must be an integer"); return;}
				if(val < 1){ alertimp("Number must be positive"); return;}
				if(type == "paramsampmax") nsampmax = val;
				else nsampevmax = val;
				break;
				
			case "termination":
				type = gettag("type"); if(type == ""){ cannotfind(); return;}
				switch(type){
				case "none": termtype = 1; break;
				case "converged":
					termtype = 0; 
					val = gettag("ess"); 
					if(val != ""){ 
						val = Number(val); if(isNaN(val)){ alertimp("ESS must be a number"); return;}
						if(val < 1){ alertimp("ESS must be positive"); return;}
						ESSmin = val;
					}
					val = gettag("R"); 
					if(val != ""){ 
						val = Number(val); if(isNaN(val)){ alertimp("R must be a number"); return;}
						if(val < 1){ alertimp("R must be positive"); return;}
						GRmax = val;
					}
					break;
				default: alertimp("'"+type+"' not recognised"); return;
				}
				break;
				
			case "description":
				te = gettag("text"); if(te == ""){ cannotfind(); return;}
				descnote = te;
				break;
				
			case "datadesc":
				te = gettag("text"); if(te == ""){ cannotfind(); return;}
				datanote = te; datanoteon = 1;
				break;
				
			case "classdesc":
				te = gettag("text"); if(te == ""){ cannotfind(); return;}
				claname = gettag("classification"); if(te == ""){ cannotfind(); return;}
				cl = 0; while(cl < ncla && cla[cl].name != claname) cl++;
				if(cl == ncla){ alertimp("Cannot find '"+claname+"'"); return;}
				cla[cl].desc = te; cla[cl].descon = 1;
				break;
				
			case "initpopulation":
				table = gettag("table"); if(table == ""){ cannotfind(); return;}
				
				if(loadsiminit(table) == 1){ alertimp(errormsg); return;}
				if(setsimind() == 1){ alertimp(errormsg); return;}
				siminit = "Load";
				showsimindinit = 0;
				break;
				
			case "data":
				type = gettag("type"); if(type == ""){ cannotfind(); return;}
				table = gettag("table"); if(table == ""){ cannotfind(); return;}
				if(loadtable(table) == 1){ alertimp(errormsg); return;}
				addingdata = 3;
			 
				switch(type){
				case "state": datatype = "state"; break;
				case "capture": datatype = "cap"; break;
				case "captureid": datatype = "capid"; break;
				case "capturepd": datatype = "cappd"; break;
				case "transition": datatype = "trans"; break;
				case "source": datatype = "trans"; break;
				case "sink": datatype = "trans"; break;
				case "move": datatype = "move"; break;
				case "sourcemove": datatype = "move"; break;
				case "sinkmove": datatype = "move"; break;
				case "population": datatype = "pop"; break;
				case "derived": datatype = "der"; break;
				case "presence": datatype = "presence"; break;
				default: alertimp("Data type '"+type+"' not recognised"); return;
				}
					
				initdatatemp();
				if(ncol != ncoldefmax){ alertimp("Data table does not have "+ncoldefmax+" columns"); return;} 
					
				switch(type){
				case "state":
					claname = gettag("classification"); if(claname == ""){ cannotfind(); return;}
					cl = 0; while(cl < ncla && cla[cl].name != claname) cl++;
					if(cl == ncla){ alertimp("Cannot find '"+claname+"'"); return;}
					datatemp.cl = cl;
					datapossetup();
					
					datatemp.type = "simple";
					for(n = 0; n < datatemp.pos.length; n++){
						val = datatemp.pos[n];
						
						states = gettag(val); 
						if(states == ""){
							for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == val) break; 
							if(c == cla[cl].ncomp){ alertimp("Require observation model for '"+val+"'"); return;}
						}
						else{
							stlist=states.split(",");
							
							for(c = 0; c < cla[cl].ncomp; c++){
								datatemp.posbinary[n][c] = 0; 
								datatemp.posexpression[n][c] = "0";
							}
							for(nn = 0; nn < stlist.length; nn++){
								for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == stlist[nn]) break;
								if(c == cla[cl].ncomp) break;
								
								datatemp.posbinary[n][c] = 1; 
								datatemp.posexpression[n][c] = "1";
								if(stlist.length == 1) datatemp.posref[n] = c;
							}
							if(nn < stlist.length){
								if(stlist.length != cla[cl].ncomp){ alertimp("Observation model '"+states+"' not understood"); return;}
								for(nn = 0; nn < stlist.length; nn++){
									var sp = stlist[nn].split(":");
									if(sp.length != 2){ alertimp("Observation model '"+states+"' not understood"); return;}
									for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == sp[0]) break;
									if(c == cla[cl].ncomp){ alertimp("Observation model '"+states+"' not understood"); return;}
									datatemp.posexpression[n][c] = sp[1];
								}
								datatemp.type = "expression";
							}
							else{
								if(stlist.length > 1 && datatemp.type != "expression") datatemp.type = "binary";
							}
						}
					}
					break;
				
				case "capture":
					var pdon = gettag("pd");
					switch(pdon){
					case "on":
						datatemp.obspd = "set"; 
						pd = gettag("detectprob");
						if(pd != ""){ 
							datatemp.pdsame = "same";
							datatemp.pd = pd;
						}
						else{
							datatemp.pdsame = "dif";
						}
						break;
						
					case "off":
						datatemp.obspd="all"; 
						break;
					
					default: alertimp("Property 'pd' must be set to 'on' or 'off'"); return;
					}
					break;
					
				case "captureid":
					break;
					
				case "capturepd":
					break;
			
				case "transition": case "move":
					fr = gettag("from"); if(fr == ""){ cannotfind(); return;}
					to = gettag("to"); if(to == ""){ cannotfind(); return;}
					for(cl = 0; cl < ncla; cl++){
						for(tr = 0; tr < cla[cl].ntra; tr++){	
							if((cla[cl].tra[tr].i >= 0 && cla[cl].comp[cla[cl].tra[tr].i].name==fr) &&
							   (cla[cl].tra[tr].f >= 0 && cla[cl].comp[cla[cl].tra[tr].f].name==to)) break; 
						}
						if(tr < cla[cl].ntra) break;
					}
					if(cl == ncla){ alertimp("Could not find transition"); return;}
					
					transcl = cl; transni = fr; transnf = to;
					transty = "trans";

					if(type == "transition"){
						whichind = "all";
						for(cl2 = 0; cl2 < ncla; cl2++){
							clagendata[cl2] = "All";
							if(cl2 != cl){
								val = gettag(cla[cl2].name);
								if(val != ""){
									for(c = 0; c < cla[cl2].ncomp; c++) if(cla[cl2].comp[c].name == val) break;
									if(c == cla[cl].ncomp){ alertimp("Cannot find '"+val+"'"); return;}
									clagendata[cl2] = val; whichind = "sub";
								}
							}
						}
					
						tmin = gettag("min"); if(tmin == ""){ cannotfind(); return;}
						tmax = gettag("max"); if(tmax == ""){ cannotfind(); return;}
						
						tmin = Number(tmin); if(isNaN(tmin)){ alertimp("'min' must be a number"); return;}
						tmax = Number(tmax); if(isNaN(tmax)){ alertimp("'max' must be a number"); return;}
						if(tmin > tmax){ alertimp("'min' must be less than 'max'"); return;}
							
						tgenmin = tmin; tgenmax = tmax;
						
						obspd = "all";
						var pdon = gettag("pd");
						switch(pdon){
						case "on":
							obspd = "set"; 
							datatemp.pd = gettag("detectprob");
							break;
						}
					}
					break;
					
				case "source": case "sourcemove":
					if(type == "source"){
						for(cl2 = 0; cl2 < ncla; cl2++){
							clagendata[cl2] = "All";
							val = gettag(cla[cl2].name);
							if(val != ""){
								for(c = 0; c < cla[cl2].ncomp; c++) if(cla[cl2].comp[c].name == val) break;
								if(c == cla[cl].ncomp){ alertimp("Cannot find '"+val+"'"); return;}
								clagendata[cl2] = val;
							}
						}
						tmin = gettag("min"); if(tmin == ""){ cannotfind(); return;}
						tmax = gettag("max"); if(tmax == ""){ cannotfind(); return;}
					
						tmin = Number(tmin); if(isNaN(tmin)){ alertimp("'min' must be a number"); return;}
						tmax = Number(tmax); if(isNaN(tmax)){ alertimp("'max' must be a number"); return;}
						if(tmin > tmax){ alertimp("'min' must be less than 'max'"); return;}
						tgenmin = tmin; tgenmax = tmax;
					}
					
					for(cl = 0; cl < ncla; cl++){
						for(tr = 0; tr < cla[cl].ntra; tr++){	
							if(cla[cl].tra[tr].i == -2){
								if(type == "source") break;
								if(clagendata[cl] == "All" || clagendata[cl] == cla[cl].comp[cla[cl].tra[tr].f].name) break; 
							}
						}
						if(tr < cla[cl].ntra) break;
					}
					if(cl == ncla){ alertimp("Could not find source transition"); return;}
					
					transcl = -1; transni = ""; transnf = ""; transty = "+";
					break;
					
				case "sink": case "sinkmove":
					if(type == "sink"){
						for(cl2 = 0; cl2 < ncla; cl2++){
							clagendata[cl2] = "All";
							val = gettag(cla[cl2].name);
							if(val != ""){
								for(c = 0; c < cla[cl2].ncomp; c++) if(cla[cl2].comp[c].name == val) break;
								if(c == cla[cl].ncomp){ alertimp("Cannot find '"+val+"'"); return;}
								clagendata[cl2] = val;
							}
						}
						tmin = gettag("min"); if(tmin == ""){ cannotfind(); return;}
						tmax = gettag("max"); if(tmax == ""){ cannotfind(); return;}
					
						tmin = Number(tmin); if(isNaN(tmin)){ alertimp("'min' must be a number"); return;}
						tmax = Number(tmax); if(isNaN(tmax)){ alertimp("'max' must be a number"); return;}
						if(tmin > tmax){ alertimp("'min' must be less than 'max'"); return;}
						tgenmin = tmin; tgenmax = tmax;
					}
					
					for(cl = 0; cl < ncla; cl++){
						for(tr = 0; tr < cla[cl].ntra; tr++){	
							if(cla[cl].tra[tr].f == -2){
								if(type == "sink") break;
								if(clagendata[cl] == "All" || clagendata[cl] == cla[cl].comp[cla[cl].tra[tr].i].name) break; 
							}
						}
						if(tr < cla[cl].ntra) break;
					}
					if(cl == ncla){ alertimp("Could not find sink transition"); return;}
					
					transcl = -1; transni = ""; transnf = ""; transty = "-";
					break;
					
				case "population":
					for(cl2 = 0; cl2 < ncla; cl2++){
						clagendata[cl2] = "All";
						val = gettag(cla[cl2].name);
						if(val != ""){
							for(c = 0; c < cla[cl2].ncomp; c++) if(cla[cl2].comp[c].name == val) break;
							if(c == cla[cl].ncomp){ alertimp("Cannot find '"+val+"'"); return;}
							clagendata[cl2] = val;
						}
					}
					break;
		
				case "derived":
					dergensel = gettag("quantity"); if(dergensel == ""){ cannotfind(); return;}
					
					d = 0; while(d < derive.length && derive[d].name != dergensel) d++; 
					if(d == derive.length){ alertimp("Cannot find derived quantity '"+quantity+"'"); return;}
					break;
				
				case "presence":
					break;
				}
				
				name = gettag("name");
				if(name != ""){
					for(d = 0; d < data.length; d++){ if(d != dataselected && data[d].name == name) break;}
					if(d < data.length){ alertimp("Data name '"+name+"' already used"); return;}
					datatemp.name = name;
				}
				if(ncol != ncoldefmax){ alertimp("The table should have '"+ncoldefmax+"' columns"); return;}

				if(adddata() == 1){ alertimp("There was a problem loading data table '"+name+"':\n"+errormsg); return;}
				break;
		
				
			default: alertimp("Command '"+type+"' not recognised"); return;
			}
			
			for(n = 0; n < tags.length; n++){
				if(tags[n].done != 1){ alertimp("Tag '"+tags[n].name+"' not used"); return;}
			}
		}
	}

	for(cl = 0; cl < ncla; cl++){    // works out any pairs
		if(cltrach[cl] == 1){
			for(tr = 0; tr < cla[cl].ntra-1; tr++){
				for(tr2 = tr+1; tr2 < cla[cl].ntra; tr2++){
					if(cla[cl].tra[tr].i == cla[cl].tra[tr2].f && cla[cl].tra[tr].f == cla[cl].tra[tr2].i){
						if(cla[cl].tra[tr].i >= 0 && cla[cl].tra[tr].f >= 0){
							cla[cl].tra[tr].pair = 1; cla[cl].tra[tr2].pair = 1; 
						}						
					}
				}
			}
		}
	}
	
	for(cl = 0; cl < ncla; cl++){   // Sets up any added classifications
		if(clcompch[cl] == 1 || cltrach[cl] == 1){
			var xmin = large, xmax = -large, ymin = large, ymax = -large, wmax = 0, wmaxname;
			for(c = 0; c < cla[cl].ncomp; c++){
				w = textwidth(cla[cl].comp[c].name,"20px georgia"); if(w > wmax){ wmax = w; wmaxname = cla[cl].comp[c].name;}
				
				x = cla[cl].comp[c].x + cla[cl].comp[c].w/2; y = cla[cl].comp[c].y + cla[cl].comp[c].h/2;
				if(x < xmin) xmin = x; if(x > xmax) xmax = x;
				if(y < ymin) ymin = y; if(y > ymax) ymax = y;
			}
			for(tr = 0; tr < cla[cl].ntra; tr++){
				for(n = 0; n < cla[cl].tra[tr].p.length; n++){
					x = cla[cl].tra[tr].p[n].x; y = cla[cl].tra[tr].p[n].y;
					if(x < xmin) xmin = x; if(x > xmax) xmax = x;
					if(y < ymin) ymin = y; if(y > ymax) ymax = y;
				}					
			}
			
			if(xmin == xmax){ xmin--; xmax++;}
			if(ymin == ymax){ ymin--; ymax++;}

			fac = 0.95*(height-120)/(ymax-ymin);
			fac2 = 0.95*(width-menux-140)/(xmax-xmin); if(fac2 < fac) fac = fac2;
			
			h = 50;
			ratio = (0.3*h+textwidth(wmaxname,Math.floor(h*0.7)+  "px georgia"))/h;
			if(ratio < 1) ratio = 1;
			
			h = large;
			for(c = 0; c < cla[cl].ncomp; c++){
				for(cc = c+1; cc < cla[cl].ncomp; cc++){
					dx = cla[cl].comp[c].x - cla[cl].comp[cc].x; if(dx < 0) dx = -dx;
					dy = cla[cl].comp[c].y - cla[cl].comp[cc].y; if(dy < 0) dy = -dy;
					if(dx/(dy+0.0001) > ratio){ if(h > dx/ratio) h = dx/ratio;}
					else{ if(h > dy) h = dy;}
				}
			}

			h *= fac*0.8;
			if(h > 50) h = 50;
			w = h*ratio;
		
			xmid = (xmax+xmin)/2; ymid = (ymax+ymin)/2;
			for(c = 0; c < cla[cl].ncomp; c++){
				cla[cl].comp[c].x = fac*(cla[cl].comp[c].x-xmid)+(width-menux)/2 - w/2;
				cla[cl].comp[c].y = fac*(cla[cl].comp[c].y-ymid)+height/2-10;
				cla[cl].comp[c].w = w;
				cla[cl].comp[c].h = h;
			}			
		
			for(tr = 0; tr < cla[cl].ntra; tr++){
				for(n = 0; n < cla[cl].tra[tr].p.length; n++){
					cla[cl].tra[tr].p[n].x = fac*(cla[cl].tra[tr].p[n].x-xmid)+(width-menux)/2;
					cla[cl].tra[tr].p[n].y = fac*(cla[cl].tra[tr].p[n].y-ymid)+height/2-10;
				}
			}
			for(tr = 0; tr < cla[cl].ntra; tr++) findpline(cl,tr);
		}
	}

	if(definemodel() ==	1){ alertimp2(errormsg); return;}
	helptype = 208;
}

function loadtable(st)                                    // Imports a table of data
{
	var j, lines=st.split('\n');
	
	row=[];
	for(j = 0; j < lines.length; j++){
		row[j] = lines[j].split('\t');
	}
	nrow = lines.length;
	if(nrow > 0){
		ncol = row[0].length; ncoldef = ncol;
		for(j = 1; j < nrow; j++){
			if(row[j].length != ncol){ errormsg="The rows do not all have the same size"; return 1;}
		}
	}
	else{ errormsg = "Table doe not contain any data"; return 1;}
	return 0;	
}

function gettag(st)                                       // Gets the value of a tag
{
	var i;
	tagst = st;
	for(i = 0; i < tags.length; i++) if(tags[i].name == st){ tags[i].done = 1; return tags[i].value;}
	return "";
}

function cannotfind()                                     // Error massage if a tag cannot be found
{
	alertimp("Cannot find the '"+tagst+"' tag for '"+typest.type+"'",typest.j);
}

