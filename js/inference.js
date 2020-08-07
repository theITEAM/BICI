function inferencebuts()                                  // Intialises the buttons on the inference page
{
	var x, y, dx = width-menux, ddx = 140, xx, yy;

	res = infres;
	
	subna = infpagename[pagesubsub[INFERENCEPAGE][2]];

	cv = resultcv;

	switch(pagesub[INFERENCEPAGE]){
	case 0: databuts(); return;
	
	case 1: 
		addbutton("Prior",menux+tab,30,0,0,-1,TITLEBUT,21,-1);
		priorbuts();
		addbutton("Next",width-105,height-45,90,30,NEXTINFAC,NEXTBUT,0,-1);	
		return;
	}
	
	x = menux+tab; y = 30;
 
	if(warning.length > 0){ addwarning(); return;}
	
	switch(subna){
	case "Start":
		if(infres.result != 0){
			cornx = menux+tab+20; corny = 80;
			addbutton("",cornx,corny,setupwidth,setupheight,CANVASBUT,CANVASBUT,-1,-1);
		
			if(advop == 0){
				addbutton("Setup",x,y,0,0,-1,TITLEBUT,23,-1);
				drawsetup();
			}
			else{
				addbutton("Advanced options",x,y,0,0,-1,TITLEBUT,24,-1);
				drawadvop();
			}
		}
		else{
			if(loading == 0) startloading();
			addbutton("",menux+(width-menux)/2,height/2+60,0,0,-1,PROGRESSBUT,0,-1);
			addbutton("Cancel",width-105,height-45,90,30,STARTAC,ADDDATABUT,0,-1);
		}
		return;
		
	case "Traces":
	case "Prob. Dist.":
		curvevar=[]; curvech=[]; curvelab=[];
		for(vs = 0; vs < res.varsel.length; vs++){
			v = res.varsel[vs];
			
			switch(subna){		
			case "Prob. Dist.":
				if(res.runsel == res.nch){ for(ch = 0; ch < res.nch; ch++){ curvevar.push(v); curvech.push(ch);}}
				else{ curvevar.push(v); curvech.push(res.runsel);}
				break;
				
			case "Traces":
				if(res.runsel3 == res.nch){ for(ch = 0; ch < res.nch; ch++){ curvevar.push(v); curvech.push(ch);}}
				else{ curvevar.push(v); curvech.push(res.runsel3); }
				break;
			}
		}
	
		de = "";
		if(res.varsel.length < 2){
			for(vs = 0; vs < res.varsel.length; vs++){ de += res.vardesc[res.varsel[vs]]; if(vs < res.varsel.length-1) de += " / ";}
		}
		else{
			for(vs = 0; vs < res.varsel.length; vs++){ de += res.varname[res.varsel[vs]]; if(vs < res.varsel.length-1) de += " / ";}
		}
		
		switch(subna){
		case "Traces": addbutton(de,x,30,0,0,-1,TITLEBUT,44,-1); showtraces();break;
		case "Prob. Dist.": addbutton(de,x,30,0,0,-1,TITLEBUT,45,-1); showdistribution(); break;
		}	
		break;

	case "Dependency":
		de = "Dependency of ";
		for(vs = 0; vs < depsel.length; vs++){ de += histoplot[depsel[vs]].name; if(vs < depsel.length-1) de += " / ";}
		de +=" against "+depopsel;
		addbutton(de,x,30,0,0,-1,TITLEBUT,44,-1);
		showhisto(); 
		break;
		
	case "Trans. Dist.":
		curvevar=[];
		curvech=[];
		curvelab=[];
		for(vs = 0; vs < transsel.length; vs++){
			v = transsel[vs];
			if(res.runsel == res.nch){
				for(ch = 0; ch < res.nch; ch++){ curvevar.push(v); curvech.push(ch);}
			}
			else{
				curvevar.push(v); curvech.push(res.runsel); 
			}
		}
	
		de = "";
		if(transsel.length < 2){
			for(vs = 0; vs < transsel.length; vs++) de += transpos[transsel[vs]].desc;
		}
		else{
			de = "Transition durations for: "
			for(vs = 0; vs < transsel.length; vs++){ de += transpos[transsel[vs]].name; if(vs < transsel.length-1) de += ", ";}
		}
		
		addbutton(de,x,30,0,0,-1,TITLEBUT,45,-1);
		showdistribution();
		break;
		
	case "Individuals":
		indplot();
		break;
		
	case "Statistics":
		addbutton("Statistics",x,30,0,0,-1,TITLEBUT,47,-1);
		showstatistics();
		break;
	
	case "Correlation":
		if(res.varselx >= 0 && res.varsely >= 0){ 
			switch(res.scattype){
			case "KDE": addbutton("Posterior distribution",x,30,0,0,-1,TITLEBUT,81,-1); break;
			case "scatter": addbutton("Scatter plot",x,30,0,0,-1,TITLEBUT,46,-1); break;
			}
			
			showscatter();
		}
		else{
			addbutton("Correlation",x,30,0,0,-1,TITLEBUT,49,-1); y += 50;
		
			if(corlist.length == 0){
				addbutton("Select variables from the right-hand menu.",x+18,y,900,0,-1,PARAGRAPHBUT,1,-1);
			}
			else{
				if(corlist.length < 2) addbutton("Select variables from the right-hand menu.",x+18,y,900,0,-1,PARAGRAPHBUT,1,-1);
				else addbutton("Click on matrix squares to show corresponding posterior.",x+18,y,900,0,-1,PARAGRAPHBUT,1,-1);
				y += 30;
				
				addbutton("",x+20,y,corwidth,corheight,CANVASBUT,CANVASBUT,-1,-1);
						
				showcorrelation();
					
				tablexfrac = corwidth/ytot;
				if(tablexfrac < 1) addbutton("",x+20,y+corheight+10,corwidth,13,SLIDEAC,XSLIDEBUT,-1,-1);
					
				tableyfrac = corheight/ytot;
				if(tableyfrac < 1) addbutton("",x+20+corwidth+10,y,13,corheight,SLIDEAC,YSLIDEBUT,-1,-1);
			}
		}
		break;
		
	case "Populations":
		showpopulations(0);
		break;
		
	case "Pos. Pred. Ch.":
		showpopulations(1);
		break;
		
	case "Derived":
		showpopulations(0); 
		break;
	}

	switch(subna){
	case "Traces": case "Prob. Dist.": case "Trans. Dist.": case "Statistics": case "Correlation": case "Dependency":
		rightmenu(width-rwid,75,subna); 
		break;
	case "Derived": 
		rightmenu(width-rwid,75,"der"); 
		break;
	case "Populations": case "Pos. Pred. Ch.": 
		rightmenu(width-rwid,75,"pop");
		break;
	}	
	
	if(infres.result == 1) addbutton("Stop!",width-115,height-45,100,30,STARTAC,UPLOADBUT,0,-1);
}

function addwarning()                                    // Warning page
{
	var x, y, w, wmax;
	
	x = menux+tab; y = 30;
	addbutton("There was a problem...",x,y,0,0,-1,TITLEBUT,-1,-1); y += 50;
		
	addbutton("",x+20,y,0,0,-1,WARNBUT,-1,-1);
	
	wmax = warning.length; if(wmax > 10) wmax = 10;
	for(w = 0; w < wmax; w++){
		addbutton("",x+75,y,0,0,-1,BULLETBUT,1,-1);
		addbutton(warning[w],x+87,y,800,0,-1,PARAGRAPHBUT,1,-1);
		y += 50;
		if(w > 10) break;
	}
	addbutton("Back",width-115,height-45,100,30,ERRORBACKAC,BACKBUT,0,-1);
}

function drawsetup()                                     // Page showing inference setup options
{
	x = 10; y = 5;
		
	addcanbutton("Number of independent runs",x,y,150,0,-1,TEXTBUT,-1,-1); y += 25;
	
	addcanbutton("Running on more than one computational core in parallel speeds up ",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 25;
	addcanbutton("MCMC convergence and facilitates convergence disgnostics.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 25;
	addcanbutton(ncpu+" processing cores are detected on your computer.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 25;
	
	gdropinfo.push({val:nchain, x:width-155, y:130, dx:100, dy:26, style:5, options:[1,2,3,4,5,6,7,8], d:d, i:i, j:j, click:"nchain"});

	y += 25;
			
	addcanbutton("Time range for inference",x,y,150,0,-1,TEXTBUT,-1,-1); y += 25;
	
	if(tpostmin != large) addcanbutton("Begin: "+tpostmin,width-340,y,135,30,MINMAXBIGBUT,MINMAXBIGBUT,-1,4);
	if(tpostmax != -large) addcanbutton("End: "+tpostmax,width-340,y+30,135,30,MINMAXBIGBUT,MINMAXBIGBUT,-1,5);

	addcanbutton("This represents the range in time over which interence is performed.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 25;

	if(datares.tmin == undefined){ addcanbutton("This can only be set after data has been loaded.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 50;}
	else{
		y += 5;
		if(data.length == 0){
			addcanbutton("Note, data must be loaded before inference can start.",x+15,y,150,0,-1,TEXTBUT2,-1,-1);
		}
		else{
			addcanbutton("Note, this range must contain the time period of the data between" 
						,x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 25;
			addcanbutton(tpre(datares.tmin,4)+" and "+tpre(datares.tmax,4)+", but may be extended for future and past posterior predition.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 25;
		}
	}
	
	y += 25;
	te = "Advanced options..."; addcanbutton(te,x,y,textwidth(te,"14px arial"),20,ADVOPAC,LINKBUT,-1,-1);
	
	addbutton("Start",width-105,height-45,90,30,STARTAC,GREENBUT,0,-1);
}

function drawadvop()                                     // Page for advanced setup options
{
	x = 10; y = 5;
		
	addcanbutton("Samples",x,y,150,0,-1,TEXTBUT,-1,-1); y += 25;
			
	addcanbutton(nsampmax,width-360,y-5,145,30,MINMAXBIGBUT,MINMAXBIGBUT,-1,18);
	addcanbutton("Maximum number of parameter samples.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 30;
	
	addcanbutton(nsampevmax,width-360,y-5,145,30,MINMAXBIGBUT,MINMAXBIGBUT,-1,19);
	addcanbutton("Maximum number of event sequence samples.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 25;
		
	y += 25;
		
	addcanbutton("Limit",x,y,150,0,-1,TEXTBUT,-1,-1); y += 30;
	
	addcanbutton("The maximum number of allowable individuals.",x+15,y,150,0,-1,TEXTBUT2,-1,-1);
	addcanbutton(indmaxnumber,width-360,y-8,129,30,MINMAXBIGBUT,MINMAXBIGBUT,0,99); y += 25;
	
	y += 25;
			
	addcanbutton("Termination",x,y,150,0,-1,TEXTBUT,-1,-1); y += 20;
	
	if(termtype == 0){
		addcanbutton("ESS min.: "+ESSmin,width-360,y,145,30,MINMAXBIGBUT,MINMAXBIGBUT,-1,20);
		if(nchain > 1) addcanbutton("R̂ max.: "+GRmax,width-360,y+30,145,30,MINMAXBIGBUT,MINMAXBIGBUT,-1,21);
	}
	
	if(termtype == 2){
		addcanbutton("Iter.: "+itermax,width-360,y,145,30,MINMAXBIGBUT,MINMAXBIGBUT,-1,22);
	}

	y += 10;
	addcanbutton("None",x+15,y,110,22,CANRADIOBUT,CANRADIOBUT,1,CANRADIOTERM);
	addcanbutton("Iterations",x+15+140,y,110,22,CANRADIOBUT,CANRADIOBUT,2,CANRADIOTERM);
	addcanbutton("Converged",x+15+290,y,110,22,CANRADIOBUT,CANRADIOBUT,0,CANRADIOTERM);
	
	y += 25;
	
	switch(termtype){
	case 0:
		addcanbutton("Execution is terminated when the analysis has converged.",x+15,y,150,0,-1,TEXTBUT2,-1,-1);  
		break;
	case 1:
		addcanbutton("Execution continues until stopped by the user.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); 
		break;
	case 2:
		addcanbutton("Execution is terminated after a defined number of iterations.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); 
		break;
	}
	
	addbutton("Back",width-105,height-45,90,30,ADVOPAC2,BACKBUT,0,-1);	
}

function startinference(sim,numch)                        // Starts inference
{
	var index=[], indinfo, sta=[], cl, cl2, tr, ch, st;
	
	loading = 0; if(sim == 1) simres.result = -1; else infres.result = -1;

	stopinference();
	
	if(sim == 0){ if(converttoobs("data") == 0) return;}
	collectvariables();	
		
	if(checkprior(1) == 1) return;
		
	if(page == INFERENCEPAGE){ if(checkprior(0) == 1) return;}
	
	if(checkbeforeinference() == 1) return;

	nchrun = numch;	
	st = '<?xml version="1.0" encoding="UTF-8" ?>\n';
	st += '<BICI version="1.0">\n';
	st += '<model>\n';
		
	if(sim == 0){
		for(p = 0; p < param.length; p++){
			if(param[p].dist == 0){
				name = getparamname(param,p); type = param[p].prior;
				st += '\t<prior parameter="'+name+'" class="'+param[p].classname+'" type="'+type+'"';
				if(typeof param[p].val[0] != 'undefined') st += ' val1="'+param[p].val[0]+'"';
				if(typeof param[p].val[1] != 'undefined')  st += ' val2="'+param[p].val[1]+'"';
				st += '/>\n';
			}
		}
	}
	else{
		for(p = 0; p < paramsim.length; p++){
			val = paramsim[p].sim; name = getparamname(paramsim,p);
			if(isNaN(val)){ changepage(SIMULATEPAGE,0,0); alertp("Parameter "+name+" must be set"); return;}
			
			st += '\t<set parameter="'+name+'" class="'+paramsim[p].classname+'" value="'+val+'"/>\n';
		}
	}

	for(p = 0; p < param.length; p++){
		if(param[p].dist == 1){
			name = getparamname(param,p); type = param[p].prior;
			st += '\t<prior parameter="'+name+'" class="'+param[p].classname+'" type="'+type+'"';
			if(typeof param[p].val[0] != 'undefined')  st += ' val1="'+param[p].val[0]+'" ';
			if(typeof param[p].val[1] != 'undefined')  st += ' val2="'+param[p].val[1]+'"';
			st += '/>\n';
		}
	}

	if(sim == 0){
		for(p = 0; p < paramagesmooth.length; p++){
			if(paramagesmooth[p].type != "None"){
				st += '\t<agesmooth parameter="'+paramagesmooth[p].name+'" type="'+paramagesmooth[p].type+'" val="'+paramagesmooth[p].val+'"/>\n';	
			}
		}
		
		for(p = 0; p < paramtimesmooth.length; p++){
			if(paramtimesmooth[p].type != "None"){
				st += '\t<timesmooth parameter="'+paramtimesmooth[p].name+'" type="'+paramtimesmooth[p].type+'" val="'+paramtimesmooth[p].val+'"/>\n';	
			}
		}
	}
	
	if(simst == 1) rest = simres; else rest = infres; 
	for(d = 0; d < derive.length; d++){  // Expands out derived
		var clst=[], ind=[];
		eqn = eqntidy(derive[d].eq);

		n = derive[d].dep.length;
		if(n > 0){
			for(j = 0; j < n; j++){
				cl = 0; while(cl < ncla && cla[cl].name != derive[d].dep[j]) cl++;
				if(cl == ncla){ alertp("Cannot find "+cla[cl].name+" in derived quantity"); return;}
				clst[j] = cl;
				ind[j] = 0;
			}
			
			do{
				eq = eqn;
				z = 0;
				do{
					while(z < eq.length && eq.substr(z,1) != "{" && eq.substr(z,1) != "[") z++;
					if(z < eq.length){
						if(eq.substr(z,1) == "{") end = "}"; else end = "]";
						
						while(z < eq.length && eq.substr(z,1) != "_" && eq.substr(z,1) != end) z++;
						
						if(eq.substr(z,1) == "_"){
							z++;
							do{
								zi = z; while(z < eq.length && eq.substr(z,1) != "," && eq.substr(z,1) != end) z++;
								if(z == eq.length){ alertp("Brackets do not match up"); return;}
							
								fra = eq.substr(zi,z-zi);
								for(j = 0; j < n; j++){
									if(fra == cla[clst[j]].name){
										rep = cla[clst[j]].comp[ind[j]].name;
										eq = eq.substr(0,zi)+rep+eq.substring(z);
										z = zi + rep.length;
										break;
									}
								}
								if(eq.substr(z,1) == end) break;
								z++;
							}while(1 == 1)
						}
					}
				}while(z < eq.length);
				
				sub = ""; for(j = 0; j < n; j++){ if(j > 0) sub +=","; sub += cla[clst[j]].comp[ind[j]].name;}
			
				st += '\t<derived name="'+ derive[d].name+'_'+sub+'" value="'+eq+'"/>\n'; 
				
				j = 0; do{ flag = 0; ind[j]++; if(ind[j] == cla[clst[j]].ncomp){ ind[j] = 0; j++; flag = 1;}}while(j < n && flag == 1);
			}while(j < n);
		}
		else st += '\t<derived name="'+ derive[d].name+'" value="'+eqn+'"/>\n'; 
	}
	
	for(cl = 0; cl < ncla; cl++){
		st += '\t<classification name="'+cla[cl].name+'" value="';
		for(j = 0; j < cla[cl].ncomp; j++){ st += cla[cl].comp[j].name; if(j < cla[cl].ncomp-1) st += ",";}
		st +=  '">\n';
		
		for(i = 0; i < cla[cl].ntra; i++){
			tr = cla[cl].tra[i];
			
			if(tr.type != "Source" && tr.type != "Sink"){
				st += '\t\t<transition name="trans'+cl+'_'+i+'" from="'+cla[cl].comp[tr.i].name+'" to="'+cla[cl].comp[tr.f].name+'" type="'+tr.type+'" ';
			
				switch(tr.type){
				case "Exponential":
					if(tr.i < 0 || tr.f < 0) alertp("Error code EC17");
	
					st += poisout(tr);
					st += '\t\t</transition>\n';
					break;

				case "Gamma":
					if(tr.i < 0 || tr.f < 0) alertp("Error code EC18");		
					st += '>\n';
					
					var ratelist = getrateforfile(tr.mean,tr.ratetime);
					for(li = 0; li < ratelist.length; li++){
						st += '\t\t\t<value ';
						for(cl2 = 0; cl2 < ncla; cl2++){
							if(ratelist[li].filt[cl2] != "All"){
								st += cla[cl2].name + '="'+ratelist[li].filt[cl2]+'" ';
							}
						}
						st += 'mean="'+ratelist[li].rate+'" shape="'+tr.shape+'"/>\n';
					}
					st += '\t\t</transition>\n';
					break;
					
				case "Weibull":
					if(tr.i < 0 || tr.f < 0) alertp("Error code EC20");
					
					st += '>\n';
					
					var ratelist = getrateforfile(tr.lam,tr.ratetime);
					for(li = 0; li < ratelist.length; li++){
						st += '\t\t\t<value ';
						for(cl2 = 0; cl2 < ncla; cl2++){
							if(ratelist[li].filt[cl2] != "All"){
								st += cla[cl2].name + '="'+ratelist[li].filt[cl2]+'" ';
							}
						}
						st += ' λ="'+ratelist[li].rate+'" k="'+tr.k+'"/>\n';
					}
					st += '\t\t</transition>\n';
					break;
					
				 case "Fixed": case "Grow": case "Set":
					switch(tr.type){
					case "Fixed": st += 'time="'+tr.time+'" '; break;
					case "Grow": st += 'age="'+tr.age+'" '; break;	 
					case "Set": st += 'time="'+tr.time+'" '; break;	 
					}
					st += '/>\n';
					break;
				}
			}
		}
		st += '\t</classification>\n';
	}

	for(cl = 0; cl < ncla; cl++){
		for(i = 0; i < cla[cl].ntra; i++){
			tr = cla[cl].tra[i];
		 
			switch(tr.type){
			case "Source":
				st += '\t<source type="Exponential" '+cla[cl].name+'="'+cla[cl].comp[tr.f].name+'"';
				st += poisout(tr,cl);
				st += '\t</source>\n';
				break;
			
			case "Sink":
				st += '\t<sink type="Exponential" '+cla[cl].name+'="'+cla[cl].comp[tr.i].name+'"';
				st += poisout(tr);
				st += '\t</sink>\n';
				break;
			}
		}
	}
	
	st += '</model>\n';
	
	st += '\n';
	 
	
	var probexp=[];
	
	st += '<data>\n';	

	if(induoflag == 1) st += "\t<unobserved/>\n";
	
	simst = sim;

	var r, rr, d, dd;
	if(sim == 0){
		for(d = 0; d < data.length; d++){
			switch(data[d].variety){
			case "cap":
				for(r = 0; r < data[d].comps.length; r++){
					st += '\t<capture name="'+data[d].capname[r]+'" ';
					if(data[d].comps[r] != "All"){
						sta = data[d].comps[r].split(",");
						for(j = 0; j < sta.length; j++){
							for(cl = 0; cl < ncla; cl++){
								for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == sta[j]) break;
								if(c < cla[cl].ncomp) break;
							}
							if(cl == ncla) alertp("Error code EC21");
							st += cla[cl].name+'="'+cla[cl].comp[c].name+'" ';
						}
					}
			
					st += 'time="'+ data[d].t[r]+'">\n';
					st += '\t\t<probability value="';
				
					if(data[d].obspd == "all") st += "1";
					else{
						if(data[d].pdsame == "same") st +=  eqntidy(data[d].pd);
						else{
							for(dd = 0; dd < data.length; dd++){
								if(data[dd].variety == "cappd"){
									for(rr = 0; rr < data[dd].capname.length; rr++){
										if(data[dd].capname[rr] == data[d].capname[r]){
											st +=  eqntidy(data[dd].pd[rr]);
											break;
										}
									}
									if(rr < data[dd].capname.length) break;									
								}
							}
							if(dd == data.length){
								changepage(INFERENCEPAGE,0,0);			
								alertp("The capture '"+data[d].capname[r]+"' needs information about the detection probability. Please upload a 'capture PD' file.");
								return 1;
							}
						}
					}
					st += '"/>\n';
					st += '\t</capture>\n';
				}
				break;
				
			case "capid":
				for(r = 0; r < data[d].capname.length; r++){  // checks that captures in capture id actual exist
					capna = data[d].capname[r];
					for(dd = 0; dd < data.length; dd++){
						if(data[dd].variety == "cap"){
							for(j = 0; j < data[dd].capname.length; j++){ if(data[dd].capname[j] == capna) break;}
							if(j < data[dd].capname.length) break;
						}
					}
					if(dd == data.length){
						alertp("In '"+data[d].name+"' individual '"+data[d].id[r]+"' is captured in '"+capna+"', but this does not exist.");
						return 1;
					}
				}
				break;
				
			case "trans":
				switch(data[d].transty){
				case "+": st += '\t<capev name="capev'+d+'" type="source" '; break;
				case "-": st += '\t<capev name="capev'+d+'" type="sink" '; break;
				case "trans": st += '\t<capev name="capev'+d+'" type="trans" from="'+data[d].transi+'" to="'+data[d].transf+'"'; break;
				}
				if(data[d].whichind != undefined && data[d].whichind != "all"){
					for(cl = 0; cl < ncla; cl++){
						if(data[d].filt[cl] != "All") st += ' '+cla[cl].name+'="'+data[d].filt[cl]+'"';
					}
				}
				st += ' tmin="'+data[d].tmin+'" tmax="'+data[d].tmax+'"';
				if(data[d].obspd != undefined && data[d].obspd != "all"){
					st += ' pd="'+ eqntidy(data[d].pd)+'"';
				}
				else st += ' pd="1"';
				st += "/>\n";
				break;
				
			case "move": break;
			
			case "pop":
				for(r = 0; r < data[d].t.length; r++){
					st += '\t<population time="'+data[d].t[r]+'" value="'+data[d].val[r]+'" sd="'+data[d].sd[r]+'" ';
					for(cl = 0; cl < ncla; cl++) if(data[d].popcl[cl] != "All") st += cla[cl].name+'="'+data[d].popcl[cl]+'" ';
					st += "/>\n";
				}
				break;
				
			case "der":
				for(r = 0; r < data[d].t.length; r++){
					st += '\t<derived name="';
					st +=  data[d].name;
					if(data[d].dep.length > 0){
						st += "_";
						for(k = 0; k < data[d].dep.length; k++){ if(k != 0) st += ","; st += data[d].depval[r][k];}
					}
					
					st += '" time="'+data[d].t[r]+'" value="'+data[d].val[r]+'" sd="'+data[d].sd[r]+'" ';
					st += "/>\n";
				}
				break;
			}
		}
	}

	if(sim == 1) indinfo = indsim; else indinfo = inddata;

	for(i = 0; i < indinfo.nindtotal; i++){
		var indi = copy(indinfo.ind[i]);

		for(cl = 0; cl < ncla; cl++){
			index[cl] = 0;
			e = 0;
			while(e < indi.cl[cl].ev.length){
				ev = indi.cl[cl].ev[e];
				if(ev.variety == "transtmin" || ev.variety == "transtmax") indi.cl[cl].ev.splice(e,1);
				else e++;				
			}
		}	
		do{
			tnext = large;
			for(cl = 0; cl < ncla; cl++){
				e = index[cl];
				if(e < indi.cl[cl].ev.length){
					ev = indi.cl[cl].ev[e];
					t = ev.t;
					if(t < tnext || (t == tnext && ev.variety != "state")){ clnext = cl; tnext = t; ty = ev.variety;}	
				}
			}
			if(tnext == large) break;
	
			switch(ty){
			case "state": case "presence":
				st += '\t<observation id="'+indi.id+'" time="'+tnext+'" ';
				if(sim == 0){
					var capna = "";
					
					for(dd = 0; dd < data.length; dd++){  // First looks for "capture ID" file
						if(data[dd].variety == "capid"){
							for(r = 0; r < data[dd].id.length; r++){
								if(data[dd].id[r] == indi.id){
									for(ddd = 0; ddd < data.length; ddd++){ 
										if(data[ddd].variety == "cap"){
											for(rr = 0; rr < data[ddd].capname.length; rr++){
												if(data[ddd].capname[rr] == data[dd].capname[r]){
													if(data[ddd].t[rr] == tnext){
														capna = data[dd].capname[r]; break;
													}
												}	
											}
										}	
										if(capna != "") break;
									}
								}	
								if(capna != "") break;								
							}
						}
						if(capna != "") break;
					}
			
					if(capna == ""){  // Looks to directly find capture
						var cappos=[];						
						for(dd = 0; dd < data.length; dd++){ 
							if(data[dd].variety == "cap"){
								for(r = 0; r < data[dd].capname.length; r++){
									if(data[dd].t[r] == tnext){
										for(j = 0; j < cappos.length; j++){
											if(cappos[j] == data[dd].capname[r]) break;
										}
										if(j == cappos.length) cappos.push(data[dd].capname[r]);
									}
								}
							}
						}
						
						switch(cappos.length){
						case 0: break;
						case 1: capna = cappos[0]; break;
						default:
							changepage(INFERENCEPAGE,0,0);				
							alertp("Individual '"+indi.id+"' at time '"+tnext+"' can be detected in '"+cappos[0]+"' or '"+cappos[1]+"'. Please upload a 'capture ID' file to clarify which individuals are detected in which captures.");
							return 1;
						}
					}
					
					if(capna != ""){
						st += 'capture="'+capna+'" ';
					}
				}
		
				if(ty == "presence"){ index[clnext]++; st +='/>\n'; break;}
				
				for(cl = 0; cl < ncla; cl++){
					st += cla[cl].name+'="';
					
					flag = 0;
					e = index[cl];
					if(e < indi.cl[cl].ev.length){
						ev = indi.cl[cl].ev[e];
						t = ev.t;
						if(t == tnext){
							flag = 1;
							if(ev.variety != "state"){ alert("Cannot have state measurements at the same time as other measurements."); return;}
							switch(ev.variety){
							case "state":
								singrs = "";
								for(j = 0; j < cla[cl].ncomp; j++){
									if(ev.probexp[j] == "1"){ if(singrs == "") singrs = cla[cl].comp[j].name; else break;}
									else{
										if(ev.probexp[j] != "0") break;
									}
								}
							
								if(j == cla[cl].ncomp && singrs != "") st += singrs;
								else{
									for(j = 0; j < cla[cl].ncomp; j++){ st += cla[cl].comp[j].name+':'+ev.probexp[j]; if(j < cla[cl].ncomp-1) st += "|";}
								}
								break;
							}
							index[cl]++;
						}
					}
					
					if(flag == 0){
						for(j = 0; j < cla[cl].ncomp; j++){ st += cla[cl].comp[j].name+":1"; if(j < cla[cl].ncomp-1) st += "|";}
					}
					st += '" ';
				}
				st += '/>\n';
				break;
			
			case "trans":
				cl = clnext;
				ev = indi.cl[cl].ev[index[cl]];
				
				st += '\t<transition id="'+indi.id+'" time="'+tnext+'" capev="capev'+ev.obsdata+'" from="'+ev.transi+'" to="'+ev.transf+'"/>\n'
	
				index[cl]++;
				break;
				
			case "move":
				cl = clnext;
				ev = indi.cl[cl].ev[index[cl]];
				switch(ev.transty){
				case "+": st += '\t<move id="'+indi.id+'" time="'+tnext+'" type="source"/>\n'; break;
				case "-": st += '\t<move id="'+indi.id+'" time="'+tnext+'" type="sink"/>\n'; break;
				case "trans": st += '\t<move id="'+indi.id+'" time="'+tnext+'" type="trans" from="'+ev.transi+'" to="'+ev.transf+'"/>\n'; break;
				}
				index[cl]++;
				break;
				
			case "transtmin": case "transtmax": alertp("Error code EC22"); break;
			}
		}while(1 == 1);
	}	
	st += '</data>\n';		

	if(sim == 0){
		st += '<inference nsamp="';
		if(termtype == 2) st += itermax;
		else st += "20000001";
		st += '" ';
	}
	else{
		st += '<simulation nsamp="';
		if(simty == 0) st += 1;
		else st += simnumber;
		st += '" ';
	}
	
	if(sim == 1) st += 'tmin="'+tsimmin+'" tmax="'+tsimmax+'"';	
	else st += 'tmin="'+tpostmin+'" tmax="'+tpostmax+'"';	
	
	st +=' indmax="'+indmaxnumber+'"/>\n';		
	
	st += '</BICI>\n';

	vcalc = 0;
	
	chooseaxis = -1;
	
	tableyfr = 0;
	
	startloading();

	var indlist=[]; for(i = 0; i < indinfo.nindtotal; i++) indlist.push(indinfo.ind[i].id);
		
	simsta = 0;
	
	if(sim == 1){
		if(simstarttext.substr(0,5) == "Rerun"){
			simres.result = 0; simres.burninev = 0; simres.evthin = 1; simres.nindmax = 0; 
			simres.burnin = 0; simres.sampthin = 1;				
		}		
		else{
			simres={ result:0, burninev:0, evthin:1, nindmax:0, indlist:indlist, 
					tmin:tsimmin, tmax:tsimmax, agecl:ncla-2, settimecl:ncla-1, 
					nch:nchrun, ch:[], DX:0, nderive:0, derive:[],
					nvar:0, vartype:[], varname:[], vardesc:[], varset:0, filter:[], filtervar:[],
					burnin:0, sampthin:1, filt:0, popfilt:[], popfiltcheck:[],
					
					clasel:0, varsel:[], dersel:[], popview:"Graph", indview:"Timeline", varselx:-1, varsely:-1, indfilt:0, sampfilt:"All", popshow:[],
					nrunlab:0, runlab:[], runsel:0, nrunlab2:0, runlab2:[], runsel2:0, nrunlab3:0, runlab3:[], runsel3:0, nrunpoplab:0, runpoplab:[], runpopsel:0
					};
			resset(simres);
		}
		
		for(ch = 0; ch < nchrun; ch++){
			simres.ch[ch] = { nsamp:0, nsampread:0, varval:[], nsampev:0, sampev:[], nsampevread:0, derivepl:[]};
		}	
	}
	else{
		infres={ result:0, globaldiagnostics:[], diagnostics:[], inddata:copy(inddata),
		        burninev:0, evthin:1, nindmax:0, indlist:indlist, tdatamin:datares.tmin, 
				tdatamax:datares.tmax,  tmin:tpostmin, tmax:tpostmax, agecl:ncla-2, 
				settimecl:ncla-1, nch:nchrun, ch:[], DX:0, nderive:0, derive:[], param:copy(param),
				nvar:0, vartype:[], varname:[], vardesc:[], varset:0, filter:[], filtervar:[], 
				nsamppl:100, burnin:0, sampthin:1, filt:0, popfilt:[], popfiltcheck:[],
				
				clasel:0, varcheck:[], varsel:[], dersel:[], popview:"Graph", indview:"Timeline", varselx:-1, varsely:-1, scattype:"KDE", scatnum:5000, indfilt:0, sampfilt:"All", popshow:[],
				nrunlab:0, runlab:[], runsel:0, nrunlab2:0, runlab2:[], runsel2:0, nrunlab3:0, runlab3:[], runsel3:0, nrunpoplab:0, runpoplab:[], runpopsel:0
				};
		 
		for(ch = 0; ch < nchrun; ch++){
			infres.ch[ch] = { nsamp:0, nsampread:0, varval:[], nsampev:0, sampev:[], nsampevread:0, derivepl:[]};
		}
		resset(infres);
	
		ppcres={ burninev:0, evthin:1, nindmax:0,  tmin:tpostmin, tmax:tpostmax, 
				agecl:ncla-2, settimecl:ncla-1, nch:nchrun, ch:[],
				filt:0, popfilt:[], popfiltcheck:[],
				
				clasel:0, varsel:[], dersel:[], popview:"Graph", indview:"Timeline", varselx:-1, varsely:-1, indfilt:0, sampfilt:"All", popshow:[],
				nrunlab:0, runlab:[], runsel:0, nrunlab2:0, runlab2:[], runsel2:0, nrunlab3:0, runlab3:[], runsel3:0, nrunpoplab:0, runpoplab:[], runpopsel:0
				};
			
		for(ch = 0; ch < nchrun; ch++){ ppcres.ch[ch] = {nsampev:0, sampev:[], nsampevread:0, derivepl:[]}};	
		resset(ppcres);
		
		transplotinit();
	}
	
	for(ch = 0; ch < nchrun; ch++){ exe[ch] = 0; leftover[ch] = "";}	
	percentrun = 0; percentload = 0;
	buttoninit();
	
	if(sim == 0) runtype = "inf"; else runtype = "sim";
 
	st = st.replace(/{all}/g, "{}");
	 
	var fs = require('fs');
	fs.writeFile("init.xml",st, function(err) {
		if(err) { alertp("There was a problem: "+err); return;}
		setTimeout(function(){ startspawn(0);}, 100);
	}); 
}

function resset(res)                                       // Sets initial properties of result objects
{
	var cl, fl;

	if(cla[0].ncomp > 10) res.popview = "Model"; 

	res.varsel.push(0);
	res.dersel.push(0);
		
	res.nrunlab = nchrun+2; for(c = 0; c < nchrun; c++) res.runlab[c] = "Run "+(c+1); res.runlab[nchrun] = "All Runs"; res.runlab[nchrun+1] = "Combine"; 
	res.nrunlab2 = nchrun+1; for(c = 0; c < nchrun; c++) res.runlab2[c] = "Run "+(c+1); res.runlab2[nchrun] = "Combine"; 
	res.nrunlab3 = nchrun+1; for(c = 0; c < nchrun; c++) res.runlab3[c] = "Run "+(c+1); res.runlab3[nchrun] = "All Runs"; 
	res.nrunpoplab = nchrun+1; for(c = 0; c < nchrun; c++) res.runpoplab[c] = "Run "+(c+1); res.runpoplab[nchrun] = "Combine"; 

	if(nchrun == 1){ res.runsel = 0; res.runsel2 = 0; res.runsel3 = 0; res.runpopsel = 0;} 
	else{ res.runsel = nchrun+1; res.runsel2 = nchrun; res.runsel3 = nchrun; res.runpopsel = nchrun;}

	for(cl = 0; cl < ncla; cl++){ 
		res.popshow[cl]=[]; for(c = 0; c < cla[cl].ncomp; c++){ if(c < 10) res.popshow[cl][c] = 1; else res.popshow[cl][c] = 0;}
	}

	for(cl = 0; cl < ncla; cl++){ res.popfilt[cl] = -1; res.popfiltcheck[cl]=[];}
}

function poisout(tr,sourcecl)                                    // Adds Poisson transition to init.xml file
{
	var st="", fac, fl, li, cl2, j, k;
	
	st += '>\n';
	var ratelist = getrateforfile(tr.rate,tr.ratetime);
	for(li = 0; li < ratelist.length; li++){
		st += '\t\t\t<value ';
		fac = 1;
		for(cl2 = 0; cl2 < ncla; cl2++){
			if(ratelist[li].filt[cl2] != "All") st += cla[cl2].name + '="'+ratelist[li].filt[cl2]+'" ';
			else{
				if(cl2 != sourcecl){
					// checks to see if classification name appears in parameters or population
					// if it does this implies that the expression is used multiple times
					// and so doesn't require the factor below
					
					fl = 0;
					if(sourcecl != undefined){
						for(j = 0; j < tr.ratevardep.length; j++) if( tr.ratevardep[j] == cla[cl2].name) fl = 1;
						if(fl == 0){
							var spl = ratelist[li].rate.split("{");
							for(j = 1; j < spl.length; j++){
								var spl2 = spl[j].split("}");
								if(spl2.length != 2) alert("pr");
								var spl3 = spl2[0].split(",");
								for(k = 0; k < spl3.length; k++){ if(spl3[k] == cla[cl2].name){ fl = 1; break;}}
							}
						}
					}	
					if(fl == 0) fac *= cla[cl2].ncomp;
				}
			}
		}
		if(sourcecl != undefined && fac != 1) st += 'rate="'+(1/fac)+'*'+ratelist[li].rate+'"/>\n';
		else st += 'rate="'+ratelist[li].rate+'"/>\n';
	}	
	return st;
}

function startspawn(ch)                                    // Starts execution of C++ code
{		
	if(runtype == "") return;

	var num = Math.floor(Math.random()*1000); 
	if(runtype == "inf") num = ch;

	switch(ver){
	case "windows": child[ch] = spawn('Execute/a.exe',["init.xml",num]); break;
	case "linux": child[ch] = spawn('Execute/a.out',["init.xml",num]); break;
	case "mac": child[ch] = spawn('Execute/a.out',["/tmp/init.xml",num]); break;
	}

	exe[ch] = 1;

	percentload = Math.floor((100*(ch+1))/nchrun);
	percent = Math.floor(0.5*percentload+0.5*percentrun);
	loadingfunc(0);

	funct(child[ch],ch);

	cv.clearRect (menux+(width-menux)/2-100,height/2+60-20,200,40);
	if(ch < nchrun-1){
		setTimeout(function(){ startspawn(ch+1);}, 100);
	}
	else{
		if(simst == 0 && termtype == 0){
			timestart = (new Date()).getTime()
			convid = setTimeout(function(){ convergencetestinit();},5000);
		}
	}	
}
	
function stopinference()                                    // Used to stop inference
{
	var ch;
	if(runtype != ""){
		runtype = "";
		for(ch = 0; ch < nchrun; ch++){ if(exe[ch] == 1){ exe[ch] = 0; child[ch].kill();}}
		running = 0;
		if(intervalid != undefined){ clearInterval(intervalid); intervalid = undefined;}
		if(convid != undefined){ clearTimeout(convid); convid = undefined;}
	}
}

function funct(chi,ch)                                      // Gathers output of C++ file
{
	chi.stdout.on('data', function (data) {
		var s, st, li, v, n, j, off, ch2, bu, min, Lmin, val, av, ans, dT, fi;
		
		if(runtype == "") return;
		
		st = leftover[ch]+data;
		var lines = st.split('\n');
		leftover[ch] = lines[lines.length-1];

		for(li = 0; li < lines.length-1; li++){
			st = lines[li];
			tabs = st.split("|");
			switch(tabs[0]){
			case '1': case 'c':                // Variable names
				switch(tabs[0]){
				case '1': rest = infres; break;
				case 'c': rest = simres; break;
				}

				if(rest.varset == 0){
					rest.nvar = Math.floor((tabs.length-1)/3);
					for(v = 0; v < rest.nvar; v++){
						if(tabs[0] == '1') rest.varcheck[v]=0;
						rest.vartype[v] = tabs[3*v+1]; rest.varname[v] = tabs[3*v+2]; rest.vardesc[v] = tabs[3*v+3];
						
						fi = 0; while(fi < rest.filter.length && rest.filter[fi] != rest.vartype[v]) fi++;
						if(fi == rest.filter.length){
							rest.filter.push(rest.vartype[v]);
							rest.filtervar[fi]=[];
						}
						
						rest.filtervar[fi].push(v);
					}
					histoinit(rest); 
					rest.varset = 1;		
				}
				for(v = 0; v < rest.nvar; v++) rest.ch[ch].varval[v]=[];
				break;
			
			case '0': case 'd':               // Variable values
				switch(tabs[0]){
				case '0': rest = infres; break;
				case 'd': rest = simres; break;
				}
				
				if(rest.ch[ch].nsampread%rest.sampthin == 0){
					var stot = 0; for(ch2 = 0; ch2 < rest.nch; ch2++) stot += rest.ch[ch2].nsamp;  
					if(stot >= nsampmax){
						for(ch2 = 0; ch2 < rest.nch; ch2++){
							nsa = Math.floor(rest.ch[ch2].nsamp/2);
							for(s = 0; s < nsa; s++){
								for(v = 0; v < rest.nvar; v++) rest.ch[ch2].varval[v][s] = rest.ch[ch2].varval[v][2*s];	
							}
							for(v = 0; v < rest.nvar; v++) rest.ch[ch2].varval[v].length = nsa;
							rest.ch[ch2].nsamp = nsa;
						}	
						rest.burnin = Math.floor(rest.burnin/2);	
						rest.nsamppl = Math.floor(rest.nsamppl/2);					
						rest.sampthin *= 2;
					}
					s = rest.ch[ch].nsamp;
					for(v = 0; v < rest.nvar; v++) rest.ch[ch].varval[v][s] = Number(tabs[v+1]);	
					
					rest.ch[ch].nsamp++;
					bu = large*large; for(ch2 = 0; ch2 < rest.nch; ch2++) if(rest.ch[ch2].nsamp < bu) bu = rest.ch[ch2].nsamp;
					bu = Math.floor(bu/5);
					if(bu > rest.burnin) rest.burnin = bu;
				}
				rest.ch[ch].nsampread++;
				break;
				
			case '3':                 // Loading percentage
				if(simsta == 1){
					per = Number(tabs[1]); if(per > percent) percent = per; 
				}
				else{
					if(ch == nchrun-1){
						per = Number(tabs[1]); if(per > percentrun) percentrun = per; percent = Math.floor(0.5*percentload+0.5*percentrun);
					}
				}
				break;
				
			case '4':                 // List of individual ids
				if(ch == 0){
					nid = tabs.length-1;
					for(i = 0; i < nid; i++) id[i] = tabs[i+1];
				}
				break;
			
			case '5': case 'a': case 'b': // Event data
				switch(tabs[0]){
				case '5': rest = infres; break;
				case 'a': rest = ppcres; break;
				case 'b': rest = simres; break;
				}
				
				if(process.memoryUsage().rss/(1024*1024) > 1400 && runtype != ""){
					rest.result = 2;
					alertp("Unfortunately SIRE has run out of available memory. Consider reducing the maximum number of samples and starting again.");
					stopinference(); 
					changepage(INFERNECEPAGE,2,0);
				}
				else{		
					if(rest.ch[ch].nsampevread%rest.evthin == 0){
						var stot = 0; for(ch2 = 0; ch2 < rest.nch; ch2++) stot += rest.ch[ch2].nsampev; 
						if(stot >= nsampevmax){  // Thins samples
							for(ch2 = 0; ch2 < rest.nch; ch2++){
								var nsa = Math.floor(rest.ch[ch2].nsampev/2);
								for(s = 0; s < nsa; s++) rest.ch[ch2].sampev[s] = rest.ch[ch2].sampev[2*s];
								rest.ch[ch2].sampev.length = nsa;
								rest.ch[ch2].nsampev = nsa;
							}
							rest.burninev = Math.floor(rest.burninev/2);	
							s = rest.ch[ch].nsampev; 
							rest.evthin *= 2;
						}
						
						z = 1;
						s = rest.ch[ch].nsampev;
						nindto = parseInt(tabs[z]); z++;
						if(nindto > rest.nindmax) rest.nindmax = nindto;

						rest.ch[ch].sampev[s]={nind:nindto, ind:[]}; 
						ressa = rest.ch[ch].sampev[s];
						for(i = 0; i < nindto; i++){
							tbir = Number(tabs[z]); z++;
							nev = parseInt(tabs[z]); z++;
							ressa.ind[i]={tbirth:tbir, nev:nev, evc:[], evt:[]};
							ressai = ressa.ind[i];
							for(k = 0; k < nev; k++){
								ressai.evc[k] = (parseInt(tabs[z])); z++;
								val = Number(tabs[z]); z++;
								ressai.evt[k] = val;
							}
						}
						rest.ch[ch].nsampev++;
						bu = large; 
						for(ch2 = 0; ch2 < rest.nch; ch2++) if(rest.ch[ch2].nsampev < bu) bu = rest.ch[ch2].nsampev;
						bu = Math.floor(bu/5);
						if(bu > rest.burninev) rest.burninev = bu;
					}
				}
				rest.ch[ch].nsampevread++;
				break;
				
			case '6':            // Chain diagnostics
				infres.diagnostics[ch] = st;
				break;
				
			case '7':            // Global diagnostics
				infres.globaldiagnostics[ch] = st;
				break;
				
			case '8':            // Start simulating
				simsta = 1; percent = 0;
				break;
			
		
			case 'e':            // Warning
				warning = st.substr(2,st.length-3).trim().split("*"); if(warning[warning.length-1] == "") warning.pop();
				if(simst == 1){ simres.result = -1; changepage(SIMULATIONPAGE,2,0);}
				else{ infres.result = -1; changepage(INFERENCEPAGE,2,0);}
				stopinference();
				loading = 0;				
				break;
			
			case 'f':            // Major warning
				alertp("Warning: "+st);
				break;
			
			case 'g':            // Derived variables
				if(simst == 1) rest = simres; else rest = infres; 

				rest.DX = Number(tabs[1]);
				rest.nderive = Number(tabs[2]);
				for(d = 0; d < rest.nderive; d++){
					rest.derive.push(tabs[3+d]);
					rest.ch[ch].derivepl[d]=[];
				}
				break;
				
			case 'h':            // Derived variables
				if(simst == 1) rest = simres; else rest = infres; 
				z=1;
				for(d = 0; d < rest.nderive; d++){
					for(i = 0; i < rest.DX; i++){
						rest.ch[ch].derivepl[d][i] = {t:tabs[z], min:tabs[z+1], av:tabs[z+2], max:tabs[z+3]};
						z += 4;
					}
				}
				break;
			
 			case 'z':            // Comment
				pr("Comment: "+st);
				break;
				
			default:              // Unknown output
				alertp("Unexpected program output:"+st);
				break;
			}
		}

		switch(runtype){
		case "inf":
			if(infres.result == 0 && infres.burnin > 0){
				for(ch2 = 0; ch2 < infres.nch; ch2++) if(infres.ch[ch2].nsampev == 0) break;
				if(ch2 == infres.nch){		
					loading = 0; infres.result = 1; startrunning();
					buttoninit(); buttoninit();
					pagesubsub[INFERENCEPAGE][2] = 1; 
					changepage(INFERENCEPAGE,2,-1);
					if(intervalid == undefined){
						intervalid = setInterval(function(){   // Ensures that plots are updated
							if(page == INFERENCEPAGE && pagesub[page] == 2 && pagesubsub[page][2] == 1) buttoninit();				
						}, 2000);
					}						
				}
			}
			break;
		
		case "sim":		
			if(simres.result == 0 && (simty == 0 && simres.ch[0].nsampev == 1) || (simty == 1 && simres.ch[0].nsampev == simnumber)){ 
				loading = 0; simres.result = 1; gentype = -1;
				switch(simstarttext){
				case "Rerun2": changepage(SIMULATEPAGE,2,2); break;
				default: changepage(SIMULATEPAGE,2,1); break;
				}
			}
		}
	});

	chi.stderr.on('data', function (data) {
		alertp('Error: ' + data);
	});

	chi.on('close', function (code) {
		if(runtype != ""){
			var fi = 1;
			if(termtype == 2){
				fi = 0;
				exe[ch] = 0;
				for(ch = 0; ch < nchrun; ch++) if(exe[ch] == 1) break;
				if(ch == nchrun) fi = 1;
			}
			
			if(fi == 1){
				switch(runtype){
				case "inf":
					if(infres.result == 1) infres.result = 2;
					if(infres.result == 0){ infres.result = -1; helptype = 109;}
					break;
					
				case "sim":
					if(simres.result == 1) simres.result = 2; 
					if(simres.result == 0){ simres.result = -1; helptype = 109;}
					break;
				}
				stopinference();
				
				buttoninit();
				loading = 0;
			}
		}
	});
}

function checkbeforeinference()                        	  // Checks model and data before inference is started
{
	var cl, i, tr, cl2, j;
	
	if(time.length > 0){
		switch(page){
		case SIMULATEPAGE:
			if(time[0] < tsimmin){
				alertp("The time transition cannot be before the simulation time.");
				changepage(MODELPAGE,ncla-1,-1); return 1;
			}
			if(time[time.length-1] > tsimmax){
				alertp("The time transition cannot be after the simulation time.");
				changepage(MODELPAGE,ncla-1,-1); return 1;
			}
			break;
		case INFERENCEPAGE:
			if(time[0] < tpostmin){
				alertp("The time transition cannot be before the inference time.");
				changepage(MODELPAGE,ncla-1,-1); return 1;
			}
			if(time[time.length-1] > tpostmax){
				alertp("The time transition cannot be after the inference time.");
				changepage(MODELPAGE,ncla-1,-1); return 1;
			}
			break;
		}
	}
	
	if(page == INFERENCEPAGE){
		for(i = 0; i < inddata.ind.length; i++){
			a = inddata.ind[i];
			for(cl = 0; cl < a.cl.length; cl++){
				b = a.cl[cl];
				for(e = 0; e < b.ev.length; e++){
					if(b.ev[e].col){
						if(b.ev[e].col.length == 0){
							alertp("Zero probability for data");
							selectdataev(i,inddata.ind.length,cl,e);
							return 1;
						}
					}						
				}
			}
		}
	}
	
	return 0;
}