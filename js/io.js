function save(withres)                                    // Saves a BICI file
{
	st = "";
	
	fileToLoadlist=[];
	
	st += "modelsetup:"+JSON.stringify(modelsetup)+"\n";

	st += "cla:"+JSON.stringify(cla)+"\n";
	st += "ncla:"+JSON.stringify(ncla)+"\n";
	
	st += "age:"+JSON.stringify(age)+"\n";
	st += "time:"+JSON.stringify(time)+"\n";
	
	st += "param:"+JSON.stringify(param)+"\n";
	st += "paramsim:"+JSON.stringify(paramsim)+"\n";

	st += "paramagesmooth:"+JSON.stringify(paramagesmooth)+"\n";
	st += "paramtimesmooth:"+JSON.stringify(paramtimesmooth)+"\n";
	
 	st += "derive:"+JSON.stringify(derive)+"\n";
 	
	st += "termtype:"+JSON.stringify(termtype)+"\n";
	
	st += "tsimmin:"+JSON.stringify(tsimmin)+"\n";
	st += "tsimmax:"+JSON.stringify(tsimmax)+"\n";

	st += "tpostmin:"+JSON.stringify(tpostmin)+"\n";
	st += "tpostmax:"+JSON.stringify(tpostmax)+"\n";

	st += "simpopinitset:"+JSON.stringify(simpopinitset)+"\n";
	st += "simpopinitlist:"+JSON.stringify(simpopinitlist)+"\n";
	st += "simty:"+JSON.stringify(simty)+"\n";
	//simindinit = [];
	st += "simindinit:"+JSON.stringify(simindinit)+"\n";
	st += "siminit:"+JSON.stringify(siminit)+"\n";

	st += "clagendata:"+JSON.stringify(clagendata)+"\n";
	st += "detectprob:"+JSON.stringify(detectprob)+"\n";
	st += "datagenttype:"+JSON.stringify(datagenttype)+"\n";
	st += "tgenmin:"+JSON.stringify(tgenmin)+"\n";
	st += "tgenmax:"+JSON.stringify(tgenmax)+"\n";
	st += "dtgen:"+JSON.stringify(dtgen)+"\n";
	st += "tgenuserdef:"+JSON.stringify(tgenuserdef)+"\n";

	st += "data:"+JSON.stringify(data)+"\n";
	
	st += "nchain:"+JSON.stringify(nchain)+"\n";

	st += "indshow:"+JSON.stringify(indshow)+"\n";
	
	st += "claselect:"+JSON.stringify(claselect)+"\n";
	st += "transpos:"+JSON.stringify(transpos)+"\n";
	st += "transmat:"+JSON.stringify(transmat)+"\n";
	st += "transsel:"+JSON.stringify(transsel)+"\n";
	st += "transindfilt:"+JSON.stringify(transindfilt)+"\n";
	st += "depsel:"+JSON.stringify(depsel)+"\n";
	st += "histoplot:"+JSON.stringify(histoplot)+"\n";
	st += "depop:"+JSON.stringify(depop)+"\n";
	st += "depopsel:"+JSON.stringify(depopsel)+"\n"; 
	st += "distCI:"+JSON.stringify(distCI)+"\n"; 
	st += "distmean:"+JSON.stringify(distmean)+"\n"; 
	
	if(withres == 1) st += "indsim:"+JSON.stringify(indsim)+"\n";
	st += "inddata:"+JSON.stringify(inddata)+"\n";

	st += "simnumber:"+JSON.stringify(simnumber)+"\n";
	st += "indmaxnumber:"+JSON.stringify(indmaxnumber)+"\n";
	
	if(withres == 1){
		st += "simres:"+JSON.stringify(simres)+"\n";
		st += "infres:"+JSON.stringify(infres)+"\n";
		st += "ppcres:"+JSON.stringify(ppcres)+"\n";
	}
	st += "datares:"+JSON.stringify(datares)+"\n";
	st += "data:"+JSON.stringify(data)+"\n";
	st += "simdata:"+JSON.stringify(simdata)+"\n";
	st += "radlab:"+JSON.stringify(radlab)+"\n";
	st += "radstyle:"+JSON.stringify(radstyle)+"\n";
	st += "radstyle2:"+JSON.stringify(radstyle2)+"\n";
	st += "radnamlab:"+JSON.stringify(radnamlab)+"\n";
	st += "nid:"+JSON.stringify(nid)+"\n";
	st += "id:"+JSON.stringify(id)+"\n";
	st += "widthst:"+JSON.stringify(width)+"\n";
	st += "heightst:"+JSON.stringify(height)+"\n";
	st += "datanote:"+JSON.stringify(datanote)+"\n";
	st += "datanoteon:"+JSON.stringify(datanoteon)+"\n";
	st += "descnote:"+JSON.stringify(descnote)+"\n";
	st += "nsampmax:"+JSON.stringify(nsampmax)+"\n";
	st += "nsampevmax:"+JSON.stringify(nsampevmax)+"\n";
	st += "GRmax:"+JSON.stringify(GRmax)+"\n";
	st += "ESSmin:"+JSON.stringify(ESSmin)+"\n";
	st += "itermax:"+JSON.stringify(itermax)+"\n";
	st += "fileToLoadlist:"+JSON.stringify(fileToLoadlist)+"\n";
	st += "induoflag:"+JSON.stringify(induoflag)+"\n";
	
	return st;
}

function load()                                           // Loads a BICI file
{		
	var lines = textFromFile.split('\n');
	var n = 0, j, st;

	simres={}; infres={}; ppcres={}; indsim=[]; paramagesmooth=[]; paramtimesmooth=[]; 
	induoflag = 0;
	
	for(n = 0; n < lines.length; n++){
		st = lines[n];
		if(st.length > 0){
			j = 0; while(j < st.length && st.substr(j,1) != ":") j++; if(j == st.length){ alertp("Problem loading");}
			code = st.substr(0,j); st = st.substr(j+1);
			if(st != "undefined"){
				switch(code){
				case "modelsetup": modelsetup = JSON.parse(st); break;
				case "cla": cla = JSON.parse(st); break;
				case "ncla": ncla = cla.length; break;
				case "age": age = JSON.parse(st); break;
				case "time": time = JSON.parse(st); break;
				case "param": param = JSON.parse(st); break;
				case "paramsim": paramsim = JSON.parse(st); break;
				case "paramagesmooth": paramagesmooth = JSON.parse(st); break;
				case "paramtimesmooth": paramtimesmooth = JSON.parse(st); break;
				case "derive": derive = JSON.parse(st); break;
				case "termtype": termtype = JSON.parse(st); break;
				case "tsimmin": tsimmin = JSON.parse(st); break;
				case "tsimmax": tsimmax = JSON.parse(st); break;
				case "tpostmin": tpostmin = JSON.parse(st); break;
				case "tpostmax": tpostmax = JSON.parse(st); break;
				case "ndatavar": ndatavar = JSON.parse(st); break;
				case "datavar": datavar = JSON.parse(st); break;
				case "simpopinitset": simpopinitset = JSON.parse(st); break;
				case "simpopinitlist": simpopinitlist = JSON.parse(st); break;
				case "simty": simty = JSON.parse(st); break;
				case "simindinit": simindinit = JSON.parse(st); break;
				case "siminit": siminit = JSON.parse(st); break;
				case "clagendata": clagendata = JSON.parse(st); break;
				case "detectprob": detectprob = JSON.parse(st); break;
				case "datagenttype": datagenttype = JSON.parse(st); break;
				case "tgenmin": tgenmin = JSON.parse(st); break;
				case "tgenmax": tgenmax = JSON.parse(st); break;
				case "dtgen": dtgen = JSON.parse(st); break;
				case "tgenuserdef": tgenuserdef = JSON.parse(st); break;
				case "data": data = JSON.parse(st); break;
				case "nchain": nchain = JSON.parse(st); break;
				case "indshow": indshow = JSON.parse(st); break;
				case "claselect": claselect = JSON.parse(st); break;
				case "transpos": transpos = JSON.parse(st); break;
				case "transmat": transmat = JSON.parse(st); break;
				case "transsel": transsel = JSON.parse(st); break;
				case "transindfilt": transindfilt = JSON.parse(st); break;
				case "depsel": depsel = JSON.parse(st); break;
				case "histoplot": histoplot = JSON.parse(st); break;
				case "depop": depop = JSON.parse(st); break;
				case "depopsel": depopsel = JSON.parse(st); break;
				case "distCI": distCI = JSON.parse(st); break;
				case "distmean": distmean = JSON.parse(st); break;
				
				case "indsim": indsim = JSON.parse(st); break;
				case "inddata": inddata = JSON.parse(st); break;
				case "simnumber": simnumber = JSON.parse(st); break;
				case "indmaxnumber": indmaxnumber = JSON.parse(st); break;
				case "simres": simres = JSON.parse(st); break;
				case "infres": infres = JSON.parse(st); break;
				case "ppcres": ppcres = JSON.parse(st); break;
				case "datares": datares = JSON.parse(st); break;		
				case "data": data = JSON.parse(st); break;
				case "simdata": simdata = JSON.parse(st); break;
				case "radlab": radlab = JSON.parse(st); break;
				case "radstyle": radstyle = JSON.parse(st); break;
				case "radstyle2": radstyle2 = JSON.parse(st); break;
				case "radnamlab": radnamlab = JSON.parse(st); break;
				case "nid": nid = JSON.parse(st); break;
				case "id": id = JSON.parse(st); break;
				case "widthst":  widthold = JSON.parse(st); break;
				case "heightst":  heightold = JSON.parse(st); break;
				case "datanote": datanote = JSON.parse(st); break;
				case "datanoteon": datanoteon = JSON.parse(st); break;
				case "descnote": descnote = JSON.parse(st); break;
				case "nsampmax": nsampmax = JSON.parse(st); break;
				case "nsampevmax": nsampevmax = JSON.parse(st); break;
				case "GRmax": GRmax = JSON.parse(st); break;
				case "ESSmin": ESSmin = JSON.parse(st); break;
				case "itermax": itermax = JSON.parse(st); break;
				case "fileToLoadlist": fileToLoadlist = JSON.parse(st); break;
				case "induoflag": induoflag = JSON.parse(st); break;
				}
			}
		}
	}

	initpages();
	
	modelstart = 1;
	setsize();	
	
	initcompmult();
	settimeon();
	
	//param = []; collectvariables(); converttoobs(ty);
	//changepage(MODELPAGE,0,-1);

	return 0;
}

function loadexamp(val)                                   // Load one of the example files
{
	stopinference();
			 
	startloading();
	
	examploaded = val;
	
	pr("Example:"+val);
	xmlhttp.open("GET","Examples/"+val+".bici",true); xmlhttp.send();
}
