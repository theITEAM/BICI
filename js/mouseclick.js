function mouseclick(x,y)                                  // Fires when the mouse is clicked 
{
	var flag;
	
	if(gdropsel >= 0 && buttype[over] != GDROPSELBUT) gdropsel = -1;

	mx = x; my = y;
	
	if(helptype >= 0 && !(over >= 0 && buttype[over] == HELPBACKBUT) && !(over >= 0 && buttype[over] == DELFILEBUT)) helptype = -1;

	if(over >= 0){ finalaction(over); over = -1; canover = -1;}
	
	buttoninit(); mousemove(mx,my);
}

function mousedblclick()                                   // Fires for a double click
{
	var dtx, dy, ttmid, ymid, xx, yy;

	if(buttype[over] == CANVASBUT){
		if(canover == -1){
			if(arrow == 1){
				if(page == MODELPAGE){
					if(addtransmode == 0 && addcompmode == 0 && addsourcemode == 0){
						zoomin(pagesub[MODELPAGE],mx-menux,my,1.3,0,0);
						buttoninit();
					}
					if(addtransmode == 1){ deletetransmode(); buttoninit();}
					addcompmode = 0; addsourcemode = 0;
				}
			}
			if(arrow == 2){
				if((page == SIMULATEPAGE && pagesub[page] == 2 && simpagename[pagesubsub[page][2]] == "Individuals") ||
				   (page == INFERENCEPAGE && pagesub[page] == 2 && infpagename[pagesubsub[page][2]] == "Individuals") ||
				   (page == INFERENCEPAGE && pagesub[page] == 0 && pagesubsub[page][0] == 1)){
					indzoom(mx-tlinexmin-cornx,0.5);
				}
			}
		}
	}
}
	
function finalaction(i)                                   // What happens when a button is clicked
{
	var x, y, dx, dy, text, val, val2, ac, linkratest, tcvtemp, namest, subnamest;

	makesure = -2;

	x = butx[i]; y = buty[i]; dx = butdx[i]; dy = butdy[i], text = buttext[i]; val = butval[i]; val2 = butval2[i]; ac = butac[i];
	if(mx < x || mx > x+dx || my < y || my > y+dy) return;
	
	switch(ac){
	case TABBUT:		
		if(modelstart == 1){
			if(modelsetup == 0 && val > 2) helptype = 99;
			else{
				if(modelsetup == -1 && val >= 2) helptype = 201;
				else{
					loading = 0;
					changepage(val,-1,-1);
				}					
			}
		}
		break;
    
	case PAGESUBBUT:
		loading = 0;
		if(page == MODELPAGE && val == pagesub[MODELPAGE] && val >= 0 && val < ncla-2){
			if(modelsetup == 0){
				if(cla[val].name == "Time" || cla[val].name == "Age"){ selectbuby = 264+val*26; selectbub = EDITAGETIMEBUB;}
				else{ selectbuby = ymodel + (val+0.5)*dysub; selectbub = EDITCLASSBUB;}
	
				selectbubx = 3; selectbubval = val; selectbubdx = 0; selectbubdy =0; 
				bubble();
			}
			else selectbub = EDITBUB; 
			bubble();
		}	
		else changepage(-1,val,-1);
		break;
	
	case PAGESUBSUBBUT:
		changepage(-1,-1,val);
		break;
	
	case IMPORTAC:
		fitype = IMPORTFILE;
		buboff();		
		ById("fileToLoad").value = "";
		ById("fileToLoad").accept=".txt,.csv"; 
		loadsta();
		over = -1; canover = -1;
		break;

	case EXPORTAC2: exporton = 0; break;
	
	case EXPORTAC: exporton = 1; break;
	case EXPORTAC2: exporton = 0; break;
	
	case UPLOADBUT:
		fitype = val
		ById("fileToLoad").value = "";
		ById("fileToLoad").accept=".txt"; 
		over = -1; canover = -1;
		loadsta();
		break;
	
	case CANVASBUT:
		if(canover != -1 && canbutac[canover] != -1) canfinalaction(canover);
		if(canover == -1 && addtransmode == 1 &&  cla[transtempcl].tra[transtempk].i >= 0){
			addtransp();
			if(addsourcemode == -1){ cla[transtempcl].tra[transtempk].f = -2; addsourcemode = 0; completetrans(transtempcl,-2);}
		}
		if(canover == -1 && addcompmode == 1) addcomp();
		if(canover == -1 && addsourcemode == 1) addsource();
		break;
		
    case RELOADAC:
		stopinference();
		infres={}; ppcres={}; simres={};
	    location.reload(true);
		break;

	case ERRORBACKAC:
		warning=[]; infres.result = -1;
		break;
	
	case LOADBUT:
	 	fitype = BICIFILE;
		ById("fileToLoad").value = "";
		ById("fileToLoad").accept=".bici";
		loadsta();
		over = -1; canover = -1;
		break;
	
	case SAVEAC:
		if(modelstart == 1) saveon = 1;
		break;
	
	case CHECKBUT:
		res.popshow[val][val2] = 1-res.popshow[val][val2];
		plotinitfl = 0;
		break;
		
	case CHECKBUT2:
		indshow[val] = 1-indshow[val];
		break;
	
	case CHECKBUT3:
		res.popfiltcheck[val][val2] = 1-res.popfiltcheck[val][val2];
		plotinitfl = 0;
		break;
		
	case CHECKBUT4:
		transindcheck[val] = 1-transindcheck[val];
		plotinitfl = 0;
		break;
	
	case CHECKBUT5:
		distCI = 1-distCI;
		plotinitfl = 0;
		break;
		
	case CHECKBUT6:
		distmean = 1-distmean;
		plotinitfl = 0;
		break;
		
	case CHECKBUT7:
		res.varcheck[val] = 1-res.varcheck[val];
		startcalccor();
		break;
		
	case RADIOBUT:
		switch(val2){
		case RADIOTERM: termtype = val; break;
		case RADIOVAR: if(!ctrlkey) res.varsel=[]; res.varsel.push(val); xaxisauto = 1; yaxisauto = 1; plotinitfl = 0; break;
		case RADIOTRANS: if(!ctrlkey) transsel=[]; transsel.push(val); xaxisauto = 1; yaxisauto = 1; plotinitfl = 0; break;
		case RADIOCHOOSE: if(chooseaxis == 0) res.varselx = val; else res.varsely = val; chooseaxis = -1; break;
		case RADIOLABEL: radlab = val; break;
		case RADIOSTYLE: radstyle = val; break;
		case RADIOSTYLE2: radstyle2 = val; break;
		case RADIONAMLAB: radnamlab = val; break;
		case RADIODER: if(!ctrlkey) res.dersel=[]; res.dersel.push(val); xaxisauto = 1; yaxisauto = 1; plotinitfl = 0; break;
		case RADIODEP: if(!ctrlkey) depsel=[]; depsel.push(val); xaxisauto = 1; yaxisauto = 1; plotinitfl = 0; break;
		case RADIOTRANSTY: transty = val; break;
		case RADIOSCATTY: res.scattype = val; plotinitfl = 0; break;
		}
		break;
		
	case STARTAC:
		switch(page){
		case SIMULATEPAGE:
			simstarttext = text+pagesubsub[SIMULATEPAGE][2];
			if(text == "Rerun") changepage(SIMULATEPAGE,2,0);
			switch(text){
			case "Start": case "Rerun":
				if(siminit == "" || (siminit == "Load" && simindinit.length == 0)){ alertp("The initial population must be defined"); changepage(-1,1,-1); return;}
				if(runtype == "inf"){ helptype = 110; return;}
			
				simres.result = 0; startloading();
				setTimeout(function(){ simsetup();},300);
				break;
			
			case "Stop!": stopinference(); simres.result = 2; break;
			case "Exit!": case "Cancel": stopinference(); simres.result = -1; changepage(-1,2,-1); break;
			}
			break;
			
		case INFERENCEPAGE:
			switch(text){
			case "Start":
				if(data.length == 0){ alertp("Data must first be loaded!"); changepage(INFERENCEPAGE,0,0); return;}
			
				if(runtype == "sim"){ alertp("Simulation must be stopped before inference can begin."); return;}
				
				if(runtype == "inf"){
					if(!confirm("Do you want to stop the current inference?")) return;
					stopinference();
				}
				infres.result = 0; startloading();
				setTimeout(function(){ startinference(0,nchain);},300);
				break;
			
			case "Stop!": stopinference(); infres.result = 2; break;
			case "Exit!": case "Cancel": stopinference(); infres.result = -1; loading = 0; changepage(-1,2,-1); break;
			}
			break;
		}
		break;

	case STOPSTARTAC:
		stopinference();
		simres.result = 0; startloading();
		setTimeout(function(){ simsetup();},300);
		break;
		
	case BAYESBUT:
		var val = prompt("Please select a value for "+res.varname[res.varsel[0]]);
		if(val){
			if(isNaN(val)) alertp("This is not a number!");
			else calcbayesfac(val);
		}
		break;
		
	case EXPORTMINIBUT:
		ById("fileToSave").value = "";
		switch(val){
		case 0: exporttype = 4; ById("fileToSave").accept=".png"; savesta(); break;
		case 1: exporttype = 5; ById("fileToSave").accept=".txt"; savesta(); break;
		case 2: exporttype = 10; ById("fileToSave").accept=".bici"; savesta(); break;
		case 3: exporttype = 0; ById("fileToSave").accept=".bici"; savesta(); break;
		case 4: exporttype = 4; exportimage(0); break;
		case 5: exporttype = 7; ById("fileToSave").accept=".png"; savesta(); break;
		case 6: exporttype = 7; exportmodel(0); break;
		}
		saveon = 0; exporton = 0; over = -1; canover = -1;
		break;
	
	case EXPORTPARAMAC:
		exportparam();
		break;
		
	case EXPORTSTATEAC:
		exportstate();
		break;
		
	case EXPORTDIAGNOSTICAC:
		exporttype = 3;
		ById("fileToSave").value = "";
		ById("fileToSave").accept=".txt";
		savesta();
		over = -1; canover = -1;
		break;
	
	case DONEAC:
		if(addingdata == 2){ addingdata = 3; return;}
		tableyfr = 0;
		
		if(adddata() == 0){
			converttoobs("data");
			collectvariables();
			addingdata = 0;
			if(newfile == 1){ fileToLoadlist.push({IDcol:IDcol, name:fileToLoad.name, text:textFromFile});}
		}
		break;
	
	case BACKEVDATAAC2: datashow = "table"; break;
	
	case BACKEVDATAAC3: datashow = "seltrans"; break;
	
	case BACKEVDATAAC: addingdata--; break;
	
	case BACKAC:
		obsloaded = 1;
		setcolumns();
		break;
		
	case BACKAC2:
		obsloaded = 2;
		break;
	
	case POPFILTBUT:
		new_win = nw.Window.open('addpopfilt.html', { show:false});
		break;
		
	case ADDCLASSBUT:
		switch(val){
		case -1: 
			var necl = addclass(newname("Name"));
			changepage(-1,necl,-1);
			buboff(1);
			selectbubx = 3; //selectbuby = 261+necl*26;
			selectbubdx = 0; selectbubdy =0; selectbub = EDITCLASSBUB; selectbubval = necl;
			break;
		case 0: setageon(); break;
		case 1: settimeon(); break;
		}
		break;
		
	case SIMTYAC:
		simty = val;
		break;
	
	case CHOOSECLASSAC:
		datatemp.cl = val;
		datapossetup();
		break;
		
	case CANCELBUT2:
		addingdata = 0;
		changepage(INFERENCEPAGE,0,0);
		break;
	
	case NEXTDESCAC:
		changepage(MODELPAGE,-1,0);
		if(modelsetup == -1){
			pagesub[MODELPAGE] = -1; 
			modelsetup = 0;
			addclassbubon();
		}
		break;
		
	case NEXTBUT:
		tableyfr = 0;
		if(datashow == "table"){
			if(checkdata(1) == 0){
				switch(datatype){
				case "state": datashow = "obsmodel"; break;
				case "cap": datashow = "pdmodel"; break;
				case "trans": case "move": datashow = "seltrans"; break;
				case "pop": datashow = "pop"; break;
				case "der": datashow = "der"; break;
				default: alert("Problem EC99"); break;
				}
			}
		}
		break;
	
	case GENSTARTAC:
		gendatastart(val);
		break;
	
	case GENCANCAC:
		gentype = -1;
		break;
	
	case NEXTSIMAC:
		changepage(SIMULATEPAGE,pagesub[SIMULATEPAGE]+1,0);
		break;
		
	case NEXTINFAC:
		if(pagesub[page] == 0 && pagesubsub[page][0] == 0) changepage(page,pagesub[page],1);
		else changepage(page,pagesub[page]+1,0);
		break;
	
	case NEXTGENAC:
		switch(drawgennum){
		case -2:
			if(datagenttype == 0 && (isNaN(tgenmin) || isNaN(tgenmax) || tgenmin > tgenmax)) alertp("Not a valid time range");
			else changegennum(-1);
			break;
			
		case -1:
			sellist=[];
			for(cl = 0; cl < ncla-1; cl++) if(selclass[cl] == 1) sellist.push(cl);
			if(sellist.length == 0) alertp("At least one classification must be selected.");
			else changegennum(0);
			break;
			
		default:
			if(gentype == 0 && transselcl == -1 && datagenttype == 3){ alert("Must set transtion"); return;}
			
			changegennum(drawgennum+1);
			if(drawgennum == sellist.length) generatestatedata();
			if(drawgennum == sellist.length+1) generatecapturedata();
			break;
		}
		break;
	
	case NEXTGENEVAC:
		generateeventdata();
		changegennum(drawgennum+1);
		break;
		
	case NEXTGENAC2:
		transcl = -1; transni = ""; transnf = "";
		if(page == SIMULATEPAGE) changegennum(1);
		else datashow = "transfilt";
			
		break;
	
	case NEXTGENAC3:
		generatepopdata();
		changegennum(drawgennum+1);
		break;
		
	case NEXTGENAC4:
		generatederdata();
		changegennum(drawgennum+1);
		break;
		
	case BACKGENAC:
		changegennum(drawgennum-1);
		if(drawgennum == sellist.length) generatestatedata();
		break;
		
	case SAVEGENAC:
		exporttype = 6;
		ById("fileToSave").value = "";
		ById("fileToSave").accept=".txt";
		savesta();	
		over = -1; canover = -1;	
		break;
	
	case BACKAC3:
		datashow = "table";
		break;
		
	case BACKBUT:
		if(addingdata == 1){
			if(datatemp.cl != -1) datatemp.cl = -1;
			else datashow = "table";
		}
		else{
			addingdata = 0;
		}
		tableyfr = 0;
		break;

	case TABBACKAC:
		n = tablehist.length-1;
		if(n >= 0){
			colname = tablehist[n].colname;
			ncol = tablehist[n].ncol;
			ncoldef = tablehist[n].ncoldef;
			row = tablehist[n].row;
			rowwidth = tablehist[n].rowwidth;
			tablehist.pop();
			setcolumns();
		}
		break;
		
	case ADDDATAAC:
		switch(val2){
		case "state":
			for(cl = 0; cl < ncla-1; cl++) if(cla[cl].ncomp > 1) break;
			if(cl == ncla-1){ helptype = 200; return;}
			break;
			
		case "trans": case "move":
			for(cl = 0; cl < ncla-1; cl++) if(cla[cl].ntra > 0) break;
			if(cl == ncla-1){ helptype = 202; return;}
			break;
			
		case "der":
			if(derive.length == 0){ helptype = 203; return;}
			break;
		}
		
		for(cl = 0; cl < ncla; cl++) clagendata[cl] = "All";

		fitype = val; datatype = val2;
		ById("fileToLoad").value = "";
		ById("fileToLoad").accept=".txt,.csv";
		over = -1; canover = -1;
		helptype = 0;
		break;

	case IMPBUT:
		helptype = -1;
		if(val == 0) loadsta();
		break;
		
	case SIMINFBUT:
		switch(val){
		case 0: changepage(SIMULATEPAGE,0,-1); break;
		case 1: changepage(INFERENCEPAGE,0,-1); break;
		}
		break;
		
	case HELPCLOSEBUT:
		helptype = -1;
		break;
		
	case GDROPBUT:
		gdropfrac = 0;
		gdropsel = val;
		break;
		
	case GDROPSELBUT:
		if(gdropnum < val.length){
			if(mx > x+dx-17){
				dddy = gdropnum*gdropdy;
				ymin = y+Math.floor(1.15*gdropdy)+gdropfrac*dddy;
				ddx = 14; ddy = Math.floor(dddy*gdropnum/val.length);
				if(my >= ymin+ddy) gdropfrac += gdropnum/val.length;
				else{
					if(my <= ymin) gdropfrac -= gdropnum/val.length;
				}
				
				if(gdropfrac < 0) gdropfrac = 0; if(gdropfrac > 1-gdropnum/val.length) gdropfrac = 1-gdropnum/val.length;
				return;
			}				
		}
		
		if(gdropselop != "&*&"){
			var gd = gdropinfo[gdropsel];
			switch(gd.click){
			case "dataname":
				var ev;
				if(page == INFERENCEPAGE){
					if(pagesub[page] == 0) ev = inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
					else  ev = infres.inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
				}
				else ev = indsim.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
				data[ev.obsdata[gd.j]].val[ev.obsdatai[gd.j]] = gdropselop
				converttoobs("data");
				break;
				
			case "depop":
				depopsel = gdropselop;
				depsel=[]; 
				var h; for(h = 0; h < histoplot.length; h++){ if(histoplot[h].classname == depopsel){ depsel.push(h); break;}}
				break;
				
			case "dataname2":
				var ev;
				if(page == INFERENCEPAGE){
					if(pagesub[page] == 0) ev = inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
					else ev = infres.inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
				}
				else ev = indsim.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
				data[ev.obsdata].val[ev.obsdatai] = gdropselop;
				converttoobs("data");
				break;
				
			case "clasel":
				res.clasel = gdropselj;
				break;
			
			case "tranclasel":
				transcl = gdropselj;
				break;
				
			case "viewpopinit":
				viewpopinitset = 0; while(viewpopinitset < ncla && cla[viewpopinitset].name != gdropselop) viewpopinitset++;
				if(viewpopinitset == ncla) alertp("Error code EC7");
				break;
			
			case "simpopinit":
				simpopinitset = 0; while(simpopinitset < ncla && cla[simpopinitset].name != gdropselop) simpopinitset++;
				if(simpopinitset == ncla) alertp("Error code EC8");
				viewpopinitset = simpopinitset;
				break;
				
			case "filtfilt":
				cl = gd.cl;
				if(gdropselop == "All") res.popfilt[cl] = -1;
				else{
					if(gdropselop == "Checkbox"){
						res.popfilt[cl] = -2;
						res.popfiltcheck[cl]=[]; for(c = 0; c < cla[cl].ncomp; c++) res.popfiltcheck[cl][c] = 1;
					}
					else{
						for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == gdropselop) break;
						if(c == cla[cl].ncomp) alertp("Error code EC10");
						else res.popfilt[cl] = c;
					}
				}
				break;
			
			case "popview":
				res.popview = gdropselop;
				break;
			
			case "indview":
				res.indview = gdropselop; indplotst=[];
				break;
			
			case "dergen":
				dergensel = gdropselop;
				break;
			
			case "gendata":
				clagendata[gd.cl] = gdropselop;
				break;
				
			case "agesmooth":
				paramagesmooth[gd.j].type = gdropselop;
				paramagesmooth[gd.j].val = "";
				break;
				
			case "timesmooth":
				paramtimesmooth[gd.j].type = gdropselop;
				paramtimesmooth[gd.j].val = "";
				break;
				
			case "prior":
				param[gd.j].prior = gdropselop;
				switch(gdropselop){
				case "Unbounded": param[gd.j].val = []; break;
				case "Flat": param[gd.j].val = ["",""]; break;
				case "Gamma": param[gd.j].val = ["",""]; break;
				case "Normal": param[gd.j].val = ["",""]; break;
				case "Log-Normal": param[gd.j].val = ["",""]; break;
				case "Exponential": param[gd.j].val = [""]; break;
				case "Beta": param[gd.j].val = ["",""]; break;
				case "Weibull": param[gd.j].val = ["",""]; break;
				case "Fix": param[gd.j].val = [""]; break;
				}
				break;
			
			case "sampfilt":
				res.sampfilt = gdropselop;
				break;
		
			case "tabledrop":
				row[gd.j][gd.i] = gdropselop;
				break;
				
			case "nchain":
				nchain = parseInt(gdropselop);
				break;
			
			case "indfilt":
				res.indfilt = gdropselj;
				break;
			
			case "transindfilt":
				transindfilt = gdropselop;
				break;
				
			case "filter":
				res.filt = gdropselj; 
				if(page == INFERENCEPAGE){
					switch(infpagename[pagesubsub[page][2]]){
					case "Statistics":
						startcalc(); 
						break;
					case "Traces": case "Prob. Dist.": case "Trans. Dist.":
						res.varsel=[]; if(res.filtervar[res.filt].length > 0) res.varsel.push(res.filtervar[res.filt][0]);
						tableyfr = 0;
						break;
					}
				}
				if(page == SIMULATEPAGE){
					switch(simpagename[pagesubsub[page][2]]){
					case "Statistics":
						startcalc(); 
						break;
					}
				}
				break;
				
			case "scatnum":
				res.scatnum = gdropselop;
				break;
			
			case "runlab":
				res.runsel = gdropselj;
				if(page == INFERENCEPAGE && pagesub[page] == 2 && infpagename[pagesubsub[page][2]] == "Statistics") startcalc();
				
				if(page == SIMULATEPAGE && pagesub[page] == 2 && simpagename[pagesubsub[page][2]] == "Statistics") startcalc();
				if(page == INFERENCEPAGE && pagesub[page] == 2 && infpagename[pagesubsub[page][2]] == "Correlation") startcalccor();
				break;
				
			case "runlab2":
				res.runsel2 = gdropselj;
				if(page == INFERENCEPAGE && infpagename[pagesubsub[page][2]] == "Statistics") startcalc();
				if(page == SIMULATEPAGE && simpagename[pagesubsub[page][2]] == "Statistics") startcalc();
				break;
					
			case "runlab3":
				res.runsel3 = gdropselj;
				if(page == INFERENCEPAGE && infpagename[pagesubsub[page][2]] == "Statistics") startcalc();
				if(page == SIMULATEPAGE && simpagename[pagesubsub[page][2]] == "Statistics") startcalc();
				break;
				
			case "runpoplab":
				res.runpopsel = gdropselj;
				res.sampfilt = "All";
				break;
				
			case "transsel":
				transselcl = gd.clop[gdropselj];
				transseltr = gd.trop[gdropselj];
				break;
			}
			xaxisauto = 1; yaxisauto = 1; 
			gdropsel = -1; plotinitfl = 0; indplotst=[];
		}
		break;
	
	case CHOOSEALLAC:
		tablehist.push({colname:copy(colname), ncol:ncol,  ncoldef:ncoldef, row:copy(row), rowwidth:copy(rowwidth)});
	
		colname.splice(1,0,"Time");
		for(r = 0; r < nrow; r++){
			row[r].splice(1,0,"All");
			rowwidth[r].splice(1,0,100);
		}
		ncol++;
		ncoldef++; tablexfr = 0;
		setcolumns();
		break;
	
	case CHOOSETIMEBUT:
		var val = prompt("Please enter the observation time");
		if(isNaN(val)) alertp("This is not a number!");
		else{
			tablehist.push({colname:copy(colname), ncol:ncol,  ncoldef:ncoldef, row:copy(row), rowwidth:copy(rowwidth)});
	
			colname.splice(1,0,"Time");
			for(r = 0; r < nrow; r++){
				row[r].splice(1,0,val);
				rowwidth[r].splice(1,0,100);
			}
			ncol++;
			ncoldef++; tablexfr = 0;
			setcolumns();
		}
		break;
		
	case ADDDERIVEAC:
		var newna = prompt("Please enter the name of the derived variable");
		if(newna){	
			if(addderive(newna) == 1) alertp("The derived name cannot contain any of the following characters: space, underscore, comma, '[', ']', '{', '}'.");
			else{
				buttoninit();
				for(i = 0; i < ncanbut; i++) if(canbuttype[i] == EQBUT && canbutval[i] == derive.length-1) break; 
				if(i < ncanbut) selbut(i);  
			}
		}
		break;
	
	case ADDDEPBUT:
		for(p = 0; p < param.length; p++) if(param[p].dist == 0) break;
		if(p == param.length){ helptype = 204; return;}
		adddepfl = 1;
		break;
	
	case CANCELDEPAC:
		adddepfl = 0;
		break;
		
	case SIMINITPOPBACKAC: 
		if(setsimind() == 0){ canover = -1; tableyfr = 0; showsimindinit = 0;}
		break;
	
	case SIMLOADINITAC:
		fitype = SIMINITFILE;
		ById("fileToLoad").value = "";
		ById("fileToLoad").accept=".txt,.csv"; 
		over = -1; canover = -1;
		loadsta();
		break;
		
	case EDITDESCAC:
		selbuttype(SPEECHBUT2); 
		break;
		
	case SIMVIEWINITAC: showsimind(); break;

	case CHOOSEEXISTAC: existflag = 1; break;
	
	case NEWMODBUT: startnewmodel(); break;
		
	case RELOADSTATBUT: if(vcalc >= res.filtervar[res.filt].length) startcalc(); break;

	case RELOADCORBUT: startcalccor(); 
		break;

	case CLEARBUT2: donemodel(); break;
	
	case ADDFILTAC: ncoldefmax++; break;
	
	case REMADDFILTAC: ncoldefmax--; break;
	
	case SLIDEAC: break;
	
	case LOADFILEBUT:
		newfile = 0;
		textFromFile = fileToLoadlist[val].text;
		loadedfile();
		if(fileToLoadlist[val].IDcol >= 0 && fileToLoadlist[val].IDcol < ncol && 
			datatype != "cap" && datatype != "pop" && datatype != "der"){ 
			selecttablecol(fileToLoadlist[val].IDcol);
		}
		helptype = -1;
		break;
		
	case HELPBACKBUT: break;

	case HELPICONBUT: helptype = val; break;
	
	case HELPBACKBUT: break;
	
	case LOADEXAC: loadexamp(exampst); break;
	
	case CHECKALLBUT: 
		if(text == "Check All") vv = 1; else vv = 0;

		switch(val){
		case 0:	for(v = 0; v < res.indlist.length; v++) transindcheck[v] = vv; plotinitfl = 0; break;
		case 1: for(v = 0; v < cla[val2].ncomp; v++) res.popfiltcheck[val2][v] = vv; plotinitfl = 0; break;
		case 2: for(v = 0; v < cla[val2].ncomp; v++) res.popshow[val2][v] = vv; plotinitfl = 0; break;
		case 3: 
			for(v = 0; v < res.filtervar[val2].length; v++) res.varcheck[res.filtervar[val2][v]] = vv; 
			startcalccor();
			break;
		}
		break;		
		
	case DELFILEBUT:
		fileToLoadlist.splice(val,1);
		break;
		
	case SIMINITCANCELAC:
		showsimindinit = 0;
		break;
	
	case WHITEBUT: break;
	
	case ADVOPAC2: advop = 0; break;
	
	case ALLAC: 
		var r;
		for(r = 0; r < nrow; r++){ row[r].splice(1,0,"All"); rowwidth[r].splice(1,0,0);}
		colname.splice(1,0,"Compartments");
		ncol++; ncoldef++;
		calcrowwidth(); setcolumns();
		break;
	
	case PDFBUT:
		alertp("Manual not yet added!");
		break;
		
	case PDFBUT2:
		alertp("Paper not yet added!");
		break;
	
	case DONOTHINGBUT: helptype = -1; break;
		
	case CORBACKAC: res.varselx = -1; res.varsely = -1; break;
	
	default: alertp("Error code EC4"); break;
	}
}

