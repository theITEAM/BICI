function canbuttonplot()                                  // Plots all the canvas buttons on the page
{
	var i, ov;

	cv = graphcv;
	cv.clearRect(0,0,graphcan.width,graphcan.height);
	for(i = 0; i < ncanbut; i++){
		ov = 0; if(i == canover) ov = 1;
		canbutplot(i,ov);
	}
	cv = maincv;
}

function canfinalaction(i)                                // Click on a canvas button
{
	val = canbutval[i]; val2 = canbutval2[i]; 
	x = Math.floor(canbutx[i]); y = Math.floor(canbuty[i]); dx = Math.floor(canbutdx[i]); dy = Math.floor(canbutdy[i]);
	text = canbuttext[i]; 

	switch(canbutac[i]){
	case TABLECOLBUT:
		if(val == 0) IDcol = val;
	
		if(selectbubst == JOINBUB){ joincols(selectbubval,val); return;}
		
		if(val >= ncoldef) selecttablecol(val);
		break;

	case COMPBUT:
		if(addtransmode == 1){
			if(cla[transtempcl].tra[transtempk].i == -1) cla[transtempcl].tra[transtempk].i = val2;
			else{
				if(val2 != cla[transtempcl].tra[transtempk].i) completetrans(val,val2);
				else deletetransmode();
			}
		}
		else{
			if(modelsetup == 0){ if(timeup-timedown < 300 || dragged == 0) selbut(i);}
			else selectbub = EDITBUB; 
		}
		break;
	
	case TESTNAMEBUT:
		selbut(i);
		break;
	
	case COMPSMALLBUT:
		if(addingdata != 2) selbut(i);
		break;
		
	case COMPSMALLPERBUT:
		if(addingdata != 2) selbut(i);
		break;
	
	case COMPSMALLBUT2:
		if(addingdata >= 0 && addingdata != 2){
			switch(datatemp.type){
			case "binary":
				datatemp.posbinary[val][val2] = 1-datatemp.posbinary[val][val2];
				break;
			case "expression":
				selbut(i);
				break;
			}
		}
		break;
	
	case SENSAC:
		datagen[clgl].sensitive[val] = 1-datagen[clgl].sensitive[val];
		break;
		
	case COMPSMALLBUT3:
		if(addingdata >= 0 && addingdata != 2){
			datatemp.sensitive[val] = 1-datatemp.sensitive[val];
		}
		break;
	 
	case COMPSMALLBUT6:
		if(addingdata != 2){
			switch(datatemp.variety){
			case "state":
				if(datatemp.type == "simple") datatemp.posref[val] = val2;
				else{
					datatemp.posbinary[val][val2] = 1-datatemp.posbinary[val][val2];
				}
				break;
			}
		}
		break;
		
	case BUBBLEOFFBUT:
		buboff(1);
		break;
	
	case DONEBUT:
		buboff();
		break;
		
	case COLBUT:
		cla[selectbubval].comp[selectbubval2].col = collist[val];
		break;
		
	case ZOOMINBUT:
	case ZOOMOUTBUT:
		buboff();
		
		switch(page){
		case MODELPAGE:
			midx = (width-menux)/2; midy = height/2;
			fac = 1.3; if(canbutac[i] == ZOOMOUTBUT) fac = 1.0/fac;
			zoomin(pagesub[MODELPAGE],midx,midy,fac);
			break;
			
		case SIMULATEPAGE: case INFERENCEPAGE:
			fac = 0.5; if(canbutac[i] == ZOOMOUTBUT) fac = 1.0/fac;
			indzoom(tablewidth/2,fac); indplotst=[];
			break;
		}
		break;
	
	case INDAC:
		if(page == SIMULATEPAGE) smin = 0; else smin = res.burninev;
		i = val;
		indzoom(tablewidth/2,1);
		tma =  -large; tmi = large;
		for(ch = 0; ch < res.nch; ch++){
			w = res.ch[ch];
			for(s = smin; s < res.ch[ch].nsampev; s++){
				ws = w.sampev[s].ind[i];
				t = ws.evt[0]; if(t < tmi) tmi = t;
				t = ws.evt[ws.nev-1]; if(t > tma) tma = t;
			}
		}
		tmi -= 0.0001*(tma-tmi); tma += 0.0001*(tma-tmi)
		xaxisfixmin = tmi;
		xaxisfixmax = tma;
		break;
		
	case BUBBLEBUT:
		break;
		
	case TRANSBUT:
		cl = pagesub[MODELPAGE];
		addtrans(cl,"Exponential");
		buboff();
		break;
	
	case TRANSPBUT:
		if(timeup-timedown < 300 || dragged == 0) removetransp(val,val2);
		break;
		
	case TRANSPLBUT:
		switch(page){
		case MODELPAGE:
			if(addtransmode == 0){
				if(modelsetup == 0) seltrans(val,val2);
				else selectbub = EDITBUB;
			}
			break;
			
		case SIMULATEPAGE:
		case INFERENCEPAGE:
			cl = val; 
			transcl = cl;
			ci = cla[cl].tra[val2].i; cf = cla[cl].tra[val2].f;
			if(ci >= 0 && cf >= 0){
				transni = cla[cl].comp[ci].name; 
				transnf = cla[cl].comp[cf].name;	
				if(page == SIMULATEPAGE) changegennum(1);
				else{
					if(datatype == "move") adddata();
					else datashow = "transfilt";
				}
			}
			break;
		}
		break;
		
	case ADDPARAMBUT:
		popdx = 160; 
		if(page == MODELPAGE && pagesub[page] >= 0){ popdy = 220+(ncla-1)*20; if(ncla <= 1) popdy -= 20;}
		else popdy = 180;
		new_win = nw.Window.open('addparam.html', { show:false});
		break;
		
	case ADDPOPAC:
		new_win = nw.Window.open('addpop.html', { show:false});	
		break;
	
	case FUNCTBUT:
		inserttext(functi[val]+"()");
		a.selectionStart--;  a.selectionEnd--;
		break;

	case OPBUT:
		inserttext(opbut[val]);
		break;
	
	case NUMBUT:
		inserttext(numbut[val]);
		break;

	case GREEKBUT:
		inserttext(greek[val]);
		break;
		
	case ALPBUT:
		inserttext(alp[val]);
		break;
		
	case DELCOMPAC:
		cl = selectbubval; i = selectbubval2;
		if(cla[cl].name == "Time"){
			if(i == 0){ buboff(1); alertp("Cannot delete");}
			else{ buboff(1); time.splice(i-1,1); settimeon();}
		}
		else{
			if(cla[cl].name == "Age"){
				if(i == 0){ buboff(1); alertp("Cannot delete");}
				else{ buboff(1); age.splice(i-1,1); setageon();}
			}
			else{
				if(confirm("Are you sure you want to delete this compartment?")) deletecomp(cl,i);
			}
		}
		break;
		
	case DELTRANSAC:
		deletetrans(selectbubval,selectbubval2);
		selectbub = -1; buboff();
		break;
		
	case DELTRANSAC2:
		if(cla[cl].name == "Time"){ buboff(1); time.splice(selectbubval2,1); settimeon();}
		else{
			if(cla[cl].name == "Age"){ buboff(1); age.splice(selectbubval2,1); setageon();}
		}
		break;
	
	case EVSOURCEBUT:
	case EVSINKBUT:
		selbut(i);
		break;
		
	case MINMAXBUT: case MINMAXBIGBUT: case SELPARAMBUT:
		selbut(i);
		break;
	
	case PDBUT:
		selbut(i);
		break;
		
	case EDITDBUT:
		selbut(i);
		break;
		
	case GENUSERBUT:
		selbut(i);
		break;
		
	case CANCELBUT:
		buboff(1);
		break;
		
	case CONVERTAC:
		convertdates();
		buboff();
		break;
		
	case TABLEBUT:
		if(page == INFERENCEPAGE && pagesub[page] == 0 && datatemp.variety == "cappd"){
			selectelement(val2,val,PDBUT);
		}		
		else selectelement(val2,val,TABLEBUT);
		break;
	
	case TABLEDROPBUT:
		if(addingdata != 2){ selectelement(val2,val,TABLEDROPBUT);}
		break;
	
	case TABLEBUT3:
		selbut(i);
		break;
		
	case TABLEHEADBUT:
		if(selectbubst == JOINBUB){ joincols(selectbubval,val); return;}
	 	selecthead(val,TABLEHEADBUT);
		break;
		
	case DELROWAC:
		buboff(1);
		row.splice(selectbubval2,1);
		rowwidth.splice(selectbubval2,1);
		nrow--;	
		setcolumns();
		break;
		
	case INSROWAC:
		buboff(1);
		var newrow=[], j, n = row[selectbubval2].length;
		for(j = 0; j < n; j++) newrow[j] = "";
		row.splice(selectbubval2,0,newrow);
		var newrowwidth=[]; newrowwidth.length = rowwidth[selectbubval2].length;
		rowwidth.splice(selectbubval2,0,newrowwidth);
		nrow++;	
		calcrowwidth(); setcolumns();
		break;
		
	case JOINCOLSAC:
		buboff(1);
		selectbub = JOINBUB;
		break;
		
	case SEARCHAC:
		inpon = 0;
		selectbub = SEARCHBUB;
		break;
		
	case DOSEARCHAC:
		buboff(0);
		break;
		
	case REPLACEAC:
		inpon = 0;
		selectbub = REPLACEBUB;
		break;
	
	case DOREPLACEAC:
		doreplace();
		break;
		
	case DELROWSAC:
		inpon = 0;
		selectbub = DELROWSBUB;
		break;
		
	case DODELROWSAC:
		dodelrows();
		break;
	
	case SORTAC:
		if(val == 0) dosort(selectbubval,"alp","swap");
		else dosort(selectbubval,"num","swap");
		buboff(1);
		buttoninit();
		break;
		
	case CANRADIOBUT:
		if(val2 != CANRADIOXAXIS && val2 != CANRADIOYAXIS) buboff(1);
		switch(val2){
		case CANRADIORATE: cla[selectbubval].tra[selectbubval2].ratetime = val; break;
		case CANRADIOCLASS: if(addingdata == 1 || addingdata == 3) datatemp.cl = val; break;
		case CANRADIOCHECK: if(addingdata == 1 || addingdata == 3) datatemp.type = val; break;
		case CANRADIOTESTRES: if(addingdata == 1 || addingdata == 3) datatemp.postestres[Math.floor(val/1000)] = val%1000; break;
		case CANRADIODEFINE: simpopinitset = val; break;
		case CANRADIOSIMTY: simty = val; break;
		case CANRADIOTERM: termtype = val; break;
		case CANRADIOXAXIS: 
			xaxisauto = val;
			if(xaxisauto == 0){
				xaxisfixmin = parseFloat(axxmin.toPrecision(4));
				xaxisfixmax = parseFloat(axxmax.toPrecision(4));
			}
			else{ indplotstclear(); inpon = 0; ById("add").innerHTML = "";}
			break;
		case CANRADIOYAXIS:
			yaxisauto = val;
			if(yaxisauto == 0){
				yaxisfixmin = parseFloat(axymin.toPrecision(4));
				yaxisfixmax = parseFloat(axymax.toPrecision(4));
			}
			else{ inpon = 0; ById("add").innerHTML = "";}
			break;
		case CANRADIOGENT: datagenttype = val; break;
		case CANRADIOGENT2: datagenttype2 = val; break;
		case CANRADIOCHECK2: datagen[sellist[drawgennum]].type = val; break;
		case CANRADIOWHICH: whichind = val; break;
		case CANRADIOPD: obspd = val; break;
		case CANRADIOPD2: datatemp.obspd = val; break;
		case CANRADIOPD3: datatemp.pdsame = val; break;
		case CANRADIOMANUAL: siminit = val; break;
		case CANRADIOERRBAR: errbar = val; break;
		case CANRADIOINDUO: induoflag = val; break;
		}
		break;
		
	case EQBUT:
		selbut(i);
		break;
		
	case ADDCOMPAC:
		buboff();
		addcompmode = 1; 
		cl = pagesub[0];
		
		dragcoli = 7; 
		do{
			dragcol = collist[dragcoli];
			for(cl2 = 0; cl2 < ncla; cl2++){
				for(k = 0; k < cla[cl2].ncomp; k++) if(cla[cl2].comp[k].col == dragcol) break;
				if(k < cla[cl2].ncomp) break;
			}
			if(cl2 == ncla) break;
			dragcoli++;
		}while(1 == 1);
		dragtext = newname("Comp");
		dragdh = cla[cl].zoomfac*50;
		dragww = Math.floor(textwidth(dragtext, Math.floor(dragdh*0.7)+"px georgia"))+dragdh;	
		break;
	
	case ADDTRANSAC:
		selbut(i);
		break;
		
	case ADDLINKAC:
		buboff();	
		if(text == "Load"){
			fitype = TRANSFILE;
			ById("fileToLoad").value = "";
			ById("fileToLoad").accept=".txt"; 
			loadsta();
			over = -1; canover = -1;
		}
		else{
			if(text == "Source"){
				for(cl = 0; cl < ncla; cl++){
					if(cl != pagesub[MODELPAGE]){
						for(tr = 0; tr < cla[cl].ntra; tr++){
							if(cla[cl].tra[tr].type == "Source"){ helptype = 206; return;}
						}
					}							
				}
				addsourcemode = 1; 		
			}
			else{
				if(text == "Sink"){
					for(cl = 0; cl < ncla; cl++){
						if(cl != pagesub[MODELPAGE]){
							for(tr = 0; tr < cla[cl].ntra; tr++){
								if(cla[cl].tra[tr].type == "Sink"){ helptype = 205; return;}
							}
						}							
					}
					addsourcemode = -1; 
				}

				cl = pagesub[MODELPAGE];
				if(text == "Markovian") text = "Exponential";
				addtrans(cl,text);
			}
		}
		break;
	
	case CROSSBUT:
		if(addtransmode == 1) deletetransmode();
		addcompmode = 0; addsourcemode = 0;
		break;
		
	case DELCLASSAC:
		buboff();
		if(confirm("Are you sure you want to delete "+cla[selectbubval].name)) delclass(selectbubval);
		break;
		
	case DONEMODELAC:
		donemodel();
		break;

	case VIEWBUT:
		datatemp = copy(data[val]);
		reloaddata();
		ncoldef = ncol; ncoldefmax = ncol;
		
		calcrowwidth();
		setcolumns();
		 
		dataselected = val;
		datashow = val2;
		addingdata = 2;
		break;
		
	case DELETEBUT:
		if(confirm("Are you sure you want to delete this data source?")){
			data.splice(val,1);
			converttoobs("data");
			collectvariables();
		}
		break;
	
	case EVBUT: case EVBUT2: case EVBUT3: case EVBUT4: case EVBUT5: case EVBUT6: case EVBUT7:
		selbut(i);
		break;
	
	case XLABELBUT: chooseaxis = 0; break;
	case YLABELBUT: chooseaxis = 1; break;
	
	case XTICKBUT:
		selbut(i);
		break;
		
	case YTICKTRBUT:
		selbut(i);
		break;
		
	case SUBBUT:
		param = priorsubbut(param,val,val2);
		break;
		
	case SUBBUT2:
		paramsim = priorsubbut(paramsim,val,val2);
		break;
		
	case ADDTIMEPOINTAC:
		newtime = prompt("Enter time:");
		if(newtime){
			if(isNaN(newtime)) alertp("Not a number!");
			else{
				buboff(1);
				newtime = parseFloat(newtime);
				j = 0; while(j < time.length && newtime > time[j]) j++;
				if(newtime == time[j]) alertp("Time point already exists");
				else{
					time.splice(j,0,newtime);
					settimeon();
				}
			}
		}
		break;
		
	case ADDAGEPOINTAC:
		newtime = prompt("Enter time:");
		if(newtime){
			if(isNaN(newtime)) alertp("Not a number!");
			else{
				buboff(1);
				newtime = parseFloat(newtime);
				if(newtime <= 0) alertp("Number must be positive!");
				else{
					j = 0; while(j < age.length && newtime > age[j]) j++;
					if(newtime == age[j]) alertp("Time point already exists");
					else{
						age.splice(j,0,newtime);
						setageon();
					}
				}
			}
		}
		break;
		
	case CANUPLOADBUT: selclass[val] = 1-selclass[val]; break;
	
	case PLAYBUT: 
		buboff(); playing = 1-playing; 
		if(playing == 1){
			if(playtime == POPX-1) playtime = 0;
			playstartt = getsec();
			playtimeinit =  playtime;
			playanim();
		}
		break;
		
	case PLAYLINEBUT: playtime = Math.floor(POPX*(mx-x-cornx)/dx); playing = 0; break;
	
	case LINKBUT: val.posrefmore[val2] = 1-val.posrefmore[val2]; break;
	
	case COMPPOPINITBUT: case COMPFRACINITBUT: case COMPLOADINITBUT: case COMPPOPOUTPUTBUT: selbut(i); break;	
		
	case EXPANDAC: selectbub = TRANSPLBIGBUB; break;
	case REDUCEAC: selectbub = TRANSPLBUT; break;

	case EXAMPBUT: 
		if(modelstart == 1){ helptype = 103; exampst = val;}
		else loadexamp(val);
		break;
			
	case SPEECHBUT:
		if(page == MODELPAGE){ cla[pagesub[MODELPAGE]].descon = 1- cla[pagesub[MODELPAGE]].descon;}
		else datanoteon = 1-datanoteon;
		if(isNaN(datanoteon)) datanoteon = 1; 
		break;
		
	case EDITSPEECHBUT: 
		selbut(val);
		break;
		
	case PARACANBUT2:
		selbut(i);
		break;
		
	case RELOADBUT:
		if(plotinitfl == 2) plotinitfl = 0;
		break;
		
	case RELOADBUT2:
		indplotst=[];
		buttoninit();
		break;
		
	case RELOADBUT3:
		if(plotinitfl == 2) plotinitfl = 0;
		break;
		
	case NEXTSEARCHAC:
		if(searchres.length == 0) searchresnum = 0;
		else{
			searchresnum++; if(searchresnum==searchres.length) searchresnum = 0;
			selectelement(searchres[searchresnum],selectbubval,SEARCHRESBUB);
		}
		break
		
	case BACKSEARCHAC:
		if(searchres.length == 0) searchresnum = 0;
		else{
			searchresnum--; if(searchresnum == -1) searchresnum = searchres.length-1;
			selectelement(searchres[searchresnum],selectbubval,SEARCHRESBUB);
		}
		break;
	
	case UNDOAC: row = rowcopy; nrow = row.length; calcrowwidth(); setcolumns(); buboff(1); break;
	
	case DELDERAC: derive.splice(val); break;
	
	case SPEECHBUT2:
		selbut(i);
		break;
		
	case HELPICONCANBUT: helptype = val; break;
	
	case SLIDERBUT: break;
	 
	case SELPARBUT: 
		adddepfl = 0; param[val].dist = 1; param[val].prior = "Unspecified"; param[val].val[0] = ""; param[val].val[1] = "";
		break;
	
	case DELDEPAC:
		param[val].dist = 0; param[val].prior = "Unspecified"; param[val].val[0] = ""; param[val].val[1] = ""; 
		collectvariables();
		break;
	
	case INACTIVEBUT:	
		switch(page){
		case MODELPAGE: selectbub = EDITBUB; break;
		case INFERENCEPAGE: 
			if(pagesub[page] == 0){
				switch(val){
				case 0: helptype = 213; break;
				case 1: helptype = 216; break;
				default: helptype = 212; break;
				}
			}
			break; 
		}
		break;
	
	case INDBUT:
		var tmin = large, tmax = -large;

		if(page == INFERENCEPAGE && pagesub[page] == 0){
			var indi = inddata.ind[val];
			for(cl = 0; cl < ncla; cl++){
				for(e = 0; e < indi.cl[cl].ev.length; e++){
					var ev = indi.cl[cl].ev[e];
					if(ev.variety != "transtmin" && ev.variety != "transtmax"){
						if(ev.t < tmin) tmin = ev.t; if(ev.t > tmax) tmax = ev.t;
					}
				}
			}		
		}
		else{
			if(page == INFERENCEPAGE) rest = infres; else rest = simres;
			if(res.sampfilt == "All"){
				for(ch2 = 0; ch2 < rest.nch; ch2++){
					for(s = 0; s < rest.ch[ch2].nsampev; s++){
						ressa = rest.ch[ch2].sampev[s].ind[val];
						if(ressa.evt[0] < tmin) tmin = ressa.evt[0];
						if(ressa.evt[ressa.nev-1] > tmax) tmax = ressa.evt[ressa.nev-1];
					}
				}
			}
			else{
				var r = getsampfilt(); ch2 = r.ch; s = r.s; 
				ressa = rest.ch[ch2].sampev[s].ind[val];
				if(ressa.evt[0] < tmin) tmin = ressa.evt[0];
				if(ressa.evt[ressa.nev-1] > tmax) tmax = ressa.evt[ressa.nev-1];
			}
		}
		
		if(tmin == large) alert("The are no observations on this individual!");
		else{
			if(xaxisauto == 1) xaxisauto = 0;
			if(tmin == tmax){ tmin -= (xaxisfixmax-xaxisfixmin)/20; tmin += (xaxisfixmax-xaxisfixmin)/20;}
			else{ dd = tmax-tmin; tmin -= dd/6; tmax += dd/6;}
			if(!(page == INFERENCEPAGE && pagesub[page] == 0)){
				if(tmin < rest.tmin) tmin = rest.tmin;
				if(tmax > rest.tmax) tmax = rest.tmax;
				indplotstclear();
			}
			xaxisfixmin = tmin;
			xaxisfixmax = tmax;
		}	
		break;
		
	case ADVOPAC: advop = 1; buboff(1); break;
	
	case CORBUT:
		res.varselx = val; res.varsely = val2; plotinitfl = 0;
		break;
		
	default: alertp("Error code EC2 "+canbutac[i]); break;
	}
}

function getsec(){ return new Date().getTime() / 1000;}    // Gets the time in seconds

function selbut(i)                                         // Selects a button
{
	if(selectbub != -1) buboff(1);
	
	selectbubx = canbutx[i]; selectbuby = canbuty[i]; selectbubdx = canbutdx[i]; selectbubdy = canbutdy[i];
	selectbub = canbuttype[i];
	selectbubval = canbutval[i];
	selectbubval2 = canbutval2[i];
	selectbubtext = canbuttext[i];
}

function selbuttype(ty)                                    // Selects a button of a given type
{
	i = 0; while(i < ncanbut && canbuttype[i] != ty) i++;
	if(i == ncanbut){ alertp("Problem EC20"); return;}
	
	if(selectbub != -1) buboff(1);
	
	selectbubx = canbutx[i]; selectbuby = canbuty[i]; selectbubdx = canbutdx[i]; selectbubdy = canbutdy[i];
	selectbub = canbuttype[i];
	selectbubval = canbutval[i];
	selectbubval2 = canbutval2[i];
	selectbubtext = canbuttext[i];
}

function canbutplot(i,ov)                                // Plots an individual canvas button
{
	if(addtransmode == 1){
		tr = cla[transtempcl].tra[transtempk];
		if(canover == -1){  if(tr.i != -1){ addtransp(); findpline(transtempcl,transtempk); drawtrans(tr,0); tr.p.pop();}}
		else{
			if(canover != -1 && canbutac[canover] == COMPBUT){
				c = canbutval2[canover];
				if(tr.i != -1 && c != tr.i){
					tr.f = c; findpline(transtempcl,transtempk); drawtrans(tr,0); tr.f = -1;
				}
			}
		}
	}
	
	x = Math.floor(canbutx[i]); y = Math.floor(canbuty[i]); dx = Math.floor(canbutdx[i]); dy = Math.floor(canbutdy[i]);
	text = canbuttext[i]; ty = canbuttype[i];
	val = canbutval[i]; val2 = canbutval2[i];
	if(selectbub == canbuttype[i] && selectbubval == val && selectbubval2 == val2) ov = 0;
	
	switch(ty){ 
	case EXAMPMODBUT:
		cv.drawImage(examppic[val],x,y);
		break;
	
	case EXAMPBUT:
		if(val == examploaded){
			col = DRED; if(ov == 1) col = LRED;
			plottext(text,x+11,y+16,examplefont2,col);
		}
		else{
			col = RED; if(ov == 1) col = LRED;
			plottext(text,x+11,y+16,examplefont,col);
		}		
		ddx = 0; ddx2 = 5; ddy = 5; dy = 20;
		drawline(x+ddx,y+dy/2-ddy,x+ddx+ddx2,y+dy/2,col,THICKLINE)
		drawline(x+ddx,y+dy/2+ddy,x+ddx+ddx2,y+dy/2,col,THICKLINE)
		break;
		
	case TABLEHEADBUT:
		gettabcol(val,dy);
		if(ov == 1 && dy > 0){
			fillrect(x, y, dx, dy, col2);
			plottext(text,x+4,y+16,tableheadfont,BLACK);
		}
		else{
			plottext(text,x+4,y+16,tableheadfont,col);
		}
		break;
	
	case PARAMBUT:
		col = BLACK;
		plottext(text,x+4,y+16,tablefont,col);
		break;
		
	case TABLEBUT: case TABLEDROPBUT: case TABLEBUT3:	
		gettabcol(val,dy);
		if(ov == 1 && dy > 0) fillrect(x, y, dx, dy, col2);
		plottext(text,x+4,y+16,tablefont,col,dx);
		break;
	
	case SUBBUT: case SUBBUT2:
		col = DGREEN; if(ov == 1) col = GREEN;
		plottext(text,x+1,y+11,subfont,col);
		break;
	
	case BUBSELBUT:
		cv.beginPath();
		cv.lineWidth = 1;
		cv.rect(x,y,dx,dy);
		cv.strokeStyle = RED;
		cv.stroke();
		break;
		
	case TESTNAMEBUT:
		col = BLACK; if(ov == 1) col = BLUE;
		plottext(text,x+4,y+18,"18px georgia",col);
		break;
		
	case CLASSNAMEBUT:
		plottext(text,x+4,y+18,"bold 20px georgia",BLACK);
		break;
		
	case CLASSVALBUT:
		col = cla[val].val[val2].col; if(ov == 1 && drag != 8) col = darkcol(col);
		drawroundrect(x,y,dx,dy,Math.floor(0.2*dy),col,darkcol(col));
		centertext(text,x+dx/2,y+17,classvalfont,WHITE);
		break;
	
	case CANUPLOADBUT:
		col = WHITE; col2 = RED; col3 = RED; if(ov == 1){ col = LLRED;}
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2,y+21,INPUTFONT,col3);
		break;
		
	case DRAGBUT:
		cv.globalAlpha = 0.5;
		col = dragcol; if(ov == 1) col = darkcol(col);
		drawroundrect(x,y,dx,dy,Math.floor(0.2*dy),col,darkcol(col));
		centertext(text,x+dx/2,y+Math.floor(0.75*dy),Math.floor(dy*0.7)+  "px georgia",WHITE);
		cv.globalAlpha = 1;	
		break;
		
	case TABLECOLBUT:
		if(ov == 1 && val >= ncoldef) fillrect(x, y, dx, dy, LLRED);
		break;
	
	case RELOADBUT:
		if(plotinitfl == 2){ col = RED; if(ov == 1) col = LRED; reloadsign(x,y,col);}
		break;
		
	case RELOADBUT2:
		col = RED; if(ov == 1) col = LRED; reloadsign(x,y,col);
		break;
		
	case RELOADBUT3:
		if(plotinitfl == 2){ col = RED; if(ov == 1) col = LRED; reloadsign(x,y,col);}
		break;
		
	case PLAYBUT:
		if(ov == 1) col = LRED; else col = LBLUE;

		fillrect(x,y,dx,dy,col);
		if(playing == 0){
			if(playtime < POPX-1){
				ddx = 5; ddy = 5;
				polypoint[0][0] = x+ddx; polypoint[0][1] = y+ddy; 
				polypoint[1][0] = x+dx-ddx; polypoint[1][1] = y+dy/2; 
				polypoint[2][0] = x+ddx; polypoint[2][1] = y+dy-ddy; 
				drawpolygon(3,WHITE,col,THICKLINE);
			}
			else{
				cv.lineWidth = 3; 
				cv.beginPath();
				cv.arc(x+dx/2, y+dy/2+1, 8, -0.1, 1.5*Math.PI);
				cv.strokeStyle = WHITE;
				cv.stroke(); 
				drawarrow(x+dx/2+7,y+6,x-100,y+6,8,WHITE);
			}
		}
		else{
			x1 = 5, y1 = 5, dx1 = 5;
			fillrect(x+x1,y+y1,dx1,dy-2*y1,WHITE);
			fillrect(x+dx-x1-dx1,y+y1,dx1,dy-2*y1,WHITE);
		}
		break;
	
	case PLAYLINEBUT:
		xx = Math.floor(dx*playtime/(POPX-1));
		fillrect(x,y,xx,dy,LGREY);
		drawline(x+xx,y,x+xx,y+dy,DGREY,NORMLINE);
		drawrect(x,y,dx,dy,BLUE,NORMLINE);
		break;
		
	case COMPBUT: case COMPPOPINITBUT: case COMPFRACINITBUT: case COMPLOADINITBUT: case COMPPOPOUTPUTBUT:
		if(ov == 1) compplot(val,val2,x,y,dx,dy,ty,1);
		else{
			if(playing == 1 && ty == COMPPOPOUTPUTBUT) compplot(val,val2,x,y,dx,dy,ty,0);
		}
		break;
		
	case COMPSMALLPERBUT2:
		col = cla[val].comp[val2].col;
		cv.globalAlpha = 0.5;
		drawroundrect(x,y,dx,dy,Math.floor(0.2*dy),col,darkcol(col));
		centertext(text,x+dx/2,y+21,"20px georgia",WHITE);
		cv.globalAlpha = 1;
		break;
		
	case COMPSMALLBUT2:
		col = cla[datatemp.cl].comp[val2].col; 
		al = 1;
		if(ov == 1 && addingdata != 2) col = darkcol(col);
		if(text.substr(text.length-3,3) == "= 0" ) al = 0.5;
		if(addingdata == 2) al *= 0.7;
		
		cv.globalAlpha = al;
		drawroundrect(x,y,dx,dy,Math.floor(0.2*dy),col,darkcol(col));
		centertext(text,x+dx/2,y+21,"20px georgia",WHITE);
		cv.globalAlpha = 1;
		break;
	
	case COMPSMALLBUT3:
		col = cla[clgl].comp[val].col; 
		al = 1;
		if(ov == 1 && addingdata != 2) col = darkcol(col);
		if(text.substr(text.length-1,1) == "✖") al = 0.5;
		if(addingdata == 2) al *= 0.7;
		cv.globalAlpha = al;
		drawroundrect(x,y,dx,dy,Math.floor(0.2*dy),col,darkcol(col));
		centertext(text,x+dx/2,y+21,"20px georgia",WHITE);
		cv.globalAlpha = 1;
		break;
		
	case COMPSMALLBUT4:
		col = cla[clgl].comp[val2].col; 
		al = 1; if(addingdata == 2) al = 0.5;
		cv.globalAlpha = al;
		drawroundrect(x,y,dx,dy,Math.floor(0.2*dy),col,darkcol(col));
		centertext(text,x+dx/2,y+16,"16px georgia",WHITE);
		cv.globalAlpha = 1;
		break;
		
	case COMPSMALLBUT5:
		col = cla[val].comp[val2].col; 
		drawroundrect(x,y,dx,dy,Math.floor(0.2*dy),col,darkcol(col));
		centertext(text,x+dx/2,y+16,"16px georgia",WHITE);
		break;
		
	case COMPSMALLBUT6:
		col = cla[clgl].comp[val2].col; if(ov == 1 && addingdata != 2) col = darkcol(col); 
		
		if(text.substr(text.length-1,1)=="✖") cv.globalAlpha = 0.5;
		drawroundrect(x,y,dx,dy,Math.floor(0.2*dy),col,darkcol(col));
		centertext(text,x+dx/2,y+21,"20px georgia",WHITE);
		if(text.substr(text.length-1,1)=="✖") cv.globalAlpha = 1;
		break;
		
	case COMPSMALLBUT7:
		col = cla[clgl].comp[val2].col;
		drawroundrect(x,y,dx,dy,Math.floor(0.2*dy),col,darkcol(col));
		centertext(text,x+dx/2,y+21,"20px georgia",WHITE);
		break;
		
	case TRANSPLBUT:
		if(ov == 1) drawtrans(cla[val].tra[val2],ov);
		break;
	 
	case BUBBLEBUT:
		drawroundrect(x,y,dx,dy,20,LLLBLUE,LBLUE);
		if(cornbub == 0){
			drawline(bubx2,buby2,bubx3,buby3,LLLBLUE,NORMLINE);
			polypoint[0][0] = bubx1; polypoint[0][1] = buby1; 
			polypoint[1][0] = bubx2; polypoint[1][1] = buby2; 
			polypoint[2][0] = bubx3; polypoint[2][1] = buby3; 
			drawpolygon(3,LLLBLUE,LLLBLUE,THICKLINE);
			drawline(bubx1,buby1,bubx2,buby2,LBLUE,NORMLINE);
			drawline(bubx1,buby1,bubx3,buby3,LBLUE,NORMLINE);
		}
		break; 
		
	case BUBBLEOFFBUT:
		col = LRED; col2 = RED; if(ov == 1){ col = RED; col2 = DRED;}
		r = dx/2-1; r2 = r-3;
		x = x+r+1; y = y+r+1; 
		drawroundrect(x-r,y-r,2*r,2*r,r,col,col2);
		drawline(x-r2,y-r2,x+r2,y+r2,WHITE,THICKLINE);
		drawline(x+r2,y-r2,x-r2,y+r2,WHITE,THICKLINE);
		break;
		
	case BUBBLETEXTBUT:
		plottext(text,x+5,y+12,"16px arial",BLUE); 
		break;
		
	case BUBBLETEXTBUT2:
		plottext(text,x+5,y+12,"14px arial",BLACK); 
		break;

	case BUBBLETEXTBUT3:
		plottext(text,x+5,y+12,"bold 16px arial",BLUE); 
		break;
		
	case COLBUT:
		col = collist[val];
		drawroundrect(x,y,dx,dy,3,col,darkcol(col));
		if(cla[selectbubval].comp[selectbubval2].col == col){
			r = 3; if(col == BLACK) col2 = WHITE; else col2 = BLACK;
			drawroundrect(x+dx/2-r,y+dy/2-r,2*r,2*r,r,col2,col2);
		}
		else{
			if(ov == 1){ r = 3; drawroundrect(x+dx/2-r,y+dy/2-r,2*r,2*r,r,GREY,GREY);}
		}
		break;
			
	case VIEWBUT:
		if(val2 == "table"){ col = RED; col2 = DRED; col3 = WHITE; if(ov == 1) col = LRED;}
		else{
			if(val2 == "transfilt" || val2 == "pdmodel"){ col = DGREEN; col2 = DDGREEN; col3 = WHITE; if(ov == 1) col = GREEN;}
			else{ col = BLUE; col2 = DBLUE; col3 = WHITE; if(ov == 1) col = LBLUE;}
		}
		
		drawroundrect(x,y,dx,dy,7,col,col2);
		centertext(text,x+dx/2,y+17,"18px arial",col3); 
		break;
		
	case DONEBUT:
		col = LBLUE; col2 = BLUE; col3 = WHITE; if(ov == 1){ col = BLUE; col2 = DBLUE;}
		drawroundrect(x,y,dx,dy,7,col,col2);
		centertext(text,x+dx/2,y+16,"bold 17px arial",col3); 
		break;
		
	case CANCELBUT:
		col = LRED; col2 = RED; col3 = WHITE; if(ov == 1){ col = RED; col2 = DRED;}
		drawroundrect(x,y,dx,dy,7,col,col2);
		centertext(text,x+dx/2,y+16,"bold 17px arial",col3); 
		break;
		
	case GREENCANBUT:
		col = "#55FF55"; col2 = GREEN; col3 = WHITE; if(ov == 1){ col = DGREEN; col2 = DDGREEN;}
		drawroundrect(x,y,dx,dy,7,col,col2);
		centertext(text,x+dx/2,y+16,"bold 17px arial",col3); 
		break;
		
	case CLEARBUT:
		if(modelsetup == 1){ col = WHITE; col2 = RED; col3 = RED; if(ov == 1) col = LLRED;}
		else{ col = WHITE; col2 = DGREEN; col3 = DGREEN; if(ov == 1) col = LLGREEN;}
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2,y+21, "bold 18px arial",col3);
		break;
	
	case SOURCEBUT:
		col = BLACK; if(ov == 1) col = DGREY;
		fillcircle(x,y,15,col,DGREY,NORMLINE);
		centertext(text,x,y+11, "bold 30px arial",WHITE);
		break;
		
	case SOURCEMINIBUT:
		drawline(x,y,x+dx,y+dy,BLACK,THICKLINE,val);
		drawarrow(x+dx,y+dy,x,y,dy*0.4,BLACK);	
		break;
		
	case ZOOMINBUT:
	case ZOOMOUTBUT:
		cv.globalAlpha = 0.95;
		drawroundrect(x+1,y+1,dx-2,dy-2,4,WHITE,WHITE);
		cv.globalAlpha = 1;
	
		r = 8; r2 = 5; r3 = 23;
		col = BLACK; if(ov == 1) col = DGREY;
		
		xx = x+dx-r; yy = y+r;
		circle(xx,yy,r,col,THICKLINE);
		drawline(xx-r2,yy,xx+r2,yy,col,THICKLINE);
		if(canbuttype[i] == ZOOMINBUT) drawline(xx,yy-r2,xx,yy+r2,col,THICKLINE);
		th = 2.2;
		drawline(xx+r*Math.cos(th),yy+r*Math.sin(th),xx+r3*Math.cos(th),yy+r3*Math.sin(th),col,VTHICKLINE);
		break;
	
	case TRANSPBUT:
		if(ov == 1 && addtransmode == 0){
			drawroundrect(x,y,dx,dy,Math.floor(0.5*dy),GREY,GREY);
			if(text == "+") centertext(text,x+dx/2,y+dy/2+11, "bold 30px arial",WHITE);
			else centertext(text,x+dx/2,y+dy/2+8, "bold 30px arial",WHITE);
		}
		break;
	
	case REQUESTBUT:
		centertext(text,x+dx/2,y+16,"bold 16px arial",BLACK); 
		break;
		
	case CROSSBUT:
		col = RED; if(ov == 1) col = LRED;
		drawline(x+2,y+2,x+dx-2,y+dy-2,col,VTHICKLINE);
		drawline(x+dx-2,y+2,x+2,y+dy-2,col,VTHICKLINE);
		break;
		
	case TRANSADDBUT:
		if(secbutsel == val){
			plottext(text,x+16,y+16,"bold 14px arial",BLACK); 
			drawline(x+1,y+9,x+6,y+15,BLACK,VTHICKLINE);
			drawline(x+11,y+9,x+6,y+15,BLACK,VTHICKLINE);
		}
		break;
		
	case ADDPARAMBUT:
		col = WHITE; col2 = DBLUE; col3 = DBLUE; if(ov == 1){ col2 = DRED; col3 = RED;}
		drawroundrect(x,y,dx,dy,7,col,col2);
		centertext(text,x+dx/2,y+16,"17px arial",col3); 
		break;
		
	case DELBUT:
		col = LRED; col2 = RED; col3 = WHITE; if(ov == 1) col3 = RED;
		drawroundrect(x,y,dx,dy,7,col,col2);
		centertext(text,x+dx/2,y+16,"17px arial",col3); 
		break;

	case DELETEBUT:
		col = RED; if(ov == 1) col = LRED;
		drawline(x+1,y+1,x+dx-1,y+dy-1,col,THICKLINE);
		drawline(x+1,y+dy-1,x+dx-1,y+1,col,THICKLINE);
		break;
	
	case FUNCTBUT:
		col = LBLUE; if(ov == 1) col = LLBLUE;
		drawroundrect(x,y,dx,dy,4,col,darkcol(col));
		centertext(functi[val],x+dx/2,y+18,"17px arial",WHITE); 
		break;

	case OPBUT:
		col = LRED; if(ov == 1) col = LLRED;
		drawroundrect(x,y,dx,dy,4,col,darkcol(col));
		centertext(opbut[val],x+dx/2,y+18,"17px arial",WHITE); 
		break;
		
	case NUMBUT:
		col = LORANGE; if(ov == 1) col = LLORANGE;
		drawroundrect(x,y,dx,dy,4,col,darkcol(col));
		centertext(numbut[val],x+dx/2,y+18,"17px arial",WHITE); 
		break;
		
	case GREEKBUT:
		col = LORANGE; if(ov == 1) col = LLORANGE;
		drawroundrect(x,y,dx,dy,4,col,darkcol(col));
		centertext(greek[val],x+dx/2,y+18,"17px arial",WHITE); 
		break;
		
	case ALPBUT:
		col = LORANGE; if(ov == 1) col = LLORANGE;
		drawroundrect(x,y,dx,dy,4,col,darkcol(col));
		centertext(alp[val],x+dx/2,y+18,"17px arial",WHITE); 
		break;
		
	case ERRMSGBUT: 
		plottext(text,x+17,y+15,"17px arial",RED);
		if(val == 0){ drawline(x,y+4,x+10,y+4+10,RED,THICKLINE); drawline(x,y+4+10,x+10,y+4,RED,THICKLINE);}
		break;
	
	case TEXTBUT:
		col = DDGREY; if(ov == 1) col = GREY;
		plottext(text,x,y+14,"bold 16px arial",col);
		break;
		
	case PRTITLEBUT:
		col = DDGREY;
		plottext(text,x,y+14,"bold 16px arial",col);
		break;
		
	case INDBUT:
		col = DDGREY; if(ov == 1) col = GREY;
		plottext(text,x,y+14,"bold 16px arial",col);
		break;
		
	case LINKBUT:
		col = DDBLUE; if(ov == 1) col = BLUE;
		plottext(text,x,y+14,"14px arial",col);
		break;
		
	case TEXTBUT2:
		plottext(text,x+4,y+18,"17px arial",DDGREY);
		break;
		
	case EDITDBUT:
		if(ov == 1) fillrect(x, y, dx, dy, LLGREEN);
		plottext(text,x+4,y+18,tablefont,DGREEN);
		break;
		
	case MINMAXBUT:
		if(ov == 1) fillrect(x, y, dx, dy, LLGREEN);
		plottext(text,x+4,y+22,tablefont,DGREEN);
		break;
	
	case NOTACTIVEBUT:
		plottext(text,x+4,y+22,tablefont,DPURPLE);
		break;
	
	case SELPARAMBUT:
		if(ov == 1) fillrect(x, y, dx, dy, LLGREEN);
		plottext(text,x+4,y+22,tablefont,DGREEN);
		break;
		
	case NOTSELPARAMBUT:
		plottext(text,x+4,y+22,tablefont,BLACK);
		break;
		
	case MINMAXBIGBUT:
		if(ov == 1) fillrect(x, y, dx, dy, LLGREEN);
		plottext(text,x+4,y+22,"bold 20px Times",BLACK);
		break;
		
	case PDBUT:
		if(ov == 1) fillrect(x, y, dx, dy, LLGREEN);
		plottext(text,x+4,y+22,tablefont,DGREEN);
		break;
		
	case GENUSERBUT:
		if(ov == 1) fillrect(x, y, dx, dy, LLGREEN);
		plottext(text,x+4,y+22,tablefont,DGREEN,dx);
		break;
		
	case EQBUT:
		if(ov == 1) fillrect(x, y, dx, dy, LLGREEN);
		plottext(text,x+4,y+22,tablefont,DGREEN);
		break;
		
	case MINMAXHEADBUT:
		plottext(text,x+4,y+16,tableheadfont,DDGREEN);
		break;
		
	case CANRADIOBUT:
		col=DGREY;
		switch(val2){
		case CANRADIORATE: selval = cla[selectbubval].tra[selectbubval2].ratetime; break;
		case CANRADIOCLASS: selval = datatemp.cl; if(addingdata == 2) ov = 0; break;
		case CANRADIOCHECK: selval = datatemp.type; if(addingdata == 2) ov = 0; break;
		case CANRADIOTESTRES: selval = datatemp.postestres[Math.floor(val/1000)]; val = val%1000; if(addingdata == 2) ov = 0; break;
		case CANRADIODEFINE: selval = simpopinitset; break;
		case CANRADIOSIMTY: selval = simty; break;
		case CANRADIOTERM: selval = termtype; break;
		case CANRADIOXAXIS: selval = xaxisauto; break;
		case CANRADIOYAXIS: selval = yaxisauto; break;
		case CANRADIOGENT: selval = datagenttype; break;
		case CANRADIOGENT2: selval = datagenttype2; break;
		case CANRADIOCHECK2: selval = datagen[sellist[drawgennum]].type; break;
		case CANRADIOWHICH: selval = whichind; break;
		case CANRADIOPD: selval = obspd; break;
		case CANRADIOPD2: selval = datatemp.obspd; break;
		case CANRADIOPD3: selval = datatemp.pdsame; break;
		case CANRADIOMANUAL: selval = siminit; break;
		case CANRADIOERRBAR: selval = errbar; break;
		case CANRADIOINDUO: selval = induoflag; break;
		}

		if(val == selval){ if(addingdata != 2)  col2 = BLACK; else col2 = GREY;}
		else{
			if(ov == 1) col2 = DGREY;
			else col2 = WHITE
		}

		r = 7;
		drawroundrect(x+1,y+2,2*r,2*r,r,WHITE,col);

		if(val == selval){ if(addingdata != 2) col3 = BLACK; else col3 = GREY;}
		else{
			col3 = GREY;  if(ov == 1) col3 = DGREY;
		}
    
		plottext(text,x+21,y+14,"bold 16px arial",col3);
 
		if(col2 != WHITE){
			x += 2; y += 4; r -= 2;
			drawroundrect(x+1,y,2*r,2*r,r,col2,col2);
		}
		break;
		
	case SELPARBUT:
		col = WHITE; if(ov == 1) col = LLBLUE; 
		drawroundrect(x,y,dx,dy,10,col,BLUE);
		break;
		
	case ADDBUT: case ADDBUT3:
		if(text == "Compartment" || text == "Transition"){
			cv.globalAlpha = 0.8;
			drawroundrect(x+1,y+1,dx-2,dy-2,4,WHITE,WHITE);
			cv.globalAlpha = 1;
		}
		
		col = BLACK; if(ov == 1) col = GREY;
		x = canbutx[i]+8; y = canbuty[i]+dy/2; r = 6; r2 = 4;
		drawroundrect(x-r,y-r,2*r,2*r,r,col,col);
		drawline(x-r2,y,x+r2,y,WHITE,THICKLINE);
		drawline(x,y-r2,x,y+r2,WHITE,THICKLINE);
		plottext(text,x+11,canbuty[i]+16,addfont,col); 
		break;
	
	case BRACKETBUT:
		drawline(x,y,x,y+dy,BLACK,THICKLINE);
		drawline(x,y,x+dx,y,BLACK,THICKLINE);
		drawline(x,y+dy,x+dx,y+dy,BLACK,THICKLINE);
		drawline(x-8,y+Math.floor(dy/2),x,y+Math.floor(dy/2),BLACK,THICKLINE);	
		break;
		
	case TIMELINEBUT:
		righttext(text,x-10,y+4,"bold 14px ariel","#cc2222"); 

		if((page == INFERENCEPAGE && pagesub[INFERENCEPAGE] == 2) || page == SIMULATEPAGE) calcindprob(val,val2,x,y,dx); 
		else drawline(x,y,x+dx,y,BLACK,THICKLINE);
		break;
		
	case EVBUT5:
		col = BLACK; if(ov == 1) col = GREY;
		drawline(x+dx/2,y,x+dx/2,y+dy,col,THICKLINE,2);
		drawline(x+dx/2,y+3,x+dx+8,y+3,col,THICKLINE);
		drawarrow(x+dx+8,y+3,x-10,y+3,7,col);
		break;
		
	case EVBUT6:
		col = BLACK; if(ov == 1) col = GREY;
		drawline(x+dx/2,y,x+dx/2,y+dy,col,THICKLINE,2);
		drawline(x+dx/2,y+3,x-8,y+3,col,THICKLINE);
		drawarrow(x-8,y+3,x+dx+10,y+3,7,col);
		break;
		
	case EVBUT:
	case EVBUT4:
		var ev;
		if(page == INFERENCEPAGE){
			if(pagesub[page] == 0) ev = inddata.ind[val].cl[val2].ev[text];
			else ev = infres.inddata.ind[val].cl[val2].ev[text];
		}
		else ev = indsim.ind[val].cl[val2].ev[text];
		
		col = ev.col;
		jmax = col.length;
		if(jmax == 0) fillcircle(x+dx/2,y+dy/2,dx/2,BLACK,BLACK,NORMLINE);
		else{
			for(j = 0; j < jmax; j++){
				cv.beginPath();
				cv.moveTo(x+dx/2,y+dy/2);
				
				cv.arc(x+dx/2,y+dy/2,dx/2,2*Math.PI*(0.25+j/jmax),2*Math.PI*(0.25+(j+1)/jmax));
				switch(col[j]){
				case -1: co = BLACK; break;
				case -2: co = WHITE; break;
				case -3: co = GREY; break;
				default: co = cla[val2].comp[col[j]].col; break;
				}
			
				if(ov == 1) co = darkcol(co); 
				cv.fillStyle = co;
				cv.fill();
			}
		}
		circle(x+dx/2,y+dy/2,dx/2,BLACK,NORMLINE);
		break;
	
	case EVBUT7:
		fillcircle(x+dx/2,y+dy/2,dx/2,WHITE,BLACK,NORMLINE);
		break;
		
	case EVSOURCEBUT:
		col = ORANGE; if(ov == 1) col = LORANGE;
		drawline(x+dx/2,y,x+dx/2,y+dy-18,col,VTHICKLINE);
		centertext("+",x+6,y+dy-5,TICKFONT,col); 
		break;
		
	case EVSINKBUT:
		col = ORANGE; if(ov == 1) col = LORANGE;
		drawline(x+dx/2,y,x+dx/2,y+dy-18,col,VTHICKLINE);
		centertext("-",x+6,y+dy-5,TICKFONT,col); 
		break;
		
	case EVBUT2:
		if(page == INFERENCEPAGE){
			if(pagesub[page] == 0) ev = inddata.ind[val].cl[val2].ev[text];
			else ev = infres.inddata.ind[val].cl[val2].ev[text];
		}
		else ev = indsim.ind[val].cl[val2].ev[text];
		
		col = ev.col;
		co = cla[val2].comp[col[0]].col; if(ov == 1) co = darkcol(co); 
		fillrect(x,y,dx/2,dy,co);
		 
	    co = cla[val2].comp[col[1]].col; if(ov == 1) co = darkcol(co); 
		fillrect(x+dx/2,y,dx/2,dy,co);
		drawrect(x,y,dx,dy,BLACK,NORMLINE);
		break;
		
	case EVBUT3:
		if(page == INFERENCEPAGE){
			if(pagesub[page] == 0) ev = inddata.ind[val].cl[val2].ev[text];
			else  ev = infres.inddata.ind[val].cl[val2].ev[text];
		}
		else ev = indsim.ind[val].cl[val2].ev[text];
		col = ev.col;
		co = cla[val2].comp[col[0]].col; if(ov == 1) co = darkcol(co); 
		fillrect(x,y,dx,dy,co);
		drawrect(x,y,dx,dy,BLACK,NORMLINE);
		break;
		
	case CANSMALLTEXTBUT:
		plottext(text,x+2,y+15,"12px arial",val);
		break;
		
	case RESULTBUT:
		cv.drawImage(resultcan,0,0,dx,dy,x,y,dx,dy);
		drawline(x,y+dy,x,y-10,BLACK,THICKLINE);
		drawarrow(x,y-20,x,y+dy,15,BLACK);
		drawline(x,y+dy,x+dx+10,y+dy,BLACK,THICKLINE);
		drawarrow(x+dx+20,y+dy,x,y+dy,15,BLACK);
		break;
	
	case MODELPICBUT:
		cv.drawImage(resultcan,0,0,dx,dy,x,y,dx,dy);
		break;
		
	case SPEECHBUT:
		if(ov == 1) cv.globalAlpha = 0.3;
		cv.drawImage(speechpic,0,0,dx,dy,x,y,dx,dy);
		if((page == MODELPAGE && cla[pagesub[MODELPAGE]].descon == 1) || (page == INFERENCEPAGE && datanoteon == 1)){
			r = 4;
			drawline(x+dx/2-r,y+dy/2-r-2,x+dx/2+r,y+dy/2+r-2,RED,THICKLINE);
			drawline(x+dx/2-r,y+dy/2+r-2,x+dx/2+r,y+dy/2-r-2,RED,THICKLINE);
		}		
		if(ov == 1) cv.globalAlpha = 1;
		break;
	
	case PARACANBUT:
		alignparagraph(text,dx-10);
		col = DDBLUE;
		var j; 
		fillcircle(x-4,y+10,3,BLACK,BLACK,NORMLINE);
		
		for(j = 0; j < nlines; j++) plottext(lines[j],x+5,y+18+linesy[j],HELPFONT,col);
		break;
		
	case PARACANBUT2:
		alignparagraph(text,dx-10);
		col = DDBLUE;
		var j; for(j = 0; j < nlines; j++) plottext(lines[j],x+5,y+18+20*j,HELPFONT,col);
		break;
	
	case PARACANBUT3:
		alignparagraph(text,dx-5,"bold 18px arial");
		col = DDBLUE;
		var j; 
		fillcircle(x-4,y+10,3,BLACK,BLACK,NORMLINE);
		
		for(j = 0; j < nlines; j++) plottext(lines[j],x+5,y+18+linesy[j],"bold 18px arial",col);
		break;
		
	
		
	case XTICKBUT:
		col = BLACK; if(ov == 1 && selectbub != XTICKBUT) col = DGREY;
		for(j = 0; j < ntickx; j++){
			xx = Math.floor(x+dx*(tickx[j]-axxmin)/(axxmax-axxmin));
			centertext(tickx[j],xx,y+dy/2+8,TICKFONT,col); 
			drawline(xx,y,xx,y-10,BLACK,NORMLINE); 
			if(val == -2) drawline(xx,y-graphdy,xx,y-graphdy+10,BLACK,NORMLINE); 
		}
		break;
		
	 case YTICKTRBUT:
		col = BLACK; if(ov == 1 && selectbub != YTICKTRBUT) col = DGREY;
		for(j = 0; j < nticky; j++){
			yy = Math.floor(y+dy-dy*(ticky[j]-axymin)/(axymax-axymin));
			righttext(ticky[j],x+dx-5,yy+6,TICKFONT,col); 
			drawline(x+dx,yy,x+dx+10,yy,BLACK,NORMLINE); 
		}
		break;
	
	case BURNINBUT:
		col = BLACK; 
		drawline(x+dx/2,y,x+dx/2,y+dy,col,NORMLINE,2); 
		if(val == -1) plotangletext(text,x+dx/2+5,y+5,-Math.PI/2,MENUFONTSMALL,col);
		else plotangletext(text,x+dx/2+5-19,y+5,-Math.PI/2,MENUFONTSMALL,col);
		break;
	
	case TTRANSBUT:
		col = BLACK; 
		drawline(x+dx/2,y+77,x+dx/2,y+dy,col,NORMLINE,3); 
		plotangletext("Transition",x+dx/20-4,y+5,-Math.PI/2,MENUFONTSMALL,col);
		break;
		
	case HISTOBUT:
		col = BLACK; 
		plotxlabel(text,x+dx/2,y+dy/2+4,"22px times",col);
		break;
		
	case XLABELBUT:
	 	col=BLACK; if(ov == 1){ col = GREY;}
		if(chooseaxis == 0) col = RED;
		plotxlabel(text,x+dx/2,y+dy/2+4,LABELFONT,col);
		break;
		
	case YLABELBUT:
		col=BLACK; if(ov == 1){ col = GREY;}
		if(chooseaxis == 1) col = RED;
		plotylabel(text,x+dx/2+5,y+dy/2,LABELFONT,col);
		break;
	
	case HOZLABBUT:
		plotxlabel(text,x+dx/2,y+dy/2+4,CORFONT,BLACK);
		break;
		
	case VERTLABBUT:
		plotylabel(text,x+dx/2+5,y+dy/2,CORFONT,BLACK);
		break;
		
	case SLIDERBACKBUT:
		plottext(text,x-60,y+16,MENUFONTSMALL,DGREY);
		drawline(x,y+dy/2,x+dx,y+dy/2,GREY,THICKLINE);
		break;

	case SLIDERBUT:
		col = WHITE; if(ov == 1) col = LGREY; if(drag == 10) col = RED;
		fillrect(x,y,slidedx,dy,col);
		drawrect(x,y,slidedx,dy,GREY,NORMLINE);
		break;
		
	case LABBUT:
		if(page == INFERENCEPAGE && pagesub[page] == 2 && pagesubsub[page][2] == 1) val2 = 0;
		
		drawline(x,y+dy/2,x+dx,y+dy/2,val,MTHICKLINE,val2);
		plottext(text,x+dx+5,y+12,KEYFONT,BLACK); 
		break;
		
	case WHITERECTBUT:
		fillrect(x,y,dx,dy,WHITE);
		break;
		
	case INACTIVEBUT:
		cv.globalAlpha = 0.5; 
		fillrect(x,y,dx,dy,WHITE);
		cv.globalAlpha = 1; 
		break;
		
	case ARROWBUT:
		drawline(x,y,x+dx+10,y,BLACK,THICKLINE);
		drawarrow(x+dx+20,y,x,y,15,BLACK);
		break;
		
	case NEXTCANBUT:
		col = WHITE; col2 = RED; col3 = RED; if(ov == 1) col = LLRED;
		drawroundrect(x,y,dx,dy,7,col,col2);
		centertext(text,x+dx/2-5,y+16,"17px arial",col3);
		drawline(x+dx-8,y+dy/2,x+dx-14,y+5,col3,THICKLINE);
		drawline(x+dx-8,y+dy/2,x+dx-14,y+dy-5,col3,THICKLINE);
		drawline(x+dx-8-6,y+dy/2,x+dx-14-6,y+5,col3,THICKLINE);
		drawline(x+dx-8-6,y+dy/2,x+dx-14-6,y+dy-5,col3,THICKLINE);
		break;
		
	case BACKCANBUT:
		col = WHITE; col2 = RED; col3 = RED; if(ov == 1) col = LLRED;
		drawroundrect(x,y,dx,dy,7,col,col2);
		centertext(text,x+dx/2+5,y+16,"17px arial",col3);
		drawline(x+8,y+dy/2,x+14,y+5,col3,THICKLINE);
		drawline(x+8,y+dy/2,x+14,y+dy-5,col3,THICKLINE);
		drawline(x+8+6,y+dy/2,x+14+6,y+5,col3,THICKLINE);
		drawline(x+8+6,y+dy/2,x+14+6,y+dy-5,col3,THICKLINE);
		break;
	
	case SPEECHBUT2: break;
		
	case HELPICONCANBUT:
		if(val2 == -1) fillrect(x,y,dx,dy,WHITE); else fillrect(x,y,dx,dy,"#ededff"); 
		col = DRED; if(ov == 1) col = LRED;
		plottext(text,x+3,y+10,"bold 10px arial",col);
		break;
		
	case EDITSPEECHBUT: 
		col = DRED; if(ov == 1) col = LRED;
		plottext("EDIT",x+3,y+10,"bold 10px arial",col);
		break;
		
	case DESCBUT:
		cv.drawImage(descpic,x,y);
		break;
		
	case CORBUT:
		plotcorbut(text,x,y,dx,dy,ov,14);
		break;
		
	default: alertp("Error code EC3"); break;
	}
}

function addcanbutton(text,x,y,dx,dy,ac,type,val,val2)   // Adds a canvas button
{
	canbuttext[ncanbut] = text;
	canbutx[ncanbut] = x;
	canbuty[ncanbut] = y;
	canbutdx[ncanbut] = dx;
	canbutdy[ncanbut] = dy;
	canbutac[ncanbut] = ac; 	
	canbuttype[ncanbut] = type;
	canbutover[ncanbut] = -1;
	canbutval[ncanbut] = val;
	canbutval2[ncanbut] = val2;
	ncanbut++;
	
	if(type == PRTITLEBUT && val >= 0) addcanbutton("[?]",x+textwidth(text,"bold 16px arial")+1,y-2,15,20,HELPICONCANBUT,HELPICONCANBUT,val,-1);
}
