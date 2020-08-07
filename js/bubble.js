function bubble()                                         // Draws a speech bubble for selected item
{
	var x, y, w, h, dx, dy, i, j, dl, dr, du, dd, orient, hh, ww;

	cornbub = 0;

	switch(selectbub){
	case COMPBUT: bubw = 230; bubh = 199; break;
	
	case EVBUT:
		bubw = 250; 
		
		var ev;
		if(page == INFERENCEPAGE){
			if(pagesub[page] == 0) ev = inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
			else ev = infres.inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
		}
		else ev = indsim.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
		
		bubh = 86 + 30*ev.obsdata.length;
		
		wmax = 0;
		for(k = 0; k < cla[selectbubval2].ncomp; k++){
			if(ev.probexp[k] != "0"){
				st = "Pr(D|"+cla[selectbubval2].comp[k].name+") = "+ev.probexp[k];
				w = Math.floor(textwidth(st,"16px georgia"));
				if(w > wmax) wmax = w;	
				bubh += 30;
			}
		}
		wmax += 10; if(wmax < 80) wmax = 80;
		if(wmax+42 > bubw) bubw = wmax+42;
		
		datwid = 0; opwid = 50; fo = Math.floor(0.85*gdropdy)+"px arial";
		for(j = 0; j < ev.obsdata.length; j++){
			d = ev.obsdata[j]; w = textwidth(data[d].name,"14px arial")+20; if(w > datwid) datwid = w;
			for(k = 0; k < data[d].pos.length; k++){ w = Math.floor(textwidth(data[d].pos[k],fo)+30); if(w > opwid) opwid = w;}
		}
		w = datwid+opwid+20; if(w > bubw) bubw = w;
			
		if(ev.col.length == 0) bubh += 30;
		break;

	case ADDBUT: bubw = 235; bubh = 100; break;
	case PARACANBUT2: bubw = 500; bubh = 280; break;
	case EMPTYBUB: bubw = 165; bubh = 60; break;
	case NANBUB: bubw = 145; bubh = 60; break;
	case PROBBUB: bubw = 145; bubh = 60; break;
	case EVBUT2: bubw = 200; bubh = 60; break;
	case EVBUT3: bubw = 250; bubh = 140; break;
	case EVBUT4: bubw = 200; bubh = 60; break;
	case EVBUT5: bubw = 200; bubh = 60; break;
	case EVBUT6: bubw = 200; bubh = 60; break;
	case EVBUT7: bubw = 150; bubh = 60; break;
	case EVSOURCEBUT: case EVSINKBUT: bubw = 200; bubh = 60; break;
	case COMPSMALLBUT: bubw = 190; bubh = 95; break;
	case COMPSMALLPERBUT: bubw = 190; bubh = 95; break;
	case COMPSMALLBUT2: bubw = 260; bubh = 102; break;
	case TESTNAMEBUT: bubw = 260; bubh = 102; break;
	case EQBUT: bubw = 355; bubh = 195+52; bubh = 115+52;  break;
	case PDBUT: bubw = 365; bubh = 195+47; break;
	case TRANSPLBUT:
		bubw = 365;
		switch(cla[selectbubval].tra[selectbubval2].type){
		case "Exponential": case "Source": case "Sink": bubh = 195+77; break;
		case "Gamma": case "Weibull": bubh = 195+77+47; break;
		case "Fixed": bubh = 185; break;
		case "Set": case "Grow": bubw = 275; bubh = 65; break;
		}
		//bubh -= 26;
		break; 
	case TRANSPLBIGBUB: bubw = canvasdx-40; bubh = canvasdy-80; cornbub = 1; break;
	case MINMAXBUT: case MINMAXBIGBUT: bubw = 205; bubh = 165; break;
	case SELPARAMBUT: bubw = 365; bubh = 195+47; break;
	case EDITDBUT: bubw = 195; bubh = 102; break;
	case CONVERTDATEBUB: bubw = 275; bubh = 102; break;
	case MULTISTATEBUB: bubw = 275; bubh = 102; break;
	case TABLEBUT: bubw = 275; bubh = 95; break;
	case TABLEDROPBUT: bubw = 200; bubh = 60; break;
	case TABLEBUT3: bubw = 275; bubh = 102; break;
	case TABLEHEADBUT: bubw = 230; bubh = 156; break;
	case SEARCHBUB: bubw = 230; bubh = 98; break;
	case SEARCHRESBUB: bubw = 230; bubh = 118; break;
	case SPEECHBUT2: bubw = 600; bubh = 275; break;
	case REPLACEBUB: bubw = 230; bubh = 153; break;
	case REPLACEDONEBUB: bubw = 230; bubh = 63; break;
	case DELROWSDONEBUB: bubw = 230; bubh = 63; break;
	case DELROWSBUB: bubw = 230; bubh = 98; break;
	case ERRBUB: bubw = textwidth(errmsg,"17px arial")+50; bubh = 60; break;
	case ADDCLASSBUB: bubw = 180; bubh = 57; selectbuby = ymodel+Math.floor(dysub*0.5); break; 
	case EDITCLASSBUB: 
		bubw = 195; bubh = 102; selectbuby = ymodel+Math.floor(dysub*(0.5+selectbubval)); 
		break; 
	case EDITAGETIMEBUB: bubw = 195; bubh = 57; break; 
	case XTICKBUT: bubw = 190; bubh = 202; if(xaxisauto == 0) bubh = 185; else bubh = 62; break; 
	case YTICKTRBUT: bubw = 190; bubh = 202; if(yaxisauto == 0) bubh = 185; else bubh = 62; break; 
	case EDITBUB:
		bubw = 140; bubh = 57;  selectbubdx = 0; selectbubdy =0;
		if(pagesub[page] >= 0){ selectbubx = width-205; selectbuby = height-45;}
		else{ selectbubx = canvaswidth-45; selectbuby = canvasheight;}
		break;
	case GENDATABUB: bubw = 200; bubh = 300; selectbubx = 300; selectbuby = 405; selectbubdx = 0; selectbubdy =0; break;
	case GENUSERBUT: bubw = 455; bubh = 115; break;
	case COMPPOPOUTPUTBUT: bubw = 195; bubh = 95; break; 
	case COMPPOPINITBUT: case COMPFRACINITBUT: bubw = 195; bubh = 95; break; 
	case COMPLOADINITBUT: bubw = 195; bubh = 95; break; 
	case JOINBUB: bubw = 165; bubh = 79; break;
	}
	
	if(errmsg != "" && selectbub != ERRBUB){
		bubh += 15;
		
		errmsgline=[];
		if(textwidth(errmsg,"17px arial")+28 < bubw) errmsgline.push(errmsg);
		else{
			var text = errmsg
			do{
				j = text.length;
				while(textwidth(text.substr(0,j),"17px arial") > bubw-28){
					j--;
					while(j > 0 && text.substr(j,1) != " ") j--;
					if(j == 0) break;
				}
				if(j == 0) j = text.length/2;
				errmsgline.push(text.substr(0,j));
				bubh += 20;
				text = text.substring(j);
			}while(text.length > 0);
			bubh -= 20;
		}
	}	
	
	x = selectbubx; y = selectbuby; dx = selectbubdx; dy = selectbubdy; w = bubw; h = bubh; 
	switch(selectbub){
	case TABLEBUT: case TABLEDROPBUT: case TABLEBUT3: case TABLEHEADBUT: case MINMAXBUT: case MINMAXBIGBUT: 
	case SELPARAMBUT: case EDITDBUT: case GENUSERBUT: case EQBUT: case PARACANBUT2: case SEARCHRESBUB: 
	case JOINBUB: case SPEECHBUT2: case EMPTYBUB: case NANBUB: case PROBBUB: case PDBUT: case SEARCHBUB:
	case REPLACEBUB: case REPLACEDONEBUB: case DELROWSBUB: case DELROWSDONEBUB:
		addcanbutton("Bubble",x,y,dx,dy,-1,BUBSELBUT,-1,-1);
		break;
	}
	
	hh = canvasdy; ww = canvasdx;

	if(cornbub == 0){
		i = x+dx/2; j = y+dy/2;      	// works out the best way to orient the bubble
		
		i = x; j = y+dy/2; dl = i*i+(j-hh/2)*(j-hh/2);
		i = x+dx; j = y+dy/2; dr = (i-ww)*(i-ww)+(j-hh/2)*(j-hh/2);
		i = x+dx/2; j = y; du = (i-ww/2)*(i-ww/2)+j*j;
		i = x+dx/2; j = y+dy; dd = (i-ww/2)*(i-ww/2)+(j-hh)*(j-hh);

		var loop, pos=[];
		for(orient = 0; orient < 4; orient++) pos[orient] = placebub(orient,x,y,dx,dy,w,h,ww,hh)
	
		orient = -1;
		for(loop = 0; loop < 2; loop++){	
			if(dr >= dl && dr >= du && dr >= dd && !(loop == 0 && pos[RIGHT] == 1)) orient = RIGHT;
			else{
				if(dl >= dr && dl >= du && dl >= dd && !(loop == 0 && pos[LEFT] == 1)) orient = LEFT;
				else{
					if(dd >= dl && dd >= du && dd >= dr && !(loop == 0 && pos[DOWN] == 1)) orient = DOWN;
					else{
						if(!(loop == 0 && pos[UP] == 1)) orient = UP;
					}
				}
			}
			if(orient >= 0) break; 
		}
		
		if(page == INFERENCEPAGE && pagesub[page] == 1){ if(du < dd) orient = DOWN; else orient = UP;}
	
		if(page == MODELPAGE && pagesub[page] == -3){ if(du < dd) orient = DOWN; else orient = UP;}
		
		if(selectbub == ADDCLASSBUB) orient = RIGHT;
		
		if(orient == DOWN){
			if(y+dy+bubtag+h > hh){
				if(dl < dr) orient = RIGHT; else orient = LEFT;
			}
		}

		if(orient == UP){
			if(y-bubtag-h < 0){
				if(dl < dr) orient = RIGHT; else orient = LEFT;
			}
		}

		if(selectbub == EDITBUB) orient = UP;
		if(selectbub == JOINBUB) orient = DOWN;
		if(selectbub == SPEECHBUT2) orient = DOWN; 
		if(selectbub == ADDBUT) orient = UP;
		if(selectbub == GENUSERBUT) orient = DOWN; 

		placebub(orient,x,y,dx,dy,w,h,ww,hh);
		
		if(buby < 1){ ddy = 1 - buby; buby += ddy; buby1 += ddy; buby2 += ddy; buby3 += ddy;}
		if(buby+h > hh-1){ ddy = hh-1-h - buby; buby += ddy; buby1 += ddy; buby2 += ddy; buby3 += ddy;}
		if(bubx < 1){ ddx = 1 - bubx; bubx += ddx; bubx1 += ddx; bubx2 += ddx; bubx3 += ddx;}
		if(bubx+w > ww-1){ ddx = ww-1-w - bubx; bubx += ddx; bubx1 += ddx; bubx2 += ddx; bubx3 += ddx;}
		
		bubw = w; bubh = h;
	}
	else{ bubx = 20; buby = 30;}
	
	addcanbutton("Bubble",bubx,buby,w,h,BUBBLEBUT,BUBBLEBUT,-1,-1);

	addcanbutton("Bubble Close",bubx+w-22,buby+6,16,16,BUBBLEOFFBUT,BUBBLEOFFBUT,-1,-1);

	cursx = bubx+bubmar; cursy = buby + bubmar;

	bubblecontent(w);
}

function placebub(orient,x,y,dx,dy,w,h,ww,hh)             // Places the speech bubble on the page
{
	var fl = 0;
	
	switch(orient){
	case LEFT: 
		bubx = x-w-bubtag;
		bubx1 = x; buby1 = y+dy/2;
		bubx2 = bubx+w; buby2 = buby1+5;
		bubx3 = bubx+w; buby3 = buby1+bubtag2;
		buby = buby1-30;
		if(buby < 5){ buby = 5; fl = 1;}
		if(buby+h > hh-5){ buby = hh-5-h; fl = 1;}
		if(buby+bubh-20 < buby3){ buby = buby3-bubh+20; fl = 1;}
		break; 
		
	case RIGHT:
		bubx = x+dx+bubtag;
		bubx1 = x+dx; buby1 = y+dy/2;
		bubx2 = bubx; buby2 = buby1+5;
		bubx3 = bubx; buby3 = buby1+bubtag2;
		buby = buby1-30;
		if(buby < 5){ buby = 5; fl = 1;}
		if(buby+h > hh-5){ buby = hh-5-h; fl = 1;}
		if(buby+bubh-20 < buby3){ buby = buby3-bubh+20; fl = 1;}
		break; 
		
	case UP:
		buby = y-bubtag-h;
		bubx1 = x+dx/2; buby1 = y;
		bubx2 = bubx1+5; buby2 = y-bubtag;
		bubx3 = bubx1+bubtag2; buby3 = y-bubtag;
		bubx = bubx1-30;
		if(bubx < 5){ bubx = 5; fl = 1;}
		if(bubx+w > ww-5){ bubx = ww-5-w; fl = 1;}
		break;
		
	case DOWN:
		buby = y+dy+bubtag;
		bubx1 = x+dx/2; buby1 = y+dy;
		bubx2 = bubx1+5; buby2 = buby;
		bubx3 = bubx1+bubtag2; buby3 = buby;
		bubx = bubx1-30;
		if(bubx < 5){ bubx = 5; fl = 1;}
		if(bubx+w > ww-5){ bubx = ww-5-w; fl = 1;}
		break;
	}
	return fl;
}

function addbubtext(text)                                // Adds text to a bubble
{
  addcanbutton(text,cursx,cursy,0,0,-1,BUBBLETEXTBUT,-1,-1);
  cursy += 20;
}

function addinput(text,w)                                // Adds input area to bubble
{
	if(inpon == 0){
		text = ""+text;
		st = "<input id='inp' autocomplete='off' style='text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+9)+"px; left:"+(canx+cursx+12)+"px; width:"+w+"px; height:25px; background=#ffffff;' value='"+text+"' onblur='setTimeout(function() {if(inpon == 1) ById(\"inp\").focus();}, 0);' onkeypress='if(event.keyCode == 13){ if(selectbub == SEARCHBUB || selectbub == SEARCHRESBUB) buboff(0); else{ if(selectbub == DELROWSBUB) dodelrows(); else{ buboff(0); buttoninit();}}}' type='text' onfocus='foc=\"inp\";'/>";
		ById("add").innerHTML = st;
		a = ById("inp"); a.focus(); a.selectionStart = 0;  a.selectionEnd = text.length;
		inpon = 1;
	}
	else{
		a = ById("inp"); a.style = "text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+9)+"px; left:"+(canx+cursx+12)+"px; width:"+w+"px; height:25px; background=#ffffff;";
	}
	 
	cursy += 27;
}

function addtextarea(text,w,h,text2,w2)                  // Adds a text area to a bubble
{
	if(text2 != -1){
		if(h == 104) addcanbutton("Expand",cursx+w-50,cursy-20,50,15,EXPANDAC,LINKBUT,-1,-1);
		else addcanbutton("Reduce",cursx+w-50,cursy-20,50,15,REDUCEAC,LINKBUT,-1,-1);
	}
	
	if(inpon == 0){
		text = ""+text;
		st = "<textarea id='inptext' style='text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+9)+"px; left:"+(canx+cursx+12)+"px; width:"+w+"px; height:"+h+"px; background=#ffffff;' onfocus='foc=\"inptext\";' ";

		st += "type='text'>"+text+"</textarea>";
		
		if(w2 > 0){
			st += "<input id='inp' autocomplete='off' style='text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+h+37)+"px; left:"+(canx+cursx+12)+"px; width:"+w2+"px; height:25px; background=#ffffff;' value='"+text2+"' type='text' onfocus='foc=\"inp\";' />";
		}
		
		ById("add").innerHTML = st;
		a = ById("inptext"); a.focus();
		inpon = 1;
	}
	else{
		a = ById("inptext"); a.style = "text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+9)+"px; left:"+(canx+cursx+12)+"px; width:"+w+"px; height:"+h+"px; background=#ffffff;";
		
		if(w2 > 0){
			a = ById("inp"); a.style = "text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+h+37)+"px; left:"+(canx+cursx+12)+"px; width:"+w2+"px; height:25px; background=#ffffff;";
		}
	}
	 
	cursy += h+2;
}

function buboff(close)                                   // Turns a bubble off
{
	var cl, i, j, st, tr;
	
	errmsg = "";
	if(close == 1){selectbubst = selectbub; selectbub = -1; ById("add").innerHTML = ""; inpon = 0; return;}
	
	cl = selectbubval; i = selectbubval2;
	switch(selectbub){
	case SPEECHBUT2:
		st = ById("inptext").value;
		if(st == "") errmsg = "Must enter text";
		else descnote = st;
		break;
	
	case EQBUT:
		setderive(cl,ById("inp").value);
		break;
		
	case EDITDBUT:
		st = ById("inp").value; 
		if(st == "") errmsg = "Must enter a value";
		else datagen[clgl].posedit[selectbubval] = st;
		break;
	
	case PARACANBUT2:
		st = ById("inptext").value;
		if(st == "") errmsg = "Must enter text";
		else{
			switch(page){
			case MODELPAGE:  cla[pagesub[MODELPAGE]].desc = st; break;
			case INFERENCEPAGE: datanote = st; break;
			}
		}
		break;
	
	case TRANSPLBUT:	// checks that rates are correctly entered
	case TRANSPLBIGBUB:
		switch(cla[cl].tra[i].type){
		case "Exponential":	case "Source":	case "Sink":
			cla[cl].tra[i].rate = ById("inptext").value; 
			break;
			
		case "Gamma":
			cla[cl].tra[i].mean = ById("inptext").value;
			cla[cl].tra[i].shape = ById("inp").value;
			break;
		
		case "Weibull":
			cla[cl].tra[i].lam = ById("inptext").value;
			cla[cl].tra[i].k = ById("inp").value;
			break;
			
		case "Fixed":
			cla[cl].tra[i].time = ById("inp").value;
			break;			
		}
		checkrate(cl,i); 
		collectvariables();
		break;

	case COMPBUT:
		st = ById("inp").value;
		
		checkname(st,"");
		if(errmsg == ""){	
			for(cl2 = 0; cl2 < ncla; cl2++){
				for(j = 0; j < cla[cl2].ncomp; j++){
					if(cla[cl2].comp[j].name == st && !(cl2 == cl && j == i)) errmsg = "Same name as another compartment";
				}
			}
		}
	
		if(errmsg == ""){
			c = cla[cl].comp[i];
			old = c.name;
			if(old != st){
				paramchname(param,old,st);
				paramchname(paramsim,old,st);
		
				c.name = st;
				adjustrate(old,st);
				adjustdata(old,st,"comp");
				gx = c.x + c.w/2;;
				dh = c.h;		
				c.w = Math.floor(textwidth(c.name, Math.floor(dh*0.7)+"px georgia"))+dh;	
				c.x = gx - c.w/2;
				for(tr = 0; tr < cla[cl].ntra; tr++) findpline(cl,tr);
			}
		}
		break;
	
	case EDITCLASSBUB:
		st = ById("inp").value;
		
		checkname(st,cla[cl].name);
		if(errmsg == ""){	
			for(cl2 = 0; cl2 < ncla; cl2++){
				for(j = 0; j < cla[cl2].ncomp; j++){
					if(cla[cl2].comp[j].name == st) errmsg = "Same name as one of the compartments";
				}
			}
		}
		
		if(errmsg == ""){
			paramchname(param,cla[cl].name,st);
			paramchname(paramsim,cla[cl].name,st);

			old = cla[cl].name;
			cla[cl].name = st;
			adjustrate(old,st);
			adjustdata(old,st,"class");   
		}
		break;
		
	case COMPPOPINITBUT:
		st = ById("inp").value;
		if(st == "" || isNaN(st)) errmsg = "Must enter a number";
		else{
			st = parseFloat(st);
			if(st != Math.floor(st)) errmsg = "Must enter an integer";
			else{
				if(Math.floor(st) < 0) errmsg = "Must be non-negative";
				else cla[cl].comp[i].simpopinit = Math.floor(st);
			}
		}
		break;
		
	case COMPFRACINITBUT:
		st = ById("inp").value;
		if(st == "" || isNaN(st)) errmsg = "Must enter a number";
		else{
			st = parseFloat(st);
			if(Math.floor(st) < 0) errmsg = "Must be non-negative";
			else{
				sum = 0; 
				for(j = 0; j < cla[cl].ncomp-1; j++){
					if(j == i) sum += st/100;
					else sum += cla[cl].comp[i].simfracinit;
				}
				if(sum < 0 || sum > 1)  errmsg = "Out of range";
				else{
					cla[cl].comp[i].simfracinit = st/100;
					cla[cl].comp[cla[cl].ncomp-1].simfracinit = 1-sum; 
				}
			}
		}
		break;
		
	case COMPSMALLBUT2:
		st = ById("inp").value;
		if(st == "") errmsg = "Must enter an expression";
		else{
			datatemp.posexpression[cl][i] = st;
		}
		break;
	
	case TESTNAMEBUT:
		st = ById("inp").value;
		if(st == "") errmsg = "Must enter an expression";
		else{
			datatemp.testname = st;
		}
		break;
		
	case GENUSERBUT:
		st = ById("inp").value; 
		val = st.split(",");
		for(j = 0; j < val.length; j++){
			if(isNaN(val[j])) errmsg="Not all numbers";
			else{
				val[j] = parseFloat(val[j]);
				if(val[j] < simres.tmin || val[j] > simres.tmax) errmsg = "Not is simulation range";
			}				
		}
		if(errmsg == "")  tgenuserdef = st;
		break;
		
	case COMPPOPINITBUT: cla[selectbubval].comp[selectbubval2].simpopinit
		break;
	
	case SELPARAMBUT:
		st = ById("inptext").value; 
		if(st.length == 0) errmsg = "Must be set";
		else{
			nvarlist = 0; varlist=[]; vardeplist=[];   
			checkeqn(st);
			
			if(errmsg == ""){ for(j = 0; j < nvarlist; j++) if(vardeplist[j] != "") errmsg = "Cannot depend on classification";}
			
			if(errmsg == ""){
				switch(selectbubval2){
				case 0: case 30: case 32: case 34: case 35: case 36: case 38: case 40: param[selectbubval].val[0] = st; break;
				case 1: case 31: case 33: case 37: case 39: param[selectbubval].val[1] = st; break;
				}
				collectvariables();
			}
		}
		break;
		
	case MINMAXBUT: case MINMAXBIGBUT:
		st = ById("inp").value; 
		if(st.length == 0) errmsg = "Must be set to a value";
		else{
			if(isNaN(st)) errmsg = "Must be a number"; 
			else{
				st = parseFloat(st); var isint = 0; if(st == Math.floor(st)) isint = 1;
				switch(selectbubval2){
				case 0: case 30: case 32: case 34: case 35: case 36: case 38: case 40: param[selectbubval].val[0] = st; break;
				case 1: case 31: case 33: case 37: case 39: param[selectbubval].val[1] = st; break;
				case 50: if(st <= 0) errmsg="Must be positive"; else paramagesmooth[selectbubval].val = st;	break;
				case 51: if(st <= 0) errmsg="Must be positive"; else paramtimesmooth[selectbubval].val = st; break;

				case 2: paramsim[selectbubval].sim = st; break;
				case 3:
					if(st > tsimmax) errmsg = "Must be smaller than the maximum simulation time";
					else tsimmin = st; 
					break;
				case 7:
					if(st <= tsimmin) errmsg = "Must be greater than the minimum simulation time";
					else tsimmax = st; 
					break;	
				case 4: if(st > datares.tmin) errmsg = "Must be before data"; else tpostmin = st; break;
				case 5: if(st < datares.tmax) errmsg = "Must be after data"; else tpostmax = st; break;
				case 6:
					if(st <= 0) errmsg = "Must be positive";
					else simnumber = st; 
					break;
				
				case 99:
					if(st <= 0) errmsg = "Must be positive";
					else indmaxnumber = st; 
					break;
							
				case 9: if((st < simres.tmin || st > simres.tmax) && page == SIMULATEPAGE) errmsg = "Must be within the simulation time"; else tgenmin = st; break;
				case 10: if((st > simres.tmax || st < simres.tmin) && page == SIMULATEPAGE) errmsg = "Must be within the simulation time"; else tgenmax = st; break;
				case 11: if(st <= 0) errmsg="Must be positive"; else dtgen = st; break;
				case 12: if(st <= 0 || st > 1) errmsg="Must be between 0 and 1"; else Segen = st; break;
				case 13: if(st <= 0 || st > 1) errmsg="Must be between 0 and 1"; else Spgen = st; break;
				case 14: if(st <= 0 || st > 1) errmsg="Must be between 0 and 1"; else detectprob = st; break;
				case 15: if(st <= 0) errmsg="Must be positive"; else errbarscale = st; break;
				case 16: if(st <= 0) errmsg="Must be positive"; else errbarfix = st; break;
				case 17: if(st <= 0 || st > 1) errmsg="Must be between 0 and 1"; else datatemp.pd = st; break;
				case 18:
					if(st < 1) errmsg = "Must be at least 1"; 
					else{
						if(isint == 0) errmsg = "Must be integer";
						else nsampmax = st;
					}
					break;
				case 19:
					if(st < 1) errmsg = "Must be at least 1"; 
					else{
						if(isint == 0) errmsg = "Must be integer";
						else nsampevmax = st;
					}
					break;
				case 20: if(st < 100) errmsg = "Must be at least 100"; else ESSmin = st; break;
				case 21: if(st < 1) errmsg = "Must be greater than 1"; else GRmax = st; break;
				case 22:
					if(st < 500) errmsg = "Must be at least 500"; 
					else{
						if(isint == 0) errmsg = "Must be integer";
						else itermax = st;
					}
					break;
				}
			}
		}
		break;

	case PDBUT:
		st = ById("inptext").value; 
		
		nvarlist = 0; varlist=[]; vardeplist=[];   
		checkeqn(st);
		if(errmsg == ""){
			switch(datatemp.variety){
			case "cap": case "trans": datatemp.pd = st; break;
			case "capid": row[selectbubval2][1] = st; break;
			}
			setpdvar();
		}
		break;
	
	case TABLEHEADBUT:
		st = ById("inp").value.trim();
		if(st == "")  errmsg = "Cannot be blank";
		else{
			colname[selectbubval] = st; 
			setcolumns();
		}
		break;
		
	case TABLEBUT:
		val = ById("inp").value.trim(); 
		if(val == "") errmsg = "Cannot be blank";
		else{
			if(page == INFERENCEPAGE && pagesub[page] == 1){
				if(selectbubval == 1){
					if(isNaN(val)) errmsg = "Must be a number";
				}
			}
		}
			
		if(errmsg == ""){
			row[selectbubval2][selectbubval] = ById("inp").value.trim(); 
			cv.font = tableheadfont; 
			rowwidth[selectbubval2][selectbubval] = Math.floor(cv.measureText(ById("inp").value).width);
			setcolumns();
		}
		break;
		
	case TABLEBUT3:
		if(ById("inp").value == "") errmsg = "Must enter a name";
		else{
			for(d = 0; d < data.length; d++){
				if(d != selectbubval && data[d].name == ById("inp").value) break;
			}
			if(d < data.length) errmsg = "This name is already used.";
			else data[selectbubval].name = ById("inp").value; 
		}
		break;
		
	case SEARCHBUB:
	case SEARCHRESBUB:
		if(ById("inp").value == "") errmsg = "Must enter a search term";
		else searchterm = ById("inp").value; 
		break;
		
	case DELROWSBUB:
		searchterm = ById("inp").value; 
		break;
		
	case XTICKBUT: case YTICKTRBUT:
		if((selectbub == XTICKBUT && xaxisauto == 0) || (selectbub == YTICKTRBUT && yaxisauto == 0)){
			if(isNaN(ById("inp").value) || isNaN(ById("inp2").value)) errmsg = "Not a number!";
			else{
				num1 = parseFloat(ById("inp").value); num2 = parseFloat(ById("inp2").value)+0.0000000001;
		
				if(num1 > num2) errmsg = "Wrong order!";
				else{
					if(selectbub == XTICKBUT){ xaxisfixmin = num1; xaxisfixmax = num2;}
					else{ yaxisfixmin = num1; yaxisfixmax = num2;}
					indplotstclear();
					plotinitfl = 0;
				}
			}
		}
		if(errmsg == "" && page == INFERENCEPAGE) indplotst=[];
		break;
	}

	if(errmsg != "") return;
	selectbubst = selectbub;
	selectbub = -1; ById("add").innerHTML = ""; inpon = 0; errmsg = "";
	
	switch(selectbubst){
	case SEARCHBUB:
	case SEARCHRESBUB:
		dosearch();
		break;
	}
}

function addbubcol(colnum)                               // Adds a pallete button to bubble
{
	addcanbutton("",cursx-2,cursy,26,15,COLBUT,COLBUT,colnum,-1);
	cursx += 31;
}

function addbubbut(text,type,ac)                         // Adds a button to a bubble
{
  var w, x, y;
  w = textwidth(text,BUTFONT)+20;
  x = bubx + bubw-w-8; y = buby+bubh-25;
  addcanbutton(text,x,y,w,20,ac,type,-1,-1);
}

function addbubbut2(text,type,ac,text2,type2,ac2)        // Adds two buttons to a bubble
{
  var w, x, y;
  w = textwidth(text2,BUTFONT)+20; if(type == BACKCANBUT) w += 15;
  x = bubx + bubw-w-8; y = buby+bubh-25;
  addcanbutton(text2,x,y,w,20,ac2,type2,-1,-1);

  w = textwidth(text,BUTFONT)+20; if(type == BACKCANBUT) w += 15;
  x = x-w-8; y = buby+bubh-25;
  addcanbutton(text,x,y,w,20,ac,type,-1,-1);
}

function addbubbut3(text,type,ac,text2,type2,ac2,text3,type3,ac3)  // Adds three buttons to a bubble
{
  var w, x, y;
  
  w = textwidth(text3,BUTFONT)+20; if(type == BACKCANBUT) w += 15;
  x = bubx + bubw-w-8; y = buby+bubh-25;
  addcanbutton(text3,x,y,w,20,ac3,type3,-1,-1);

  w = textwidth(text2,BUTFONT)+20; if(type == BACKCANBUT) w += 15;
  x = x-w-8; y = buby+bubh-25;
  addcanbutton(text2,x,y,w,20,ac2,type2,-1,-1);

  w = textwidth(text,BUTFONT)+20; if(type == BACKCANBUT) w += 15;
  x = x-w-8; y = buby+bubh-25;
  addcanbutton(text,x,y,w,20,ac,type,-1,-1);
}

function inserttext(text)                                 // Adds text to a bubble 
{
	a = ById(foc);
	if(!a) return;

	startPos = a.selectionStart; endPos = a.selectionEnd;
	st = a.value;
	a.value = st.substr(0,startPos)+text+st.substr(endPos,st.length);

	a.selectionStart = startPos+text.length;  a.selectionEnd = a.selectionStart;
	a.focus();
}

function bubblecontent(w)                                 // The content of bubbles
{
	var i, sel;

	switch(selectbub){
	case EVBUT5:
		addbubtext(selectbubtext+" start time:");
		addbubtext(selectbubval);
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;
	
	case EVBUT6:
		addbubtext(selectbubtext+" end time:");
		addbubtext(selectbubval);
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;
		
	case EVBUT7:
		addbubtext("Individual");
		addbubtext("present");
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;
		
	case EVBUT4:
		ev = indsim.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
		if(ev.obsdata.length != 1) alertp("Error code EC11");
		d = ev.obsdata[0]; i = ev.obsdatai[0];
		addbubtext("Initial state:");
		addbubtext(simdata[d].val[i]);
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;
		
	case EVBUT:
		addbubtext("Observation model:");
		var ev;
		if(page == INFERENCEPAGE){
			if(pagesub[page] == 0) ev = inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
			else ev = infres.inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
		}
		else ev = indsim.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
		
		if(ev.col.length == 0){
			addcanbutton("Problem: Zero probability!",cursx,cursy,80,20,-1,BUBBLETEXTBUT2,-1,-1);
			cursy += 30;
		}		
		else{
			for(k = 0; k < cla[selectbubval2].ncomp; k++){
				if(ev.probexp[k] != "0"){
					st = "Pr(D|"+cla[selectbubval2].comp[k].name+") = "+ev.probexp[k];
					addcanbutton(st,cursx+10,cursy,wmax,20,-1,COMPSMALLBUT5,selectbubval2,k);
					cursy += 30;
				}
			}
		}
		
		addbubtext("Data:");
		for(j = 0; j < ev.obsdata.length; j++){
			d = ev.obsdata[j]; i = ev.obsdatai[j];	
			addcanbutton(data[d].name,cursx,cursy+30*j+4,datwid,20,-1,BUBBLETEXTBUT2,-1,-1);

			gdropinfo.push({val:data[d].val[i], x:cornx+cursx+datwid, y:corny+cursy+30*j, dx:opwid, dy:20, style:1, options:data[d].pos, d:d, i:i, j:j, click:"dataname"});
		}
		cursy += 30*ev.obsdata.length;
		
		addbubtext("Time: "+ev.t); cursy += 20;
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Done",DONEBUT,DONEBUT);	
		break;
	
	case EVBUT2: case EVSOURCEBUT: case EVSINKBUT:
		var ev;
		if(page == INFERENCEPAGE){
			if(pagesub[page] == 0) ev = inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
			else ev = infres.inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
		}
		else ev = indsim.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
		
		addbubtext("Data: "+data[ev.obsdata].name); 
		addbubtext("Time: "+tpre(ev.t,4));
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;	
		
	case EVBUT3:
		addbubtext("Observation model:");
		var ev;
		if(page == INFERENCEPAGE){
			if(pagesub[page] == 0) ev = inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
			else ev = infres.inddata.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
		}
		else ev = indsim.ind[selectbubval].cl[selectbubval2].ev[selectbubtext];
		
		k = ev.col[0];
		st = cla[selectbubval2].comp[k].name;
		w = textwidth(st,"16px georgia");
		addcanbutton(st,cursx+10,cursy,160,20,-1,COMPSMALLBUT5,selectbubval2,k);
		cursy += 30;
			
		addbubtext("Data:");
		d = ev.obsdata; i = ev.obsdatai;	
		addcanbutton(data[d].name,cursx,cursy+4,80,20,-1,BUBBLETEXTBUT2,-1,-1);
			
		gdropinfo.push({val:data[d].val[i], x:cornx+cursx+80, y:corny+cursy, dx:145, dy:20, style:1, options:data[d].pos, d:d, i:i, click:"dataname2"});
	
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Done",DONEBUT,DONEBUT);	
		break;
		
	case COMPBUT:
		addbubtext("Compartment name:");
		addinput(cla[selectbubval].comp[selectbubval2].name,w-20);

		cursy += 10;
		adderrmsg();
	
		addbubtext("Colour:");
		cursx = bubx+bubmar+2; 
		for(i = 0; i < 7; i++) addbubcol(i); 
		cursx = bubx+bubmar+2; cursy += 22;
		for(i = 7; i < 14; i++) addbubcol(i); 
		cursx = bubx+bubmar+2; cursy += 22;
		for(i = 14; i < 21; i++) addbubcol(i); 
		cursx = bubx+bubmar+2; cursy += 22;
		for(i = 21; i < 28; i++) addbubcol(i); 
	
		cursx = bubx+bubmar; cursy += 20;
		if(selectbubval >= ncla-2) addbubbut("Done",DONEBUT,DONEBUT);
		else addbubbut2("Delete",DELBUT,DELCOMPAC,"Done",DONEBUT,DONEBUT);
		break;
	
	case COMPSMALLBUT:
		addbubtext("Value:");
		addinput(cla[selectbubval].comp[selectbubval2].simpopinit,w-20);
		cursy += 10;
		adderrmsg();
	
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Done",DONEBUT,DONEBUT);	
		break;
		
	case COMPSMALLPERBUT:
		addbubtext("Percentage:");
		addinput(cla[selectbubval].comp[selectbubval2].simfracinit*100,w-20);
		cursy += 10;
		adderrmsg();
	
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Done",DONEBUT,DONEBUT);	
		break;
	
	case COMPSMALLBUT2:
		addbubtext("Observation probability:");
		addinput(datatemp.posexpression[selectbubval][selectbubval2],w-20);
		cursy += 10;
		adderrmsg();
	
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;
	
	case TESTNAMEBUT:
		addbubtext("Test name:");
		addinput(datatemp.testname,w-20);
		cursy += 10;
		adderrmsg();
		
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;

	case EQBUT:
		d = selectbubval;
		addbubtext("Expression:");
		addinput(derive[d].eq,w-20); cursy += 10;
		
		cursy += 3;
		adderrmsg();
		cursy += 3;
		
		for(ii = 0; ii < nfuncti; ii++){ addcanbutton("",cursx,cursy,37,26,FUNCTBUT,FUNCTBUT,ii,-1); cursx += 43;}
		cursx = bubx+bubmar+173; for(ii = 0; ii < nopbut; ii++){ addcanbutton("",cursx,cursy,19,26,OPBUT,OPBUT,ii,-1); cursx += 25;}
		
		cursy += 33;
	
		cursx = bubx+bubmar;
		for(i = 0; i < nnumbut; i++){
			addcanbutton("",cursx,cursy,22,26,NUMBUT,NUMBUT,i,-1);
			cursx += 31;
		}
		cursy += 33;
	
		cursx = bubx+bubmar+3;
		addcanbutton("Parameter",cursx,cursy,90,20,ADDPARAMBUT,ADDPARAMBUT,-1,-1);
		addcanbutton("Population",cursx+95,cursy,90,20,ADDPOPAC,ADDPARAMBUT,-1,-1);
	
		addbubbut("Done",DONEBUT,DONEBUT);
		break;
	
	case PARACANBUT2:
		addcanbutton("Description:",cursx,cursy,0,0,-1,BUBBLETEXTBUT3,-1,-1); cursy += 23;
		switch(page){
		case MODELPAGE: addtextarea( cla[pagesub[MODELPAGE]].desc,w-20,203,-1); break;
		case INFERENCEPAGE: addtextarea(datanote,w-20,203,-1); break;
		}
		cursy += 10;
		adderrmsg();
		
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Done",DONEBUT,DONEBUT);
		break;
		
	case TRANSPLBUT:
	case TRANSPLBIGBUB: 
		cl = selectbubval;
		i = selectbubval2;
		tr = cla[cl].tra[i];
		
		if(cl >= ncla-2){
			if(cl == ncla-2){
				addcanbutton("Delete age transition at "+age[i]+"?",cursx,cursy,0,0,-1,BUBBLETEXTBUT3,-1,-1); cursy += 23;
				addbubbut2("Cancel",DONEBUT,CANCELBUT,"Delete",DELBUT,DELTRANSAC2);
			}
			else{
				addcanbutton("Delete time transition at "+time[i]+"?",cursx,cursy,0,0,-1,BUBBLETEXTBUT3,-1,-1); cursy += 23;
				addbubbut2("Cancel",DONEBUT,CANCELBUT,"Delete",DELBUT,DELTRANSAC2);
			}
		}
		else{
			ty = tr.type; if(ty == "Exponential") ty = "Markovian";
			
			addcanbutton("Type: "+ty,cursx,cursy,0,0,-1,BUBBLETEXTBUT3,-1,-1); cursy += 23;
			
			if(selectbub == TRANSPLBIGBUB){
				hh = bubh-135;
				if(errmsg != ""){ if(textwidth(errmsg,"17px arial")+28 > bubw) hh -= 35; else hh -= 15;}
			}
			else hh = 104;
				
			switch(tr.type){
			case "Exponential":	case "Source": case "Sink":	
				addcanbutton("Rate",cursx,cursy,70,20,CANRADIOBUT,CANRADIOBUT,0,CANRADIORATE);
				addcanbutton("Time",cursx+70,cursy,70,20,CANRADIOBUT,CANRADIOBUT,1,CANRADIORATE);
				cursy += 22;
				
				addtextarea(tr.rate,w-20,hh);
				cursy += 7;
				break;
			
			case "Gamma":
				addbubtext("Mean");
				if(selectbub == TRANSPLBIGBUB) addtextarea(tr.mean,w-20,hh-45,tr.shape,w-20);
				else addtextarea(tr.mean,w-20,104,tr.shape,w-20);
				cursy += 7;
			
				addbubtext("Shape parameter");
				cursy += 36;
				break;
				
			case "Weibull":
				addbubtext("λ parameter");
				if(selectbub == TRANSPLBIGBUB) addtextarea(tr.mean,w-20,hh-45,tr.shape,w-20);
				else addtextarea(tr.lam,w-20,104,tr.k,w-20);
				cursy += 7;
			
				addbubtext("k parameter");
				cursy += 36;
				break;
				
			case "Fixed":
				addbubtext("Time");
				addinput(tr.time,w-20);
				cursy += 6;
				break;
			}
			
			cursy += 3;
			adderrmsg();
			
			for(ii = 0; ii < nfuncti; ii++){ addcanbutton("",cursx,cursy,37,26,FUNCTBUT,FUNCTBUT,ii,-1); cursx += 43;}
			cursx = bubx+bubmar+173; for(ii = 0; ii < nopbut; ii++){ addcanbutton("",cursx,cursy,19,26,OPBUT,OPBUT,ii,-1); cursx += 25;}
		
			if(selectbub != TRANSPLBIGBUB){	cursy += 33; cursx = bubx+bubmar;}
			for(ii = 0; ii < nnumbut; ii++){ addcanbutton("",cursx,cursy,22,26,NUMBUT,NUMBUT,ii,-1); cursx += 32;}
			cursy += 33;
			
			cursx = bubx+bubmar+3;
			addcanbutton("Parameter",cursx,cursy,90,20,ADDPARAMBUT,ADDPARAMBUT,-1,-1);
			switch(cla[cl].tra[i].type){
			case "Exponential": case "Source": case "Sink":
				addcanbutton("Population",cursx+95,cursy,90,20,ADDPOPAC,ADDPARAMBUT,-1,-1);
				break;
			}
			
			addbubbut2("Delete",DELBUT,DELTRANSAC,"Done",DONEBUT,DONEBUT);
		}
		break;
	
	case PDBUT:
		addcanbutton("Detection Probability:",cursx,cursy,0,0,-1,BUBBLETEXTBUT3,-1,-1); cursy += 23;
		switch(datatemp.variety){
		case "cap": case "trans": addtextarea(datatemp.pd,w-20,103,-1); break;
		case "cappd": addtextarea(row[selectbubval2][1],w-20,103,-1); break;
		}
		
		cursy += 7;
		cursy += 3;
		adderrmsg();
		
		for(ii = 0; ii < nfuncti; ii++){ addcanbutton("",cursx,cursy,37,26,FUNCTBUT,FUNCTBUT,ii,-1); cursx += 43;}
		cursx = bubx+bubmar+173; for(ii = 0; ii < nopbut; ii++){ addcanbutton("",cursx,cursy,19,26,OPBUT,OPBUT,ii,-1); cursx += 25;}
		
		cursy += 33; cursx = bubx+bubmar;
		for(ii = 0; ii < nnumbut; ii++){ addcanbutton("",cursx,cursy,22,26,NUMBUT,NUMBUT,ii,-1); cursx += 31;}
		cursy += 33;
		
		cursx = bubx+bubmar+3;
		addcanbutton("Parameter",cursx,cursy,90,20,ADDPARAMBUT,ADDPARAMBUT,-1,-1);
		
		
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Done",DONEBUT,DONEBUT);
		break;
	
	case GENUSERBUT:
		addbubtext("User defined times:");
		addinput(tgenuserdef,w-20);
		cursy += 10;
		adderrmsg();
		cursy += 5;
	
		cursx = bubx+bubmar+3;
		for(i = 0; i < nnumbut; i++){
			addcanbutton("",cursx,cursy,22,26,NUMBUT,NUMBUT,i,-1);
			cursx += 31;
		}
		
		addbubbut("Done",DONEBUT,DONEBUT);
		break;
	
	case SPEECHBUT2:
		addcanbutton("Description:",cursx,cursy,0,0,-1,BUBBLETEXTBUT3,-1,-1); cursy += 23;
		addtextarea(descnote,w-20,203,-1);
		cursy += 10;
		adderrmsg();
		
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Done",DONEBUT,DONEBUT);
		break;
		
		
	case MINMAXBUT: case MINMAXBIGBUT:
		switch(selectbubval2){
		case 0:
			addbubtext("Minimum:");
			addinput(param[selectbubval].val[0],w-20);
			break;
			
		case 1:
			addbubtext("Maximum:");
			addinput(param[selectbubval].val[1],w-20);
			break;
	
		case 30:
			addbubtext("Mean:");
			addinput(param[selectbubval].val[0],w-20);
			break;
			
		case 31:
			addbubtext("Standard deviation:");
			addinput(param[selectbubval].val[1],w-20);
			break;
		
		case 32:
			addbubtext("Mean (logscale):");
			addinput(param[selectbubval].val[0],w-20);
			break;
			
		case 33:
			addbubtext("Standard deviation (logscale):");
			addinput(param[selectbubval].val[1],w-20);
			break;	

		case 34:
			addbubtext("Rate:");
			addinput(param[selectbubval].val[0],w-20);
			break;
			
		case 35:
			addbubtext("Value:");
			addinput(param[selectbubval].val[0],w-20);
			break;
			
		case 36: case 40:
			addbubtext("α:");
			addinput(param[selectbubval].val[0],w-20);
			break;
			
		case 37:
			addbubtext("β:");
			addinput(param[selectbubval].val[0],w-20);
			break;
			
		case 38:
			addbubtext("λ:");
			addinput(param[selectbubval].val[0],w-20);
			break;
			
		case 39:
			addbubtext("k:");
			addinput(param[selectbubval].val[0],w-20);
			break;
		
		case 50: 
			addbubtext("Value:");
			addinput(paramagesmooth[selectbubval].val,w-20);
			break;
		
		case 51: 
			addbubtext("Value:");
			addinput(paramtimesmooth[selectbubval].val,w-20);
			break;
			
		case 2: 
			addbubtext("Value:");
			addinput(paramsim[selectbubval].sim,w-20);
			break;
		case 3: 
			addbubtext("Value:");
			addinput(tsimmin,w-20);
			break;
		case 7: 
			addbubtext("Value:");
			addinput(tsimmax,w-20);
			break;
		case 4:
			addbubtext("Minimum time:");
			addinput(tpostmin,w-20);
			break;
		case 5:
			addbubtext("Maximum time:");
			addinput(tpostmax,w-20);
			break;
		case 6: 
			addbubtext("Value:");
			addinput(simnumber,w-20);
			break;
			
		case 99: 
			addbubtext("Value:");
			addinput(indmaxnumber,w-20);
			break;
			
		case 8:
			addbubtext(":");
			addinput(simnumber,w-20);
			break;
		
		case 9:
			addbubtext("Minimum time:");
			addinput(tgenmin,w-20);
			break;
			
		case 10:
			addbubtext("Maximum time:");
			addinput(tgenmax,w-20);
			break;
			
		case 11:
			addbubtext("Step size:");
			addinput(dtgen,w-20);
			break;
			
		case 12:
			addbubtext("Sensitivity:");
			addinput(Segen,w-20);
			break;
			
		case 13:
			addbubtext("Specificity:");
			addinput(Spgen,w-20);
			break;
			
		case 14:
			addbubtext("Detection probability:");
			addinput(detectprob,w-20);
			break;
			
		case 15:
			addbubtext("Error bar scale:");
			addinput(errbarscale,w-20);
			break;	
			
		case 16:
			addbubtext("Error bar SD:");
			addinput(errbarfix,w-20);
			break;	
		
		case 17:
			addbubtext("Detection probability:");
			addinput(datatemp.pd,w-20);
			break;
			
		case 18:
			addbubtext("Maximum:");
			addinput(nsampmax,w-20);
			break;
			
		case 19:
			addbubtext("Maximum:");
			addinput(nsampevmax,w-20);
			break;
			
		case 20:
			addbubtext("Minimum:");
			addinput(ESSmin,w-20);
			break;
			
		case 21:
			addbubtext("Maximum:");
			addinput(GRmax,w-20);
			break;
			
		case 22:
			addbubtext("Iteration number:");
			addinput(itermax,w-20);
			break;
		}
		
		cursy += 15;
		adderrmsg();
	
		cursx = bubx+bubmar+3;
		
		
		for(i = 0; i <6; i++){
			addcanbutton("",cursx,cursy,22,26,NUMBUT,NUMBUT,i,-1);
			cursx += 31;
		}
		cursx = bubx+bubmar+3+15; cursy += 30;
		for(i = 6; i < nnumbut; i++){
			addcanbutton("",cursx,cursy,22,26,NUMBUT,NUMBUT,i,-1);
			cursx += 31;
		}
		
		addbubbut("Done",DONEBUT,DONEBUT);
		break;

	case CONVERTDATEBUB:
		addbubtext("This looks like a date which");
		addbubtext("needs to be converted to decimal.");
		addbubtext("Would you like to convert?");
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Convert",DONEBUT,CONVERTAC);
		break;
		
	case MULTISTATEBUB:
		addbubtext("This individual has multiple");
		addbubtext("states associated with it.");
		break;
		
	case TABLEBUT:
		addbubtext("Value:");
		addinput(row[selectbubval2][selectbubval],w-20); cursy += 10;
		adderrmsg();
		addbubbut3("Ins. row",GREENCANBUT,INSROWAC,"Del. row",CANCELBUT,DELROWAC,"Done",DONEBUT,DONEBUT);	
		break;
		
	case TABLEDROPBUT:
		addbubtext("Value:");
		
		gdropinfo.push({val:row[selectbubval2][selectbubval], x:cornx+cursx+58, y:corny+cursy-25, dx:105, dy:20, style:2, options:datatemp.pos, click:"tabledrop", j:selectbubval2, i:selectbubval});
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;
	
	case TABLEBUT3:
		addbubtext("Name:");
		addinput(data[selectbubval].name,w-20); cursy += 13;
		adderrmsg();
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;	
		
	case TABLEHEADBUT:
		addbubtext("Name:");
		addinput(colname[selectbubval],w-20); cursy += 13;
		adderrmsg();
	
		addcanbutton("Search",cursx,cursy,100,20,SEARCHAC,DELBUT,-1,-1);
		addcanbutton("Replace",cursx+110,cursy,100,20,REPLACEAC,DELBUT,-1,-1);
		cursy += 28;
		addcanbutton("Sort (A-Z)",cursx,cursy,100,20,SORTAC,DELBUT,0,-1);
		addcanbutton("Sort (0-9)",cursx+110,cursy,100,20,SORTAC,DELBUT,1,-1);
		cursy += 28;	
		addcanbutton("Delete rows",cursx,cursy,100,20,DELROWSAC,DELBUT,-1,-1);
		if(addingdata != 2) addcanbutton("Join cols",cursx+110,cursy,100,20,JOINCOLSAC,DELBUT,-1,-1);
		break;
		
	case SEARCHBUB:
		addbubtext("Search for:");
		addinput(searchterm,w-20); cursy += 13;
		adderrmsg();
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Search",DONEBUT,DOSEARCHAC);	
		break;
		
	case SEARCHRESBUB:
		addbubtext("Search for:");
		addinput(searchterm,w-20); cursy += 13;
		adderrmsg();
		addbubtext((searchresnum+1)+" out of "+searchres.length);
		addbubbut2("Back",BACKCANBUT,BACKSEARCHAC,"Next",NEXTCANBUT,NEXTSEARCHAC);	
		break;
		
	case REPLACEBUB:
		addbubtext("Replace:");
		
		st = "<input id='inp' autocomplete='off' style='text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+9)+"px; left:"
		     +(canx+cursx+12)+"px; width:"+(w-20)+"px; height:25px; background=#ffffff;' value='"+searchterm
			 +"' onkeypress='if(event.keyCode == 13) ById(\"inp2\").focus();' oninput='searchterm=ById(\"inp\").value;' type='text'/>";
		cursy += 37;
		
		addbubtext("With:");
	 
		st += "<input id='inp2' autocomplete='off' style='text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+9)+"px; left:"
		      +(canx+cursx+12)+"px; width:"+(w-20)+"px; height:25px; background=#ffffff;' value='"+replaceterm
			  +"' onkeypress='if(event.keyCode == 13) doreplace();' oninput='replaceterm=ById(\"inp2\").value;' type='text'/>";
		cursy += 37;
		adderrmsg();

		ById("add").innerHTML = st;
		a = ById("inp"); a.focus(); a.selectionStart = 0;  a.selectionEnd = text.length;
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Replace",DONEBUT,DOREPLACEAC);	
		break;
		
	case REPLACEDONEBUB:
		addbubtext(nreplace +" entries replaced");
		addbubbut2("Undo",CANCELBUT,UNDOAC,"Done",DONEBUT,CANCELBUT);	
		break;
	
	case DELROWSBUB:
		addbubtext("Delete all rows with:");
		addinput(searchterm,w-20); cursy += 13;
		adderrmsg();
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Delete",DONEBUT,DODELROWSAC);	
		break;
		
	case DELROWSDONEBUB:
		addbubtext(ndel +" entries deleted");
		addbubbut2("Undo",CANCELBUT,UNDOAC,"Done",DONEBUT,CANCELBUT);	
		break;
		
	case ERRBUB:
		addbubtext(errmsg);
		addbubbut("OK",DONEBUT,DONEBUT);
		break;
		
	case ADDCLASSBUB:
		addbubtext("First a classification"); y += 30;
		addbubtext("must be added");
		break;
		
	case EDITCLASSBUB:
	    addbubtext("Name:");
		addinput(cla[selectbubval].name,w-20);
		cursy += 10;
		adderrmsg();
	
		addbubbut2("Delete",CANCELBUT,DELCLASSAC,"Done",DONEBUT,DONEBUT);	
		break;
	
	case EDITDBUT:
	    addbubtext("Data D:");
		addinput(datagen[clgl].posedit[selectbubval],w-20);
		cursy += 10;
		adderrmsg();
	
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;
		
	case EDITAGETIMEBUB:
	    addbubtext("Delete?");
		addbubbut2("Delete",CANCELBUT,DELCLASSAC,"Done",DONEBUT,DONEBUT);	
		break;
		
	case XTICKBUT:
	case YTICKTRBUT:
		if(selectbub == XTICKBUT){ mi = xaxisfixmin; ma = xaxisfixmax; rad = CANRADIOXAXIS;}
		else{ mi = yaxisfixmin; ma = yaxisfixmax; rad = CANRADIOYAXIS;}
		 
		addcanbutton("Auto",cursx,cursy,70,20,CANRADIOBUT,CANRADIOBUT,1,rad);
		addcanbutton("Define",cursx+70,cursy,70,20,CANRADIOBUT,CANRADIOBUT,0,rad);
		
		if((xaxisauto == 0 && selectbub == XTICKBUT) || (yaxisauto == 0 && selectbub == YTICKTRBUT)){
			mi = rn(parseFloat(mi)); ma = rn(parseFloat(ma));	
			
			cursy += 30;
			addbubtext("Minimum:");
		
			if(inpon == 0){
				st = "<input id='inp' autocomplete='off' style='text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+9)+"px; left:"+(canx+cursx+12)+"px; width:"+(bubw-20)+"px; height:25px; background=#ffffff;' value='"+mi+"' onkeypress='if(event.keyCode == 13) ById(\"inp2\").focus();' type='text'/>";
				cursy += 35;
		
				addbubtext("Maximum:");
				st += "<input id='inp2' autocomplete='off' style='text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+9)+"px; left:"+(canx+cursx+12)+"px; width:"+(bubw-20)+"px; height:25px; background=#ffffff;' value='"+ma+"' onkeypress='if(event.keyCode == 13){ buboff(0); buttoninit();}' type='text'/>";
		
				ById("add").innerHTML = st;
				inpon = 1;
			}
			else{
				a = ById("inp"); a.style = "text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+9)+"px; left:"+(canx+cursx+12)+"px; width:"+(bubw-20)+"px;";
				cursy += 35;
				addbubtext("Maximum:");
				a = ById("inp2"); a.style = "text-align:left; position:absolute;  font-size:22px; top:"+(cany+cursy+9)+"px; left:"+(canx+cursx+12)+"px; width:"+(bubw-20)+"px;";
			}
			cursy += 40;
			adderrmsg();
		}
		
		addbubbut("Done",DONEBUT,DONEBUT);
		break;
		
	case EDITBUB:
		addbubtext("Click here to");
		addbubtext("edit the model!");
		break;
	
	case ADDBUT:
		addbubtext("Type:");
		
		te = "Markovian";
		var dx = textwidth(te,addfont)+19;
		addcanbutton(te,cursx,cursy,dx,20,ADDLINKAC,ADDBUT3,-1,-1);
		addcanbutton("[?]",cursx+dx,cursy,15,20,HELPICONCANBUT,HELPICONCANBUT,39,0);
		cursy += 20;
			
		te = "Gamma";
		var dx = textwidth(te,addfont)+19;
		addcanbutton(te,cursx,cursy,dx,20,ADDLINKAC,ADDBUT3,-1,-1);
		addcanbutton("[?]",cursx+dx,cursy,15,20,HELPICONCANBUT,HELPICONCANBUT,40,0);
		cursy += 20;
		
		te = "Weibull";
		var dx = textwidth(te,addfont)+19;
		addcanbutton(te,cursx,cursy,dx,20,ADDLINKAC,ADDBUT3,-1,-1);
		addcanbutton("[?]",cursx+dx,cursy,15,20,HELPICONCANBUT,HELPICONCANBUT,41,0);
		cursy += 20;
		cursy -= 60; cursx +=130;
		
		te = "Source";
		var dx = textwidth(te,addfont)+19;
		addcanbutton(te,cursx,cursy,dx,20,ADDLINKAC,ADDBUT3,-1,-1);
		addcanbutton("[?]",cursx+dx,cursy,15,20,HELPICONCANBUT,HELPICONCANBUT,42,0);
		cursy += 20;
		
		te = "Sink";
		var dx = textwidth(te,addfont)+19;
		addcanbutton(te,cursx,cursy,dx,20,ADDLINKAC,ADDBUT3,-1,-1);
		addcanbutton("[?]",cursx+dx,cursy,15,20,HELPICONCANBUT,HELPICONCANBUT,43,0);
		cursy += 20;
		break;
		
	case COMPPOPOUTPUTBUT:
		addbubtext("Mean:"+tpre(quantmean[playtime][selectbubval2],5));
		addbubtext("95% CI (min):"+tpre(quantCImin[playtime][selectbubval2],5));
		addbubtext("95% CI (max):"+tpre(quantCImax[playtime][selectbubval2],5));
		addbubtext("Time:"+tpre((res.tmin+(playtime*(res.tmax-res.tmin)/(POPX-1))),5));
		break;
		
	case COMPFRACINITBUT:
		addbubtext("Initial fraction:");
		addinput(cla[selectbubval].comp[selectbubval2].simfracinit*100,w-20);
		cursy += 10;
		adderrmsg();
	
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;
		
	case COMPPOPINITBUT:
		addbubtext("Initial population:");
		addinput(cla[selectbubval].comp[selectbubval2].simpopinit,w-20);
		cursy += 10;
		adderrmsg();
	
		addbubbut("Done",DONEBUT,DONEBUT);	
		break;
		
	case COMPLOADINITBUT: 
		addbubtext("Initial population:");
		addinput(cla[selectbubval].comp[selectbubval2].loadinit.length,w-20);
		break;
		
	case JOINBUB: 
		addbubtext("Select column to ");
		addbubtext("join to this one");
		addbubbut("Cancel",CANCELBUT,CANCELBUT);	
		break;
		
	case EMPTYBUB:
		addbubtext("This element ");
		addbubtext("cannot be empty!");
		break;
		
	case NANBUB:
		addbubtext("This must ");
		addbubtext("be a number!");
		break;
		
	case PROBBUB:
		addbubtext("There is a ");
		addbubtext("problem here!");
		break;
		
	case SELPARAMBUT:
		var quant = "value";
		switch(selectbubval2){
		case 0: quant = "Minimum"; break;
		case 1: quant = "Maximum"; break;
		case 30: quant = "Mean"; break;
		case 31: quant = "Standard deviation"; break;
		case 32: quant = "Mean (logscale)"; break;
		case 33: quant = "Standard deviation (logscale)"; break;
		case 34: quant = "Rate"; break;
		case 35: quant = "Value"; break;
		case 36: case 38: quant = "α:"; break;
		case 37: quant = "β"; break;
		case 38: quant = "λ"; break;
		case 39: quant = "k"; break;
		}
		
		addcanbutton(quant+":",cursx,cursy,0,0,-1,BUBBLETEXTBUT3,-1,-1); cursy += 23;
		
		switch(selectbubval2){
		case 0: case 30: case 32: case 34: case 35: case 36: case 38: case 40: namenow = param[selectbubval].val[0]; break;
		case 1: case 31: case 33: case 37: case 39: namenow = param[selectbubval].val[1]; break;
		}
		addtextarea(namenow,w-20,103,-1);
	
		cursy += 7;
		cursy += 3;
		adderrmsg();
		
		for(ii = 0; ii < nfuncti; ii++){ addcanbutton("",cursx,cursy,37,26,FUNCTBUT,FUNCTBUT,ii,-1); cursx += 43;}
		cursx = bubx+bubmar+173; for(ii = 0; ii < nopbut; ii++){ addcanbutton("",cursx,cursy,19,26,OPBUT,OPBUT,ii,-1); cursx += 25;}
		
		cursy += 33; cursx = bubx+bubmar;
		for(ii = 0; ii < nnumbut; ii++){ addcanbutton("",cursx,cursy,22,26,NUMBUT,NUMBUT,ii,-1); cursx += 31;}
		cursy += 33;
		
		cursx = bubx+bubmar+3;
		addcanbutton("Parameter",cursx,cursy,90,20,ADDPARAMBUT,ADDPARAMBUT,-1,-1);
		
		addbubbut2("Cancel",CANCELBUT,CANCELBUT,"Done",DONEBUT,DONEBUT);
		break;
	}
}

function adderrmsg(st)                                  // Adds an error message to a bubble
{
	var j;
	if(errmsg != ""){
		for(j = 0; j < errmsgline.length; j++){
			addcanbutton(errmsgline[j],cursx,cursy-5,bubw,20,-1,ERRMSGBUT,j,-1);
			cursy += 20;
		}
	}
}
		