function modelbuts()                                      // Draws the model page
{
	cl = pagesub[MODELPAGE];

	switch(cl){
	case -3:		// Depend
		x = menux+tab; y = 30;
		xx = width-115; yy = height-45;		
		if(adddepfl == 1){
			addbutton("Choose parameter to place distribution on:",x,y,0,0,-1,TITLEBUT,-1,-1); y += 20;
			priorbuts();
			addbutton("Cancel",xx,yy,100,30,CANCELDEPAC,CANCELBUT2,-1,-1);
		}
		else{
			addbutton("Impose distributions on parameters",x,y,0,0,-1,TITLEBUT,29,-1); y += 50;
			priorbuts();
			
			if(modelsetup == 0){
				if(ncla > 2) addbutton("Done!",xx,yy,100,30,CLEARBUT2,CLEARBUT2,-1,-1);

				te = "Distribution";
				var dx = textwidth(te,addfont)+19;
				addbutton(te,menux+tab+20,height-40,dx,20,ADDDEPBUT,ADDDEPBUT,-1,-1);
				addbutton("[?]",menux+tab+20+dx,height-40,15,20,HELPICONBUT,HELPICONBUT,96,-1);
			}
			else addbutton("Edit",xx,yy,100,30,CLEARBUT2,CLEARBUT2,-1,-1);
		}
		return;
	
	case -2:       // Derived
		x = menux+tab; y = 30;
		addbutton("Quantities derived from the model",x,y,0,0,-1,TITLEBUT,32,-1); y += 50;
	
		if(derive.length == 0) addbutton("There are currently no derived quantities.",x+18,y,900,0,-1,PARAGRAPHBUT,1,-1);
		else{
			y = 80;
			addbutton("Derived",x+20,y,100,0,-1,PRHEADBUT,-1,-1);
			addbutton("Expression",x+180,y,100,0,-1,PRHEADBUT,-1,-1);

			y = 110;
			cornx = x+20, corny = y;
			addbutton("",cornx,corny,derwidth,derheight,CANVASBUT,CANVASBUT,-1,-1);
			
			drawderived();
			
			tableyfrac = derheight/ytot;
			if(tableyfrac < 1) addbutton("",cornx+derwidth+10,corny,13,derheight-22,SLIDEAC,YSLIDEBUT,-1,-1);
		}
			
		xx = width-115; yy = height-45;		
		if(modelsetup == 0){
			if(ncla > 2) addbutton("Done!",xx,yy,100,30,CLEARBUT2,CLEARBUT2,-1,-1);
			
			textwidth(te,addfont)+19
			te = "Derived";
			var dx = textwidth(te,addfont)+19;
			addbutton(te,menux+tab+20,height-40,dx,20,ADDDERIVEAC,ADDBUT2,-1,-1);
			addbutton("[?]",menux+tab+20+dx,height-40,15,20,HELPICONBUT,HELPICONBUT,34,-1);
		}
		else addbutton("Edit",xx,yy,100,30,CLEARBUT2,CLEARBUT2,-1,-1);
	
		return;
	}
	
	x = menux+tab; y = 30;
	addbutton("",menux,0,width-menux,height,CANVASBUT,CANVASBUT,-1,-1);
	
	if(cl == -1) return;
	
	drawmodel(pagesub[MODELPAGE],0,0,width-menux,height,"noframe","normal");

	xx = width-315-90; addcanbutton("",xx,3,25,22,SPEECHBUT,SPEECHBUT,-1,-1);
	if(cla[pagesub[MODELPAGE]].descon == 1){
		dx = width-345-90;
		alignparagraph(cla[pagesub[MODELPAGE]].desc,dx);
		addcanbutton("",xx-2,3+25,28,18,EDITSPEECHBUT,EDITSPEECHBUT,ncanbut+1,-1);
		addcanbutton(cla[pagesub[MODELPAGE]].desc,20,10,dx,nlines*20+20,PARACANBUT2,PARACANBUT2,-1,-1);
	}
	
	if(ncla > 2){
		xx = width-menux-115; yy = height-45;	
		if(modelsetup == 0) addcanbutton("Done!",xx,yy,100,30,DONEMODELAC,CLEARBUT,-1,-1);
		else addcanbutton("Edit",xx,yy,100,30,DONEMODELAC,CLEARBUT,-1,-1);
		
		xx = width-menux-150; yy = height-25;	
		dx = 22; dy = 26; dd = 3;
		addcanbutton("",xx+dd,yy-dy/2,dx,dy,ZOOMINBUT,ZOOMINBUT,UP,-1);
		addcanbutton("",xx-dd-dx,yy-dy/2,dx,dy,ZOOMOUTBUT,ZOOMOUTBUT,UP,-1); 
	
		x = 0; y = height-40;
	
		if(modelsetup == 1) return;
	
		if(cla[cl].name == "Time"){ 
			te = "Time transition";
			var dx = textwidth(te,addfont)+19;
			addcanbutton(te,20,height-40,dx,20,ADDTIMEPOINTAC,ADDBUT,-1,-1);
			addcanbutton("[?]",20+dx,height-40,15,20,HELPICONCANBUT,HELPICONCANBUT,37,-1);
			return;
		}
		if(cla[cl].name == "Age"){
			te = "Age transition";
			var dx = textwidth(te,addfont)+19;
			addcanbutton(te,20,height-40,dx,20,ADDAGEPOINTAC,ADDBUT,-1,-1);
			addcanbutton("[?]",20+dx,height-40,15,20,HELPICONCANBUT,HELPICONCANBUT,38,-1);
			return;
		}

		if(addtransmode == 1){
			if(cla[transtempcl].tra[transtempk].i == -1){
				if(addsourcemode == -1) addrequest("Select compartment");		
				else addrequest("Select initial compartment");	
			}				
			else{
				if(addsourcemode == -1) addrequest("Place sink");	
				else addrequest("Select final compartment (along with any intermediary points)");	
			}
		}		
		else{
			if(addcompmode == 1) addrequest("Place compartment");
			else{
				if(addsourcemode == 1) addrequest("Places source"); 
				else{
					te = "Compartment";
					var dx = textwidth(te,addfont)+19;
					addcanbutton(te,x+20,height-40,dx,20,ADDCOMPAC,ADDBUT,-1,-1);
					addcanbutton("[?]",x+20+dx,y,15,20,HELPICONCANBUT,HELPICONCANBUT,35,-1);
			
					te = "Transition";
					var dx = textwidth(te,addfont)+19;
					addcanbutton(te,x+220,height-40,dx,20,ADDTRANSAC,ADDBUT,-1,-1);
					addcanbutton("[?]",x+220+dx,y,15,20,HELPICONCANBUT,HELPICONCANBUT,36,-1);
				}
			}
		}
		
		if(addcompmode == 1 && canover == -1){ 
			addcanbutton(dragtext,canmx-dragww/2,canmy-dragdh/2,dragww,dragdh,-1,DRAGBUT,-1,-1);
		}
		
		if(addsourcemode == 1 && canover == -1) addcanbutton("+",canmx,canmy,0,0,-1,SOURCEBUT,-1,-1);
	}
}

function modframe(cl,dx,dy)                              // Calculates the bounding box for a model  
{
	var k, c, tr, p, xmin = large, xmax=-large, ymin=large, ymax=-large, wmax = 0;

	for(k = 0; k < cla[cl].ncomp; k++){
		c = cla[cl].comp[k];
		if(c.x < xmin) xmin = c.x; if(c.x+c.w > xmax) xmax = c.x+c.w;
		if(c.y < ymin) ymin = c.y; if(c.y+c.h > ymax) ymax = c.y+c.h;
		if(c.w > wmax) wmax = c.w;
	}

	for(tr = 0; tr < cla[cl].ntra; tr++){
		for(p = 0; p < cla[cl].tra[tr].p.length; p++){
			if(cla[cl].tra[tr].p[p].x < xmin) xmin = cla[cl].tra[tr].p[p].x;
			if(cla[cl].tra[tr].p[p].x > xmax) xmax = cla[cl].tra[tr].p[p].x;
			if(cla[cl].tra[tr].p[p].y < ymin) ymin = cla[cl].tra[tr].p[p].y;
			if(cla[cl].tra[tr].p[p].y > ymax) ymax = cla[cl].tra[tr].p[p].y;
		}
	}
	
	fac = (dx/(xmax-xmin))/(1.03+1.0/cla[cl].ncomp);

	fac2 = (dy/(ymax-ymin))/(1.03+1.0/cla[cl].ncomp); if(fac2 < fac) fac = fac2;
	dyideal = Math.floor(fac*(ymax-ymin)*(1.03+1.0/cla[cl].ncomp));
	
	fac3 = 0.2*dx/wmax; if(fac3 < fac) fac = fac3;
	
	shx = 0.5*(xmax+xmin) - 0.5*dx/fac; shy = 0.5*(ymax+ymin) - 0.5*dy/fac; shfac = fac;	

	return dyideal;
}

function drawmodel(cl,xco,yco,dx,dy,frame,type)            // Draws all the elements of the model
{
	var k;

	modelty = type;
	
	if(frame == "noframe"){ shx = 0; shy = 0; shfac = 1;}
	else modframe(cl,dx,dy);	  
	
	cv = resultcv;
	fillrect(0,0,dx,dy,WHITE); 
	addcanbutton("",xco,yco,dx,dy,-1,MODELPICBUT,-1,-1);

	if(ncla == 0) return;

	for(k = 0; k < cla[cl].ntra; k++){
		tr = cla[cl].tra[k];
		if(tr.i != -1 && tr.f != -1) drawtrans(tr,0);
		if(type == "normal" || type == "seltrans") addcanbutton("links",0,0,0,0,TRANSPLBUT,TRANSPLBUT,cl,k);
	}

	if(type == "normal"){
		for(k = 0; k < cla[cl].ntra; k++){
			tr = cla[cl].tra[k];
			for(j = 0; j < tr.p.length; j++){
				x = (tr.p[j].x-shx)*shfac; y = (tr.p[j].y-shy)*shfac;
				if(tr.i == -2 && j == 0){ st = "+"; r = 15;}
				else{
					if(tr.f == -2 && j == tr.p.length-1){ st = "-"; r = 15;}
					else{ st = ""; r = 8;}
				}	
				addcanbutton(st,x-r,y-r,2*r,2*r,TRANSPBUT,TRANSPBUT,k,j);
			}
		}	
	}

	switch(type){
	case "seltrans": bu = COMPBUT; break;
	
	case "normal": case "indoutput": bu = COMPBUT; break;
	
	case "popinit":
		if(cl == simpopinitset){
			bu = COMPPOPINITBUT;
			popmax = 1; for(k = 0; k < cla[cl].ncomp; k++){ c = cla[cl].comp[k]; if(c.simpopinit > popmax) popmax = c.simpopinit;}
		}
		else{
			bu = COMPFRACINITBUT;
		}
		break;
		
	case "loadinit":
		bu = COMPLOADINITBUT;
		popmax = 1; for(k = 0; k < cla[cl].ncomp; k++){ c = cla[cl].comp[k]; if(c.loadinit.length > popmax) popmax = c.loadinit.length;}
		break;
		
	case "popoutput":
		bu = COMPPOPOUTPUTBUT;
		break;
	}
		
	for(k = 0; k < cla[cl].ncomp; k++){
		c = cla[cl].comp[k];
		x = Math.floor((c.x-shx)*shfac); y = Math.floor((c.y-shy)*shfac);
		w = Math.floor(c.w*shfac); h = Math.floor(c.h*shfac); 
		
		if(type == "indoutput") compplot(cl,k,x,y,w,h,bu,2);
		else{
			if(type != "seltrans" && type != "loadinit") addcanbutton(c.name,x+xco,y+yco,w,h,bu,bu,cl,k);
			if(playing == 0) compplot(cl,k,x,y,w,h,bu,0);
		}
		
		// works out if any sources or sinks
		sofl = 0; sifl = 0;
		for(cl2 = 0; cl2 < ncla; cl2++){
			if(cl2 != cl){ 
				for(t = 0; t < cla[cl2].ntra; t++){
					tr = cla[cl2].tra[t]; 
					if(tr.i == -2) sofl = 1;
					if(tr.f == -2) sifl = 1;
				}
			}
		}		
	}
	
	if(type == "indoutput"){
		cv.globalAlpha = 0.8;
		fillrect(0,0,dx,dy,WHITE);
		getminmax();
		 
		i = res.indfilt; cl = res.clasel;
		for(ch = chmin; ch < chmax; ch++){
			for(s = smin; s < smax; s++){
				w = res.ch[ch].sampev[s].ind[i];
				var p=[];
				for(e = 0; e < w.nev; e++){
					cc = w.evc[e];
					if(cc != NOTALIVE){
						c = cla[cl].comp[compval[cc][cl]];
						p.push({x:xco+(c.x+c.w/2-shx)*shfac, y:yco+(c.y+c.h/2-shy)*shfac, t:w.evt[e]});
					}
				}
			
				dt = p[p.length-1].t - p[0].t;
				for(e = 0; e < p.length-1; e++){
					drawgradientline(p[e].x,p[e].y, p[e+1].x,p[e+1].y,(p[e].t- p[0].t)/dt,(p[e+1].t- p[0].t)/dt);
				}
				
				r = 5;
				drawroundrect(p[0].x-r,p[0].y-r,2*r,2*r,r,BLACK,BLACK);
				centertext("Start",p[0].x,p[0].y-9,"bold 16px Times",BLUE); 
				drawroundrect(p[p.length-1].x-r,p[p.length-1].y-r,2*r,2*r,r,WHITE,BLACK);
				centertext("End",p[p.length-1].x,p[p.length-1].y+19,"bold 16px Times",BLUE); 
			}
		}
	}
}
	
function compplot(cl,k,x,y,dx,dy,bu,ov)                   // Plots a compartment button
{
	var v, text, col, col2, col3, c = cla[cl].comp[k];
	
	switch(bu){
	case COMPBUT:
		text = c.name; col = c.col; col2 = darkcol(col); col3 = WHITE;
		if(ov == 1) col = darkcol(col);
		break;

	case COMPPOPINITBUT:
		if(radnamlab == 0) text =  c.name+":"+c.simpopinit; else text = c.simpopinit; 
		if(radstyle == 0){ col = c.col; col2 = darkcol(col); col3 = WHITE;}
		else{
			v = Math.floor(255*(1-parseInt(c.simpopinit)/popmax));
			col = "rgb("+v+","+v+","+v+")"; col2 = BLACK; if(v > 128) col3 = DBLUE; else col3 = WHITE;
		}
		break;
		
	case COMPFRACINITBUT:
		if(radnamlab == 0) text =  c.name+":"+tpre((100*c.simfracinit),3)+"%"; 
		else text = tpre((100*c.simfracinit),3)+"%";
		
		if(radstyle == 0){ col = c.col; col2 = darkcol(col); col3 = WHITE;}
		else{
			v = Math.floor(255*(1-c.simfracinit));
			col = "rgb("+v+","+v+","+v+")"; col2 = BLACK; if(v > 128) col3 = DBLUE; else col3 = WHITE;
		}
		break;
	
	case COMPLOADINITBUT:
		if(radnamlab == 0)  text =  c.name+":"+c.loadinit.length;
		else text = c.loadinit.length;
		
		if(radstyle == 0){ col = c.col; col2 = darkcol(col); col3 = WHITE;}
		else{
			v = Math.floor(255*(1-c.loadinit.length/popmax));
			col = "rgb("+v+","+v+","+v+")"; col2 = BLACK; if(v > 128) col3 = DBLUE; else col3 = WHITE;
		}
		break;
		
	case COMPPOPOUTPUTBUT:
		switch(radstyle2){
		case 0:
			col = c.col; col2 = darkcol(col); col3 = WHITE;
			break;
		
		case 1:
			v = Math.floor(255*(1-quantmean[playtime][k]/quantmaxitot));
			col = "rgb("+v+","+v+","+v+")"; col2 = BLACK; if(v > 128) col3 = DBLUE; else col3 = WHITE;
			break;
			
		case 2:	
			v = Math.floor(255*(1-quantmean[playtime][k]/quantmaxi[playtime]));
			col = "rgb("+v+","+v+","+v+")"; col2 = BLACK; if(v > 128) col3 = DBLUE; else col3 = WHITE;
			break;
		}
		
		text = "";
		if(radlab != 0){
			if(radlab == 1) text = cla[cl].comp[k].name; else text = tpre(quantmean[playtime][k],3);
		}
		break;
	}
	
	dy = Math.floor(dy);
	drawroundrect(Math.floor(x),Math.floor(y),Math.floor(dx),dy,Math.floor(0.2*dy),col,col2);
	
	if(ov == 1 && bu != COMPBUT){
		cv.lineWidth = 3;
		cv.strokeStyle = RED;
		cv.stroke();
	}
	
	if(bu != COMPPOPOUTPUTBUT || text != ""){
		si = dy*0.7; 
		if(bu == COMPBUT) fo = Math.floor(si)+  "px georgia";
		else{
			do{
				fo = Math.floor(si)+  "px georgia";
				if(textwidth(text,fo) < dx*0.95) break;
				si *= 0.8;
			}while(1 == 1);	
		}
		centertext(text,x+dx/2,y+Math.floor(0.5*dy+0.3*si),fo,col3);
	}
}

function addrequest(st)                                   // Adds a request on the model page
{
	var dx = textwidth(st,"bold 16px arial");
	
	addcanbutton(st,(width-menux)/2-10,height-30,0,0,-1,REQUESTBUT,-1,-1);
	addcanbutton(st,(width-menux)/2-10 + dx/2 + 4,height-27,15,15,CROSSBUT,CROSSBUT,-1,-1);
}

function newname(na)                                      // Generates a new compartment name
{
	var flag, num, st, j, cl;
	
	num = 1; 
	do{
		flag = 0;
		st = na; if(num > 1) st += num;
		for(cl = 0; cl < ncla; cl++){
			if(cla[cl].name == st) flag = 1;
			for(j = 0; j < cla[cl].ncomp; j++) if(cla[cl].comp[j].name == st) flag = 1;
		}
		if(flag == 0) break;
		num++;
	}while(1 == 1);
	return st;	
}

function addclass(na)                                     // Adds a new classification
{	
	var cl, clin, t, tr, j;

	if(na == "Time" || na == "Age") clin = ncla;
	else clin = ncla-2;

	for(i = 0; i < simindinit.length; i++){  // Initial state for simulation
		simindinit[i].state.splice(clin,0,"");
	}
		
	for(d = 0; d < data.length; d++){  // Adjusts fitlers
		switch(data[d].variety){
		case "trans": data[d].filt.splice(clin,0,"All"); break;
		case "pop": data[d].popcl.splice(clin,0,"All"); break;
		}
	}
					
	cla.splice(clin,0,{name:na, zoomfac:1, ncomp:0, comp:[], ntra:0, tra:[], desc:"Click to add a description...", descon:0});	
	ncla++;
	return clin;
}

function delclass(cldel)                                 // Deletes an existing classification
{
	var j, fr=[], cl, t, tr;

	for(i = 0; i < simindinit.length; i++){  // Initial state for simulation
		simindinit[i].state.splice(cldel,1);
	}
	
	d = 0;
	while(d < data.length){
		switch(data[d].variety){
		case "state": case "trans": case "move":
			if(data[d].cl == cldel) data.splice(d,1);
			else{ if(data[d].cl > cldel) data[d].cl--; d++;}
			break;
		default: d++; break;
		}
	}
	
	for(d = 0; d < data.length; d++){
		switch(data[d].variety){
		case "trans": data[d].filt.splice(cldel,1); break;
		case "pop": data[d].popcl.splice(cldel,1); break;
		}
	}

	fr.push(cla[cldel].name);
	for(j = 0; j < cla[cldel].ncomp; j++) fr.push(cla[cldel].comp[j].name);

	cla.splice(cldel,1); ncla--;
	
	for(j = 0; j < fr.length; j++) adjustrate(fr[j],"");
	pagesub[MODELPAGE]--; if(pagesub[MODELPAGE] == -1) addclassbubon();
}

function addclassbubon()                                   // Starts a popup message to add a classification
{
	setTimeout(function(){
		if(ncla == 2){
			selectbubx = 3; selectbuby = yclass + dysub/2;
			selectbubdx = 0; selectbubdy = 0; selectbub = ADDCLASSBUB; buttoninit();}
		},500);
}

function addcomp()                                         // Adds a compartment and selects
{
	cl = pagesub[MODELPAGE];
	
	togrid(cl,canmx,canmy);
		
	addcompartment(cl,dragtext,dragcol,gx-dragww/2,gy-dragdh/2,dragww,dragdh,0);
	
	addcompmode = 0;
	buttoninit();
	
	cl = pagesub[MODELPAGE];
	i = 0; while(i < ncanbut && !(canbuttype[i] == COMPBUT && canbutval[i] == cl && canbutval2[i] == cla[cl].ncomp-1)) i++;
	if(i == ncanbut) alertp("find");
	selbut(i);
}

function addcompartment(cl,name,col,x,y,w,h,simpopinit)   // Adds compartment to classification
{	
	cla[cl].comp.push({name:name, col:col, x:x, y:y, w:w, h:h, loadinit:[], simpopinit:simpopinit, simfracinit:undefined});
	cla[cl].ncomp++;
}

function deletecomp(cl,i)                                 // Deletes an existing compartment
{
	var tr, d, dat, j, fl;
	
	buboff(1);
	
	for(tr = cla[cl].ntra-1; tr >= 0; tr--){
		if(cla[cl].tra[tr].i == i || cla[cl].tra[tr].f == i) deletetrans(cl,tr);
		else{
			if(cla[cl].tra[tr].i > i) cla[cl].tra[tr].i--;
			if(cla[cl].tra[tr].f > i) cla[cl].tra[tr].f--;
		}
	}

	d = 0;  // adjusts the data
	while(d < data.length){
		fl = 0;
		dat = data[d];
		if(dat.variety == "state" && dat.cl == cl){
			for(j = 0; j < dat.pos.length; j++){
				if(dat.posref[j] == i) fl = 1;
				if(dat.posref[j] > i) dat.posref[j]--;	
				dat.posbinary[j].splice(i,1);
				dat.posexpression[j].splice(i,1);
			}
			dat.sensitive.splice(i,1);
		}
		if(fl == 1){ 
			points = []; points.push("Data source '"+dat.name+"' has been removed because compartment no longer exists.");
			helptype = 207;
			data.splice(d,1);
		}
		else d++;
	}
	
	cla[cl].comp.splice(i,1);
	cla[cl].ncomp--;
}

function addtransp()                                     // Adds intermediary point on a transition
{
	togrid(cl,mx-menux,my)
	tr = cla[transtempcl].tra[transtempk];
	tr.p.push({x:gx, y:gy});
	findpline(transtempcl,transtempk);
}

function addsource()                                     // Adds a source
{
	addsourcemode = 0;
	cl = pagesub[MODELPAGE];
	addtrans(cl,"Source");
	tr = cla[cl].tra[cla[cl].ntra-1];
	tr.i = -2;
	
	togrid(cl,canmx,canmy);
	tr.p=[];
	tr.p[0]={x:gx, y:gy};
	cl = pagesub[MODELPAGE];
}

function addtrans(cl,ty)                                  // Adds a transition
{
	addtransmode = 1;
	ty = ty.toUpperCase().substr(0,1) +	ty.toLowerCase().substring(1);	

	transtempcl = cl;
	transtempk = cla[cl].ntra;

	switch(ty){
	case "Exponential": case "Source": case "Sink":
		cla[cl].tra.push({type:ty, i:-1, f:-1, l:0, midx:0, midy:0, middx:0, middy:0, pline:[], p:[], ratetime:0, rate:"", nratevar:0, ratevar:[], ratevardep:[], pair:0});
		break;
		
	case "Gamma":
		cla[cl].tra.push({type:ty, i:-1, f:-1, l:0, midx:0, midy:0, middx:0, middy:0, pline:[], p:[], mean:"", shape:"", nratevar:0, ratevar:[], ratevardep:[], pair:0});
		break;
		
	case "Weibull":
		cla[cl].tra.push({type:ty, i:-1, f:-1, l:0, midx:0, midy:0, middx:0, middy:0, pline:[], p:[], lam:"", k:"", nratevar:0, ratevar:[], ratevardep:[], pair:0});
		break;
		
	case "Fixed":
		cla[cl].tra.push({type:ty, i:-1, f:-1, l:0, midx:0, midy:0, middx:0, middy:0, pline:[], p:[], time:"", nratevar:0, ratevar:[], ratevardep:[], pair:0});
		break;
		
	default: alertp("Error code EC40"); break;
	}
	
	cla[cl].ntra++;
}

function deletetrans(cl,i)                               // Deletes an existing transition
{
	var tr;
	if(cla[cl].tra[i].pair == 1){
		for(tr = 0; tr < cla[cl].ntra; tr++){
			if(cla[cl].tra[tr].i == cla[cl].tra[i].f && cla[cl].tra[tr].f == cla[cl].tra[i].i){
				cla[cl].tra[tr].pair = 0; findpline(cl,tr);
			}
		}			
	}
	
	cla[cl].tra.splice(i,1);
	cla[cl].ntra--;
}

function deletetransmode()                                // Stops the "add transtion" mode
{
	cla[transtempcl].tra.pop();
	cla[transtempcl].ntra--;
	addtransmode = 0;	
}

function completetrans(cl,k)                              // Completes a transition
{
	var tr;
	
	for(tr = 0; tr < cla[transtempcl].ntra-1; tr++){
		if(cla[cl].tra[tr].i == cla[cl].tra[transtempk].i && cla[cl].tra[tr].f == k){
			alertp("Cannot have more than one transition between the same compartments");
			deletetransmode(); return;
		}
		else{
			if(cla[cl].tra[tr].i == k && cla[cl].tra[tr].f == cla[cl].tra[transtempk].i){
				cla[cl].tra[transtempk].pair = 1; 
				cla[cl].tra[tr].pair = 1; 
				findpline(cl,tr);
			}
		}
	}
	
	cla[cl].tra[transtempk].f = k;
	findpline(cl,transtempk);
	seltrans(cl,transtempk);
	addtransmode = 0;
}

function seltrans(cl,i)                                  // Selects a transtion
{
	var tr;
	
	if(selectbub != -1) buboff(); if(selectbub != -1) return;

	tr = cla[cl].tra[i];
	selectbubx = tr.midx; selectbuby = tr.midy; selectbubdx = 0; selectbubdy = 0;
	switch(page){
	case MODELPAGE: selectbub = TRANSPLBUT; break;
	}
	selectbubval = cl;
	selectbubval2 = i;
}

function zoomin(cl,midx,midy,fac,dx,dy)                   // Zoom in
{
	togrid(cl,midx,midy);
	midx = gx; midy = gy;

	if(dx == undefined) dx = 0; if(dy == undefined) dy = 0;
	
	cla[cl].zoomfac *= fac;
	for(k = 0; k < cla[cl].ncomp; k++){
		c = cla[cl].comp[k];
		c.x = midx + fac*(c.x - midx) + dx;
		c.y = midy + fac*(c.y - midy) + dy;
		c.w *= fac; c.h *= fac;
	}
	for(k = 0; k < cla[cl].ntra; k++){
		tr = cla[cl].tra[k];
		for(j = 0; j < tr.p.length; j++){
			tr.p[j].x = midx + fac*(tr.p[j].x - midx) + dx;
			tr.p[j].y = midy + fac*(tr.p[j].y - midy) + dy;
		}
		for(j = 0; j < tr.pline.length; j++){
			tr.pline[j].x = midx + fac*(tr.pline[j].x - midx) + dx;
			tr.pline[j].y = midy + fac*(tr.pline[j].y - midy) + dy;
		}
		tr.midx = midx + fac*(tr.midx - midx) + dx;
		tr.midy = midy + fac*(tr.midy - midy) + dy;
		tr.l *= fac;
	}
}

function transover(cl,i,mmx,mmy)                         // Determines if mouse is over a transition
{	
	var tr, p, np, x, y, xx, yy, l;

			
	if(addtransmode == 1) return 0;
	
	tr = cla[cl].tra[i];
	
	if((tr.i < 0 || tr.f < 0) && modelty == "seltrans") return 0;
	
	if(tr.i < 0 && tr.f < 0) return 0;
	
	np = tr.pline.length;
	for(p = 0; p < np-1; p++){
		x = (tr.pline[p].x-shx)*shfac; y = (tr.pline[p].y-shy)*shfac;
		xx = (tr.pline[p+1].x-shx)*shfac; yy = (tr.pline[p+1].y-shy)*shfac;
		
		dx = xx-x; dy = yy-y;
		dmx = mmx-x; dmy = mmy-y;
		l = (dx*dmx + dy*dmy)/(dx*dx+dy*dy);
		if(l > 0 && l < 1){
			l = (dx*dmy - dy*dmx)/Math.sqrt(dx*dx+dy*dy);
			if(l < 5 && l > -5){ transoverj = p; return 1;}
		}
	}

	return 0;
}

function findpline(cl,i)                                  // Converts points in a transition to points on the canvas
{
	var tr, l, l2, np, p;

	tr = cla[cl].tra[i];
	tr.pline=[];
	if(tr.i != -1){
		k = tr.i;
		if(k >= 0){
			xbeg = cla[cl].comp[k].x + cla[cl].comp[k].w/2;
			ybeg = cla[cl].comp[k].y + cla[cl].comp[k].h/2;
		}
		
		k = tr.f;
		if(k >= 0){
			xend = cla[cl].comp[k].x + cla[cl].comp[k].w/2;
			yend = cla[cl].comp[k].y + cla[cl].comp[k].h/2;
		}
		
		np = tr.p.length;
		if(np == 0 && tr.pair == 1){
			dx = yend - ybeg; dy = -(xend - xbeg); r = Math.sqrt(dx*dx+dy*dy); dx /= r; dy /= r;
			fac = 0.05*(cla[cl].comp[tr.i].h + cla[cl].comp[tr.f].h);
			xbeg += fac*dx; ybeg += fac*dy; xend += fac*dx; yend += fac*dy; 
		}
		
		for(loop = 0; loop < 2; loop++){
			if(loop == 0){
				k = tr.i;
				if(k >= 0){
					xc = xbeg; yc = ybeg; 
					if(np > 0){ x = tr.p[0].x; y = tr.p[0].y;}else{ if(tr.f >= 0){ x = xend; y = yend;}}
				}
			}
			else{ 
				k = tr.f;
				if(k >= 0){				
					xc = xend; yc = yend; 
					if(np > 0){ x = tr.p[np-1].x; y = tr.p[np-1].y;}else{ x = xbeg; y = ybeg;}
				}
			}
			
			if(k >= 0){
				var lx=[], ly=[], jmax = 4;
				
				xx = cla[cl].comp[k].x; yy = cla[cl].comp[k].y;
				dx = cla[cl].comp[k].w; dy = cla[cl].comp[k].h;
				
				nl = 0;
				r = Math.floor(0.2*dy);
				for(j = 0; j <= jmax; j++){ th = j*Math.PI/(2*jmax); lx[nl] = xx+dx-r+r*Math.sin(th); ly[nl] = yy+r-r*Math.cos(th); nl++;}
				for(j = 0; j <= jmax; j++){ th = j*Math.PI/(2*jmax); lx[nl] = xx+dx-r+r*Math.cos(th); ly[nl] = yy+dy-r+r*Math.sin(th); nl++;}
				for(j = 0; j <= jmax; j++){ th = j*Math.PI/(2*jmax); lx[nl] = xx+r-r*Math.sin(th); ly[nl] = yy+dy-r+r*Math.cos(th); nl++;}
				for(j = 0; j <= jmax; j++){ th = j*Math.PI/(2*jmax); lx[nl] = xx+r-r*Math.cos(th); ly[nl] = yy+r-r*Math.sin(th); nl++;}

				nx =xc - x; ny = yc - y;

				for(ii = 0; ii < nl; ii++){
					x2 = lx[ii]; y2 = ly[ii];
					nx2 = lx[(ii+1)%nl]-x2; ny2 = ly[(ii+1)%nl]-y2;

					r = nx*ny2-ny*nx2;
					if(r != 0){
						al = ((x2-x)*ny2-(y2-y)*nx2)/r;
						be = ((x2-x)*ny-(y2-y)*nx)/r;
						if(al >= 0 && al <= 1 && be >= 0 && be <= 1){
							nx *= al; ny *= al;
						}
					}
				}
				tr.pline.push({x:x+nx, y:y+ny});
			}
			if(loop == 0){
				for(p = 0; p < np; p++){ tr.pline.push({x:tr.p[p].x, y:tr.p[p].y});}
			}
		}
	}

	np = tr.pline.length;
	l = 0;
	for(p = 0; p < np-1; p++){
		dx = tr.pline[p+1].x - tr.pline[p].x; dy = tr.pline[p+1].y - tr.pline[p].y;
		l += Math.sqrt(dx*dx+dy*dy);
	}

	l2 = 0; 
	for(p = 0; p < np-1; p++){
		dx = tr.pline[p+1].x - tr.pline[p].x; dy = tr.pline[p+1].y - tr.pline[p].y;
		dl = Math.sqrt(dx*dx+dy*dy);
		if(l2+dl > l/2){
			f = (l/2-l2)/dl;
			tr.middx = dx/Math.sqrt(dx*dx+dy*dy); tr.middy = dy/Math.sqrt(dx*dx+dy*dy);
			tr.midx = tr.pline[p].x+f*dx; tr.midy = tr.pline[p].y+f*dy;
			break;
		}
		l2 += dl;
	}
	tr.l = l;
}

function drawtrans(tr,ov,exp)                            // Draws a transition
{
	var fosize = "20px Times", fody=0, fosubsize = "12px Times", fosubdy= 4, loop, loopmin;
	
	if(exp == 1){
		fosize = "40px Times", fody=-8, fosubsize = "24px Times", fosubdy=16;
	}
	
	if(tr.i == -1 && tr.f == -1) return;
	
	switch(tr.type){
	case "Exponential": case "Source": case "Sink": dash = 0; break;
	case "Gamma": case "Weibull": dash = 1; break;
	case "Fixed": dash = 2; break;
	case "Grow": dash = 2; break;
	case "Set": dash = 2; break;
	}

	if(addtransmode == 1 && i == transtempk) col = GREY;
	else{ col = BLACK; if(ov == 1) col = DGREY;}

	if(page == SIMULATEPAGE){
		col = LGREY;
		if(pagesub[page] == 2 && simpagename[pagesubsub[page][2]] == "Gen. Data"){
			col = BLACK; if(ov == 1) col = RED;
		}
	}	
	
	if(page == INFERENCEPAGE && pagesub[page] == 0){ col = BLACK; if(ov == 1) col = RED;}
	
	if((tr.i < 0 || tr.f < 0) && modelty == "seltrans") col = LLGREY;

	if(tr.l*shfac > 25) li = THICKLINE; else li = NORMLINE;
	var arrowsi = 10;
	if(exp == 1){ li = VTHICKLINE; arrowsi = 20;}
	
	np = tr.pline.length;

	for(p = 0; p < np-1; p++){
		drawline((tr.pline[p].x-shx)*shfac,(tr.pline[p].y-shy)*shfac,(tr.pline[p+1].x-shx)*shfac,(tr.pline[p+1].y-shy)*shfac,col,li,dash);
	}
	if(tr.f >= 0 && tr.l*shfac > 25){
		drawarrow((tr.pline[np-1].x-shx)*shfac,(tr.pline[np-1].y-shy)*shfac,(tr.pline[np-2].x-shx)*shfac,(tr.pline[np-2].y-shy)*shfac,arrowsi,col);
	}
	
	if(tr.i != -1 && tr.f != -1 && tr.l*shfac > 40){	
		switch(tr.type){
		case "Exponential": case "Source": case "Sink": stbig = tr.rate; stsmall = "f_"; break;		
		case "Gamma": stbig = "Γ("+tr.mean+","+tr.shape+")";  stsmall = "Γ_"; break;
		case "Weibull": stbig = "W("+tr.lam+","+tr.k+")";  stsmall = "W_"; break;
		case "Fixed": stbig = "Fix("+tr.time+")"; stsmall = "f_"; break;
		case "Grow": stbig = "";  stsmall = ""; break;
		case "Set": stbig = ""; stsmall = ""; break;
		}
		
		switch(tr.type){
		case "Exponential": case "Source": case "Sink":case "Gamma": case "Weibull": case "Fixed":
			for(j = 0; j < tr.ratevar.length; j++){ stsmall += tr.ratevar[j]; if(j < tr.ratevar.length-1) stsmall += ",";}
			break;
		}
		
		if(stbig != ""){
			stbig = stbig.replace(/[\n\r]+/g, '');
			stbig = stbig.replace(/×/g,"");
			
			loopmin = 0; if(stbig.length > 30) loopmin = 1;
			for(loop = loopmin; loop < 2; loop++){
				switch(loop){
				case 0: st = stbig; break;
				case 1: st = stsmall; break;
				}
				
				var lab=[], j, jst;
				
				w = 0; 
				j = 0; 
				do{
					jst = j;
					while(j < st.length && st.substr(j,1) != "_") j++;
					frag = st.substr(jst,j-jst);
					frag = frag.replace(/\[/g,"");
					frag = frag.replace(/\]/g,"");
					lab.push({frag:frag, size:"big", w:w, font:fosize, dy:fody});
					w += textwidth(frag,fosize);
						
					if(j < st.length){
						j++; 
						jst = j;
						while(j < st.length && st.substr(j,1) != "]") j++;
						frag = st.substr(jst,j-jst);
						frag = frag.replace(/\[/g,"");
						frag = frag.replace(/\]/g,"");
						if(frag != ""){
							lab.push({frag:frag, size:"small", w:w, font:fosubsize, dy:fosubdy});
							w += textwidth(frag,fosubsize);
						}
						j++;
					}
				}while(j < st.length);
				
				if(w < 0.9*tr.l*shfac) break;
			}
			
			if(w < 0.9*tr.l*shfac){
				p = -8;
				midx = (tr.midx -shx)*shfac; midy = (tr.midy -shy)*shfac;
				
				nxmid = -tr.middy; nymid = tr.middx; 
				th = Math.acos(nymid);
				if(nxmid < 0) th *= -1;
				for(j = 0; j < lab.length; j++){
					if(tr.middx > 0){
						plotangletext(lab[j].frag,midx+(lab[j].w-w/2)*tr.middx+(-8+lab[j].dy)*nxmid,
											      midy+(lab[j].w-w/2)*tr.middy+(-8+lab[j].dy)*nymid,th,lab[j].font,col);
					}
					else{
						plotangletext(lab[j].frag,midx+(-lab[j].w+w/2)*tr.middx+(-19-lab[j].dy)*nxmid,
											      midy+(-lab[j].w+w/2)*tr.middy+(-19-lab[j].dy)*nymid,th+Math.PI,lab[j].font,col);
					}
				}
			}
		}
	}
	
	if(col == LLGREY) col2 = LGREY; else col2 = DGREY;
	r = 15; if(exp == 1) r = 30;
	if(tr.i == -2){
		x = (tr.pline[0].x-shx)*shfac; y = (tr.pline[0].y-shy)*shfac;
		fillcircle(x,y,r,col);
		circle(x,y,r,col2,NORMLINE);
		if(exp == 1) centertext("+",x,y+17, "bold 50px arial",WHITE);	
		else centertext("+",x,y+11, "bold 30px arial",WHITE);	
	}
	
	if(tr.f == -2){
		x = (tr.pline[np-1].x-shx)*shfac; y = (tr.pline[np-1].y-shy)*shfac;	
		fillcircle(x,y,r,col);
		circle(x,y,r,col2,NORMLINE);
		if(exp == 1) centertext("-",x,y+16,"bold 50px arial",WHITE);	
		else centertext("-",x,y+8,"bold 30px arial",WHITE);	
	}
}

/*
function labtrunc(stin,tr)                               // Truncates a label (so it isn't too long)
{
	var li = stin.split("\n"), j, te, st, stsmall; 
	st = "";for(j = 0; j < li.length; j++){ te = li[j].trim(); if(te != ""){ if(st != "") st += ","; st += te;}}

	//if(st.length > 30){
	stsmall = "f_";
	for(j = 0; j < tr.ratevar.length; j++){ stsmall += tr.ratevar[j]; if(j < tr.ratevar.length-1) st += ",";}
	//if(st.length > 30) st = st.substr(0,30)+"...";
	stsmall += ")";
	//}
	return {st:st, stsmall:stsmall};
}
*/

function togrid(cl,x,y)                                   // Shifts canvas points onto a grid
{
	var f = gridsize*cla[cl].zoomfac;
	gx = Math.round(x/f)*f;
	gy = Math.round(y/f)*f;
}

function breaktrans()                                     // Add intermidiary point in a transition
{
	cl = dragval;
	k = dragval2;
			
	tr = cla[cl].tra[k];
	
	np = tr.p.length;
	tr.p[np] = {x:0,y:0};
	if(np > 0){
		for(j = np; j > transoverj; j--){ tr.p[j].x = tr.p[j-1].x; tr.p[j].y = tr.p[j-1].y;}
	}
	togrid(cl,mx-menux,my)
	tr.p[transoverj].x = gx; tr.p[transoverj].y = gy;
	
	dragval = dragval2;
	dragval2 = transoverj;
	mxst = mx; myst = my; mxlast = mx; mylast = my;
	drag = 5;
}

function removetransp(k,j)                               // Removes intermidiary point in a transition 
{
	cl = pagesub[MODELPAGE];
	tr = cla[cl].tra[k];
	
	if((tr.i == -2 && j == 0) || (tr.f == -2 && j == tr.p.length-1)) return;
	tr.p.splice(j,1);
	findpline(cl,k);
}

function adjustrate(fr,to)                           	 // Changes rates when model changes 
{
	var cl, tr, cl2, j;
	
	for(cl = 0; cl < ncla-2; cl++){
		for(tr = 0; tr < cla[cl].ntra; tr++){
			var tra =  cla[cl].tra[tr];
			
			for(j = 0; j < tra.ratevardep.length; j++) tra.ratevardep[j] = adjusttext(tra.ratevardep[j],fr,to);
		
			switch(tra.type){
			case "Exponential": case "Source": case "Sink": tra.rate = adjusttext2(tra.rate,fr,to); break;
			case "Gamma": tra.mean = adjusttext2(tra.mean,fr,to); tra.shape = adjusttext2(tra.shape,fr,to); break;
			case "Weibull": tra.lam = adjusttext2(tra.lam,fr,to); tra.k = adjusttext2(tra.k,fr,to); break;
			}	
		}
	}
}

function adjusttext(te,fr,to)                           // Adjusts text to reflect change in model 
{
	var list = te.split(","), k;
	for(k = 0; k < list.length; k++){ if(list[k] == fr) list[k] = to;}
	te = ""; for(k = 0; k < list.length; k++){ te += list[k]; if(k < list.length-1) te += ",";}
	return te;
}
				
function adjusttext2(te,fr,to)                          // Adjusts text to reflect change in model 
{
	var i, ty, fl;

	i = 0;
	while(i < te.length){
		while(i < te.length && te.substr(i,1) != "[" && te.substr(i,1) != "{") i++;
		if(i < te.length){
			ty = te.substr(i,1);
			fl = 0;
			if(ty == "["){
				while(i < te.length && te.substr(i,1) != "]" && te.substr(i,1) != "_") i++;
				if(te.substr(i,1) != "_") fl = 1;
			}
			i++;
			if(fl == 0){
				ist = i;
				while(i < te.length && !(ty == "[" && te.substr(i,1) == "]")
							  	    && !(ty == "{" && te.substr(i,1) == "}")) i++;
				if(i < te.length){
					te = te.substr(0,ist) + adjusttext(te.substr(ist,i-ist),fr,to) + te.substring(i);
					i = ist+1;
				}
			}
		}		
	}
	return te;
}

function getclfromcomp(cnam)                           // Gets the classification number from compartment name
{
	for(cl = 0; cl < ncla; cl++){
		for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == cnam) break;
		if(c < cla[cl].ncomp) break;
	}
	if(cl == ncla) alert("Cannot get cl from comp");
	return cl;
}
	
function checkrate(cl,i)                               // Checks the rate expression
{
	var tr = cla[cl].tra[i];
	
	nvarlist = 0; varlist=[]; vardeplist=[]; // Makes a list of variables in transition
	
	switch(tr.type){
	case "Exponential": case "Source": case "Sink": checkeqn(tr.rate); break;
	case "Gamma": checkeqn(tr.mean); checkeqn(tr.shape); break;
	case "Weibull": checkeqn(tr.lam); checkeqn(tr.k); break;
	case "Fixed": checkeqn(tr.time); break;
	}
	
	if(errmsg != ""){ inpon = 0; return;}

	tr.nratevar = nvarlist;
	tr.ratevar=[]; tr.ratevardep=[];
	for(n = 0; n < nvarlist; n++){ tr.ratevar[n] = varlist[n]; tr.ratevardep[n] = vardeplist[n];}
}

function getrateforfile(st,ratetime)                     // Converts rate expression into xml format
{
	var lines=[], ratelist=[], filt=[];
	
	lines = getratelines(st);

	for(li = 0; li < lines.length; li++){
		temp = lines[li].split(":");
		switch(temp.length){
		case 1:
			if(ratetime == 1) ra = "1/("+lines[0]+")"; else ra = lines[0];
			for(cl = 0; cl < ncla; cl++) filt[cl] = "All";
			ratelist.push({rate:ra, filt:filt});
			break;
			
		case 2:
			if(temp.length != 2) errmsg("not two");
			temp2 = temp[0].split(",");
			var filt=[]; for(cl = 0; cl < ncla; cl++) filt[cl] = "All";
			for(j = 0; j < temp2.length; j++){
				for(cl = 0; cl < ncla; cl++){
					for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == temp2[j]) break;
					if(c < cla[cl].ncomp) break;
				}
				if(cl == ncla) alertp("Error code EC41");
				filt[cl] = temp2[j];
			}
			if(ratetime == 1) ra = "1/("+temp[1]+")"; else ra = temp[1];
		   	ratelist.push({rate:ra, filt:filt});
			break;
		
		default: errmsg("Error code EX98"); break;
		}
	}
	return ratelist;
}

function eqntidy(st)
{
	st = st.replace(/[\r]/g, '');
	st = st.replace(/[\t]/g, '');
	st = st.replace(/ /g, '');
	st = st.replace(/\u00d7/g, '*');
	st = st.replace(/\u2215/g, '/');	
	return st;
}

function getratelines(st)                               // Converts rate expression into lines
{
	var lines=[];

	st = eqntidy(st);
	
	lines = st.split('\n');
	li = 0; 
	while(li < lines.length){
		lines[li] = lines[li].trim();
		if(lines[li] == "") lines.splice(li,1);
		else li++;
	}

	li = 1;
	while(li < lines.length){
		temp = lines[li].split(":");
		if(temp.length > 2){ errmsg = "Two colons on one line";}
		if(temp.length == 1){ lines[li-1] += lines[li]; lines.splice(li,1);}
		else li++;
	}
	return lines;
}

function checkeqn(st)                              // Checks the validity of an equation
{
	var an=[], j, jj, ii, val, nnn, lines=[], def=[];

	lines = getratelines(st);
	if(lines.length == 0){ errmsg = "Please enter an expression"; return;}

	for(li = 0; li < lines.length; li++){
		st = lines[li];
		for(cl2 = 0; cl2 < ncla; cl2++) def[cl2] = "";
			
		//if(lines.length > 1){ // looks for the prefix
		temp = st.split(":");
		if(temp.length > 2){ errmsg = "Line should only contain one colon"; return;}
		
		if(temp.length == 2){	
			an = temp[0].split(",");  // makes sure the bit before the colon is correct
			for(k = 0; k < an.length; k++){
				
				for(cl2 = 0; cl2 < ncla; cl2++){
					for(c = 0; c < cla[cl2].ncomp; c++) if(cla[cl2].comp[c].name == an[k]) break;
					if(c < cla[cl2].ncomp){ 
						if(def[cl2] != ""){ errmsg =cla[cl2].name+" is defined more than once."; return;}
						def[cl2] = an;		
						break;
					}
				}
				if(cl2 == ncla){ errmsg = "'"+an[k]+"' is not a compartment name"; return;}
			}	
			st = temp[1];
		}
		stst = st;

		ii = 0;                                // removes variables
		while(ii < st.length){
			if(st.substr(ii,1) == "["){
				j = ii+1; while(j < st.length && st.substr(j,1) != "]") j++;
				if(j == st.length){ errmsg = "Brackets do not match up"; return;}
				
				for(jj = ii+1; jj < j; jj++){
					if(st.substr(jj,1) == " "){ errmsg = "Cannot have spaces in variables"; return;}
				}				
			
				jj = ii+1; while(jj < j && st.substr(jj,1) != "_") jj++;
				val = st.substr(ii+1,jj-(ii+1));
				if(val == ""){ errmsg = "Set variable name"; return;} 
				nnn = 0; while(nnn < nvarlist && varlist[nnn] != val) nnn++;
				if(nnn == nvarlist){
					varlist[nvarlist] = val;
					if(jj == j) vardeplist[nvarlist] = "";
					else vardeplist[nvarlist] = st.substr(jj+1,j-(jj+1)); 
					nvarlist++;
				}
				if(jj < j){  // makes sure the class names are correct
					an = st.substr(jj+1,j-(jj+1)).split(',');
					for(k = 0; k < an.length; k++){
						for(cl2 = 0; cl2 < ncla; cl2++){
							if(cla[cl2].name == an[k]){ cladepflag[cl2] = 1; break;}
						}
						if(cl2 == ncla){ errmsg = "'"+an[k]+"' is not a classifcation name"; return;}
					}
				}
				st = st.substr(0,ii)+"[&q]"+st.substr(j+1);
				ii += 3;
			}
			else ii++;
		}

		ii = 0;                                       // removes variables
		while(ii < st.length-1){
			if(st.substr(ii,1) == "{"){
				popflag = 1;
				 
				j = ii+1; while(j < st.length && st.substr(j,1) != "}") j++;
				if(j == st.length){ errmsg = "Brackets do not match up"; return;}
				
				st2 = st.substr(ii+1,j-(ii+1));
				if(st2 == "all" || st2 == ""){ st = st.substr(0,ii)+"[&q]"+st.substr(j+1); ii += 4;}
				else{						
					an = st2.split(',');
					for(k = 0; k < an.length; k++){
						for(cl2 = 0; cl2 < ncla; cl2++){
							if(cla[cl2].name == an[k]){ cladepflag[cl2] = 1; break;}
							for(kk = 0; kk < cla[cl2].ncomp; kk++) if(cla[cl2].comp[kk].name == an[k]) break;
							if(kk < cla[cl2].ncomp) break;
						}
						if(cl2 == ncla){ errmsg = "'"+an[k]+"' is not recognised"; return;}
					}
					if(k == an.length){
						st = st.substr(0,ii)+"[&q]"+st.substr(j+1);
						ii += 3;
					}
					else break;
				}
			}
			else ii++;
		}
	                            
		ii = 0;      // removes numbers
		while(ii < st.length){
			if(isnum(st.substr(ii,1)) == 1){
				j = ii+1; while(j < st.length && isnum(st.substr(j,1)) == 1) j++;
				if(j < st.length && st.substr(j,1) == "."){
					j++; while(j < st.length && isnum(st.substr(j,1)) == 1) j++;
				}
				st = st.substr(0,ii)+"[&q]"+st.substr(j);
				ii += 4;
			}
			else ii++;
		}
	
		st = st.replace(/\+/g,"[&o]");
		st = st.replace(/-/g,"[&m]");
		st = st.replace(/\//g,"[&o]");
		st = st.replace(/\u2215/g,"[&o]");
		st = st.replace(/\*/g,"[&o]");
		st = st.replace(/\^/g,"[&o]");
		st = st.replace(/\u00d7/g,"[&o]");
		
		do{   // tries to reduce
			lenbeg = st.length;	
			st = st.replace("[&q][&o][&q]","[&q]");
			st = st.replace("[&q][&m][&q]","[&q]");
			st = st.replace("[&q][&q]","[&q]");
			st = st.replace("([&q])","[&q]");
			st = st.replace("cos[&q]","[&q]");
			st = st.replace("sin[&q]","[&q]");
			st = st.replace("exp[&q]","[&q]");
			st = st.replace("log[&q]","[&q]");
			
			if(st.length == lenbeg){
				st = st.replace("[&m][&q]","[&q]");
			}
		}while(st.length > 4 && st.length < lenbeg);
		
		if(st != "[&q]"){
			errmsg = "There is a problem with this expression!"; return;
		}
	}
}

function isnum(st)                                      // Checks for integer number
{
	if(st == "0" || st == "1" || st == "2" || st == "3" || st == "4" || st == "5" || st == "6" || st == "7" || st == "8" || st == "9") return 1; else return 0;  
}

function settimeon()                                    // Adds the "set time" timeline on
{
	var cl, j, na=[], x, xx=[], w, h;
	
	cl = 0; while(cl < ncla && cla[cl].name != "Time") cl++;
	if(cl == ncla) cl = addclass("Time");
	
	cla[cl].comp=[];
	cla[cl].ncomp=0;

	x = 0;
	for(j = 0; j < time.length+1; j++){
		if(j == 0){
			if(j == time.length) st = "All";
			else st = "T<"+time[j];
		}
		else{
			if(j == time.length) st = "T>"+time[j-1];
			else st = "T"+ time[j-1]+"-"+time[j];
		}
		na[j] = st;
	
		xx[j] = x;
		x += 1; if(j < time.length) x += 0.6;
	}
	
	widx = width-menux; widy = height;
	sca = (widx-100)/x; if(sca > 200) sca = 200;
	w = sca; h = w/3;
	
	for(j = 0; j < na.length; j++){
		if(na.length == 1) co = PURPLE;
		else{ f = j/(na.length-1); co = "#" + hex(255*f+120*(1-f)) + hex(150*f+0*(1-f)) + hex(255*f+120*(1-f));}
		cla[cl].comp.push({name:na[j], col:co, x:Math.floor(widx/2-sca*x/2+xx[j]*sca), y:widy/2-h/2, w:w, h:h});
		cla[cl].ncomp++;
	}
	
	cla[cl].tra=[];
	cla[cl].ntra=0;

	for(j = 0; j < na.length-1; j++){
		cla[cl].tra.push({type:"Set", time:time[j], i:j, f:j+1, l:0, midx:0, midy:0, middx:0, middy:0, p:[], pline:[]});
		cla[cl].ntra++;
	}
	for(tr = 0; tr < cla[cl].ntra; tr++) findpline(cl,tr);
}
	
function setageon()                                     // Adds the age timeline
{
	var cl, j, na=[], x, xx=[], w,h;
	
	cl = 0; while(cl < ncla && cla[cl].name != "Age") cl++;
	if(cl == ncla) cl = addclass("Age");
	cla[cl].comp=[];
	cla[cl].ncomp=0;

	x = 0;
	for(j = 0; j < age.length+1; j++){
		if(j == 0){
			if(j == age.length) st = "A0+";
			else st = "A0-"+age[j];
		}
		else{
			if(j == age.length) st = "A"+age[j-1]+"+";
			else st = "A"+ age[j-1]+"-"+age[j];
		}
		na[j] = st;
	
		xx[j] = x;
		x += 1; if(j < age.length) x += 0.6;
	}
	
	widx = width-menux; widy = height;
	
	sca = (widx-100)/x; if(sca > 200) sca = 200;
	w = sca; h = w/3;
	
	for(j = 0; j < na.length; j++){
		if(na.length == 1) co = GREY;
		else{ f = j/(na.length-1); co = "#" + hex(230*f+60*(1-f)) + hex(230*f+60*(1-f)) + hex(230*f+60*(1-f));}
		cla[cl].comp.push({name:na[j], col:co, x:Math.floor(widx/2-sca*x/2+xx[j]*sca), loadinit:[], simpopinit:0, simfracinit:1.0/na.length, y:widy/2-h/2, w:w, h:h});
		cla[cl].ncomp++;
	}
	
	cla[cl].tra=[];
	cla[cl].ntra=0;

	for(j = 0; j < na.length-1; j++){
		cla[cl].tra.push({type:"Grow", age:age[j], i:j, f:j+1, l:0, midx:0, midy:0, middx:0, middy:0, p:[], pline:[]});
		cla[cl].ntra++;
	}
	for(tr = 0; tr < cla[cl].ntra; tr++) findpline(cl,tr);
}

function checkname(st,except)                           // Check a new name for validity
{
	var cl2, j, k, sch="_[]{},*=|'\"";
	 
	if(st == "") errmsg = "Must enter name";
	else{
		if(st == "Age") errmsg = "Cannot use 'Age'";
		if(st == "Time") errmsg = "Cannot use 'Time'";
		
		for(cl2 = 0; cl2 < ncla; cl2++){
			if(st == cla[cl2].name && st != except) errmsg = "Same name as a classification";
		}
		
		for(j = 0; j < sch.length; j++){
			for(k = 0; k < st.length; k++) if(st.substr(k,1) == sch.substr(j,1)) errmsg = "Cannot use '"+sch[j]+"' in the name";
		}
	}
}

function startnewmodel()                                // Initialises a new model
{
	if(modelstart == 1){
		if(confirm("Are you sure you want to discard the currently loaded model?")); else return; 
	}
	initvar();
	
	nchrun=0;
	derive=[];
	data=[]; simdata=[];
	param=[]; paramsim=[];
	 
	ncla=0; cla=[];
	age=[];
	time=[];
	setageon(); settimeon();
	
	simres={}; infres={};
	modelstart = 1; modelsetup = -1;
	induoflag = 0;
	fileToLoadlist=[];
	examploaded = "New";
	
	changepage(DESCPAGE,-1,0);
}

function initvar()                                     // Initialises variables for a new model
{
	var i;
	siminit="";
	
	datanote = "Click here to place a comment about the data."; datanoteon = 0;
	descnote = "Place a comment about the data/analysis (optional).";
 
 	nsampevmax = 1000; nsampmax = 10000; GRmax = 1.1; ESSmin = 400;
	indmaxnumber = 100000;
	
	nchain = 3;
	tpostmin=large; tpostmax=-large;
	termtype = 0;
	
	simpopinitset = 0;
	tsimmin = 0; tsimmax = 100; simty = 0;
	detectprob=1; datagenttype=0; datagenttype2=0; tgenmin=0; tgenmax=100; dtgen=10; tgenuserdef="";
	Segen=1; Spgen=1;
	addtransmode = 0; addcompmode = 0; addsourcemode = 0;
	kde = 10; kdemin = 1; kdemax = 20;   
	distCI = 1; distmean = 1;
	initpages();
}

function addderive(newna)                                // Adds a derived quantity
{
	var i;
	
	for(i = 0; i < newna.length; i++){
		var c = newna.substr(i,1);
		if(c == " " || c == "_" || c == "," || c == "[" || c == "]" || c == "{" || c == "}"){
			return 1;
		}				
	}
			
	derive.push({name:newna, eq:"", popdep:0, dep:[], ndervar:0, dervar:[]});
	return 0;
}		
	
function setderive(d,val)                                // Sets a derived quantity
{
	var cl2, i;
	
	nvarlist = 0; varlist=[]; vardeplist=[]; popflag = 0; cladepflag=[]; for(cl2 = 0; cl2 < ncla; cl2++) cladepflag[cl2] = 0;
	checkeqn(val);
	if(errmsg != "") return 1;
	else{
		for(i = 0; i < nvarlist; i++){
			if(vardeplist[i].length > 1){ errmsg = "Variable "+vardep[i]+ " cannot have any dependencies"; return 1;}
		}
		
		derive[d].eq = val;
		derive[d].popdep = popflag;
		derive[d].dep = []; for(cl2 = 0; cl2 < ncla; cl2++) if(cladepflag[cl2] == 1) derive[d].dep.push(cla[cl2].name);
		
		pr(derive);
		pr(d+"d");
		derive[d].ndervar = nvarlist; for(i = 0; i < nvarlist; i++) derive[d].dervar[i] = varlist[i]; 
		return 0;
	}
}	

function drawderived()                                  // Draws the derived buttons
{
	var x, y;
	
	x = 0; y = 8; nob = 0;
	
	for(d = 0; d < derive.length; d++){ addob(x,y,OBDERIVE,d); y += priordy;}
	y += 40;
		
	ytot = y;
	
	placeob();	
}

function donemodel()                                     // When "done is clicked on model page
{	
	buboff();
	if(modelsetup == 1){
		modelsetup = 0;
		if(runtype != ""){ stopinference();}
		infres={}; simres={};
		pagesubsub[SIMULATEPAGE][2] = 0; pagesubsub[INFERENCEPAGE][2] = 0; 
	}
	else definemodel();	
}

function definemodel()                                   // Called when the model is defined
{
	var cl, tr, d, j, ci, cf;

	errmsg="";
	for(cl = 0; cl < ncla; cl++){
		for(tr = 0; tr < cla[cl].ntra; tr++){
			checkrate(cl,tr);
			if(errmsg != ""){
				changepage(MODELPAGE,0,-1); seltrans(cl,tr); errormsg = errmsg; return 1;
			}
		}
	}

	emsg = "";
	for(cl = 0; cl < ncla; cl++){
		if(cla[cl].ncomp == 0){ emsg = "Classifciation '"+cla[cl].name+"' must have at least one compartment!"; break;}
	}
	if(emsg != ""){ alertp(emsg); errormsg = emsg; return 1;}
	
	if(checkprior(1) == 1){ errormsg = errmsg; return 1;}
	
	
	collectvariables();
	for(cl = 0; cl < ncla; cl++){
		if(!popinit[cl]) popinit[cl]=[];
		for(j = 0; j < cla[cl].ncomp; j++){
			if(cl == 0){ if(!popinit[cl][j]) popinit[cl][j] = 10;}
			else{
				popinit[cl][j] = Math.floor(100/cla[cl].ncomp);
			}
		}
	}
	
	simpopinitlist=[];
	for(cl = 0; cl < ncla-1; cl++) simpopinitlist.push(cla[cl].name);
	if(simpopinitlist.length == 0) simpopinitlist.push(cla[0].name);

	for(cl = 0; cl < ncla-2; cl++){     // looks at whether simfracinit needs to be defined
		sum = 0;
		for(j = 0; j < cla[cl].ncomp; j++){
			sum += cla[cl].comp[j].simfracinit;
			if(cla[cl].comp[j].simfracinit == undefined) break;
		}

		if(j < cla[cl].ncomp || sum > 1.000001 || sum < 0.9999999){
			for(j = 0; j < cla[cl].ncomp; j++) cla[cl].comp[j].simfracinit = 1.0/cla[cl].ncomp;
		}
	}

	claselect=[]; for(cl = 0; cl < ncla; cl++) claselect.push(cla[cl].name);
		
	for(cl = 0; cl < ncla; cl++){ if(cl == 0 || cla[cl].ncomp > 1) indshow[cl] = 1; else indshow[cl] = 0;}
	
	initcompmult();
		
	points=[];
	d = 0;
	while(d < data.length){ // changes need to the data and observation models
		dat = data[d];
		fl = 0;
		switch(dat.variety){
		case "trans": case "move":
			switch(dat.transty){
			case "trans":
				cl = dat.cl;
				if(cl >= ncla){ fl = 2; points.push("Data source '"+dat.name+"' has been removed because classification no longer exists."); break;}
				
				for(tr = 0; tr < cla[cl].ntra; tr++){
					ci = cla[cl].tra[tr].i;
					cf = cla[cl].tra[tr].f; 
					if(dat.transi == cla[cl].comp[ci].name && dat.transf == cla[cl].comp[cf].name) break;
				}
				if(tr == cla[cl].ntra){ fl = 2; points.push("Data source '"+dat.name+"' has been removed because transition no longer exists."); break;}
				break;
			
			case "+":
				for(cl = 0; cl < ncla; cl++){
					for(tr = 0; tr < cla[cl].ntra; tr++){
						ci = cla[cl].tra[tr].i;
						cf = cla[cl].tra[tr].f;
						if(ci == -2){
							if(dat.variety == "move") break;
							if(dat.filt[cl] == "All" || dat.filt[cl] == cla[cl].comp[cf].name) break; 
						}
					}
					if(tr < cla[cl].ntra) break;
				}
				if(cl == ncla){ fl = 2; points.push("Data source '"+dat.name+"' has been removed because source transition no longer exists."); break;}
				break;
				
			case "-":
				for(cl = 0; cl < ncla; cl++){
					for(tr = 0; tr < cla[cl].ntra; tr++){
						ci = cla[cl].tra[tr].i;
						cf = cla[cl].tra[tr].f;
						if(cf == -2){
							if(dat.variety == "move") break;
							if(dat.filt[cl] == "All" || dat.filt[cl] == cla[cl].comp[ci].name) break; 
						}
					}
					if(tr < cla[cl].ntra) break;
				}
				if(cl == ncla){ fl = 2; points.push("Data source '"+dat.name+"' has been removed because source transition no longer exists."); break;}
				break;
			}
			
			if(dat.variety == "trans"){
				if(dat.filt.length != ncla) alert("wrong size");
				for(cl = 0; cl < ncla; cl++){
					na = dat.filt[cl];
					if(na != "All"){ 
						for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == na) break; 
						if(c == cla[cl].ncomp) alert("transfilt problem");
					}
				}
			}
			break;
			
		case "pop":
			if(dat.popcl.length != ncla) alert("wrong size");
			for(cl = 0; cl < ncla; cl++){
				na = dat.popcl[cl];
				if(na != "All"){ 
					for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == na) break; 
					if(c == cla[cl].ncomp) alert("popfilt problem");
				}
			}
			break;
			
		case "state":
			cl = dat.cl; 
			if(cl >= ncla){fl = 2; points.push("Data source '"+dat.name+"' has been removed because classification no longer exists.");}
			else{
				nc = cla[cl].ncomp;
				len = dat.sensitive.length;
				if(len < nc){ fl = 1; for(k = len; k < nc; k++) dat.sensitive[k] = 0;}
				else{ if(len > nc){ fl = 1; dat.sensitive.length = nc;}}
				
				for(j = 0; j < dat.pos.length; j++){
					len = dat.posbinary[j].length;
					if(len < nc){ fl = 1; for(k = len; k < nc; k++) dat.posbinary[j][k] = 0;}
					else{ if(len > nc){ fl = 1; dat.posbinary[j].length = nc;}}
				
					len = dat.posexpression[j].length;
					if(len < nc){ fl = 1; for(k = len; k < nc; k++) dat.posexpression[j][k] = "";}
					else{ if(len > nc){ fl = 1; dat.posexpression[j].length = nc;}}	
				}
			}
			break;
		}
		
		if(fl != 2){
			datatemp = data[d];
			reloaddata();
			if(checkdata() == 1){ alertp(errormsg); return 1;}
		}
		
		if(fl == 2) data.splice(d,1);
		else{
			if(fl == 1) points.push("As a result of the alteration to the model, changes have been made to the observation model for data '"+dat.name+"'. Please check!");
			d++;
		}
	}

	modelsetup = 1; 
	if(points.length > 0) helptype = 207;
	else helptype = 209;
	
	converttoobs(ty);
	
	return 0;
}
