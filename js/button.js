function buttoninit()                                      // Initialises the buttons used on pages
{
	var dx = 120, dx2 = 100,  i, overnew = -1, xx, ddy, dddy;
	//var dy = 32, ddy = 42, dby = 13;
	var ddy, dby;

	if(pagesize == "small"){ ddy = 34; dby = 9, dysub = 23, dysubsub = 17;}
	else{ ddy = 38; dby = 12, dysub = 23, dysubsub = 19;}
	
	over = -1; canover = -1;
	nbut = 0; ncanbut = 0;
	gdropinfo=[];
	
	ninput2 = 0;

	if(addingdata > 0) addingdatabuts();
	else{
		switch(page){
		case HOMEPAGE: homebuts(); break;
		case DESCPAGE: descbuts(); break;
		case MODELPAGE: modelbuts(); break;
		case SIMULATEPAGE: simulatebuts(); break;
		case INFERENCEPAGE: inferencebuts(); break;
		}
	}
 
	addbutton("",0,0,menux,canh,-1,MENUBACKBUT,-1,-1);

	y = 90; 
 
	addbutton("Home",0,y,menux,ddy,TABBUT,TABBUT,HOMEPAGE,-1); y += ddy; y += dby;

	addbutton("Description",0,y,menux,ddy,TABBUT,TABBUT,DESCPAGE,-1); y += ddy; y += dby;

	addbutton("Model",0,y,menux,ddy,TABBUT,TABBUT,MODELPAGE,-1); y += ddy;
	
	if(page == MODELPAGE){
		ymodel = y;
		nout = 0; for(cl = 0; cl < ncla-2; cl++){ out[nout] = cla[cl].name; nout++;}
		addsub();
		xx = 15;
			
		if(modelsetup <= 0){
			te = "Class.";
			yclass = y;
			addbutton("",xx,y,menux-15,dysub,-1,PAGESUBBACKBUT,-1,-1);
			addbutton(te,xx,y,menux-xx,dysub,ADDCLASSBUT,ADDCLASSBUT,-1,-1);
			if(pagesub[MODELPAGE] == -1) fo = MENUFONTSMALL; else fo = MENUFONTSMALL2;
			addbutton("[?]",xx+37+textwidth(te,fo),y+3,15,15,HELPICONBUT,HELPICONBUT,33,0);
			y += dysub;
		}
		
		//y -= 2;
		if(ncla > 2){
			te = "Age";
			addbutton("",xx,y,menux-15,dysub,-1,PAGESUBBACKBUT,-1,-1);
			addbutton(te,xx,y,menux-15,dysub,PAGESUBBUT,PAGESUBBUT,ncla-2,-1);
			if(pagesub[MODELPAGE] == ncla-2) fo = MENUFONTSMALL; else fo = MENUFONTSMALL2;
			addbutton("[?]",xx+24+textwidth(te,fo),y+3,15,15,HELPICONBUT,HELPICONBUT,30,-2);
			y += dysub;
			te = "Time";
			addbutton("",xx,y,menux-15,dysub,-1,PAGESUBBACKBUT,-1,-1);
			addbutton(te,xx,y,menux-15,dysub,PAGESUBBUT,PAGESUBBUT,ncla-1,-1);
			if(pagesub[MODELPAGE] == ncla-1) fo = MENUFONTSMALL; else fo = MENUFONTSMALL2;
			addbutton("[?]",xx+24+ textwidth(te,fo),y+3,15,15,HELPICONBUT,HELPICONBUT,31,-2);
			y += dysub;
		}
		
		te = "Derived";
		addbutton("",xx,y,menux-15,dysub,-1,PAGESUBBACKBUT,-1,-1);
		addbutton(te,xx,y,menux-15,dysub,PAGESUBBUT,PAGESUBBUT,-2,-1);
		if(pagesub[MODELPAGE] == -2) fo = MENUFONTSMALL; else fo = MENUFONTSMALL2;
		y += dysub;
		
		te = "Distributions";
		addbutton("",xx,y,menux-15,dysub,-1,PAGESUBBACKBUT,-1,-1);
		addbutton(te,xx,y,menux-15,dysub,PAGESUBBUT,PAGESUBBUT,-3,-1);
		if(pagesub[MODELPAGE] == -3) fo = MENUFONTSMALL; else fo = MENUFONTSMALL2;
		y += dysub;			
	}
	y += dby;
	
	addbutton("Simulation",0,y,menux,ddy,TABBUT,TABBUT,SIMULATEPAGE,-1); y += ddy;
  
	if(page == SIMULATEPAGE){
		nout = 0;
		out[nout]="Parameters"; nout++;
		out[nout]="Init. Pop."; nout++;
		out[nout]="Run"; nout++;
	
		addsub();
	
		if(pagesub[SIMULATEPAGE] == 2){
			nout = 0;
			out[nout]="Start"; nout++;
			if(simres.result >= 1){
				out[nout]="Populations"; nout++;
				if(simres.nderive > 0){ out[nout]="Derived"; nout++;}
				out[nout]="Individuals"; nout++;
				out[nout]="Statistics"; nout++;
				out[nout]="Gen. Data"; nout++;
			}
			addsubsub();
			for(j = 0; j < nout; j++) simpagename[j] = out[j];
		}
	}
	y += dby;
	
	addbutton("Inference",0,y,menux,ddy,TABBUT,TABBUT,INFERENCEPAGE,-1); y += ddy;

	if(page == INFERENCEPAGE){
		nout = 0;
		out[nout]="Data"; nout++; 
		addsub();
	
		if(pagesub[INFERENCEPAGE] == 0 && data.length > 0){
			nout = 0;
			out[nout]="Sources"; nout++;
			out[nout]="Individuals"; nout++;
			addsubsub();
			for(j = 0; j < nout; j++) infpagename[j] = out[j];
		}
		
		nout = 0;
		out[nout]="Prior"; nout++; 
		out[nout]="Posterior"; nout++; 
		addsub(1);
		
		if(pagesub[INFERENCEPAGE] == 2){
			nout = 0;
			out[nout]="Start"; nout++; 
			if(infres.result >= 1){
				out[nout]="Traces"; nout++; 
				out[nout]="Prob. Dist."; nout++;
				out[nout]="Correlation"; nout++;
				if(histoplot.length > 0){ out[nout]="Dependency"; nout++;}
				out[nout]="Trans. Dist."; nout++;
				out[nout]="Individuals"; nout++;
				out[nout]="Statistics"; nout++;
				out[nout]="Populations"; nout++;
				out[nout]="Pos. Pred. Ch."; nout++;
				if(infres.nderive > 0){ out[nout]="Derived"; nout++;}
			}
			addsubsub();
			for(j = 0; j < nout; j++) infpagename[j] = out[j];
		}
	}
	y += dby;

    addbutton("",5,20,0,0,-1,LOGOBUT,-1,-1);

	if(selectbub >= 0) bubble();
	
	lastbut = nbut;

	if(fin == 0) addbutton("Reload",width-205-4*90,0,75,22,RELOADAC,LOADBUT,-1,0); 
	
	if(page == HOMEPAGE){
		x = saveoptions(width-135);
		addbutton("Load",x,0,75,22,LOADBUT,LOADBUT,-1,0); 
		addbutton("[?]",width-55,0,15,20,HELPICONBUT,HELPICONBUT,10,-1);
	}
	else{
		if((page == INFERENCEPAGE && pagesub[page] == 2)|| (page == SIMULATEPAGE && pagesub[page] == 2)){
			if(pagesubsub[page][2] > 0){
				x = saveoptions(width-135);
				x = exportoptions(x);
				addbutton("[?]",width-55,0,15,20,HELPICONBUT,HELPICONBUT,11,-1);
			}
		}
		else{
			if(page == MODELPAGE && pagesub[page] >= 0) exportoptions(width-135-90); 
			addbutton("Import",width-135,0,75,22,IMPORTAC,LOADBUT,-1,4); 
			addbutton("[?]",width-55,0,15,20,HELPICONBUT,HELPICONBUT,12,-1);
		}
	}
	
	for(i = 0; i < gdropinfo.length; i++){
		var gd = gdropinfo[i];
		addbutton(gd.val,gd.x,gd.y,gd.dx,gd.dy,GDROPBUT,GDROPBUT,i,gd.style);
	}
	
	if(gdropsel >= 0){
		gdropbut = nbut;
		var gd = gdropinfo[gdropsel];

		gdropdy = gd.dy;
		gdropnum = gd.options.length;
		dy = (1+gdropnum)*gdropdy+8;
		if(gd.y+dy > height-10){
			gdropnum = Math.floor((height-10-gd.y-10)/gdropdy)-1; if(gdropnum < 1) gdropnum = 1; 
			dy = (1+gdropnum)*gdropdy+8;
		}
		addbutton(gd.val,gd.x,gd.y,gd.dx,dy,GDROPSELBUT,GDROPSELBUT,gd.options,gd.style);
	}
	
	if(helptype >= 0) helpbuts();

	for(i = nbut-1; i >= 0; i--){
		if(butac[i] >= 0 && mx >= butx[i] && mx <= butx[i]+butdx[i] && my >= buty[i] && my <= buty[i]+butdy[i]){ over = i; break;}
	}
	
	buttonplot();

	relo = 1;
}

function addsub(add)                                      // Used to add a sub-menu
{
	var i, xx = 15, dddy = dysub;
	
	for(i = 0; i < nout; i++){
		ii = i; if(add) ii += add;
		addbutton("",xx,y,menux-xx,dddy,-1,PAGESUBBACKBUT,-1,-1);
		addbutton(out[i],xx,y,menux-xx,dddy,PAGESUBBUT,PAGESUBBUT,ii,-1);
		y += dddy;
	}
}

function addsubsub()                                      // Used to add a sub-sub-menu
{
	var i, xx = 28; dddy = dysubsub;

	for(i = 0; i < nout; i++){
		addbutton(out[i],xx,y,menux-xx,dddy,PAGESUBSUBBUT,PAGESUBSUBBUT,i,-1); y += dddy;
	}
	y += 5;
}

function changepage(pagenew, s, ss)                        // Changes the page
{
	var te, pagest=page, addst, sss;
	
	clearTimeout(calcto);
	if(addtransmode == 1) deletetransmode();

	adddepfl = 0;
	ytotright = 1; advop = 0;
	playtime = 0; playing = 0;
	addingdata = 0; gentype = -1;			
	
	if(selectbub != -1){ buboff(); if(selectbub != -1) return;}
	plotinitfl = 0;
	
	addingdata = 0;
	if(pagenew != page) ById("add").innerHTML = "";
	xaxisauto = 1; yaxisauto = 1; 
	
	sss = pagesub[page]; if(sss < 0) sss += 10; pageyst[page][sss][pagesubsub[page][sss]] = tableyfr;
	
	if(pagenew != -1) page = pagenew;
	if(s != -1) pagesub[page] = s;
	if(ss != -1) pagesubsub[page][pagesub[page]] = ss;

	sss = pagesub[page]; if(sss < 0) sss += 10; tableyfr = pageyst[page][sss][pagesubsub[page][sss]];

	if(!(page == SIMULATEPAGE && pagesub[page] == 2 && pagesubsub[page][2] == 0) &&
	   !(page == INFERENCEPAGE && pagesub[page] == 2 && pagesubsub[page][2] == 0)) warning=[];
	
	if((page == SIMULATEPAGE && pagesub[page] == 2 && simpagename[pagesubsub[page][2]] == "Populations") || 
	   (page == INFERENCEPAGE && pagesub[page] == 2 && infpagename[pagesubsub[page][2]] == "Populations") ||
	   (page == INFERENCEPAGE && pagesub[page] == 2 && infpagename[pagesubsub[page][2]] == "Pos. Pred. Ch.")) plotinitfl = 0;
	
	if(page == INFERENCEPAGE && pagesub[page] == 2 && infpagename[pagesubsub[page][2]] == "Statistics"){ res = infres; startcalc();}
	if(page == SIMULATEPAGE && pagesub[page] == 2 && simpagename[pagesubsub[page][2]] == "Statistics"){ res = simres; startcalc();}
	
	if(page == INFERENCEPAGE && pagesub[page] == 2 && infpagename[pagesubsub[page][2]] == "Correlation"){ res = infres; res.varselx = -1; startcalccor();}
	
	if(page == SIMULATEPAGE && pagesub[page] == 2 && simpagename[pagesubsub[page][2]] == "Individuals") indplotst=[];
	if(page == INFERENCEPAGE && pagesub[page] == 2 && infpagename[pagesubsub[page][2]] == "Individuals") indplotst=[];
	
	buttoninit();
}

function hexToRgb(hex) {                                  // Changes colour format
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
        r: parseInt(result[1], 16),
        g: parseInt(result[2], 16),
        b: parseInt(result[3], 16)
    } : null;
}

function hex(c) {                                         // Converts decimal to hexidecimal
	var hex = (Math.floor(c)).toString(16);
	return hex.length == 1 ? "0" + hex : hex;
}

function darkcol(col)                                     // Darkens a colour
{ 
	if(col == BLACK) return GREY;
	var bigint, r, g, b, frac = 0.7;
	bigint = parseInt(col.substring(1), 16);	
	r = (bigint >> 16) & 255;
	g = (bigint >> 8) & 255;
	b = bigint & 255;
	return "#" + hex(frac*r) + hex(frac*g) + hex(frac*b);
}

function lightcol(col)                                    // Lightens a colour
{ 
	var bigint, r, g, b, frac = 0.13;
	bigint = parseInt(col.substring(1), 16);
	r = 255-(255-((bigint >> 16) & 255))*frac;
	g = 255-(255-((bigint >> 8) & 255))*frac;
	b = 255-(255-(bigint & 255))*frac;
	return "#" + hex(r) + hex(g) + hex(b);
}

function addbuttonlim(text,x,y,dx,dy,ac,type,val,val2)    // Limits the position of a button on right menu
{	
	if(y+dy > ybeg+rightheight || y < ybeg) return;
	addbutton(text,x,y,dx,dy,ac,type,val,val2);
}

function addbutton(text,x,y,dx,dy,ac,type,val,val2)       // Adds a button onto a page
{
	var i;
	if(type == CANVASBUT){ canx = x; cany = y;}
 
	if(y > height || y+dy < 0) return;
	if(type == TEXTWITHLINKBUT){
		var sp = text.split("here"), w1, w2;
		w1 = textwidth(sp[0],HELPFONT);
		w2 = textwidth("here",HELPFONT);
		addbutton(sp[0],x+5,y,0,0,0,TEXTWITHLINKBUT2,-1,-1);
		addbutton(sp[1],x+5+w1+w2,y,0,0,0,TEXTWITHLINKBUT2,-1,-1);
		text = "here"; x = x+5+w1; dx = w2;
	}
 
	if(type == CANVASBUT){ canvaswidth = dx; canvasheight = dy;}
 
	buttext[nbut] = text;
	butx[nbut] = x;
	buty[nbut] = y;
	butdx[nbut] = dx;
	butdy[nbut] = dy;
	butac[nbut] = ac; 	
	buttype[nbut] = type;
	butover[nbut] = -1;
	butval[nbut] = val;
	butval2[nbut] = val2;
	nbut++;
	
	if(type == TITLEBUT && val >= 0) addbutton("[?]",x+textwidth(text,HUGEFONT)+15,y,15,20,HELPICONBUT,HELPICONBUT,val,-1);
	
	if(type == SMALLTEXTBUT && val2 >= 0){
		addbutton("[?]",parseFloat(x)+textwidth(text,"12px arial")+2,y+3,15,15,HELPICONBUT,HELPICONBUT,val2,-1);		
	}
}

function splitintosub(st,font)                             // Splits a string into substrings to accout for subscripts
{
	var lab=[], j, font2, pos2, fsi, si;

	pos2 = font.search("px");
	si = parseInt(font.substr(pos2-2,2));
	fsi = Math.floor(0.7*si);
	font2 = font.substr(0,pos2-2)+fsi+font.substring(pos2);
		
	j = 0; w = 0;
	do{
		jst = j;
		while(j < st.length && st.substr(j,1) != "_") j++;
		frag = st.substr(jst,j-jst);
		dw = textwidth(frag,font)+1;
		lab.push({frag:frag, size:"big", w:w, dw:dw, font:font, dy:0});
		w += dw;
			
		if(j < st.length){
			j++; 
			jst = j;
			while(j < st.length && st.substr(j,1) != " " && st.substr(j,1) != "=" &&  st.substr(j,1) != "=" && st.substr(j,1) != "]" && st.substr(j,1) != "}") j++;
			frag = st.substr(jst,j-jst);
			if(frag != ""){
				dw = textwidth(frag,font2);
				lab.push({frag:frag, size:"small", w:w, dw:dw, font:font2, dy:0.2*si});
				w += dw;
			}
		}
	}while(j < st.length);
	
	return lab;
}

function plottext(text,x,y,font,col,width)                // Print left-aligned text
{
	if(text == undefined){ text = "undefined"; pr("undefined text");}
	if(text.length == 0) return;

	text = text.toString();
		
	if(text.search("_") == -1){   // looks for subscripts
		cv.font = font;
		if(width){
			if(cv.measureText(text).width > width-5){
				while(cv.measureText(text+"...").width > width-5) text = text.substr(0,text.length-1);
				text += "...";
			}	
		}
		cv.textAlign = 'left';
		cv.fillStyle = col;
		cv.fillText(text, x, y);
	}
	else{
		lab = splitintosub(text,font);
		
		cv.textAlign = 'left';
		cv.fillStyle = col;
		for(j = 0; j < lab.length; j++){
			cv.font = lab[j].font;
			if(lab[j].w+lab[j].dw > width-5){
				te = lab[j].frag;
				while(x+lab[j].w+cv.measureText(te+"...").width > width-5 && te.length > 0) te = te.substr(0,te.length-1);
				te += "...";
				cv.fillText(te,x+lab[j].w,y+lab[j].dy);
				break;
			}
			cv.fillText(lab[j].frag,x+lab[j].w,y+lab[j].dy);
		}
	}
}

function plotangletext(text,x,y,th,font,col,width)         // Draws text at an angle
{
	cv.font = font;
	if(width){
		if(cv.measureText(text).width > width-5){
			while(cv.measureText(text+"...").width > width-5) text = text.substr(0,text.length-1);
			text += "...";
		}	
	}
	
	cv.save();
	cv.translate(x, y);
	cv.rotate(-th);
	cv.textAlign = 'left';
	cv.fillStyle = col;
	cv.fillText(text, 0, 0);
	cv.restore();
}

function centerplotangletext(text,x,y,th,font,col,width)   // Draws centered text at an angle
{
	text = text.toString();
	if(text.search("_") == -1){   // looks for subscripts
		cv.font = font;
		if(width){
			if(cv.measureText(text).width > width-5){
				while(cv.measureText(text+"...").width > width-5) text = text.substr(0,text.length-1);
				text += "...";
			}	
		}
		
		cv.save();
		cv.translate(x, y);
		cv.rotate(-th);
		cv.textAlign = 'center';
		cv.fillStyle = col;
		cv.fillText(text, 0, 0);
		cv.restore();
	}
	else{
		lab = splitintosub(text,font);
		w = lab[lab.length-1].w+lab[lab.length-1].dw;
		cv.textAlign = 'left';
		cv.fillStyle = col;
		cv.save();
		cv.translate(x, y);
		cv.rotate(-th);	
		for(j = 0; j < lab.length; j++){
			cv.font = lab[j].font;
			cv.fillText(lab[j].frag, lab[j].w-w/2,lab[j].dy);
		}
		cv.restore();
	}
}

function centertext(text,x,y,font,col,width)               // Centers some text in an area
{
	if(text == undefined) text = "undefined";
	
	text = text.toString();
	if(text.search("_") == -1){   // looks for subscripts
		cv.font = font;

		if(width){
			if(cv.measureText(text).width > width-10){
				while(cv.measureText(text+"...").width > width-10) text = text.substr(0,text.length-1);
				text += "...";
			}	
		}
		
		cv.textAlign = 'center';
		cv.fillStyle = col;
		cv.fillText(text, x, y);
	}
	else{
		lab = splitintosub(text,font);
		w = lab[lab.length-1].w+lab[lab.length-1].dw;
		cv.textAlign = 'left';
		cv.fillStyle = col;
		for(j = 0; j < lab.length; j++){
			cv.font = lab[j].font;
			cv.fillText(lab[j].frag,x+lab[j].w-w/2,y+lab[j].dy);
		}
	}
}

function righttext(text,x,y,font,col,width)              // Right-aligns text
{
	cv.font = font;
	if(width){
		if(cv.measureText(text).width > width-5){
			while(cv.measureText(text+"...").width > width-5) text = text.substr(0,text.length-1);
			text += "...";
		}	
	}
	cv.textAlign = 'right';
	cv.fillStyle = col;
	cv.fillText(text, x, y);
}

function textwidth(text,font)                             // Returns the width of some text
{ 
	text = text.toString();
	if(text.search("_") == -1){ 
		cv.font = font;
		return cv.measureText(text).width;
	}
	else{
		var lab = splitintosub(text,font);
		var j = lab.length-1;
		return lab[j].w+lab[j].dw;
	}
}

function drawrect(x1, y1, x2, y2, col, style)             // Draws a rectangle
{
	cv.beginPath();
	cv.lineWidth = style;
	cv.rect(x1,y1,x2,y2);
	cv.strokeStyle = col;
	cv.stroke();
}

function reloadsign(x,y,col)                              // Draws reload symbol
{
	var dx = 20, dy = 20;
	cv.lineWidth = 3; 
	cv.beginPath();
	cv.arc(x+dx/2, y+dy/2+1, 8, -0.1, 1.5*Math.PI);
	cv.strokeStyle = col;
	cv.stroke(); 
	drawarrow(x+dx/2+7,y+3,x-100,y+6,8,col);
	centertext("Reload",x+dx/2,y+dy+9,"bold 10px georgia",col);
}

function fillrect(x1, y1, x2, y2, col)                    // Draws a filled rectangle
{
	cv.beginPath();
	cv.rect(x1,y1,x2,y2);
	cv.fillStyle = col;
	cv.fill();
}

function drawgradientline(x1, y1, x2, y2, f1, f2)         // Draws a line with a colour gradient
{
	var grad= cv.createLinearGradient(x1,y1,x2,y2);
	v = Math.floor(f1*225); grad.addColorStop(0, "rgb("+v+","+v+","+v+")");
	v = Math.floor(f2*225); grad.addColorStop(1, "rgb("+v+","+v+","+v+")");

	cv.lineWidth = 4; 
	cv.beginPath();
	cv.moveTo(x1,y1);
	cv.lineTo(x2,y2);
	cv.strokeStyle = grad;
	cv.stroke();
}

function setdash(dash)                                    // Sets the dash type
{	
	var sm = 3*scale, md = 9*scale, bi = 14*scale;
	
	switch(dash%10){
	case 0: cv.setLineDash([]); break;
	case 1: cv.setLineDash([md, sm]); break;
	case 2: cv.setLineDash([sm, sm]); break;
	case 3: cv.setLineDash([md, sm, sm, sm]); break;
	case 4: cv.setLineDash([bi, md]); break;
	case 5: cv.setLineDash([bi, sm, sm, sm]); break;
	case 6: cv.setLineDash([bi, sm, md, sm]); break;
	case 7: cv.setLineDash([bi, sm, sm, sm, sm, sm]); break;
	case 8: cv.setLineDash([bi, sm, md, sm, sm, sm]); break;
	case 9: cv.setLineDash([bi, sm, md, sm, md, sm]); break;
	}
}

function drawline(x1, y1, x2, y2, col, style, dash)       // Draws a line
{
	cv.lineWidth = style;

	if(dash) setdash(dash);

	cv.beginPath();
	cv.moveTo(x1,y1);
	cv.lineTo(x2,y2);
	cv.strokeStyle = col;
	cv.stroke();
	
	if(dash) cv.setLineDash([]);
}

function drawpolygon(npoint,col,col2, style)             // Draws a polygon
{
	var i;

	cv.lineWidth = style;

	cv.beginPath();
	cv.moveTo(polypoint[0][0],polypoint[0][1]);
	for(i = 1; i < npoint; i++) cv.lineTo(polypoint[i][0],polypoint[i][1]);
	cv.closePath();
	cv.fillStyle = col;
	cv.fill();
	cv.strokeStyle = col2;
	cv.stroke();
}

function drawroundrect(x,y,dx,dy,r,col,col2)            // Draws a rectangle with rounded corners
{
	var th, i, nth = Math.floor(r/2);
	if(nth < 1){
		cv.lineWidth = 1;	
		cv.beginPath();
		cv.beginPath();
		cv.moveTo(x+r,y);
		cv.lineTo(x+dx-r,y);
		cv.lineTo(x+dx,y+r);
		cv.lineTo(x+dx,y+dy-r);
		cv.lineTo(x+dx-r,y+dy);
		cv.lineTo(x+r,y+dy);
		cv.lineTo(x,y+dy-r);
		cv.lineTo(x,y+r);
		cv.closePath();
		cv.fillStyle = col; 
		cv.fill();
		cv.strokeStyle = col2;
		cv.stroke();
		return;
	}
	
	cv.lineWidth = 1;
	cv.beginPath();
 
	cv.moveTo(x+r,y);

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	}

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}

	cv.closePath();
	cv.fillStyle = col; 
	cv.fill();
	cv.strokeStyle = col2;
	cv.stroke();
}

function drawloadrect(x,y,dx,dy,r,col,col2)              // Draws rectangle rounded at bottom
{
	var th, i, nth = Math.floor(r/2);
	if(nth < 1) nth = 1;

	cv.lineWidth = 1;
	cv.beginPath();
 
	cv.moveTo(x,y);
	cv.lineTo(x+dx,y);

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}
	cv.lineTo(x,y);

	cv.closePath();
	cv.fillStyle = col; 
	cv.fill();
	cv.strokeStyle = col2;
	cv.stroke();
}

function drawroundmenutop(x,y,dx,dy,r,col,col2)          // Draws rounded rectangle with title
{
	var th, i, nth = Math.floor(r/2);
	if(nth < 1) nth = 1;

	cv.lineWidth = 1;
	cv.beginPath();
 
	cv.moveTo(x+r,y);

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	}
	cv.lineTo(x+dx,y+dy);
	cv.lineTo(x,y+dy);
	
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}

	cv.closePath();
	cv.fillStyle = col; 
	cv.fill();
	cv.strokeStyle = col2;
	cv.stroke();
}

function circle(x,y,r,col,style)                         // Draws circle
{
	cv.lineWidth = style;
	cv.beginPath();
	cv.arc(x,y,r,0,2*Math.PI);
	cv.strokeStyle = col;
	cv.stroke();
}

function fillcircle(x,y,r,col,col2,style)                // Draws filled circle
{
	cv.lineWidth = style;
	cv.beginPath();
	cv.arc(x,y,r,0,2*Math.PI);
	cv.fillStyle = col;
	cv.fill();
  
	cv.strokeStyle = col2;
	cv.stroke();
}

function drawcorners(x,y,dx,dy,r,col,col2)               // Draws the corners of the main page
{
	var th, i, nth = Math.floor(r/3);
	if(nth < 1) nth = 1;

	cv.lineWidth = 2;
	cv.beginPath();
	cv.fillStyle = col; 

	cv.moveTo(x+dx,y);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	} 
	cv.closePath();
	cv.fill();

	cv.moveTo(x+dx,y+dy);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();

	cv.moveTo(x,y+dy);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}
	cv.closePath();
	cv.fill();

	cv.moveTo(x,y);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();
}

function alignparagraph(text,dx,fontset)                  // Aligns lines of a paragraph within a certain width
{
	var i, ist, len, yy, font, di;

	dx -= 8;
	if(fontset){ font = fontset; lineheight = 25;} else{ font = HELPFONT; lineheight = 23;}

	text = ""+text;
	nlines = 0;
	i = 0; len = text.length;

	di = Math.floor(len*dx/textwidth(text,font));
	
	yy = 0;
	while(i < len){
		ist = i;
		i += di; if(i > len) i = len; while(i < len && text.substr(i,1) != " ") i++;
		 
		if(textwidth(text.substr(ist,i-ist),font) > dx){
			do{
				i--; while(text.substr(i,1) != " " && i > ist) i--;
				
			}while(textwidth(text.substr(ist,i-ist),font) > dx);
			i++;
		}
		else{
			do{
				ist2 = i;
				i++;
				while(i < len && text.substr(i,1) != " ") i++;
			}while(i < len && textwidth(text.substr(ist,i-ist),font) < dx);
			i = ist2;
		}
		
		lines[nlines] = text.substr(ist,i-ist);
		while(i < len && text.substr(i,1) == " ") i++;
		linesy[nlines] = yy;
		yy += lineheight;
		nlines++;
	}

	hsto = yy;
}

function drawarrow(x,y,x2,y2,size,col)                    // Draws an arrow
{
	var nx, ny, px, py, r;

	nx = x2-x; ny = y2-y;
	r = Math.sqrt(nx*nx+ny*ny);
	if(size > r/5) size = r/5;
	nx *= size/r; ny *= size/r; 
	px = 0.5*ny; py = -0.5*nx;

	polypoint[0][0] = x; polypoint[0][1] = y; 
	polypoint[1][0] = x+Math.round(nx*0.8); polypoint[1][1] = y+Math.round(ny*0.8);
	polypoint[2][0] = x+Math.round(nx+px); polypoint[2][1] = y+Math.round(ny+py);
	drawpolygon(3,col,col,NORMLINE);
	polypoint[2][0] = x+Math.round(nx-px); polypoint[2][1] = y+Math.round(ny-py);
	drawpolygon(3,col,col,NORMLINE);
}

function plotlabel(text,x,y,font,col,width)               // Plots the variables on right menu
{
	plottext(text,x,y,font,col,width);
}

function plotxlabel(text,x,y,font,col)                    // Plots the x label on a graph
{
	centertext(text,x,y,font,col); 
}

function plotylabel(text,x,y,font,col)                    // Plots the y label on a graph
{
	centerplotangletext(text,x,y,Math.PI/2,font,col);
}

function overlaytext(text,x,y)                            // Plots text on the overlay layer
{
	cvover.font = "bold 12px georgia";

	overx = x; overy = y;
	overdx = cvover.measureText(text).width+10;
	overdy = 18;
	cvover.beginPath();
	cvover.rect(overx,overy,overdx,overdy);
	cvover.fillStyle = DDBLUE;
	cvover.fill();
	cvover.strokeStyle = BLACK;
	cvover.stroke();
	
	cvover.textAlign = 'left';
	cvover.fillStyle = WHITE;
	cvover.fillText(text,x+5,y+13);
}

function startloading()                                    // Starts the loading symbol
{
	if(loading == 0){ loading = 1; percent = 0; loadoff = 0; loadingfunc(1);}
}

function loadingfunc(ite)                                  // Animates the loading symbol
{
	var x, y, r, rr, i;
	x = 0.5*(menux+width); y = height/2; r = 20; rr = 4;
	cvover.clearRect(x-(r+rr+1), y-(r+rr+1),2*(r+rr+1),2*(r+rr+1));
	if(loading == 1){
		if(ById("pw").offsetLeft != ById("overlay").offsetLeft) ById("overlay").style.left = ById("pw").offsetLeft;
	
		loadoff++;
		
		cvover.lineWidth = 1;
		cvover.fillStyle = RED; 
		for(i = 0; i < 12; i++){
			th = (i-loadoff)*2*Math.PI/12.0;
			cvover.globalAlpha = Math.exp(-i*0.5);
			cvover.beginPath();
			cvover.arc(Math.floor(x+r*Math.sin(th)),Math.floor(y+r*Math.cos(th)),rr,0,2*Math.PI);
			cvover.fill();
		}
		cvover.globalAlpha = 1;

        cvover.font = "14px arial";
		cvover.textAlign = 'center';
		cvover.fillStyle = RED;
		cvover.fillText(percent,x,y+4);
		if(ite == 1) setTimeout(function(){ loadingfunc(1)},110);
	}
}

function startrunning()                                   // Starts the running animation around logo
{
	if(running == 0){ running = 1; runningfunc(1);}
}

function runningfunc(ite)                                 // Animates the running animation around logo
{
	var x, y, dx, dy, fr, rr;
	
	x = 3; y = 29; dx = 128; dy = 48; l = 2*(dx+dy); rr = 3;
	cvover.clearRect(x-5, y-5,dx+10,dy+10);
	if(running == 1){
		if(ById("pw").offsetLeft != ById("overlay").offsetLeft) ById("overlay").style.left = ById("pw").offsetLeft+"px";
	
		loadoff++;
		
		cvover.lineWidth = 1;
		cvover.fillStyle = BLUE; 
		for(i = 0; i < 12; i++){
			fr = ((12*10000-i+loadoff)/12.0)/3; fr = fr-Math.floor(fr)
			fr *= l;
			if(fr < dx){ xx = x+fr; yy = y;}
			else{
				fr -= dx;
				if(fr < dy){ xx = x+dx; yy = y+fr;}
				else{
					fr -= dy;
					if(fr < dx){ xx = x+dx-fr; yy = y+dy;}
					else{
						fr -= dx;
						xx = x; yy = y+dy-fr;
					}
				}
			}
			cvover.globalAlpha = Math.exp(-i*0.5);
			cvover.beginPath();
			cvover.arc(Math.floor(xx),Math.floor(yy),rr,0,2*Math.PI);
			cvover.fill();
		}
		cvover.globalAlpha = 1;

     	setTimeout(function(){ runningfunc(1)},80);
	}
}

function indbuttonplot(i)                                 // Plots an individual button
{
	cv.clearRect (butx[i],buty[i],butdx[i],butdy[i]);
	butplot(i,0);
}

function gettabcol(val,dy)                                // Gets the colour for a column in a table
{
	if(val == -1){ col = BLACK; return;}
	
	if(page == INFERENCEPAGE && pagesub[page] == 0 && addingdata == 0){ // Data sources page
		if(dy == 0) col = BLACK; else col = BLUE;  
		col2 = LLBLUE; 
		return;
	}

	if(val < ncoldef){
		switch(val){
		case 0: col = BLUE; col2 = LLBLUE; break;
		case 1: col = DBLUE; col2 = LLBLUE; break;
		default: col = DDBLUE; col2 = LLBLUE; break;
		}
	}
	else{ col = LLGREY; col2 = LRED;}
}
