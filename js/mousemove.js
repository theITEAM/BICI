function mousemove(x,y)                                   // Fires when the mouse moves
{
	var overnew = -1;

	dragged = 1;
	mx = x; my = y;

	switch(drag){    // Options for dragging different elements
	case 1:  // xslide
		tablexfr = tablexfrst + (mx-mxst)/canvaswidth;
		if(tablexfr > 1-tablexfrac) tablexfr = 1-tablexfrac; if(tablexfr < 0) tablexfr = 0;
		break;
	
	case 2:  // yslide
		tableyfr = tableyfrst + (my-myst)/canvasheight;
		if(tableyfr > 1-tableyfrac) tableyfr = 1-tableyfrac; if(tableyfr < 0) tableyfr = 0;
		break;
		
	case 3:
		cl = dragval;
		togrid(cl,mxlast-mxst,mylast-myst); gxi = gx; gyi = gy;
		togrid(cl,mx-mxst,my-myst); 
		mxlast = mx; mylast = my;
	
		cla[cl].comp[dragval2].x += gx-gxi;
		cla[cl].comp[dragval2].y += gy-gyi;
		var t;
		for(t = 0; t < cla[cl].ntra; t++){
			tr = cla[cl].tra[t];
			if(tr.i == dragval2 || tr.f == dragval2) findpline(cl,t);
		}
		break;
		
	case 4:
		cl = pagesub[MODELPAGE];
		togrid(cl,mxlast-mxst,mylast-myst); gxi = gx; gyi = gy;
		togrid(cl,mx-mxst,my-myst); 
		mxlast = mx; mylast = my;
		zoomin(cl,0,0,1,gx-gxi,gy-gyi);
		break;
	
	case 11: // drags indidual time lines
		tableyfr -= (my-mylast)/ytot;
		if(tableyfr > 1-tableyfrac) tableyfr = 1-tableyfrac; if(tableyfr < 0) tableyfr = 0;
		mylast = my;
		break;
		
	case 5:
		cl = pagesub[MODELPAGE];
		togrid(cl,mxlast-mxst,mylast-myst); gxi = gx; gyi = gy;
		togrid(cl,mx-mxst,my-myst); 
		mxlast = mx; mylast = my;
		
		cla[cl].tra[dragval].p[dragval2].x += gx-gxi;
		cla[cl].tra[dragval].p[dragval2].y += gy-gyi;
		findpline(cl,dragval);
		break;
		
	case 6:
		if((mx-mxst)*(mx-mxst) + (my-myst)*(my-myst) > 25 && page == MODELPAGE && modelsetup == 0) breaktrans();
		break;
		
	case 9:
		if(mx-mxst < colx[dragval] && dragval > 0){  // switches columns
			movecol(dragval,dragval-1);
			setcolumns();
			dragval--;
		}
		else{
			if(obsloaded == 1) imax = ncol; else imax = ncoldef;	
			
			if(mx-mxst > colx[dragval] + colw[dragval] + 10 && dragval < imax-1){  // switches columns
				movecol(dragval,dragval+1);
				setcolumns();
				dragval++;
			}
		}
		break;
	
	case 10:
		kdeold = kde;
		kde = kdest + (mx-mxst)*(kdemax-kdemin)/(slidex-slidedx);
		if(kde > kdemax) kde = kdemax; if(kde < kdemin) kde = kdemin;
		if(kde != kdeold) plotinitfl = 0;
		break;
	}
 
	if(drag >= 1 && drag != 101) buttoninit();
	
	if(exporton == 1){
		i = 0; while(i < nbut && !(buttype[i] == LOADBUT && butval2[i] == 3)) i++; if(i == nbut) alertp("Error code EC42");
		if(mx < butx[i] || mx > butx[i]+butdx[i] || my > buty[i]+butdy[i]){ exporton = 0; buttoninit();}
	}
	
	if(saveon == 1){
		i = 0; while(i < nbut && !(buttype[i] == LOADBUT && butval2[i] == 1)) i++; if(i == nbut) alertp("Error code EC43");
		if(mx < butx[i] || mx > butx[i]+butdx[i] || my > buty[i]+butdy[i]){ saveon = 0; buttoninit();}
	}
	
	canbut = -1;
	for(i = nbut-1; i >= 0; i--){
		if(buttype[i] == CANVASBUT) canbut = i;
		if(butac[i] >= 0 && mx >= butx[i] && mx <= butx[i]+butdx[i] && my >= buty[i] && my <= buty[i]+butdy[i]){
			overnew = i; break;
		}
	}

	if(canbut != -1){
		canovernew = -1;
		if(overnew == over && over == canbut){
			xx =  mx-butx[over]; yy = my-buty[over];
			canmx = xx; canmy = yy;
			for(i = ncanbut-1; i >= 0; i--){
				if(canbuttype[i] == TRANSPLBUT){
					if(transover(canbutval[i],canbutval2[i],canmx,canmy) == 1){ canovernew = i; break;}
				}
				else{
					if(canbuttype[i] != MODELPICBUT && canbutac[i] >= 0 && xx >= canbutx[i] && xx < canbutx[i]+canbutdx[i] &&
										   yy >= canbuty[i] && yy < canbuty[i]+canbutdy[i]){
						canovernew = i;
						break;
					}
				}
			}
		}
	}	
	
	if(canbut != -1){
		if((canovernew != canover && over != -1 && buttype[over] == CANVASBUT) || addtransmode == 1){
			canover = canovernew;
			butplot(canbut,-1);
		}

		if(addcompmode == 1 || addsourcemode == 1){ buttoninit();}
		
		if(overnew != canbut && addtransmode == 1){ deletetransmode(); buttoninit();}
	}

	if(gdropsel >= 0){
		butplot(gdropbut,1);
		if(buttype[overnew] != GDROPSELBUT && buttype[overnew] != GDROPBUT && drag == 0){ gdropsel = -1; buttoninit(); }
	}
			
	if(canover >= 0 && canbuttype[canover] == SLIDERBUT) butplot(canover,-1);
	
	if(overnew != over || buttype[over] == XSLIDEBUT || buttype[over] == YSLIDEBUT){
		if(overx != -1){ cvover.clearRect(overx-1,overy-1,overdx+2,overdy+2); overx = -1;}
	
		overold = over; over = overnew;
		
		if(overold >= 0 && !(buttype[overold] == LOADBUT && exporton > 0) && buttype[overold] != HELPBACKBUT){
			butplot(overold,-1);
		}
		if(overnew >= 0 && !(buttype[overnew] == LOADBUT && exporton > 0) && buttype[overnew] != HELPBACKBUT){
			butplot(overnew,-1);
		}
	}

	arrownew = 0;
	if(addingdata == 0 && page == MODELPAGE && pagesub[page] >= 0 && pagesub[page] < ncla-2 && over == canbut && canover == -1 && addtransmode == 0 &&
       addsourcemode == 0 && addcompmode == 0 && ncla > 0){
		arrownew = 1;
	}
	
	if(over == canbut && canover == -1 && 
		((page == SIMULATEPAGE && pagesub[page] == 2 && simpagename[pagesubsub[page][2]] == "Individuals" && res.indview == "Timeline")
				|| (page == INFERENCEPAGE && pagesub[page] == 0 && pagesubsub[page][0] == 1)
				|| (page == INFERENCEPAGE && pagesub[page] == 2 && infpagename[pagesubsub[page][2]] == "Individuals" && res.indview == "Timeline"))){
		xq = mx-menux-tlinexmin-40; yq = my-40;
		if(xq > 0 && xq < tablewidth && yq >0 && yq < indtableheight) arrownew = 2;
	}
	
	if(drag == 8 || drag == 9){
		arrownew = 1;
	}	
	
	if(arrownew != arrow){           // Sets the cursor (i.e. a hand is used for grabbing)
		arrow = arrownew;
		if(arrow > 0) document.getElementById("bod").style.cursor = "-webkit-grab"; 
		else document.getElementById("bod").style.cursor = "";
	}
	
	if(exporton == 1){ if(over < 0 || (buttext[over] != "Export" && buttype[over] != EXPORTMINIBUT)){ exporton = 0; buttoninit();}}
}

