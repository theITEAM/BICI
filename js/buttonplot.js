function buttonplot()                                     // Plots all the buttons on the current page
{
	var i, ov; 
	allbutplot = 1;

	cv = maincv;
	cv.clearRect(0,0,canw,canh);
	for(i = 0; i < nbut; i++){
		ov = 0; if(i == over) ov = 1;
		butplot(i,ov);
	}
	allbutplot = 0;
}

function butplot(i,ov)                                    // Plots individual button
{ 
	var x, y, dx, dy, col, col2, text, off, w, val, val2, xx, yy, dd;

	x = Math.floor(butx[i]); y = Math.floor(buty[i]); dx = Math.floor(butdx[i]); dy = Math.floor(butdy[i]); text = buttext[i];
	val = butval[i]; val2 = butval2[i];

	switch(buttype[i]){ 
	case TABBUT:
		arrowdown = 0;
	    if(butval[i] == page){
			if(page != HOMEPAGE && page != DESCPAGE) arrowdown = 1;
			
			fillrect(x,y,dx,dy,WHITE); 
			col = "#003366";
			font = MENUFONT;
			off = 0;
		}
		else{
			fillrect(x,y,dx,dy,"#ddddff"); 
			font= MENUFONT2;
			off = -1;
			col = "#3366ff";
			if(over==i) col = "#ff6666";
			if(modelstart == 0)  col = "#aaaaff";
			else{
				if((modelsetup == -1 && val > 1) || (modelsetup == 0 && val > 2)) col = "#aaaaff";
			}
		} 
		if(pagesize == "big") off += 3;
		
		plottext(text,x+28,y+24+off,font,col);
		
		if(arrowdown == 0){
			ddx = 15; ddx2 = 5; ddy = 5;
			drawline(x+ddx,y+dy/2-ddy,x+ddx+ddx2,y+dy/2,col,THICKLINE)
			drawline(x+ddx,y+dy/2+ddy,x+ddx+ddx2,y+dy/2,col,THICKLINE)
		}
		else{
			ddx = 15; ddx2 = 4; ddy = 4;
			drawline(x+ddx-ddx2,y+dy/2-ddy+1,x+ddx,y+dy/2+ddy+1,col,THICKLINE);
			drawline(x+ddx+ddx2,y+dy/2-ddy+1,x+ddx,y+dy/2+ddy+1,col,THICKLINE);
		}
		break;
		
	case PAGESUBBACKBUT:
		fillrect(x, y, dx, dy, "#ededff");
		break;
		
	case PAGESUBBUT:
		arrowdown = 0;
		if(butval[i] == pagesub[page]){
			if(!(page == SIMULATEPAGE && val == 2) && !(page == INFERENCEPAGE && val == 2) && !(page == INFERENCEPAGE && val == 0 && data.length > 0)) fillrect(x+dx-8,y+dy/2-6,8,12,WHITE); 
			else arrowdown = 1;
			
			if(text == "Time" || text == "Age") col = "#2222cc";
			else{
				if(text == "Derived" || text == "Distributions") col = "#444444";
				else{
					col = "#cc2222"; if(page == MODELPAGE && over==i) col = "#dd8888";
				}
			}
		
			font = MENUFONTSMALL;
			off = 0;
		}
		else{
			font = MENUFONTSMALL2;
			off = -1;
			if(text == "Time" || text == "Age"){ col = "#7777ff"; if(over==i) col = "#333377";}
			else{
				if(text == "Derived" || text == "Distributions"){ col = "#444444"; if(over==i) col = "#999999";}
				else{ col = "#ff2222"; if(over==i) col = "#dd8888";}
			}
		} 

		fillrect(x,y,26+textwidth(text,font), dy, "#ededff");
		//if(pagesize == "big")
		off += 1;
		plottext(text,x+24,y+20+off-4,font,col);
		if(arrowdown == 0){
			ddx = 11; ddx2 = 4; ddy = 4;
			drawline(x+ddx,y+dy/2-ddy,x+ddx+ddx2,y+dy/2+1,col,THICKLINE);
			drawline(x+ddx,y+dy/2+ddy,x+ddx+ddx2,y+dy/2+1,col,THICKLINE);
		}
		else{
			ddx = 11; ddx2 = 4; ddy = 3;
			drawline(x+ddx-ddx2,y+dy/2-ddy,x+ddx,y+dy/2+ddy+2,col,THICKLINE);
			drawline(x+ddx+ddx2,y+dy/2-ddy,x+ddx,y+dy/2+ddy+2,col,THICKLINE);
		}
		break;

	case PAGESUBSUBBUT:
		fillrect(x, y, dx, dy, "#ddddff");
		if(butval[i] == pagesubsub[page][pagesub[page]]){
			fillrect(x+dx-8,y+dy/2-6,8,12,WHITE); 
			col = "#222222"; if(over==i && modelsetup == 0) col = "#888888";
			font = MENUFONTSMALL;
			off = 0;
		}
		else{
			font = MENUFONTSMALL2;
			off = -1;
			col = "#222222";
			if(over==i) col = "#888888";
		} 
		if(pagesize == "big") off += 2;
		
		plottext(text,x+8,y+16+off-2,font,col);
		break;
		
	case ADDCLASSBUT:
		if(val == -1){ col = "#ff2222"; if(over==i) col = "#dd8888";}
		else{ col = "#2222cc"; if(over==i) col = "#8888dd";} 
		
		fillrect(x,y,26+textwidth(text,font), dy, "#ededff");	
		x = butx[i]+18; y = buty[i]+dy/2; r = 6; r2 = 4;
		drawroundrect(x-r,y-r,2*r,2*r,r,col,col);
		drawline(x-r2,y,x+r2,y,WHITE,THICKLINE);
		drawline(x,y-r2,x,y+r2,WHITE,THICKLINE);
		off = 0; if(pagesize == "big") off += 2;
		plottext(text,x+11,buty[i]+16+off,MENUFONTSMALL,col);
		break;
		
	case UPLOADBUT:
		col = WHITE; col2 = RED; col3 = RED;  if(i == over){ col = LLRED;}
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2,y+21,INPUTFONT,col3);
		break;
		
	case GREENBUT:
		col = WHITE; col2 = DGREEN; col3 = DGREEN;  if(i == over){ col = LLGREEN;}
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2,y+21,INPUTFONT,col3);
		break;
		
	case ADDDATABUT:
		col = WHITE; col2 = BLUE; col3 = BLUE;  if(i == over){ col = LLBLUE;}
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2,y+21,INPUTFONT,col3);
		break;
		
	case CANCELBUT2:
		col = WHITE; col2 = DGREEN; col3 = DGREEN; if(i == over){ col = LLGREEN;}
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2,y+21,INPUTFONT,col3);
		break;
	
	case NEXTBUT:
		col = WHITE; col2 = RED; col3 = RED;  if(i == over){ col = LLRED;}
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2-5,y+21,INPUTFONT,col3);
		drawline(x+dx-8,y+dy/2,x+dx-14,y+5,col3,THICKLINE);
		drawline(x+dx-8,y+dy/2,x+dx-14,y+dy-5,col3,THICKLINE);
		drawline(x+dx-8-6,y+dy/2,x+dx-14-6,y+5,col3,THICKLINE);
		drawline(x+dx-8-6,y+dy/2,x+dx-14-6,y+dy-5,col3,THICKLINE);
		break;
		
	case BACKBUT:
		col = WHITE; col2 = RED; col3 = RED;  if(i == over){ col = LLRED;}
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2+5,y+21,INPUTFONT,col3);
		drawline(x+8,y+dy/2,x+14,y+5,col3,THICKLINE);
		drawline(x+8,y+dy/2,x+14,y+dy-5,col3,THICKLINE);
		drawline(x+8+6,y+dy/2,x+14+6,y+5,col3,THICKLINE);
		drawline(x+8+6,y+dy/2,x+14+6,y+dy-5,col3,THICKLINE);
		break;
		
	case SMBACKBUT:
		col = WHITE; col2 = BLUE; col3 = BLUE;  if(i == over){ col = LLBLUE;}
		drawroundrect(x,y,dx,dy,8,col,col2);
		centertext(text,x+dx/2+5,y+19,INPUTFONT,col3);
		drawline(x+8,y+dy/2,x+14,y+5,col3,THICKLINE);
		drawline(x+8,y+dy/2,x+14,y+dy-5,col3,THICKLINE);
		drawline(x+8+6,y+dy/2,x+14+6,y+5,col3,THICKLINE);
		drawline(x+8+6,y+dy/2,x+14+6,y+dy-5,col3,THICKLINE);
		break;
	
	case NTEXTBUT:
		plottext(text,x+4,y+18,"17px arial",DDGREY);
		break;
	
	case BULLETBUT:
		drawroundrect(x+5,y+8,6,8,3,BLACK,BLACK);
		break;
		
	case TITLEBUT:
		plottext(text,x+15,y+19,HUGEFONT,DBLUE);
		dx = width-(x-menux)-10; 
		drawline(x+12,y+25,dx,y+25,LBLUE,NORMLINE);
		break; 
	
	case BOLDMEDBUT:
		plottext(text,x+2,y+18,"bold 17px arial",val);
		break;
		
	case SMALLTEXTBUT:
		plottext(text,x+2,y+15,"12px arial",val);
		break;
		
	case TEXTSUBTITLEBUT:
		col = DDGREY; if(over == i) col = GREY;
		plottext(text,x,y+14,"bold 16px arial",col);
		break;
		
	case TEXTWITHLINKBUT:
		fillrect(x,y,dx,dy,WHITE);
		col = BLUE; if(over == i) col = LBLUE;
		plottext("here",x,y+18,HELPFONT,col);
		break;
	
	case TEXTWITHLINKBUT2:
		plottext(text,x,y+18,HELPFONT,DDBLUE);
		break;
		
	case SMTEXTBUT:
		plottext(text,x+5,y+17,"14px arial",DDBLUE);
		break;
		 
	case PARAGRAPHBUT:
		alignparagraph(text,dx);
		for(j = 0; j < nlines; j++) plottext(lines[j],x+5,y+18+linesy[j],HELPFONT,DDBLUE);
		break;
		
	case MENUBACKBUT:
		fillrect(x,y,dx,dy,"#ddddff");
		drawcorners(menux,0,canw-menux,dy,30,"#ddddff","#9999ff")
		break;
		
	case LOGOBUT:
		cv.drawImage(logopic,x,y);
		break;
		
	case NEWMODBUT:
		cv.clearRect(x,y,dx,dy);
		if(over == i) cv.globalAlpha = 0.3;
		cv.drawImage(newmodpic,x,y,dx,dy);
		if(over == i) cv.globalAlpha = 1;
		break;
	
	case PDFBUT:
		cv.clearRect(x,y,dx,dy);
		if(over == i) cv.globalAlpha = 0.3;
		cv.drawImage(manualpic,x,y,dx,dy);
		if(over == i) cv.globalAlpha = 1;
		break;
	
	case PDFBUT2:
		cv.clearRect(x,y,dx,dy);
		if(over == i) cv.globalAlpha = 0.3;
		cv.drawImage(paperpic,x,y,dx,dy);
		if(over == i) cv.globalAlpha = 1;
		break;
	
	case WARNBUT:
		cv.drawImage(warnpic,x,y);
		break;
	
	case CANVASBUT:
		canbuttonplot();
		fillrect(x,y,dx,dy,WHITE); 
		//fillrect(x,y,dx,dy,LRED); 
	
		cv.drawImage(graphcan,0,0,dx,dy,x,y,dx,dy);
		if(dy == height) drawcorners(menux,0,canw-menux,dy,30,"#ddddff","#9999ff");
		canvasdx = dx; canvasdy = dy;
		
		for(bb = lastbut; bb < nbut; bb++){
			ov = 0; if(bb == over) ov = 1;
			butplot(bb,ov);
		}
		break;
	 
	case XSLIDEBUT:
		if(tablexfrac < 1 && tablexfr+tablexfrac > 1){ tablexfr = 1- tablexfrac; buttoninit();}
		
		fillrect(x-1,y-1,dx+2,dy+2,WHITE);
		x1 = x+tablexfr*dx; x2 = x+dx*(tablexfr+tablexfrac);
		col = GREY; if((over == i && mx >= x1 && mx <= x2) || drag == 1) col = DDGREY;
		fillrect(x1,y,x2-x1,dy,col);
		drawrect(x,y,dx,dy,DDGREY,NORMLINE);
		break;
	 
	 case YSLIDEBUT:
		if(tableyfrac < 1 && tableyfr+tableyfrac > 1){ tableyfr = 1-tableyfrac; buttoninit();}
		
		fillrect(x-1,y-1,dx+2,dy+2,WHITE);
		var dry = dy*tableyfrac, pad = 0;
		if(dry < 50){ pad = (50-dry)/2; dry = 50;}
		sliy1 = y + pad + (tableyfr+tableyfrac/2)*(dy-2*pad) - dry/2;
		sliy2 = y + pad + (tableyfr+tableyfrac/2)*(dy-2*pad) + dry/2;
		
		coheight = dy-2*pad;
		col = GREY; if((over == i && my >= sliy1 && my <= sliy2) || drag == 2) col = DDGREY;
		drawrect(x,y,dx,dy,DDGREY,NORMLINE);
		fillrect(x,sliy1,dx,sliy2-sliy1,col);
		break;
	
	case WHITEBUT:
		fillrect(x,y,dx,dy,WHITE);
		break;
		
	case LOADBUT:
		switch(val2){
		case 0: col = RED; col2 = DRED; if(over == i){ col = DRED; col2 = DDRED;}; break;
		case 1: col = BLUE; col2 = DBLUE; if(over == i && modelstart == 1 && saveon == 0){ col = DBLUE; cl2 = DDBLUE;}; break;
		case 2: col = BLUE; col2 = DBLUE; if(over == i){ col = DBLUE; cl2 = DDBLUE;}; break;
		case 3: col = LORANGE; col2 = ORANGE; if(over == i && exporton == 0){ col = ORANGE; cl2 = DORANGE;}; break;
		case 4: col = LPURPLE; col2 = PURPLE; if(over == i){ col = PURPLE; col2 = DPURPLE;} break;  	
		}
		drawloadrect(x,y,dx,dy,7,col,col2);
		centertext(text,x+dx/2,y+16,MENUFONT,WHITE);
		if(val2 == 1 && modelstart == 0){
			cv.globalAlpha = 0.9;
			fillrect(x-1,y-1,dx+2,dy+2,WHITE);
			cv.globalAlpha = 1;
		}
		break;
		
	case CHECKBUT:
		fillrect(x,y,dx,dy,WHITE);
		drawrect(x+2,y+2,dy-4,dy-4,BLACK,NORMLINE);

		if(res.popshow[val][val2] == 1){
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,GREY);
			else fillrect(x+4,y+4,dy-8,dy-8,BLACK);
		}
		else{
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,LLGREY);
		}
		
		plottext(text,x+dy+4,y+14,"bold 14px arial",BLACK);
		break;
		
	case CHECKBUT2:
		fillrect(x,y,dx,dy,WHITE);
		drawrect(x+2,y+2,dy-4,dy-4,BLACK,NORMLINE);

		if(indshow[val] == 1){
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,GREY);
			else fillrect(x+4,y+4,dy-8,dy-8,BLACK);
		}
		else{
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,LLGREY);
		}
		
		plottext(text,x+dy+4,y+14,"bold 14px arial",BLACK);
		break;

	case CHECKBUT3:
		fillrect(x,y,dx,dy,WHITE);
		drawrect(x+2,y+2,dy-4,dy-4,BLACK,NORMLINE);

		if(res.popfiltcheck[val][val2] == 1){
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,GREY);
			else fillrect(x+4,y+4,dy-8,dy-8,BLACK);
		}
		else{
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,LLGREY);
		}
		
		plottext(text,x+dy+4,y+14,"bold 14px arial",BLACK);
		break;
		
	case CHECKBUT4:
		fillrect(x,y,dx,dy,WHITE);
		drawrect(x+2,y+2,dy-4,dy-4,BLACK,NORMLINE);

		if(transindcheck[val] == 1){
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,GREY);
			else fillrect(x+4,y+4,dy-8,dy-8,BLACK);
		}
		else{
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,LLGREY);
		}
		
		plottext(text,x+dy+4,y+14,"bold 14px arial",BLACK);
		break;
		
	case CHECKBUT5:
		fillrect(x,y,dx,dy,WHITE);
		drawrect(x+2,y+2,dy-4,dy-4,BLACK,NORMLINE);

		if(distCI == 1){
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,GREY);
			else fillrect(x+4,y+4,dy-8,dy-8,BLACK);
		}
		else{
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,LLGREY);
		}
		
		plottext(text,x+dy+4,y+14,"bold 14px arial",BLACK);
		break;
		
	case CHECKBUT6:
		fillrect(x,y,dx,dy,WHITE);
		drawrect(x+2,y+2,dy-4,dy-4,BLACK,NORMLINE);

		if(distmean == 1){
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,GREY);
			else fillrect(x+4,y+4,dy-8,dy-8,BLACK);
		}
		else{
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,LLGREY);
		}
		
		plottext(text,x+dy+4,y+14,"bold 14px arial",BLACK);
		break;
		
	case CHECKBUT7:
		fillrect(x,y,dx,dy,WHITE);
		drawrect(x+2,y+2,dy-4,dy-4,BLACK,NORMLINE);

		if(res.varcheck[val] == 1){
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,GREY);
			else fillrect(x+4,y+4,dy-8,dy-8,BLACK);
		}
		else{
			if(over == i) fillrect(x+4,y+4,dy-8,dy-8,LLGREY);
		}
		
		plottext(text,x+dy+4,y+14,"bold 14px arial",BLACK);
		break;
		
	case CLEARBUT2:
		if(modelsetup == 1){ col = WHITE; col2 = RED; col3 = RED; if(over == i) col = LLRED;}
		else{ col = WHITE; col2 = DGREEN; col3 = DGREEN; if(over == i) col = LLGREEN;}
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2,y+21, "bold 18px arial",col3);
		break;
		
	case RADIOBUT:
		fillrect(x,y,dx,dy+4,WHITE);
	
		col=DGREY;
		switch(val2){
		case RADIOTERM: selval = termtype; break;
		case RADIOVAR: selval = -1; for(vs = 0; vs < res.varsel.length; vs++){ if(res.varsel[vs] == val) selval = val;} break;
		case RADIOTRANS: selval = -1; for(vs = 0; vs < transsel.length; vs++){ if(transsel[vs] == val) selval = val;} break;
		case RADIOCHOOSE: selval = -1; break;
		case RADIOLABEL: selval = radlab; break;
		case RADIOSTYLE: selval = radstyle; break;
		case RADIOSTYLE2: selval = radstyle2; break;
		case RADIONAMLAB: selval = radnamlab; break;
		case RADIODER: selval = -1; for(vs = 0; vs < res.dersel.length; vs++){ if(res.dersel[vs] == val) selval = val;} break;
		case RADIODEP: selval = -1; for(vs = 0; vs < depsel.length; vs++){ if(depsel[vs] == val) selval = val;}  break;
		case RADIOTRANSTY: selval = transty; break;
		case RADIOSCATTY:  selval = res.scattype; break;
		}
	
		if(val == selval) col2 = BLACK;
		else{
			if(over == i) col2 = DGREY;
			else col2 = WHITE
		}

		r = 7;
		drawroundrect(x+1,y+2,2*r,2*r,r,WHITE,col);

		if(val == selval) col3 = BLACK; 
		else{
			col3 = GREY;  if(over == i) col3 = DGREY;
		}
    
		if(val2 == RADIOTRANSTY) plotlabel(buttext[i],x+20,y+14,"bold 16px arial",col3,dx-20);
		else plotlabel(buttext[i],x+20,y+14,RADIOFONT,col3,dx-20);
		
		if(col2 != WHITE){
			x += 3; y += 5; r -= 3;
			drawroundrect(x+1,y,2*r,2*r,r,col2,col2);
		}
		break;
		
	case BAYESBUT:
		cv.clearRect(x-1,y-1,dx+2,dy+2); r = 6;
		col = BLACK; if(over == i) col = DDGREY;
		drawroundrect(x,y,dx,dy,r,col,DDPURPLE);
		centertext(buttext[i],butx[i]+butdx[i]/2,buty[i]+15,RADIOFONT,WHITE);
		break;
	
	case GDROPSELBUT:
		r = Math.floor(0.4*gdropdy);
		switch(val2){
		case 1:
			fillrect(x,y,dx,dy,LLLBLUE);
			col3=BLACK;
			drawroundrect(x,y,dx,dy,r,WHITE,LBLUE); 
			drawline(x,y+gdropdy,x+dx,y+gdropdy,LBLUE,NORMLINE);
			break;
		case 2:
			fillrect(x,y,dx,dy,WHITE);
			col3=WHITE;
			drawroundrect(x,y,dx,dy,r,WHITE,DDBLUE); 
			drawroundmenutop(x,y,dx,gdropdy,r,DBLUE,DDBLUE)	 
			break;
		case 3:
			fillrect(x,y,dx,dy,LLLBLUE);
			col3=WHITE;
			drawroundrect(x,y,dx,dy,r,WHITE,DRED); 
			drawroundmenutop(x,y,dx,gdropdy,r,RED,DRED)	 
			break;
		case 4:
			fillrect(x,y,dx,dy,WHITE);
			col3=BLACK;
			drawroundrect(x,y,dx,dy,r,WHITE,LBLUE); 
			drawline(x,y+gdropdy,x+dx,y+gdropdy,LBLUE,NORMLINE); 
			break;
		case 5:
			fillrect(x,y,dx,dy,WHITE);
			col3=BLACK;
			drawroundrect(x,y,dx,dy,r,WHITE,LRED); 
			drawline(x,y+gdropdy,x+dx,y+gdropdy,LRED,NORMLINE); 
			break;
		case 6:
			fillrect(x,y,dx,dy,WHITE);
			col3=BLACK;
			drawroundrect(x,y,dx,dy,r,WHITE,LGREEN); 
			drawline(x,y+gdropdy,x+dx,y+gdropdy,LGREEN,NORMLINE); 
			break;
		case 7:
			fillrect(x,y,dx,dy,WHITE);
			col3=BLACK;
			drawroundrect(x,y,dx,dy,r,WHITE,LORANGE); 
			drawline(x,y+gdropdy,x+dx,y+gdropdy,LGREEN,NORMLINE); 
			break;
		}
		
		x1 = x+dx-17; y1 = y+Math.floor(0.3*gdropdy); ddx = 5; ddy = 6;
		drawline(x1,y1,x1+ddx,y1+ddy,col3,THICKLINE); drawline(x1+2*ddx,y1,x1+ddx,y1+ddy,col3,THICKLINE);
		centertext(text,x+dx/2-8,y+Math.floor(0.8*gdropdy),Math.floor(0.85*gdropdy)+"px arial",col3,dx);
		
		gdropslider = 0;
		if(gdropnum < val.length){
			dddy = gdropnum*gdropdy;
			
			if(drag == 101){
				gdropfrac = gdropfracst + (my-gdropmyst)/(gdropnum*gdropdy);
				if(gdropfrac < 0) gdropfrac = 0; if(gdropfrac > 1-gdropnum/val.length) gdropfrac = 1-gdropnum/val.length;
			}
			dx -= 20;
			xmin = x+dx+3; ymin = y+Math.floor(1.15*gdropdy)+gdropfrac*dddy; ddx = 14; ddy = Math.floor(dddy*gdropnum/val.length);
			
			drawrect(xmin,y+23,14,dddy,GREY,NORMLINE);
			if(mx >= xmin && mx < xmin+ddx && my >= ymin && my < ymin+ddy){ col = LGREY; gdropslider = 1;}
			else col = GREY;
			fillrect(xmin,ymin,14,ddy,col);
			
		}		
		
		gdropselop = "&*&"; gdropselj = -1;
		jmin = Math.floor(val.length*gdropfrac);
		for(j = jmin; j < jmin+gdropnum; j++){
			yy = y+Math.floor(1.2*gdropdy) + (j-jmin)*gdropdy;
			if(my > yy && my < yy + gdropdy && mx < x+dx){
				fillrect(x+3,yy,dx-6,gdropdy,LLRED);
				gdropselop = val[j]; gdropselj = j;
			}				
			centertext(val[j],x+dx/2,yy+Math.floor(0.8*gdropdy),Math.floor(0.85*gdropdy)+"px arial",BLACK,dx);		
		}
		break;
		
	case GDROPBUT:	
		r = Math.floor(0.4*dy);
		if(over == i) cv.globalAlpha = 0.7;
		switch(val2){
		case 1: fillrect(x,y,dx,dy,LLLBLUE); col1 = WHITE; col2 = LBLUE; col3 = BLACK; break;
		case 2: fillrect(x,y,dx,dy,WHITE); col1 = DBLUE; col2 = DDBLUE; col3 = WHITE; break;
		case 3: fillrect(x,y,dx,dy,LLLBLUE); col1 = RED; col2 = DRED; col3 = WHITE; break;
		case 4: fillrect(x,y,dx,dy,WHITE); col1 = WHITE; col2 = LBLUE; col3 = BLACK; if(over == i) col1 = LLBLUE; break;
		case 5: fillrect(x,y,dx,dy,WHITE); col1 = WHITE; col2 = LRED; col3 = BLACK; if(over == i) col1 = LLRED; break;
		case 6: fillrect(x,y,dx,dy,WHITE); col1 = WHITE; col2 = LGREEN; col3 = BLACK; if(over == i) col1 = LLGREEN; break;
		case 7: fillrect(x,y,dx,dy,WHITE); col1 = WHITE; col2 = LORANGE; col3 = BLACK; if(over == i) col1 = LLORANGE; break;
		}
		drawroundrect(x,y,dx,dy,r,col1,col2);
		x1 = x+dx-17; y1 = y+Math.floor(0.3*dy); ddx = 5; ddy = 6;
		drawline(x1,y1,x1+ddx,y1+ddy,col3,THICKLINE); drawline(x1+2*ddx,y1,x1+ddx,y1+ddy,col3,THICKLINE);
		
		if(val2 == 1) centertext(text,x+dx/2-8,y+Math.floor(0.8*dy),Math.floor(0.85*dy)+"px arial",col3,dx-10);
		if(over == i) cv.globalAlpha = 1;
		if(val2 != 1) centertext(text,x+dx/2-8,y+Math.floor(0.8*dy),Math.floor(0.85*dy)+"px arial",col3,dx-10);
		break;
	
	case STATBUT:
		plotlabel(text,x,y,STATFONT,BLACK);
		break;
	
	case STATBUT2:
		plottext(text,x,y,"bold 18px Times",BLACK);
		break;
		
	case STATBUT3:
		plottext(text,x,y,STATFONT,BLACK);
		break;	
 
	case PLOTTIMEBUT:
		drawline(x,y,x+dx+20,y,BLACK,THICKLINE);
		drawarrow(x+dx+30,y,x,y,15,BLACK);
		break;
		
	case EXPORTBUT:
		col = ORANGE; col2 = DORANGE; col3 = WHITE;  if(i == over){ col = LORANGE;}
		drawroundrect(butx[i],buty[i],butdx[i],butdy[i],12,col,col2);
		centertext(buttext[i],butx[i]+butdx[i]/2,buty[i]+21,INPUTFONT,col3);
		break;
	
	case PROGRESSBUT:
		var ch;
		for(ch = 0; ch < nchrun; ch++) if(exe[ch] == 0) break;
		if(ch < nchrun){
			if(nchrun > 1) centertext("Initialising run "+(ch+1),menux+(width-menux)/2,height/2+60,TILEFONT,RED);
			else centertext("Initialising",menux+(width-menux)/2,height/2+60,TILEFONT,RED);
		}
		else{
			if(page == SIMULATEPAGE){
				if(simsta == 0) centertext("Setup",x,y,TILEFONT,RED);
				else centertext("Simulating",x,y,TILEFONT,RED);
			}
		}
		break;
		
	case EXPORTMINIBUT:
		if(val2 == -1){ col = ORANGE; col2 = DORANGE; col3 = WHITE;  if(i == over) col = LORANGE;}
		else{ col = DBLUE; col2 = DDBLUE; col3 = WHITE; if(i == over) col = LBLUE;}
		drawroundrect(butx[i],buty[i],butdx[i],butdy[i],5,col,col2);
		centertext(buttext[i],butx[i]+butdx[i]/2,buty[i]+13,RADIOFONT,col3);
		break;
		
	case ADDBUT2: case ADDDEPBUT:
		fillrect(x,y,dx,dy,"#ffffff");
	
		col = BLACK; if(over == i) col = GREY;
		x = butx[i]+8; y = buty[i]+dy/2; r = 6; r2 = 4;
		drawroundrect(x-r,y-r,2*r,2*r,r,col,col);
		drawline(x-r2,y,x+r2,y,WHITE,THICKLINE);
		drawline(x,y-r2,x,y+r2,WHITE,THICKLINE);
		plottext(text,x+11,buty[i]+16,ADDFONT,col); 
		break;
	
	case CHOOSETIMEBUT:
		cv.clearRect(x,y,dx,dy);
		col = BLUE; if(over == i) col = LBLUE;
		plottext(text,x+4,y+19,HELPFONT,col);
		break;
	
	case RELOADSTATBUT:
		cv.clearRect(x,y,dx,dy);
		if(vcalc >= res.filtervar[res.filt].length){ col = RED; if(over==i) col = LRED;} else col = LLRED;
		reloadsign(x,y,col);
		break;
		
	case RELOADCORBUT:
		cv.clearRect(x,y,dx,dy);
		cv.clearRect(x-10,y-10,dx+20,dy+20);
		col = RED; if(over==i)  col = LRED;
		reloadsign(x,y,col);
		break;
		
	case PRHEADBUT:
		plottext(text,x+4,y+16,tableheadfont,DDBLUE);
		break;
		
	case IMPBUT:
		fillrect(x,y,dx,dy,"#000044");
		col = WHITE; col2 = WHITE; col3 = DDBLUE;  if(i == over){ col = LRED; col2 = col;}
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2,y+21,INPUTFONT,col3);
		break;
		
	case SIMINFBUT:
		fillrect(x,y,dx,dy,"#000044");
		if(val == 0){ col = LBLUE; if(i == over) col = LLBLUE;}
		else{ col = LRED; if(i == over) col = LLRED;}
		col2 = col; col3 = DDBLUE; 
		drawroundrect(x,y,dx,dy,12,col,col2);
		centertext(text,x+dx/2,y+21,INPUTFONT,col3);
		break;
		
	case LOADFILEBUT:
		fillrect(x,y,dx,dy,"#000044");
		col = LRED; if(i == over) col = DRED;
		ddx = 3; ddx2 = 6; ddy = 6;
		drawline(x+ddx,y+dy/2-ddy,x+ddx+ddx2,y+dy/2,col,THICKLINE)
		drawline(x+ddx,y+dy/2+ddy,x+ddx+ddx2,y+dy/2,col,THICKLINE)
		plottext(text,x+15,y+18,HELPBUTFONT,col);
		break;
		
	case HELPBACKBUT:
		drawroundrect(x,y,dx,dy,20,"#000044","#000000");  
		plottext(val,x+5+20,y+38,HELPBUTFONT,WHITE);
		y += 50;
		
		var d;
		for(d = 0; d < text.length; d++){
			col = WHITE;
			dd = 0; te = text[d]; 
			if(te.substr(0,5) == "<bul>"){ dd = 20; te = te.substring(5); col = LLBLUE; fillcircle(x+5+27,y+12,4,col,col,NORMLINE);} 
			alignparagraph(te,dx-35-dd,HELPBUTFONT2);
			var j; for(j = 0; j < nlines; j++){ plottext(lines[j],x+5+20+dd,y+18,HELPBUTFONT2,col); y += 25;}
			y += 10;
		}
		break;
		
	case HELPCLOSEBUT:
		fillrect(x,y,dx,dy,"#000044");
		col = WHITE; if(i == over) col = RED;
		drawline(x+3,y+3,x+dx-3,y+dy-3,col,VTHICKLINE);
		drawline(x+3,y+dy-3,x+dx-3,y+3,col,VTHICKLINE);
		break;
		
	case HELPICONBUT:
		if(val2 == -1) fillrect(x,y,dx,dy,WHITE); else fillrect(x,y,dx,dy,"#ededff"); 
		
		if(val2 == -2){ col = DBLUE; if(i == over) col = LBLUE;}
		else{ col = DRED; if(i == over) col = LRED;}
		plottext(text,x+3,y+10,"bold 10px arial",col);
		break;
		
	case CHECKALLBUT:
		fillrect(x,y,dx,dy,WHITE); 
		if(text == "Check All"){ col = DBLUE; if(i == over) col = LBLUE;}
		else{ col = DRED; if(i == over) col = LRED;}
		
		plottext(text,x+3,y+10,"bold 10px arial",col);
		break;
		
	case DELFILEBUT:
		fillrect(x,y,dx,dy,DDBLUE); 
		col = RED; if(i == over) col = LLRED;
		drawline(x+1,y+1,x+dx-1,y+dy-1,col,VTHICKLINE);
		drawline(x+1,y+dy-1,x+dx-1,y+1,col,VTHICKLINE);
		break;
	
	case DONOTHINGBUT: break;
		
	default: alertp("Error code EC1"); break;
  }
}
