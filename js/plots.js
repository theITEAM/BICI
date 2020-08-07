function showtraces()                                     // Sets up trace plot
{
	var i, w, ww, x, y, yy, c, cc = 0, tx, ty, mm, val, vv;
 
	if(xaxisauto == 1){
		axxmin = 0; 
		for(cu = 0; cu < curvevar.length; cu++){
			ch = curvech[cu];
			while(res.ch[ch].nsamp > res.nsamppl) res.nsamppl =  Math.floor(1.5*res.nsamppl);
		}
	
		axxmax = res.nsamppl;
		axxmax *= 1.001;
	}
	else{ axxmin  = xaxisfixmin; axxmax = xaxisfixmax;}
	
	if(yaxisauto == 1){
		axymin = large; axymax = -large;
		for(cu = 0; cu < curvevar.length; cu++){
			v = curvevar[cu]; ch = curvech[cu];
		
			if(ch == res.nch+1){
				for(ch = 0; ch < res.nch; ch++){
					for(i = res.burnin; i < res.ch[ch].nsamp; i++){
						val = res.ch[ch].varval[v][i];
						if(val < axymin) axymin = val; 
						if(val > axymax) axymax = val; 
					}
				}
			}
			else{
				for(i = res.burnin; i < res.ch[ch].nsamp; i++){
					val = res.ch[ch].varval[v][i];
					if(val < axymin) axymin = val; 
					if(val > axymax) axymax = val; 
				}
			}
		}

		if(axymin >= axymax-0.00001){ axymin--; axymax++;} 
		if(axymax > 0 && axymin > 0){ if(axymin < 0.3*axymax) axymin = 0;}
		axymax += 0.001*(axymax-axymin);
	}
	else{ axymin = yaxisfixmin; axymax = yaxisfixmax;}
	
	setxtics();
	w = setytics();   
	
	if(res.varsel.length == 1) laby = res.varname[res.varsel[0]]; else laby = "Parameter Value";
		
	graphframe(menux+10,80,40+w,60,20,40,w,"Samples",laby,"tr");
	xx = Math.floor(graphdx*(res.burnin-axxmin)/(axxmax-axxmin));
	if(xx >= 0 && xx < graphdx) addcanbutton("",40+w+xx-cc,40,2*cc,graphdy,BURNINBUT,BURNINBUT,-1,-1); 
	
	drawtrace();
	loading = 0;
}

function graphframe(x,y,x1,y1,x2,y2,w,labx,laby,arg)      // Draws the frame for graphs
{	
	graphdx = (width-360)-w; graphdy = height-185+40-y2;
	if(arg == "genpop" || arg == "gender"){ graphdx = (width-270)-w; graphdy = height-185-y2;}
	
	addbutton("",x,y,x1+x2+graphdx,y1+y2+graphdy,CANVASBUT,CANVASBUT,-1,-1);
	addcanbutton("",x1,y2,graphdx,graphdy,-1,RESULTBUT,-1,-1);		
	if(arg == "scat"){
		addcanbutton(labx,x1,y2+graphdy+30,graphdx,30,XLABELBUT,XLABELBUT,-1,-1);
		addcanbutton(laby,0,y2,30,graphdy,YLABELBUT,YLABELBUT,-1,-1);
	}
	else{
		addcanbutton(labx,x1,y2+graphdy+30,graphdx,30,-1,XLABELBUT,-1,-1);
		addcanbutton(laby,0,y2,30,graphdy,-1,YLABELBUT,-1,-1);
	}
	
	if(ntickx > 0) addcanbutton("X ticks",x1,y2+graphdy,graphdx,30,XTICKBUT,XTICKBUT,-1,-1);
	if(nticky > 0) addcanbutton("Y ticks",x1-w,y2,w,graphdy,YTICKTRBUT,YTICKTRBUT,-1,-1);
	
	lablist=[]; labcollist=[];
	
	switch(arg){
	case "tr":
		for(cu = 0; cu < curvevar.length; cu++){
			v = curvevar[cu]; ch = curvech[cu];
			if(ch < infres.nch){ 
				if(res.varsel.length > 1){ de = res.varname[v]; if(infres.nch > 1) de += " run "+(ch+1);}
				else{ de = ""; if(infres.nch > 1) de += "Run "+(ch+1);}
			}
			else de = res.varname[v];
			
			lablist[cu] = de;
		}
		labcollist = chcol;
		break;

	case "trandi":
		for(cu = 0; cu < curvevar.length; cu++){
			v = curvevar[cu]; ch = curvech[cu];
			if(ch < infres.nch){ 
				if(depsel.length > 1){ de = transpos[v].name; if(infres.nch > 1) de += " run "+(ch+1);}
				else{ de = ""; if(infres.nch > 1) de += "Run "+(ch+1);}
			}
			else de = transpos[v].name;
			
			lablist[cu] = de;
		}
		labcollist = chcol;
		break;
	
	case "pop":
		for(cu = 0; cu < nquant; cu++){ lablist.push(quantname[cu]); labcollist.push(quantcol[cu]);}
		break;
		
	case "histo":
		for(cu = 0; cu < depsel.length; cu++) lablist.push(histoplot[depsel[cu]].name);
		labcollist = chcol;
		break;
	}
	
	if(arg == "scat"){
		if(res.runsel == res.nch && res.varselx >= 0 && res.varsely >= 0){
			for(ch = 0; ch < res.nch; ch++) lablist[ch] = "Run "+(ch+1);
			labcollist = chcol;
		}
	}

	if(lablist.length > 1){
		cy = 10;
		var i = 0, ist;
		do{
			var xli=[];
			ist = i;
			x = 0;
			while(i < lablist.length){
				dx = 70 + textwidth(lablist[i],KEYFONT);
				if(x + dx > graphdx-20 && xli.length > 0) break;
				xli.push(x);
				x += dx;
				i++;
			}
		
			for(ii = ist; ii < i; ii++) addcanbutton(lablist[ii],x1+x2+graphdx-20-x+xli[ii-ist],cy,40,13,-1,LABBUT,labcollist[ii],ii);	
			cy += 20;
		}while(i < lablist.length);
	}
}

function drawtrace(wid)                                   // Draws trace plot
{
	var step, i, burn;
	step = Math.floor(axxmax/1000); if(step < 1) step = 1;

	cv.clearRect(0,0,graphdx,graphdy);

	if(wid) cv.lineWidth = wid;
	else cv.lineWidth = 1;
	
	for(cu = 0; cu < curvevar.length; cu++){
		v = curvevar[cu]; ch = curvech[cu];
	
		cv.beginPath();
		burn = 1;
		cv.strokeStyle = chcol[cu]; 
		for(i = 0; i < res.ch[ch].nsamp; i+=step){
			if(burn == 1 && i >= res.burnin){ cv.globalAlpha = 0.2; cv.stroke(); cv.globalAlpha = 1; cv.beginPath(); burn = 0;}
			x = graphdx*(i-axxmin)/(axxmax-axxmin);
			y = graphdy*(1-(res.ch[ch].varval[v][i]-axymin)/(axymax-axymin));
			if(i == 0) cv.moveTo(x,y); else cv.lineTo(x,y);
		}
		cv.stroke();
	}
}

function rn(n)                                            // Rounds a number when put in interface
{
	var x;
	n = parseFloat(n);

	x = 0;
	do{
		num = n.toFixed(x);
		if((num-n)*(num-n) < 0.000000000000001) break;
		x++;
	}while(1 == 1);
	return num;
}

function showcorrelation()
{
	var v, vv;

	nob = 0;
	if(corlist.length == 0){
		addcanbutton(val,x,y,cordy,cordy,CORBUT,CORBUT,val2,val3);
	}
	for(v = 0; v < corlist.length; v++){
		addob(0,cordx+v*cordy,OBHOZLAB,res.varname[corlist[v]]);
		addob(cordx+v*cordy,0,OBVERTLAB,res.varname[corlist[v]]);
		for(vv = 0; vv < corlist.length; vv++){
			addob(cordx+v*cordy,cordx+vv*cordy,OBCOR,cormat[v][vv],corlist[v],corlist[vv]);
		}
	}		
	ytot = cordx+corlist.length*cordy;
	graphdx = ytot; graphdy = ytot;
	
	placeob();
	
	addbutton("Reload",menux+27,height-48,20,25,RELOADCORBUT,RELOADCORBUT,-1,-1);
}

function drawcortable()
{
	var v, vv;
	
	cordx *= 2; cordy *= 2; 
	fillrect(0,0,cordx+cordy*corlist.length,cordx+cordy*corlist.length,WHITE);
	for(v = 0; v < corlist.length; v++){
		plotxlabel(res.varname[corlist[v]],cordx/2,cordx+v*cordy+cordy/2+8,"40px times",BLACK);
		plotylabel(res.varname[corlist[v]],cordx+v*cordy+cordy/2+8,cordx/2,"40px times",BLACK);
		for(vv = 0; vv < corlist.length; vv++){
			 plotcorbut(cormat[v][vv],cordx+v*cordy,cordx+vv*cordy,cordy,cordy,0,28);
		}
	}	
	cordx /= 2; cordy /= 2; 
}

function plotcorbut(text,x,y,dx,dy,ov,fosi)
{
	var num;
	if(text == "---"){
		fillrect(x,y,dx,dy,WHITE);
		centertext(text,x+dx/2,y+dy/2+5*fosi/14,"bold "+fosi+"px arial",BLACK,dx);
	}
	else{
		if(text > 0){
			num = Math.floor(255-text*255);
			fillrect(x,y,dx,dy,"rgb(255,"+num+","+num+")");
		}
		else{
			num = Math.floor(255+text*255);
			fillrect(x,y,dx,dy,"rgb("+num+","+num+",255)");
		}
		if(fosi == 14) centertext(text.toFixed(2),x+dx/2,y+dy/2+5,"bold "+fosi+"px arial",BLACK,dx);
		else  centertext(text.toFixed(2),x+dx/2,y+dy/2+10," "+fosi+"px arial",BLACK,dx);
	}
	if(ov == 1){
		cv.globalAlpha = 0.2;
		fillrect(x,y,dx,dy,WHITE);
		fillrect(x,y,dx,dy,BLACK);
		cv.globalAlpha = 1;
		drawrect(x+2,y+2,dx-4,dy-4,BLACK,2);
	}
}

function showscatter()                                   // Sets up a scatter plot
{ 
	var x, y, k, ch, dx, dx2, w, cc = 6, xst, yst;
	
	setscatteraxis();

	setxtics();
	w = setytics();

	switch(res.scattype){
	case "KDE":	if(plotinitfl == 0) calc2dkernal(); break;
	case "scatter":
		if(plotinitfl == 0){ 	
			var totsamp=[];
			for(ch = chmin; ch < chmax; ch++)
				for(i = res.burnin; i < res.ch[ch].nsamp; i++){	totsamp.push({ch:ch,i:i});
			}
		
			scatsamp=[];
			for(k = 0; k < res.scatnum; k++){
				if(totsamp.length == 0) break;
				j = Math.floor(Math.random()*totsamp.length);
				scatsamp.push({ch:totsamp[j].ch, i:totsamp[j].i});
				totsamp.splice(j,1);
			}
			plotinitfl = 2;
		}
		break;
	}
	
	labx = res.varname[res.varselx];
	laby = res.varname[res.varsely];
	
	if(res.scattype == "scatter") graphframe(menux+10,80,40+w,60,20,40,w,labx,laby,"scat");
	else graphframe(menux+10,80,40+w,60,20,40,w,labx,laby,"");
		
	cv.clearRect(0,0,graphdx,graphdy);
	drawscatterplot();
	
	addcanbutton("Reload",18,30+graphdy+34,20,20,RELOADBUT3,RELOADBUT3,-1,-1); 
	
	if(res.scattype == "KDE"){
		addcanbutton("Smooth:",graphdx+w-slidex+45,graphdy+75,slidex,slidedy,-1,SLIDERBACKBUT,-1,-1); 
		addcanbutton("",graphdx+w-slidex+45+Math.floor(slidex*(kde-kdemin)/(kdemax-kdemin) - slidedx/2),graphdy+75,slidedx,slidedy,SLIDERBUT,SLIDERBUT,-1,-1); 
	}
	
	loading = 0;
}

function setscatteraxis()                                 // Sets axes for scatter plot
{
	var v, i;
	
	if(res.runsel < res.nch){ chmin = res.runsel; chmax = res.runsel+1;}
	else{ chmin = 0; chmax = res.nch;}

	if(xaxisauto == 1){
		axxmin = large; axxmax = -large;
		for(ch = chmin; ch < chmax; ch++){
			for(i = res.burnin; i < res.ch[ch].nsamp; i++){
				val = res.ch[ch].varval[res.varselx][i];
				if(val < axxmin) axxmin = val; 
				if(val > axxmax) axxmax = val; 
			}
		}
		if(Math.sqrt(axxmin*axxmin) < 0.2*Math.sqrt(axxmax*axxmax)) axxmin = 0;
	}
	else{ axxmin = xaxisfixmin; axxmax = xaxisfixmax;}

	if(yaxisauto == 1){
		axymin = large; axymax = -large;
		for(ch = chmin; ch < chmax; ch++){
			for(i = res.burnin; i < res.ch[ch].nsamp; i++){
				val = res.ch[ch].varval[res.varsely][i];
				if(val < axymin) axymin = val; 
				if(val > axymax) axymax = val; 
			}
		}
		if(Math.sqrt(axymin*axymin) < 0.2*Math.sqrt(axymax*axymax)) axymin = 0;
	}
	else{ axymin = yaxisfixmin; axymax = yaxisfixmax;}
}

function calc2dkernal()
{
	var x, y;

	if(kernid != undefined){ clearInterval(kernid); kernid = undefined;}
	
	Kbin=[]; for(x = 0; x < KX; x++){ Kbin[x]=[]; for(y = 0; y < KY; y++) Kbin[x][y] = 0;}
	
	startloading();
	kernid = setTimeout(function(){ calc2dkernal2();}, 100);

	plotinitfl = 1;
}

function calc2dkernal2()
{
	var ch, i, x, y, xx, yy, d, max;
	
	for(ch = chmin; ch < chmax; ch++){
		for(i = res.burnin; i < res.ch[ch].nsamp; i++){			
			x = Math.floor(KX*(res.ch[ch].varval[res.varselx][i]-axxmin)/(axxmax-axxmin));
			y = Math.floor(KY*(res.ch[ch].varval[res.varsely][i]-axymin)/(axymax-axymin));
		
			for(xx = Math.floor(x-kde); xx <= Math.floor(x+kde); xx++){
				for(yy = Math.floor(y-kde); yy <=  Math.floor(y+kde); yy++){
					if(xx >= 0 && xx < KX && yy >= 0 && yy < KY){
						d = Math.sqrt((x-(xx+0.5))*(x-(xx+0.5)) + (y-(yy+0.5))*(y-(yy+0.5)));
						if(d < kde) Kbin[xx][yy] += (1-d/kde);
					}
				}
			}			
		}
	}
	
	max = 0;
	for(x = 0; x < KX; x++){ for(y = 0; y < KY; y++){ if(Kbin[x][y] > max) max = Kbin[x][y];}}
	for(x = 0; x < KX; x++){ for(y = 0; y < KY; y++) Kbin[x][y] /= max;}

	plotinitfl = 2;
	buttoninit();
}

function drawscatterplot(wid)                             // Draws a scatter plot
{
	var i, j, k, max, val, x, y, ch, x1,x2, y1, y2;
	
	switch(res.scattype){
	case "KDE":
		for(x = 0; x < KX; x++){ 
			for(y = 0; y < KY; y++){
				col = Math.floor(255-255*Kbin[x][y]);
				if(col != 255){
					x1 = Math.floor(x*graphdx/KX); x2 = Math.floor((x+1)*graphdx/KX);
					y1 = graphdy - Math.floor((y+1)*graphdy/KY); y2 = graphdy - Math.floor(y*graphdy/KY); 
					fillrect(x1,y1,x2-x1,y2-y1,"rgb("+col+","+col+","+255+")");
				}
			}
		}
		break;
	
	case "scatter":
		for(k = 0; k < scatsamp.length; k++){
			ch = scatsamp[k].ch; i = scatsamp[k].i;
			if(res.runsel == res.nch) col = chcol[ch]; else col = BLACK;
		
			cv.fillStyle = col; 
			x = Math.floor(graphdx*(res.ch[ch].varval[res.varselx][i]-axxmin)/(axxmax-axxmin));
			y = Math.floor(graphdy - graphdy*(res.ch[ch].varval[res.varsely][i]-axymin)/(axymax-axymin));
			if(wid) fillcircle(x,y,wid,col,col,NORMLINE);
			else cv.fillRect(x-1,y-1,3,3);
		}
		break;
	}
}

function transinit()                                       // Initialises the graph showing transition distributions 
{
	var cl, ch, v, cl2, cu, c, chmax, chmin;
	
	startloading();

	res = infres;
	
	transcfilt=[]; calclist=[]; culim=[];
	for(cu = 0; cu < curvevar.length; cu++){
		culim.push({mean:large, min:-large, max:large});
			
		transcfilt[cu]=[]; calclist[cu]=[];
		
		v = curvevar[cu]; ch = curvech[cu];
		cl = transpos[v].cl;
		for(c = 0; c < ncomp; c++){
			for(cl2 = 0; cl2 < ncla; cl2++){
				if(cl2 != cl){
					if(res.popfilt[cl2] >= 0){ if(compval[c][cl2] != res.popfilt[cl2]) break;}
					else{
						if(res.popfilt[cl2] == -2){
							if(res.popfiltcheck[cl2][compval[c][cl2]] == 0) break;
						}							
					}	
				}
			}
			if(cl2 == ncla) transcfilt[cu][c] = 1; else transcfilt[cu][c] = 0; 
		}
		
		
		if(ch == infres.nch+1){ chmin = 0; chmax = infres.nch;}
		else{ chmin = ch; chmax = ch+1;}
	
		for(ch = chmin; ch < chmax; ch++){
			for(s = infres.burninev; s < infres.ch[ch].nsampev; s++){
				calclist[cu].push({ch:ch, s:s});
			}
		}
	}

	labx = "Duration";
	ord=[];
	plotinitfl = 1;
	calcto = setTimeout(function(){ transcalc(0,0);}, 20);
}

function transcalc(cu,n)                                  // Calculates the transition distribution graph
{		
	var ch, s, ww, www, e, cc, tt, tim;

	tim = (new Date()).getTime();
	do{
		ch = calclist[cu][n].ch; s = calclist[cu][n].s; v = curvevar[cu]; cl = transpos[v].cl;
	
		ww = infres.ch[ch].sampev[s];
		for(i = 0; i < ww.nind; i++){
			if(transindfilt == "All" || transindcheck[i] == 1){
				www = ww.ind[i];
				c = www.evc[0]; t = www.evt[0];
				if(transmat[ncomp][c] == v){ ord.push({dt:t, type:0});} // source
				
				for(e = 1; e < www.nev; e++){
					cc = www.evc[e];
						
					if(transcfilt[cu][c] == 1){ 
						tt = www.evt[e];
						
						if(transmat[c][cc] == v){
							if(c != ncomp && transpos[v].i != compval[c][cl]) alertp("pr1");
							if(cc != ncomp && transpos[v].f != compval[cc][cl]) alertp("pr2");
							ord.push({dt:tt-t, type:0});
						}
						else{
							if(transpos[v].i >= 0){
								if(transpos[v].i == compval[c][cl]){ // leaves initial state without going to final state
									if(cc == ncomp || transpos[v].i != compval[cc][cl]){
										ord.push({dt:tt-t, type:1});
									}
								}
								else{
									if(cc != ncomp){
										if(transpos[v].i == compval[cc][cl]){ // enters initial state
											t = tt;
										}
									}
								}
							}
						}
					}
					c = cc;
				}
			}
		}
		n++;
	}while(n < calclist[cu].length && (new Date()).getTime()-tim < 100);
	
	if(n == calclist[cu].length){				
		ord.sort(function(a, b){return a.dt - b.dt}); 
			
		cusamp[cu]=[];
		numtot = ord.length; numrem = 0;
		
		for(i = 0; i < ord.length; i++){
			if(ord[i].type == 0){
				w = numtot/(numtot-numrem); if(w > 10) w = 10;
				cusamp[cu].push({val:ord[i].dt, w:w});
				numtot--;
			}
			else{
				numrem++;
			}
		}
		ord=[]; n = 0; cu++;
	}
		
	if(cu == curvevar.length){ plotinitfl = 2; distinit2(); loading = 0; buttoninit();}
	else{
		percent = Math.floor(100*(cu + n/calclist[cu].length)/curvevar.length);
		calcto = setTimeout(function(){ transcalc(cu,n);},20);
	}	
}

function distinit()                                       // Initialises the probability distribution graph
{
	var cu, v, ch, i;
	culim=[];
	for(cu = 0; cu < curvevar.length; cu++){
		cusamp[cu]=[];
		
		v = curvevar[cu]; ch = curvech[cu];
		if(ch == infres.nch+1){
			for(ch = 0; ch < infres.nch; ch++){
				for(i = res.burnin; i < res.ch[ch].nsamp; i++){
					cusamp[cu].push({val:res.ch[ch].varval[v][i],w:1});
				}
			}
		}
		else{
			for(i = res.burnin; i < res.ch[ch].nsamp; i++){
				cusamp[cu].push({val:res.ch[ch].varval[v][i],w:1});
			}
		}
		
		for(n = 0; n < cusamp[cu].length; n++) tempCI[n] = cusamp[cu][n].val;
		
		var me=large, mi=-large, ma =large;
		if(distCI == 1 || distmean == 1){ 
			calcCI(cusamp[cu].length);
			if(distCI == 1){ mi = CImin; ma = CImax;}
			if(distmean == 1) me = mean;
		}
		culim.push({mean:me, min:mi, max:ma});
	}
	plotinitfl = 2; loading = 0;
	distinit2();
	
	if(res.varsel.length == 1) labx = res.varname[res.varsel[0]]; else labx = "Parameter Value";
}

function distinit2()                                      // Initialises the probability distribution graph
{
	var x, y, dx, dx2, w, cc = 6, xst, yst, i, j, tx, val, bin=[], wei=[], range = 10, novari, ww, s, ord=[], num, num2, cl, cl2;
	
	if(xaxisauto == 1){
		axxmin = large; axxmax = -large;
		for(cu = 0; cu < curvevar.length; cu++){
			for(ii = 0; ii < cusamp[cu].length; ii++){
				val = cusamp[cu][ii].val;
				if(val < axxmin) axxmin = val; 
				if(val > axxmax) axxmax = val; 
			}
		}

		novari = 0; if(axxmin >= axxmax-0.00001){ axxmin--; axxmax++; novari = 1;} 
		if(axxmax > 0 && axxmin > 0){ if(axxmin < 0.3*axxmax) axxmin = 0;}

		if(subna == "Prob. Dist."){
			for(cu = 0; cu < curvevar.length; cu++){  // looks for the boundaries in the prior
				v = curvevar[cu];
				vpri = findparam(res.varname[v]);
				if(vpri < infres.param.length && infres.param[vpri].prior == "Flat"){
					valmin = parseFloat(infres.param[vpri].val[0]);
					if(valmin < axxmin && axxmin-valmin < 0.2*(axxmax-axxmin)) axxmin = valmin;
					valmax = parseFloat(infres.param[vpri].val[1]);
					if(valmax > axxmax && valmax-axxmax < 0.2*(axxmax-axxmin)) axxmax = valmax;
				}
			}
		}
		
		if(subna == "Trans. Dist."){
			for(cu = 0; cu < curvevar.length; cu++){  // looks for the boundaries in the prior
				na = transpos[curvevar[cu]].name;
				if(na.substr(0,2) == "+ "){ valmin = infres.tmin; valmax = infres.tmax;}
				else{ valmin = 0; valmax = axxmax*1.1;}
				if(valmin < axxmin) axxmin = valmin;
				if(valmax > axxmax) axxmax = valmax;
			}
		}
		
		axxmax += 0.0001*(axxmax-axxmin);
	}
	else{ axxmin  = xaxisfixmin; axxmax = xaxisfixmax;}
	
	for(cu = 0; cu < curvevar.length; cu++){  // looks for the boundaries in the prior
		Jbinpriorfac[cu]=[];
		for(j = 0; j < JX; j++) Jbinpriorfac[cu][j] = 1;
		
		fl = 0;
		if(subna == "Prob. Dist."){
			v = curvevar[cu];
			vpri = findparam(res.varname[v]);
				
			if(vpri < infres.param.length && infres.param[vpri].prior == "Flat"){
				valmin = parseFloat(infres.param[vpri].val[0]);
				valmax = parseFloat(infres.param[vpri].val[1]);
				fl = 1;
			}
		}
		
		if(subna == "Trans. Dist."){
			v = curvevar[cu];
			na = transpos[curvevar[cu]].name;
			if(na.substr(0,2) == "+ "){ valmin = infres.tmin; valmax = infres.tmax;}
			else{ valmin = 0; valmax = large;}
			fl = 1;
		}
	
		if(fl == 1){
			valmi = JX*(valmin-axxmin)/(axxmax-axxmin);
			valma = JX*(valmax-axxmin)/(axxmax-axxmin);
			for(j = 0; j < JX; j++){
				jf = j+0.5;
				if(jf < valmi || jf > valma) Jbinpriorfac[cu][j] = large;
				else{
					d = jf-valmi;
					if(d < kde) Jbinpriorfac[cu][j] -= 0.5*(kde-d)*(kde-d)/(kde*kde);
					d = valma-jf;
					if(d < kde) Jbinpriorfac[cu][j] -= 0.5*(kde-d)*(kde-d)/(kde*kde);
				}
			}
		}
	}
		
	smoothon = 0;
	max = 0; 
	for(cu = 0; cu < curvevar.length; cu++){  // KDE
		v = curvevar[cu]; ch = curvech[cu];
	
		inton = 0;
		if(subna == "Prob. Dist."){
			if(res.vartype[v] == "Trans." || res.varname[v] == "# Ind." || res.varname[v] == "Nerr" || res.varname[v] == "# Transitions"){
				if(axxmax-axxmin < 200) inton = 1;
			}
		}		
		if(inton == 0) smoothon = 1;
		
		Jbin[cu]=[];
		for(j = 0; j < JX; j++) Jbin[cu][j] = 0;

		wsum = 0;
		for(ii = 0; ii < cusamp[cu].length; ii++){
			if(inton == 1){  // if an integer value
				valmin = JX*(cusamp[cu][ii].val-0.5-axxmin)/(axxmax-axxmin);
				valmax = JX*(cusamp[cu][ii].val+0.5-axxmin)/(axxmax-axxmin);
				if(valmax >= 0 && valmin < JX){
					w = cusamp[cu][ii].w; wsum += w;
					if(valmax-valmin > 1){
						for(j = Math.floor(valmin+1); j <= Math.floor(valmax); j++){
							if(j >= 0 && j < JX) Jbin[cu][j] += w;
						}
					}
					else{
						j = Math.floor((valmax+valmin)/2 + 0.5);
						if(j >= 0 && j < JX) Jbin[cu][j] += w;
					}
				}
			}
			else{
				val = JX*(cusamp[cu][ii].val-axxmin)/(axxmax-axxmin);
				w = cusamp[cu][ii].w; wsum += w;
				for(j = Math.floor(val-kde); j <= Math.floor(val+kde+1); j++){
					if(j >= 0 && j < JX){
						d = val-(j+0.5);
						if(d > 0){
							if(d < kde) Jbin[cu][j] += w*(1-d/kde);  // triangular
						}
						else{
							if(d > -kde) Jbin[cu][j] += w*(1+d/kde);
						}
					}
				}
			}
		}
		
		for(j = 0; j < JX; j++) Jbin[cu][j] /= wsum;
		for(j = 0; j < JX; j++) Jbin[cu][j] /= Jbinpriorfac[cu][j];
		for(j = 0; j < JX; j++) if(Jbin[cu][j] > max) max = Jbin[cu][j];
	}
	
	for(cu = 0; cu < curvevar.length; cu++){
		for(j = 0; j < JX; j++) Jbin[cu][j] /= max;
	}
	
	setxtics();
	nticky = 0;	
}

function showdistribution()                               // Shows the distribution graph
{ 
	var arg;
	
	if(subna == "Trans. Dist."){ arg = "trandi"; if(plotinitfl == 0) transinit();}	
	else{ arg = "tr"; if(plotinitfl == 0) distinit();}	
	if(plotinitfl != 2) return;
	
	w = 10;
	graphframe(menux+10,80,40+w,60,20,40,w,labx,"Probability",arg);
		
	drawdistribution();

	addcanbutton("Reload",18,30+graphdy+34,20,20,RELOADBUT3,RELOADBUT3,-1,-1); 
		
	if(smoothon == 1){
		addcanbutton("Smooth:",graphdx+w-slidex+45,graphdy+75,slidex,slidedy,-1,SLIDERBACKBUT,-1,-1); 
		addcanbutton("",graphdx+w-slidex+45+Math.floor(slidex*(kde-kdemin)/(kdemax-kdemin) - slidedx/2),graphdy+75,slidedx,slidedy,SLIDERBUT,SLIDERBUT,-1,-1); 
	}
}

function setytics()                                       // Sets up distribution y axis
{
	i = nticksize-1; while(i >= 0 && Math.floor((axymax-axymin)/ticksize[i]) < 3) i--;
    axticky = ticksize[i];

    mm = Math.floor(axymin/axticky + 0.9999);
    nticky = 0;
    while(mm*axticky < axymax){
		ticky[nticky] = rn(mm*axticky); nticky++;
		mm++;
    }

	w = 0;
	for(i = 0; i < nticky; i++){
		ww = textwidth(ticky[i],TICKFONT);
		if(ww > w) w = ww;
	}
	
	return w;
}

function setxtics()                                      // Sets up ticks on x axis
{
	i = nticksize-1; while(i >= 0 && Math.floor((axxmax-axxmin)/ticksize[i]) < 4) i--;
	axtickx = ticksize[i];

	tx = Math.floor(axxmin/axtickx + 0.999999)*axtickx;
	ntickx = 0;
	while(tx < axxmax){
		tickx[ntickx] = rn(tx); ntickx++;	
		tx += axtickx;
	}
}

function drawdistribution(wid)                           // Draws distribution
{
	var i, j, max, val, x, y, CImin, CImax, imin, imax;
  
 	max = 0; for(cu = 0; cu < curvevar.length; cu++){ for(j = 0; j < JX; j++){ val = Jbin[cu][j]; if(val > max) max = val;}}
	max *= 1.05;

    cv.clearRect(0,0,graphdx,graphdy);
	for(cu = 0; cu < curvevar.length; cu++){
		if(wid) cv.lineWidth = wid; else cv.lineWidth = 3;
	
		v = curvevar[cu]; ch = curvech[cu];
		
		cv.beginPath(); 
		for(i = 0; i < JX; i++){
			val = Jbin[cu][i];
			x = Math.floor(graphdx*(i+0.5)/JX);
			y = Math.floor(graphdy-graphdy*val/max);
			if(i == 0) cv.moveTo(x,y); else cv.lineTo(x,y);
		}
		
		setdash(cu);
		 
		cv.strokeStyle = chcol[cu];
		cv.stroke();

		cv.fillStyle = chcol[cu];
	
		var vertline=[];
		if(1 == 0){
			y = graphdy;
			cv.lineTo(x,y);
			x = Math.floor(graphdx*(0.5)/JX);
			cv.lineTo(x,y);
		}
		else{   // Used to plot 95 CI region
			cv.beginPath(); 
		
			CImin = JX*(culim[cu].min-axxmin)/(axxmax-axxmin);
			CImax = JX*(culim[cu].max-axxmin)/(axxmax-axxmin);
			
			if(CImax > 0.5){
				imin = 0;
				if(CImin > 0.5 && CImin < JX-0.5){
					imin = Math.floor(CImin+0.5);
					f = CImin - (imin+0.5-1);
					x = Math.floor(graphdx*CImin/JX);
					y = Math.floor(graphdy-graphdy*(Jbin[cu][imin-1]*(1-f) + Jbin[cu][imin]*f)/max);
					
					cv.moveTo(x,graphdy);
					cv.lineTo(x,y);
					vertline.push({x1:x, y1:y, x2:x, y2:graphdy});
				}
				else{
					x = Math.floor(graphdx*(imin+0.5)/JX);
					cv.moveTo(x,graphdy);
				}
			
				if(CImax < JX-0.5) imax = Math.floor(CImax+0.5);
				else imax = JX;
				
				for(i = imin; i < imax; i++){
					val = Jbin[cu][i];
					x = Math.floor(graphdx*(i+0.5)/JX);
					y = Math.floor(graphdy-graphdy*val/max);
					cv.lineTo(x,y);
				}
				
				if(CImax < JX-0.5){
					f = CImax - (imax+0.5-1);
					x = Math.floor(graphdx*CImax/JX);
					y = Math.floor(graphdy-graphdy*(Jbin[cu][imax-1]*(1-f) + Jbin[cu][imax]*f)/max);
					cv.lineTo(x,y);
					vertline.push({x1:x, y1:y, x2:x, y2:graphdy});
				}
				cv.lineTo(x,graphdy);
			}
		}
		
		cv.globalAlpha = 0.1;
		cv.fill();
		
		for(k = 0; k < vertline.length; k++){
			if(wid) cv.lineWidth = wid/2; else cv.lineWidth = 1;
	
			cv.globalAlpha = 0.4;
			cv.beginPath(); 
			cv.moveTo(vertline[k].x1,vertline[k].y1);
			cv.lineTo(vertline[k].x2,vertline[k].y2);
			cv.stroke();
		}
		cv.globalAlpha = 1;
		
		x = graphdx*(culim[cu].mean-axxmin)/(axxmax-axxmin);
		if(x > 0 && x < graphdx){
			cv.lineWidth = scale;
			cv.beginPath(); 
			cv.moveTo(x,0);
			cv.lineTo(x,graphdy);
			cv.stroke();
		}
	}
	setdash(0);
}

function transplotinit()                                  // Prepares to plot transition distribution graph
{
	var c, cc, cl, cl2, tr, i, f, na;
	
	res = infres;
	
	transpos=[];
	transmat=[]; for(c = 0; c <= ncomp; c++){ transmat[c]=[]; for(cc = 0; cc <= ncomp; cc++) transmat[c][cc] = -1;}
					
	for(cl = 0; cl < ncla; cl++){
		for(tr = 0; tr < cla[cl].ntra; tr++){
			i = cla[cl].tra[tr].i; f = cla[cl].tra[tr].f;
			if(i < 0) na = "+"; else na = cla[cl].comp[i].name;
			na += " → ";
			if(f < 0) na += "-"; else na += cla[cl].comp[f].name;
			
			if(i < 0){
				for(cc = 0; cc < ncomp; cc++){
					if(compval[cc][cl] == f) transmat[ncomp][cc] = transpos.length;
				}
			}
			else{
				if(f < 0){
					for(c = 0; c < ncomp; c++){
						if(compval[c][cl] == i) transmat[c][ncomp] = transpos.length;
					}
				}
				else{
					for(c = 0; c < ncomp; c++){ 
						if(compval[c][cl] == i){
							for(cc = 0; cc < ncomp; cc++){
								if(compval[cc][cl] == f){
									for(cl2 = 0; cl2 < ncla; cl2++){
										if(cl2 != cl && compval[cc][cl2] != compval[c][cl2]) break;
									}
									if(cl2 == ncla) transmat[c][cc] = transpos.length;
								}
							}
						}
					}
				}
			}
			
			transpos.push({desc:"Transistion duration for "+na, name:na, cl:cl, i:i, f:f}); 
		}
	}
	transsel=[]; if(transpos.length > 0) transsel.push(0);
	transindfilt = "All"; for(i = 0; i < res.indlist.length; i++) transindcheck[i] = 1;
}


function histoinit(rest)                                  // Initialises age and time histograms 
{
	var v, fl = 0;
	
	depop=[];
	histoplot=[];
	depsel=[];
	for(cl = 0; cl < ncla; cl++){
		fl = 0;
		for(v = 0; v < rest.nvar; v++){
			var sp = rest.varname[v].split("_");
			if(sp.length == 2){
				var dep = sp[1].split(",");
				for(d = 0; d < dep.length; d++){
					for(j = 0; j < cla[cl].ncomp; j++){
						if(dep[d] == cla[cl].comp[j].name) break;
					}
				
					if(j < cla[cl].ncomp) break;
				}
				if(d < dep.length){
					na = sp[0]; fl = 0;
					for(dd = 0; dd < dep.length; dd++){
						if(dd != d){
							if(fl == 0) na += "_"; else na += ",";
							fl = 1;
							na += dep[dd];
						}
					}
					e = 0; while(e < histoplot.length && histoplot[e].name != na) e++;
					
					if(e == histoplot.length) histoplot.push({name:na, classname:cla[cl].name, va:[], stat:[], x:[]});
					
					histoplot[e].va[j] = v;
					fl = 1;
					if(depsel.length == 0) depsel.push(e);
				}
			}
		}
		if(fl == 1) depop.push(cla[cl].name);
	}	
	
	if(depop.length > 0) depopsel = depop[0];
}

function histoplotinit()                                  // Initialises histograms
{
	var d, h, v, i, stat;
	
	startloading();
	plotinitfl = 1;
	calcto = setTimeout(function(){ histoplotcalc(0,0);}, 20);
}

function histoplotcalc(d,i)                               // Calculates histograms
{
	var tim, h, v;

	tim = (new Date()).getTime();
	do{
		h = depsel[d]; v = histoplot[h].va[i];
			
		histoplot[h].stat[i] = calcstatistics2(v,1);
		i++; if(i == histoplot[h].va.length){ d++; i = 0;}
	}while(d < depsel.length && (new Date()).getTime()-tim < 1);
	
	if(d < depsel.length){
		percent = Math.floor(100*((d + i/ histoplot[h].va.length)/depsel.length));
		calcto = setTimeout(function(){ histoplotcalc(d,i);}, 20); return;
	}

	switch(depopsel){
	case "Age": axxmin = 0; axxmax = age[age.length-1]*(age.length+1)/age.length; break;
	case "Time": axxmin = infres.tmin; axxmax = infres.tmax+0.01*(infres.tmax-infres.tmin); break;
	default: axxmin = 0; axxmax = histoplot[depsel[0]].va.length; break;
	}

	axymin = large; axymax = -large;
	for(d = 0; d < depsel.length; d++){
		h = depsel[d];
		imax = histoplot[h].va.length;
		for(i = 0; i < imax; i++){
			val = histoplot[h].stat[i].mean; if(val > axymax) axymax = val; if(val < axymin) axymin = val;
			val = histoplot[h].stat[i].CImin; if(val > axymax) axymax = val; if(val < axymin) axymin = val;
			val = histoplot[h].stat[i].CImax; if(val > axymax) axymax = val; if(val < axymin) axymin = val;
			
			switch(depopsel){
			case "Age": 
				if(i == 0) histoplot[h].x[i] = 0.5*age[0]; 
				else{
					if(i == age.length) histoplot[h].x[i] = 0.5*(axxmax+age[i-1]); 
					else histoplot[h].x[i] = 0.5*(age[i]+age[i-1]); 
				}
				break;
				
			case "Time":
				if(i == 0) histoplot[h].x[i] = 0.5*(infres.tmin + time[0]); 
				else{
					if(i == time.length) histoplot[h].x[i] = 0.5*(time[i-1] + infres.tmax); 
					else histoplot[h].x[i] = 0.5*(time[i-1] + time[i]); 
				}
				break;
				
			default: histoplot[h].x[i] = i + 0.5; break; 
			}
		}
	}
	
	if(yaxisauto == 1){
		if(Math.sqrt(axymin*axymin) < 0.2*Math.sqrt(axymax*axymax)) axymin = 0;
		dd = axymax-axymin; axymax += 0.1*(axymax-axymin); if(axymin != 0) axymin -= 0.1*(axymax-axymin);
	}
	else{ axymin = yaxisfixmin; axymax = yaxisfixmax;}
	loading = 0; plotinitfl = 2; buttoninit();
}

function drawhisto()                                      // Draws is histogram
{
	var d, h, i, x, y, ymin, ymax, r, sty;
	
	if(scale == 1){ r = 5; r2 = 5; sty = THICKLINE;}
	else{ r = 10; r2 = 10; sty = VTHICKLINE;}

	cv.clearRect(0,0,graphdx,graphdy);

	for(d = 0; d < depsel.length; d++){
		col = chcol[d];
		h = depsel[d];
		
		for(i = 0; i < histoplot[h].va.length; i++){
			x = Math.floor(graphdx*( histoplot[h].x[i]-axxmin)/(axxmax-axxmin));
			y = Math.floor(graphdy - graphdy*(histoplot[h].stat[i].mean-axymin)/(axymax-axymin));
			if(i > 0 && (depopsel == "Age" || depopsel == "Time")){
				drawline(xlast,ylast,x,y,col,sty,d);
			}
			ymin = Math.floor(graphdy - graphdy*(histoplot[h].stat[i].CImin-axymin)/(axymax-axymin));
			ymax = Math.floor(graphdy - graphdy*(histoplot[h].stat[i].CImax-axymin)/(axymax-axymin));
			
			drawline(x-r,y-r,x+r,y+r,col,sty);
			drawline(x+r,y-r,x-r,y+r,col,sty);
			drawline(x,ymin,x,ymax,col,sty);
			drawline(x-r2,ymin,x+r2,ymin,col,sty);
			drawline(x-r2,ymax,x+r2,ymax,col,sty);
			xlast = x; ylast = y;
		}
	}
	
	switch(depopsel){
	case "Age":
		for(i = 0; i < age.length; i++){
			x = Math.floor(graphdx*(age[i]-axxmin)/(axxmax-axxmin));
			drawline(x,0,x,graphdy,BLACK,sty,2);
		}
		break;
		
	case "Time":
		for(i = 0; i < time.length; i++){
			x = Math.floor(graphdx*(time[i]-axxmin)/(axxmax-axxmin));
			drawline(x,0,x,graphdy,BLACK,sty,2);
		}
		break;
		
	default:
		imax = histoplot[depsel[0]].va.length;
		for(i = 1; i < imax; i++){
			x = Math.floor(graphdx*(i-axxmin)/(axxmax-axxmin));
			drawline(x,0,x,graphdy,BLACK,sty,2);	
		}
		break;
	}
}

function showhisto()                                      // Shows a histogram
{ 
	if(plotinitfl == 0) histoplotinit();	
	if(plotinitfl != 2) return;

	w = setytics();
	if(depopsel == "Age" || depopsel == "Time") setxtics();
	else{
		ntickx = 0;
		cl = 0; while(cl < ncla && cla[cl].name != depopsel) cl++;
		imax = cla[cl].ncomp;
		for(i = 0; i < imax; i++){			
			addcanbutton(cla[cl].comp[i].name,40+w+i*graphdx/imax,40+graphdy+10,graphdx/imax,30,-1,HISTOBUT,-1,-1);
		}		
	}
	
	if(depopsel == "Age" || depopsel == "Time") graphframe(menux+10,80,40+w,60,20,40,w,depopsel,"Parameter value","histo");
	else graphframe(menux+10,80,40+w,60,20,40,w,"","Parameter value","histo");
		
	drawhisto();
	
	addcanbutton("Reload",18,30+graphdy+34,20,20,RELOADBUT3,RELOADBUT3,-1,-1); 
}

	