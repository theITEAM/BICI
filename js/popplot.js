function showpopulations(ppc)                             // Plots graphs showing variation in population number 
{
	var i, j, tx, val, bin=[], wei=[], range = 10, novari, cc = 5, c, stat=[];

	switch(ppc){
	case 0: res = infres; na = infpagename[pagesubsub[INFERENCEPAGE][2]]; break;
	case 1: res = ppcres; na = ""; break;
	case 2: res = simres; na = simpagename[pagesubsub[SIMULATEPAGE][2]]; break;
	}
	
	x = menux+tab; y = 30;

	if(na == "Derived"){
		if(res.dersel.length == 1) st = "Derived quantity "+res.derive[res.dersel[0]];
		else st = "Derived quantities";
	}
	else{
		st = "Populations in "+claselect[res.clasel];
		for(cl = 0; cl < ncla; cl++) if(res.popfilt[cl] != -1) break;
	}
	
	switch(ppc){
	case 0: addbutton(st,x,y,0,0,-1,TITLEBUT,60,-1); break;
	case 1: addbutton(st,x,y,0,0,-1,TITLEBUT,61,-1); break;
	case 2: addbutton(st,x,y,0,0,-1,TITLEBUT,62,-1); break;
	}
	
	poptmin = res.tmin; poptmax = res.tmax; 

	if(xaxisauto == 1){
		axxmin = poptmin; axxmax = poptmax;
		if(axxmax > 0 && axxmin > 0){ if(axxmin < 0.3*axxmax) axxmin = 0;}
		axxmax += 0.001*(axxmax-axxmin);
	}
	else{ axxmin = xaxisfixmin; axxmax = xaxisfixmax;}

	setxtics();
	
	if(plotinitfl == 0 && na == "Derived") popderiveinit();
	
	if(plotinitfl == 0){ popplotinit(ppc); return;}
	if(plotinitfl != 2) return; 	

	switch(res.popview){
	case "Graph":
		if(yaxisauto == 1){ axymin = plotmi; axymax = plotma;}
		else{ axymin = yaxisfixmin; axymax = yaxisfixmax;}	
		w = setytics();

		if(na == "Derived"){ if(res.dersel.length == 1) yna = res.derive[res.dersel[0]]; else yna = "Value";}
		else yna = "Sub-populations";
		
		graphframe(menux+10,80,40+w,60,20,40,w,"Time",yna,"pop");
		
		addcanbutton("Reload",26,30+graphdy+34,20,20,RELOADBUT,RELOADBUT,-1,-1); 
		
		if(page == INFERENCEPAGE && infres.tmax > infres.tdatamax+0.01*(infres.tdatamax-infres.tdatamin)){
			addcanbutton("Forward prediction",40+w+Math.floor(graphdx*(infres.tdatamax-axxmin)/(axxmax-axxmin)),30,0,graphdy,-1,BURNINBUT,-1,-1); 
		}
		
		if(page == INFERENCEPAGE && infres.tmin < infres.tdatamin -0.01*(infres.tdatamax-infres.tdatamin)){
			addcanbutton("Backward prediction",40+w+Math.floor(graphdx*(infres.tdatamin-axxmin)/(axxmax-axxmin)),30,0,graphdy,-1,BURNINBUT,0,-1 ); 
		}
		
		cl = ncla-1; 
		for(j = 0; j < time.length; j++){
			t = time[j];
			if(t >= axxmin && t < axxmax){
				addcanbutton("Timeline transition",40+w+Math.floor(graphdx*(t-axxmin)/(axxmax-axxmin)),35,0,graphdy,-1,TTRANSBUT,0,-1 ); 
			}
		}
		
		drawpopplot();
		break;
	
	case "Model":
		cornx = menux+20;
		addbutton("",cornx,80,modeldx,modeldy+15+modelplaybar,CANVASBUT,CANVASBUT,-1,-1);
		drawmodel(res.clasel,0,0,modeldx,modeldy,"frame","popoutput");
				
		addcanbutton("",modelplaybar,modeldy+15,modelplaydx-modelplaybar,modelplaybar-2,PLAYLINEBUT,PLAYLINEBUT,-1,-1);
		addcanbutton("",0,modeldy+14,modelplaybar,modelplaybar,PLAYBUT,PLAYBUT,-1,-1);
		break;
	}
}

function popderiveinit()                                  // Initialises plots for derived quantities
{
	var ti, num;

	POPX = Math.floor(100*(poptmax-poptmin)/(axxmax-axxmin)); if(POPX < 100) POPX = 100;

	quantname=[]; quantcol=[];
	nquant = res.dersel.length;
	
	quantmean=[]; quantCImin=[]; quantCImax=[];
	for(ti = 0; ti < POPX; ti++){ quantmean[ti]=[]; quantCImin[ti]=[]; quantCImax[ti]=[];}
	
	plotmi = large; plotma = -large;
	for(k = 0; k < nquant; k++){
		quantname[k] = res.derive[res.dersel[k]];
		quantcol[k] = BLACK;
	
		if(res.runpopsel < res.nch){ chmin = res.runpopsel; chmax = res.runpopsel+1; num = 1;}
		else{ chmin = 0; chmax = res.nch; num = res.nch;}
	
		d = res.dersel[k];
		
		i = -1;
		for(ti = 0; ti < POPX; ti++){
			quantmean[ti][k]=0;
			quantCImin[ti][k]=0; 
			quantCImax[ti][k]=0;
			
			//t = poptmin + (poptmax-poptmin)*ti/(POPX-1);
			t = poptmin + (poptmax-poptmin)*(ti+0.5)/POPX;
			while(i < res.DX-1 && res.ch[chmin].derivepl[d][i+1].t <= t) i++;
			
			if(i >= 0 && i < res.DX-1){
				f = t - res.ch[chmin].derivepl[d][i].t; ff = res.ch[chmin].derivepl[d][i+1].t - t;		
				for(ch = chmin; ch < chmax; ch++){	
					quantmean[ti][k] += (res.ch[ch].derivepl[d][i].av*ff + res.ch[ch].derivepl[d][i+1].av*f)/(num*(f+ff));
					quantCImin[ti][k] += (res.ch[ch].derivepl[d][i].min*ff + res.ch[ch].derivepl[d][i+1].min*f)/(num*(f+ff));
					quantCImax[ti][k] += (res.ch[ch].derivepl[d][i].max*ff + res.ch[ch].derivepl[d][i+1].max*f)/(num*(f+ff));
				}
				if(quantCImin[ti][k] < plotmi) plotmi = quantCImin[ti][k];
				if(quantCImax[ti][k] > plotma) plotma = quantCImax[ti][k];
			}
		}
	}
	
	for(d = 0; d < data.length; d++){
		dat = data[d];
		if(dat.variety == "der" && plotdataornot(d) >= 0){
			for(r = 0; r < dat.t.length; r++){
				if(dat.errbar[r].min < plotmi) plotmi = dat.errbar[r].min;
				if(dat.errbar[r].max > plotma) plotma = dat.errbar[r].max;
			}
		}
	}
		
	dplot = 0.05*(plotma-plotmi);
	plotmi -= dplot; plotma += dplot;
	
	plotinitfl = 2;
}

function popplotinit()                                    // Initialises population plot
{
	var k, pop=[];
	
	startloading();

	poppage = page; poppagesub = pagesub[page]; poppagesubsub = pagesubsub[page][pagesub[page]];

	cl = res.clasel;
	
	npop = 0;
	for(k = 0; k < cla[cl].ncomp; k++){
		if(res.popshow[cl][k] == 1 || res.popview == "Model"){
			popfilter[npop]=[];
	
			for(c = 0; c < ncomp; c++){
				fl = 1;
				for(cl2 = 0; cl2 < ncla; cl2++){
					if(cl2 == cl){
						if(compval[c][cl2] != k) fl = 0; 
					}
					else{
						if(res.popfilt[cl2] >= 0){
							if(res.popfilt[cl2] != compval[c][cl2]) fl = 0;
						}
						else{
							if(res.popfilt[cl2] == -2 && res.popfiltcheck[cl2][compval[c][cl2]] == 0) fl = 0;
						}
					}
				}
				popfilter[npop][c] = fl;
			}
			pop.push(k);
			npop++;
		}
	}

	nquant = npop;
	quantname=[];
	for(q = 0; q < nquant; q++){
		k = pop[q];
		quanteqn[q] = k;
		quantname[q] = cla[cl].comp[k].name;
		quantcol[q] = cla[cl].comp[k].col;
	}

	POPX = Math.floor(100*(poptmax-poptmin)/(axxmax-axxmin)); if(POPX < 100) POPX = 100;

	for(ti = 0; ti < POPX; ti++){
		popch[ti]=[];
		quant[ti]=[]; for(k = 0; k < nquant; k++) quant[ti][k]=[];
	}	

	for(c = 0; c <= ncomp; c++){
		cpop[c]=[];
		for(k = 0; k < npop; k++){
			cpop[c][k] = 0;
			if(c != NOTALIVE){
				cpop[c][k] = popfilter[k][c];
			}
		}
	}

	getminmax();
	
	sstart = smin;
	plotinitfl = 1;
    calcto = setTimeout(function(){ popplotcalc();}, 20);
}

function initcompmult()                                   // Initialises compval and compmult
{
	var c, cl;
	ncomp = 1; for(cl = 0; cl < ncla; cl++){ compmult[cl] = ncomp; ncomp *= cla[cl].ncomp;}
	NOTALIVE = ncomp;

	for(c = 0; c < ncomp; c++){
		compval[c]=[]; for(cl = 0; cl < ncla; cl++) compval[c][cl] = Math.floor(c/compmult[cl])%(cla[cl].ncomp);
	}	
}

function getminmax()                                      // Gets the minimum and maximum sample number
{
	var ch;
	if(res.runpopsel < res.nch){ chmin = res.runpopsel; chmax = res.runpopsel+1;}
	else{ chmin = 0; chmax = res.nch;}

	if(res.sampfilt == "All"){
		smin = res.burninev; 
		smax = large;
		for(ch = chmin; ch < chmax; ch++){
			num = res.ch[ch].nsampev;
			if(num < smax) smax = num;
		}
	}
	else{
		var r = getsampfilt();   
		chmin = r.ch; chmax = r.ch+1; 
		smin = r.s; smax = r.s+1; 
	}	
}

function popplotcalc()                                   // Calculates population variation
{
	var ch, s, ti, tim, reschsa, reschsai, c, cf, k;

	if(page != poppage || pagesub[page] != poppagesub || pagesubsub[page][pagesub[page]] != poppagesubsub){ loading = 0; return;}

	tim = (new Date()).getTime();
	do{
		s = sstart; 
		for(ch = chmin; ch < chmax; ch++){	
			if(s >= res.ch[ch].nsampev) alert("samp event out of range");
			reschsa = res.ch[ch].sampev[s];
			for(k = 0; k < npop; k++) popinit[k] = 0;
			for(ti = 0; ti < POPX; ti++){ for(k = 0; k < npop; k++) popch[ti][k] = 0;}

			imax = reschsa.nind;
			for(i = 0; i < imax; i++){
				reschsai = reschsa.ind[i];
			
				c = NOTALIVE; 
				for(e = 0; e < reschsai.nev; e++){
					cf = reschsai.evc[e];
					for(k = 0; k < npop; k++){
						d = cpop[cf][k]-cpop[c][k];
						if(d != 0){
							ti = Math.floor((POPX-1)*(reschsai.evt[e]-poptmin)/(poptmax-poptmin));
							if(ti >= POPX) break;
							if(reschsai.evt[e] <= poptmin) popinit[k] += d;
							else popch[ti][k] += d;
						}
					}
					c = cf;
				}
			}
			
			for(ti = 0; ti < POPX; ti++){
				for(k = 0; k < nquant; k++) quant[ti][k].push(popinit[k]);
				for(k = 0; k < nquant; k++) popinit[k] += popch[ti][k];
			}
		}
		sstart++;
	}while(sstart < smax && (new Date()).getTime() -tim < 100);
	percent = Math.floor(100*((sstart-smin)/(smax-smin)));
	
	if(sstart == smax){
		plotmi = large; plotma = -large; quantmaxitot = -large;
		for(ti = 0; ti < POPX; ti++){
			quantmean[ti]=[]; quantCImin[ti]=[]; quantCImax[ti]=[];
			maxi = 0;
			for(k = 0; k < npop; k++){
				tempCI = quant[ti][k];
				calcCI(quant[ti][k].length);
				quantmean[ti][k] = mean; if(mean > maxi) maxi = mean;
				quantCImin[ti][k] = CImin;
				quantCImax[ti][k] = CImax;
				t = poptmin + (poptmax-poptmin)*ti/(POPX-1);
				if(t > axxmin && t < axxmax){
					if(quantCImin[ti][k] < plotmi) plotmi = quantCImin[ti][k];
					if(quantCImax[ti][k] > plotma) plotma = quantCImax[ti][k];
				}
			}	
			quantmaxi[ti] = maxi; if(maxi > quantmaxitot) quantmaxitot = maxi;
		}
		
		for(d = 0; d < data.length; d++){
			dat = data[d];
			if(dat.variety == "pop" && plotdataornot(d) >= 0){
				for(r = 0; r < dat.t.length; r++){
					if(dat.errbar[r].max > plotma) plotma = dat.errbar[r].max;
				}
			}
		}
		
		plotmi = 0; plotma *= 1.05;
		loading = 0; plotinitfl = 2; buttoninit();
	}
	else{ calcto = setTimeout(function(){ popplotcalc();}, 20);}
}

function playanim()                                      // Plays the model result animation
{
	if(playing == 1){
		playtime = Math.floor(playtimeinit + (getsec()-playstartt)*POPX/4);
		
		if(playtime >= POPX-1){ playtime = POPX-1; playing = 0; buttoninit();}
		else buttonplot();
		
		setTimeout(function(){ playanim();}, 50);
	}
}

function drawpopplot(wid)                                // Draws the population plot
{
	var qst, k, t,x, y, ti;
	
	cv.clearRect(0,0,graphdx,graphdy);
	
	if(wid) cv.lineWidth = wid;
	else cv.lineWidth = 3;
	
	for(k = 0; k < nquant; k++){
		cv.beginPath(); 
		for(ti = 0; ti < POPX; ti++){
			t = poptmin + (poptmax-poptmin)*ti/(POPX-1);
			x = Math.floor(graphdx*(t-axxmin)/(axxmax-axxmin));
			y = Math.floor(graphdy-graphdy*(quantmean[ti][k]-axymin)/(axymax-axymin));
			if(ti == 0) cv.moveTo(x,y);
			else cv.lineTo(x,y);
		}
		
		cv.strokeStyle = quantcol[k];
		
		setdash(k)
		cv.stroke();

		cv.beginPath(); 
		for(ti = 0; ti < POPX; ti++){
			t = poptmin + (poptmax-poptmin)*ti/(POPX-1);
			x = Math.floor(graphdx*(t-axxmin)/(axxmax-axxmin));
			y = Math.floor(graphdy-graphdy*(quantCImax[ti][k]-axymin)/(axymax-axymin));
					
			if(ti == 0) cv.moveTo(x,y);
			else cv.lineTo(x,y);
		}
		
		for(ti = POPX-1; ti >= 0; ti--){
			t = poptmin + (poptmax-poptmin)*ti/(POPX-1);
			x = Math.floor(graphdx*(t-axxmin)/(axxmax-axxmin));
			y = Math.floor(graphdy-graphdy*(quantCImin[ti][k]-axymin)/(axymax-axymin));
			cv.lineTo(x,y);
		}
		cv.globalAlpha = 0.2;
		cv.fillStyle = quantcol[k];
		cv.fill();
		cv.globalAlpha = 1;
	} 	

	setdash(0);

	if(page == INFERENCEPAGE){
		for(d = 0; d < data.length; d++){
			qst = plotdataornot(d);
			if(qst >= 0){
				col = quantcol[qst];
				
				rr = 5;
				dat = data[d];
				switch(dat.variety){
				case "pop": case "der":
					for(r = 0; r < dat.t.length; r++){
						x = Math.floor(graphdx*(dat.t[r]-axxmin)/(axxmax-axxmin));
						y = Math.floor(graphdy-graphdy*(dat.val[r] -axymin)/(axymax-axymin));

						drawline(x-rr,y-rr,x+rr,y+rr,col,THICKLINE);
						drawline(x-rr,y+rr,x+rr,y-rr,col,THICKLINE);

						ymi = Math.floor(graphdy-graphdy*(dat.errbar[r].min-axymin)/(axymax-axymin));
						yma = Math.floor(graphdy-graphdy*(dat.errbar[r].max-axymin)/(axymax-axymin));
						drawline(x,ymi,x,yma,col,THICKLINE);
						drawline(x-rr,ymi,x+rr,ymi,col,THICKLINE);
						drawline(x-rr,yma,x+rr,yma,col,THICKLINE);
					}	
					break;
				}
			}
		}
	}
}

function plotdataornot(d)                                 // Detemines if data should be plotted or not
{
	var dat = data[d], cl, c, q, qst;
			
	switch(dat.variety){
	case "pop":
		if(na != "Populations" && na != "") return -1;

		if(ncla != dat.popcl.length) return -1;
		
		for(cl = 0; cl < ncla; cl++){
			if(res.clasel == cl){
				q = 0;
				if(data[d].popcl[cl] == "All") c = 0;
				else{
					for(c = 0; c < cla[cl].ncomp; c++){
						if(res.popshow[cl][c] == 1){
							if(data[d].popcl[cl] == cla[cl].comp[c].name) break;
							else q++;
						}
					}
				}
				if(c < cla[cl].ncomp) qst = q; else break; 
			}
			else{
				if(res.popfilt[cl] == -1){ if(data[d].popcl[cl] != "All") break;}
				else{ 
					if(res.popfilt[cl] == -2){
						for(c = 0; c < infcla[cl].ncomp; c++) if(res.popfiltcheck[cl][c] != 1) break;
						if(c < infcla[cl].ncomp) break;
					}
					else{
						if(data[d].popcl[cl] != infcla.comp[res.popfilt[cl]].name) break;
					}
				}
			}
		}

		if(cl == ncla) return qst;		
		break;
	
	case "der":
		if(na != "Derived") return -1;
		for(j = 0; j < res.dersel.length; j++){ if(derive[res.dersel[j]].name == dat.name) return j;}
		break;
	}
	return -1;
}

function geterrorbar(mu,sd)                              // Gets the 95% CI for gamma distribution
{
	var al, be, dx = 0.001, sumst=[], loop, loopmax = 500, xmin, xmax, xmi, xma, x, f;
	
	mu = parseFloat(mu); sd = parseFloat(sd);
	
	be = mu/(sd*sd);
	al = mu*be;
	
	xmin = mu-3*sd; if(xmin < 0.01*sd) xmin =  0.01*sd;
	xmax = mu+3*sd;

	sum = 0;
	dx = (xmax-xmin)/loopmax;
	for(loop = 0; loop < loopmax; loop++){
		sumst[loop] = sum;
		x = xmin+dx*loop;
		if(al < 80) sum += Math.pow(x,al-1)*Math.exp(-be*x)*dx;
		else sum += Math.exp(-(x-mu)*(x-mu)/(2*sd*sd));
	}
	loop = 0; while(sumst[loop+1] < 0.05*sum) loop++;
	f = (0.05*sum-sumst[loop])/(sumst[loop+1]-sumst[loop]);
	xmi = xmin+(loop+f)*dx;
	
	while(sumst[loop+1] < 0.95*sum) loop++;
	f = (0.95*sum-sumst[loop])/(sumst[loop+1]-sumst[loop]);
	xma = xmin+(loop+f)*dx;
	return {min:xmi, max:xma};
}

function getpop(w,t)                                     // Gives the population at a particular point in time
{
	var i, cl;

	num = 0;
	for(i = 0; i < w.nind; i++){
		wi = w.ind[i];
		
		c = NOTALIVE;
		e = 0; while(e < wi.nev && wi.evt[e] < t){ c = wi.evc[e]; e++;}
		
		if(c != NOTALIVE){
			for(cl = 0; cl < ncla; cl++){ if(clagendata[cl] != "All" && clagendata[cl] != cla[cl].comp[compval[c][cl]].name) break;}
			if(cl == ncla) num++;
		}
	}
	return num;	
}

function splitname(na)                                    // Splits up a variable name into name and dependencies
{
	j = 0; while(j < na.length && na.substr(j,1) != "_") j++;
	
	if(j == na.length) return {name:na,dep:[]};
	else return {name:na.substr(0,j),dep:na.substring(j+1).split(",")};
}

