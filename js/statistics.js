function startcalccor()                                   // Starts calculation of correlation
{
	var vv, vvv, num;
	
	if(corid != undefined){ clearInterval(corid); corid = undefined;}
			
	corlist=[]; for(vv = 0; vv < res.nvar; vv++){ if(res.varcheck[vv] == 1) corlist.push(vv);}
	cormat=[];
	
	loading = 0;
	if(corlist.length == 0) return;
	
	startloading();
	
	cordx = 0;
	for(vv = 0; vv < corlist.length; vv++){
		num = Math.floor(textwidth(res.varname[corlist[vv]],CORFONT));
		if(num > cordx) cordx = num;
	
		cormat[vv]=[];
		for(vvv = 0; vvv < corlist.length; vvv++) cormat[vv][vvv] = "---";
	}
	cordx += 30;
	corx = 0; cory = 0;
	
	corid = setTimeout(function(){ calccor();}, 10);
}

function calccor()
{
	var v1, v2, v1st=[], v2st=[], v1av, v1av2, v2av, v2av2, v12av, de;
	
	v1 = corlist[corx]; v2 = corlist[cory];
	
	if(res.runsel < res.nch){ chmin = res.runsel; chmax = res.runsel+1;}
	else{ chmin = 0; chmax = res.nch;}

	nvst = 0;
	for(ch = chmin; ch < chmax; ch++){
		for(i = res.burnin; i < res.ch[ch].nsamp; i++){
			v1st[nvst] = res.ch[ch].varval[v1][i];
			v2st[nvst] = res.ch[ch].varval[v2][i];
			nvst++;
		}
	}
	
	v1av = 0; v1av2 = 0; v2av = 0; v2av2 = 0; v12av = 0;
	for(i = 0; i < nvst; i++){
		v1av += v1st[i]; v1av2 += v1st[i]*v1st[i];
		v2av += v2st[i]; v2av2 += v2st[i]*v2st[i];
		v12av += v1st[i]*v2st[i];
	}
	v1av /= nvst; v1av2 /= nvst; v2av /= nvst; v2av2 /= nvst; v12av /= nvst;
	
	de = (v1av2-v1av*v1av)*(v2av2-v2av*v2av);
	if(de > 0.000000000000000001){ 	
		cormat[corx][cory] = (v12av - v1av*v2av)/Math.sqrt(de);
		cormat[cory][corx] = cormat[corx][cory];
	}
	
	corx++; if(corx == corlist.length){ cory++; corx = cory;}
	if(cory < corlist.length){
		corid = setTimeout(function(){ calccor();}, 10);
	}
	else loading = 0;
	
	buttoninit();
}

function startcalc()                                      // Starts calculation of statistics
{
	var vv;
	
	startloading();
	vcalc = 0; for(vv = 0; vv < res.filtervar[res.filt].length; vv++) vcalced[vv] = 0;
	setTimeout(function(){ calcstatistics();}, 10);
}

function calcstatistics()                                 // Calculates statistics
{
	var  v;
	
	if(vcalc >= res.filtervar[res.filt].length) return;
	v = res.filtervar[res.filt][vcalc];
	percent = Math.floor((100*vcalc)/res.filtervar[res.filt].length);

	var stat = calcstatistics2(v,1);
	
	varGR[v] = stat.GR; varESS[v] = stat.ESS; 
	varCImin[v] = stat.CImin; varCImax[v] = stat.CImax;
	varmean[v] = stat.mean;
	
	vcalced[vcalc] = 1;
	vcalc++;

	if(vcalc < res.filtervar[res.filt].length && 
	   ((page == INFERENCEPAGE && pagesub[page] == 2 && infpagename[pagesubsub[page][2]] == "Statistics") ||
	   (page == SIMULATEPAGE && pagesub[page] == 2 && simpagename[pagesubsub[page][2]] == "Statistics"))){
		setTimeout(function(){ calcstatistics();},80);
	}
	else loading = 0;
	buttoninit();
}

function calcstatistics2(v,CIon)                          // Calculates statistics
{
	var ch, nsampst=[], n, i, j, min, max,nmin, nmax, stat={ CImin:0, CImax:0, mean:0, GR:"---", ESS:"---"}, fl, burn;
	
	for(ch = 0; ch < res.nch; ch++) nsampst[ch] = res.ch[ch].nsamp;
	if(page == INFERENCEPAGE || CIon == -1){ fl = 1; burn = res.burnin;} else{ fl = 0; burn = 0;}

	n = 0;
	if(res.runsel2 < res.nch){
		ch = res.runsel2;
		for(i = burn; i < nsampst[ch]; i++){ tempCI[n] = res.ch[ch].varval[v][i]; n++;} 
	}
	else{
		for(ch = 0; ch < res.nch; ch++){
			for(i = burn; i < nsampst[ch]; i++){ tempCI[n] = res.ch[ch].varval[v][i]; n++;} 
		}

		if(res.nch > 1 && fl == 1){     // calculates the Gelman-Rubin diagnostics
			var ma = large, N, M, W, B;
			for(ch = 0; ch < res.nch; ch++) if(nsampst[ch] < ma) ma = nsampst[ch];
			
			if(ma > burn+1){
				N = ma-burn; M = infres.nch;
				var mu=[], vari=[], muav = 0;
				for(ch = 0; ch < M; ch++){ 
					valav = 0; for(i = burn; i < ma; i++) valav += res.ch[ch].varval[v][i]/N;
					varr = 0; for(i = burn; i < ma; i++) varr += (res.ch[ch].varval[v][i]-valav)*(res.ch[ch].varval[v][i]-valav)/(N-1);
					mu[ch] = valav;
					vari[ch] = varr;
					muav += mu[ch]/M;
				}
				W = 0; for(ch = 0; ch < M; ch++) W += vari[ch]/M;
				B = 0; for(ch = 0; ch < M; ch++) B += (mu[ch]-muav)*(mu[ch]-muav)*N/(M-1);
				if(W > 0.0000000000001) stat.GR = tpre(Math.sqrt(((1-1.0/N)*W+B/N)/W),4);
			}
		}
	}

	if(CIon == 1) calcCI(n); else findmeansd(n);
	stat.CImin = CImin; stat.CImax = CImax;

	stat.mean = mean;

	if(fl == 1){
		if(n == 0) stat.ESS = 0;
		if(sd > 0){
			for(i = 0; i < n; i++) tempCI[i] = (tempCI[i]-mean)/sd;
			
			sum = 1;
			for(d = 0; d < n/2; d++){             // calculates the effective sample size
				a = 0; for(i = 0; i < n-d; i++) a += tempCI[i]*tempCI[i+d]; 
				cor = a/(n-d);
				if(cor < 0) break;
				sum += 0.5*cor;			
			}
			stat.ESS = Math.floor(n/sum);
		}
	}
	
	return stat;
}

function calcCI(n)                                        // Calculates confidence intervals 
{
	var i, f;

	if(n == 0){ alertp("Error code EC45");}
	
	if(n == 1){
		CImin = tempCI[0]; CImax = tempCI[0];
	}
	else{
		var list=[];
		for(i = 0; i < n; i++) list[i] = tempCI[i];
		list.sort((a, b) => a - b);

		i = Math.floor((n-1)*0.025); f =(n-1)*0.025 - i;
		CImin = list[i]*(1-f) + list[i+1]*f;
		
		i = Math.floor((n-1)*0.975); f = (n-1)*0.975 - i;
		CImax = list[i]*(1-f) + list[i+1]*f;
	}

	findmeansd(n);	
}

function findmeansd(n)                                    // Calculates the mean and standard deviation
{
	var i, val, mean2, va;
	mean = 0; mean2 = 0;
	for(i = 0; i < n; i++){ val = tempCI[i]; mean += val; mean2 += val*val;}
	mean /= n;
	va = mean2/n - mean*mean; if(va < 0.000000001) va = 0;
	sd = Math.sqrt(va);
}
		
function showstatistics()                                 // Shows statistics
{
	var x1, x2, x3, x4, x5, dd, hei, rowmax, dy = 28;

	hei = height-185;
	rowmax = Math.floor(hei/dy);  
	hei = dy*rowmax;
	
	tableyfrac = rowmax/res.filtervar[res.filt].length;
	
	y = 95;
	vmin = Math.floor(1.001*tableyfr*res.filtervar[res.filt].length); 
	vmax = vmin+rowmax; if(vmax > res.filtervar[res.filt].length) vmax = res.filtervar[res.filt].length;

	if(width < 1000){ x1 = 40; x2 = 150; x3 = 270; x4 = 470; x5 = 560;}
	else{ x1 = 40; x2 = 160; x3 = 290; x4 = 500; x5 = 620;}

	dd = 3;
	addbutton("Variable",menux+x1,y,50,15,-1,STATBUT2,-1,-1);            
	if(page == INFERENCEPAGE){	
		addbutton("Mean",menux+x2,y,50,15,-1,STATBUT2,-1,-1);            
		addbutton("95% CI",menux+x3,y,50,15,-1,STATBUT2,-1,-1);  
		addbutton("ESS",menux+x4,y,50,15,-1,STATBUT2,-1,-1);            
		addbutton("R̂",menux+x5,y,50,15,-1,STATBUT2,-1,-1);            
	}
	else{
		if(res.ch[0].nsampread == 1) addbutton("Value",menux+x2,y,50,15,-1,STATBUT2,-1,-1); 
		else{
			addbutton("Mean",menux+x2,y,50,15,-1,STATBUT2,-1,-1);            
			addbutton("95% Credible Interval",menux+x3,y,50,15,-1,STATBUT2,-1,-1);  
		}
	}
	
	y += 40;
	for(vv = vmin; vv < vmax; vv++){
		v = res.filtervar[res.filt][vv];

		addbutton(res.varname[v],menux+x1+dd,y,50,15,-1,STATBUT,0,-1);
		if(vcalced[vv] == 0){
			addbutton("---",menux+x2+dd,y,50,15,-1,STATBUT3,-1,-1);
			if(page==INFERENCEPAGE || res.ch[0].nsampread > 1) addbutton("---",menux+x3+dd,y,50,15,-1,STATBUT3,-1,-1);
			if(page==INFERENCEPAGE) addbutton("---",menux+x4+dd,y,50,15,-1,STATBUT3,-1,-1);
			if(page==INFERENCEPAGE) addbutton("---",menux+x5+dd,y,50,15,-1,STATBUT3,-1,-1);
		}
		else{
			addbutton(tpre(varmean[v],3),menux+x2+dd,y,50,15,-1,STATBUT3,-1,-1);
			if(page==INFERENCEPAGE || res.ch[0].nsampread > 1){
				addbutton(tpre(varCImin[v],3)+"  \u2014  "+tpre(varCImax[v],3),menux+x3+dd,y,50,15,-1,STATBUT3,-1,-1);
			}
			if(page==INFERENCEPAGE) addbutton(varESS[v],menux+x4+dd,y,50,15,-1,STATBUT3,-1,-1);
			if(page==INFERENCEPAGE) addbutton(varGR[v],menux+x5+dd,y,50,15,-1,STATBUT3,-1,-1);
		}
		y += dy;
	}
	
	if(tableyfrac < 1) addbutton("",menux+680,120,13,hei,SLIDEAC,YSLIDEBUT,-1,-1);		
	
	if(page==INFERENCEPAGE){
		addbutton("ESS gives the number of independent posterior samples (ESS > 200 for good mixing).",menux+65,height-60,680,0,-1,SMTEXTBUT,1,-1); 
		
		addbutton("R̂ is the Gelman-Rubin diagnostic (0.9 < R̂ < 1.1 for good mixing).",menux+65,height-35,680,0,-1,SMTEXTBUT,1,-1); 
	
		addbutton("Reload",menux+27,height-48,20,20,RELOADSTATBUT,RELOADSTATBUT,-1,-1);
	}		
}

function calcbayesfac(va)                                 // Calculates the Bayes Factor
{
	var dj=0.1;
	
	v = curvevar[0]; ch = curvech[0];

	if(ch < nchrun){ chmin = ch; chmax = ch+1;}
	else{ chmin = 0; chmax = nchrun;}
	
	vpri = findparam(res.varname[v]);
	if(vpri == infres.param.length){ alertp("This parameter does not have a prior ditribution"); return;}
	
	if(infres.param[vpri].prior == "Flat"){
		valmin = parseFloat(infres.param[vpri].val[0]);
		valmax = parseFloat(infres.param[vpri].val[1]);
	}
	else{
		valmin = -large; valmax = large;
	}
	if(va < valmin || va > valmax){ alertp("The value must be within the range of the prior"); return;}
	
	j = JX*(va-axxmin)/(axxmax-axxmin);

	num = 0;
	for(ch = chmin; ch < chmax; ch++){ 
		for(i = res.burnin; i < res.ch[ch].nsamp; i++){
			val = JX*(res.ch[ch].varval[v][i]-axxmin)/(axxmax-axxmin);
			d = val-j;
			if(d > 0){ if(d < kde) num += 1-d/kde;}
			else{ if(d > -kde) num += 1+d/kde;}
		}
	}
	
	fac = 1;     
	val = JX*(valmin-axxmin)/(axxmax-axxmin);
	d = j-val; if(d < kde) fac -= 0.5*(kde-d)*(kde-d)/(kde*kde);
	val = JX*(valmax-axxmin)/(axxmax-axxmin);
	d = val-j; if(d < kde) fac -= 0.5*(kde-d)*(kde-d)/(kde*kde);		
	num /= fac;
	
	area = 0;
	for(j = 0; j < JX; j += dj){
		nu = 0;
		for(ch = chmin; ch < chmax; ch++){ 
			for(i = res.burnin; i < res.ch[ch].nsamp; i++){
				val = JX*(res.ch[ch].varval[v][i]-axxmin)/(axxmax-axxmin);
				d = val-j;
				if(d > 0){ if(d < kde) nu += 1-d/kde;}
				else{ if(d > -kde) nu += 1+d/kde;}
			}
		}
		fac = 1;
		val = JX*(valmin-axxmin)/(axxmax-axxmin);
		if(j < val) nu = 0;
		d = j-val; if(d < kde) fac -= 0.5*(kde-d)*(kde-d)/(kde*kde);
		val = JX*(valmax-axxmin)/(axxmax-axxmin);
		if(j > val) nu = 0;
		d = val-j; if(d < kde) fac -= 0.5*(kde-d)*(kde-d)/(kde*kde);		
		nu /= fac;
		
		area += nu*dj;
	}
	BF = num*JX/area;

	if(BF < 1){
		if(BF == 0){
			BFtext = "A very large Bayes Factor in favour of model with variable "+res.varname[res.varsel[0]]+" over a model with "+res.varname[res.varsel[0]]+" = "+va;
		}
		else BFtext = "Bayes Factor of "+(1.0/BF).toFixed(2)+" in favour of model with variable "+res.varname[res.varsel[0]]+" over a model with "+res.varname[res.varsel[0]]+" = "+va;
	}
	else BFtext = "Bayes Factor of "+BF.toFixed(2)+" in favour of model with "+res.varname[res.varsel[0]]+" = "+va+" over a model with variable "+res.varname[res.varsel[0]];
	helptype = 77;
}

function indplotstclear()                                 // Triggers when all the individual plots need to be recalculated
{
	var i, cl;
	for(i = 0; i < indplotst.length; i++){ 
		if(indplotst[i] != undefined){ for(cl = 0; cl < ncla; cl++) indplotst[i][cl] = undefined;}
	}
}

function calcindprob(i,cl,x,y,dx)                         // Calculates timelines for individuals 
{	
	var indsto;
	
	if(indplotst[i] == undefined || indplotst[i][cl] == undefined){
		var tl=[], timax, ti, prob=[], r=[], g=[], b=[], tlast, div, ndiv;
		var ev, evallt=[];
		
		if(page == INFERENCEPAGE){
			if(pagesub[page] == 0) indi = inddata.ind[i]; 
			else indi = infres.inddata.ind[i];
		}			
		else indi = indsim.ind[i];

		var small = 0.00001;
		evallt.push(axxmin);
		if(indi){
			for(cl2 = 0; cl2 < ncla; cl2++){
				ev = indi.cl[cl2].ev;
				for(e = 0; e < ev.length; e++){
					t = ev[e].t;
					if(t > axxmin && t < axxmax){
						if(cl2 == res.agecl && (e == 0 || e == ev.length-1)){ 
							tt = t-small; if(tt > axxmin && tt < axxmax) evallt.push(tt); 
							tt = t+small; if(tt > axxmin && tt < axxmax) evallt.push(tt);
						}
						else evallt.push(t);
					}
				}
			}
		}

		evallt.push(axxmax);
		evallt.sort(function(a, b){return a - b});
		
		tlast = evallt[0]; 
		tl.push(tlast);
		for(e = 0; e < evallt.length; e++){
			tnew = evallt[e];
			if(tnew != tlast){
				ndiv = Math.floor(1+400*(tnew-tlast)/(axxmax-axxmin));
				for(div = 0; div < ndiv; div++) tl.push(tlast + (tnew-tlast)*(div+1)/ndiv);
				tlast = tnew;
			}
		}
		
		timax = tl.length;
		
		kmax = cla[cl].ncomp;
		for(ti = 0; ti < timax; ti++){
			prob[ti]=[];
			for(k = 0; k <= kmax; k++) prob[ti][k] = 0;
		}
		
		getminmax();

		fl = 0;
		for(ch = chmin; ch < chmax; ch++){
			for(s = smin; s < smax; s++){
				if(s < res.ch[ch].nsampev){
					if(i < res.ch[ch].sampev[s].nind){
						ress = res.ch[ch].sampev[s].ind[i];
				
						c = NOTALIVE;
						e = 0;
						for(ti = 0; ti < timax; ti++){
							t = tl[ti];
							while(e < ress.nev && ress.evt[e] < t){ c = ress.evc[e]; e++;}
							if(c == NOTALIVE) prob[ti][kmax]++; else prob[ti][compval[c][cl]]++;
							fl = 1;
						}
					}
				}
			}
		}	
		if(fl == 0) return;
		
		for(k = 0; k < kmax; k++){
			co = cla[cl].comp[k].col;
			coconv = hexToRgb(co);
			r[k] = coconv.r; g[k] = coconv.g; b[k] = coconv.b;
		}
		r[kmax] = 255; g[kmax] = 255; b[kmax] = 255;

		gr = cv.createLinearGradient(x,0,x+dx+2,0);
	 
		for(ti = 0; ti < timax; ti++){
			rsum = 0; gsum = 0; bsum = 0; sum = 0;
			for(k = 0; k <= kmax; k++){
				rsum += prob[ti][k]*r[k]; gsum += prob[ti][k]*g[k]; bsum += prob[ti][k]*b[k];
				sum += prob[ti][k];
			}
			if(sum == 0) sum += 0.00001;

			gr.addColorStop((tl[ti]-axxmin)/(axxmax-axxmin),"rgb("+Math.floor(rsum/sum)+","+Math.floor(gsum/sum)+","+Math.floor(bsum/sum) +")"); 
		}
		x1 = Math.floor(dx*(tl[0]-axxmin)/(axxmax-axxmin)); x2 = Math.floor(dx*(tl[timax-1]-axxmin)/(axxmax-axxmin));
		indsto={x1:x1, x2:x2, gr:gr, lab:[], line:[]};
		
		// works out labels
		for(ti = 0; ti < timax; ti++){
			ma = 0; for(k = 0; k <= kmax; k++) if(prob[ti][k] >= ma){ ma = prob[ti][k]; mak = k;}
			if(ti == 0){ tist = ti; makst = mak;}
			else{
				if(mak != makst || ti == timax-1){
					if(makst < kmax && makst < cla[cl].ncomp){
						st = cla[cl].comp[makst].name;
						w = textwidth(st,"15px Times");
						xx = Math.floor(x+dx*(tl[ti]-axxmin)/(axxmax-axxmin)); 
						if(xx != x+dx) indsto.line.push(xx);
						if(w < dx*(tl[ti]-tl[tist])/(axxmax-axxmin)){
							indsto.lab.push({text:st, x:dx*((tl[tist]+tl[ti])/2-axxmin)/(axxmax-axxmin)});
						}
					}
					tist = ti; makst = mak;
				}
			}
		}
		if(indplotst[i] == undefined) indplotst[i]=[];
		indplotst[i][cl] = indsto;
	}
	else indsto = indplotst[i][cl];

	cv.fillStyle = indsto.gr; cv.fillRect(x+indsto.x1,y-10,indsto.x2-indsto.x1+2,20);

	for(j = 0; j < indsto.line.length; j++){ xx = indsto.line[j]; drawline(xx,y-10,xx,y+10,BLACK,NORMLINE);}
	for(j = 0; j < indsto.lab.length; j++) centertext(indsto.lab[j].text,x+indsto.lab[j].x,y+5,"15px Times",BLACK,0);
}

function convergencetestinit()                           // Initialises a convergence test
{
	vconvtest = 0;
	convergencetest();
}

function convergencetest()                               // Performs a convergence test
{
	var timebeg = (new Date()).getTime();
	if(infres.result != 1) return;

	while(vconvtest < infres.nvar && (new Date()).getTime() < timebeg + 80){	
		restemp = res; 
		res = infres; 
		runsel2temp = res.runsel2; 
		res.runsel2 = res.nch;
		var stat = calcstatistics2(vconvtest,-1);
		res.runsel2 = runsel2temp;
		res = restemp;
		if((stat.ESS != "---" && stat.ESS < ESSmin) || (stat.GR != "---" && stat.GR > GRmax)){
			convid = setTimeout(function(){ convergencetestinit();}, Math.floor(0.5*((new Date()).getTime()-timestart)));
			return;
		}
		vconvtest++;
	}
	if(vconvtest < infres.nvar && runtype == "inf") setTimeout(function(){ convergencetest();},80);
	else{ stopinference(); infres.result = 2;  helptype = 102; buttoninit();}
}
