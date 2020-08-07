function exportimage(filename)                            // Exports an image (e.g. a graph)
{
	var x1, y1 = 110, x2 = 70, y2 = 60, wid, hei;
	
	var na;
	switch(page){
	case INFERENCEPAGE: na = infpagename[pagesubsub[INFERENCEPAGE][2]]; break;
	case SIMULATEPAGE: na = simpagename[pagesubsub[SIMULATEPAGE][2]]; break;
	}
	
	if(na == "Correlation" && res.varselx == -1){
		outcan = document.createElement('canvas');
		
		outcan.width = 2*ytot;
		outcan.height = 2*ytot;
		outcv = outcan.getContext('2d');
		
		cv = outcv;
		
		drawcortable();
	}	
	else{
		scale = 2;
		
		graphdx *= scale; graphdy *= scale;
		
		w = 0;
		i = 0; while(i < ncanbut && canbuttype[i] != YTICKTRBUT) i++;
		if(i < ncanbut){
			for(i = 0; i < nticky; i++){
				ww = textwidth(ticky[i],"40px times");
				if(ww > w) w = ww;
			}
			x1 = 90+w;
		}	
		else x1 = 90;
		
		wid = x1+x2+graphdx; hei = y1+y2+graphdy;
		
		outcan = document.createElement('canvas');
		outcan.width = wid;
		outcan.height = hei;
		outcv = outcan.getContext('2d');
		
		graphcan.width = graphdx;
		graphcan.height = graphdy;
		
		cv = graphcv;
		cv.globalAlpha = 1;

		switch(na){
		case "Traces": drawtrace(3); break;
		case "Prob. Dist.": case "Trans. Dist.": drawdistribution(5); break;
		case "Correlation": if(res.varselx >= 0) drawscatterplot(3); break;
		case "Populations": case "Pos. Pred. Ch.": case "Derived": drawpopplot(5); break;
		case "Dependency": drawhisto(); break;
		default: alertp("Error code EC44"); break;
		}
		
		cv = outcv;
		fillrect(0,0,wid,hei,WHITE);
		 
		cv.drawImage(graphcan,0,0,graphdx,graphdy,x1,y2,graphdx,graphdy);
		
		x = x1; y = y2; dx = graphdx; dy = graphdy;

		drawline(x,y+dy,x,y-20,BLACK,VTHICKLINE);
		drawarrow(x,y-40,x,y+dy,25,BLACK);
		drawline(x,y+dy,x+dx+40,y+dy,BLACK,VTHICKLINE);
		drawarrow(x+dx+60,y+dy,x,y+dy,25,BLACK);		
		
		i = 0; while(i < ncanbut && canbuttype[i] != XLABELBUT) i++;
		if(i < ncanbut) centertext(canbuttext[i],x1+graphdx/2,hei-20,"55px times",BLACK); 
		
		i = 0; while(i < ncanbut && canbuttype[i] != YLABELBUT) i++;
		if(i < ncanbut) centerplotangletext(canbuttext[i],45,y2+dy/2,Math.PI/2,"55px times",BLACK); 
		//if(i < ncanbut) centerplotangletext(canbuttext[i],45+12,y2+dy/2,Math.PI/2,"70px times",BLACK); 
		
		i = 0; while(i < ncanbut && canbuttype[i] != XTICKBUT) i++;
		if(i < ncanbut){
			for(j = 0; j < ntickx; j++){
				xx = Math.floor(x+dx*(tickx[j]-axxmin)/(axxmax-axxmin));
				centertext(tickx[j],xx,y2+dy+40,"40px times",BLACK); 
				drawline(xx,y2+dy,xx,y2+dy-20,BLACK,THICKLINE); 
			}
		}
		
		i = 0; while(i < ncanbut && canbuttype[i] != YTICKTRBUT) i++;
		if(i < ncanbut){
			for(j = 0; j < nticky; j++){
				yy = Math.floor(y+dy-dy*(ticky[j]-axymin)/(axymax-axymin));
				righttext(ticky[j],x1-5,yy+6,"40px times",BLACK); 
				drawline(x1,yy,x1+20,yy,BLACK,THICKLINE); 
			}
		}
		
		if(na == "Dependency" && depopsel != "Age" && depopsel != "Time"){
			cl = 0; while(cl < ncla && cla[cl].name != depopsel) cl++;
			imax = cla[cl].ncomp;
			for(i = 0; i < imax; i++) centertext(cla[cl].comp[i].name,x+(i+0.5)*dx/imax,y2+dy+60,"44px times",BLACK); 
		}
		
		if(lablist.length > 1){
			//fo = "40px times";
			fo = "60px times";
			cy = 15;
			var i = 0, ist;
			do{
				var xli=[];
				ist = i;
				x = 0;
				while(i < lablist.length){
					dx = 170 + textwidth(lablist[i],fo);
					if(x + dx > graphdx-20 && xli.length > 0) break;
					xli.push(x);
					x += dx;
					i++;
				}
			
				for(ii = ist; ii < i; ii++){
					cx = x1+x2+graphdx-20-x+xli[ii-ist];
					drawline(cx,cy+12,cx+100,cy+12,labcollist[ii],VVTHICKLINE,ii);
					//plottext(lablist[ii],cx+100+10,cy+24,fo,BLACK); 
					plottext(lablist[ii],cx+100+10,cy+24+10,fo,BLACK); 
				}
				cy += 50;
			}while(i < lablist.length);
		}
		
		graphdx /= scale; graphdy /= scale;
		graphcan.width = width;
		graphcan.height = height;
	}
	
	cv = maincv;

	var dataURL = outcan.toDataURL();
	
	if(filename == 0){
		printimage(dataURL);
	}
	else{
		dataURL = dataURL.replace(/^data:image\/\w+;base64,/, "");
	
		var buf = new Buffer(dataURL, 'base64');
		var fs = require('fs');
		fs.writeFile(filename, buf, function (err) {
		   if (err) {
			   console.info("There was an error attempting to save your data.");
			   console.warn(err.message);
			   return;
		   } else if (callback) { callback();}});
	}
		
	loading = 0;
	scale = 1;
}

function printimage(dataURL)
{
	var html  = '<html><head><title></title></head>';
	html += '<body style="width: 100%; padding: 0; margin: 0;">';
	html += '<img onload="window.print(); window.close();" ';
	if(graphdx > Math.sqrt(2)*graphdy) html += 'style="width:95%"';
	else html += 'style="height:95%"';
	html += ' src="' + dataURL + '"  /></body></html>';
	html += '</body></html>';
	
	var printWindow = window.open('', 'to_print', 'left=100, top=100, height=500,width=800');

	printWindow.document.open();
	printWindow.document.write(html);
}

function exportmodel(filename)                             // Exports a text file
{
	var cl, y, dy, yst=[], xsize = 1500, ysize, dyy = 100;
	
	y = 0;
	for(cl = 0; cl < ncla; cl++){
		yst.push(y); 
		if(!(cl >= ncla-2 && cla[cl].ncomp == 1)){
			dy = modframe(cl,xsize,large);
			y += dy+dyy;
		}
	}
	yst.push(y); 
	ysize = y;
	
	outcan = document.createElement('canvas');
    outcan.width = xsize;
    outcan.height = ysize;
    outcv = outcan.getContext('2d');

	cv = outcv;
	cv.clearRect(0,0,xsize,ysize);
	fillrect(0,0,xsize,ysize,WHITE);

	cv = outcv;
	for(cl = 0; cl < ncla; cl++){
		dy = yst[cl+1] - yst[cl] - dyy;
		if(dy > 0){
			dy = modframe(cl,xsize,dy); shy -= (yst[cl]+60)/shfac;
		
			plottext(cla[cl].name+":",20,yst[cl]+60,"bold 60px ariel",BLACK); 
			for(k = 0; k < cla[cl].ntra; k++){
				tr = cla[cl].tra[k];
				if(tr.i != -1 && tr.f != -1) drawtrans(tr,0,1);
			}
	
			for(k = 0; k < cla[cl].ncomp; k++){
				c = cla[cl].comp[k];
				x =(c.x-shx)*shfac; y = 0*yst[cl] + (c.y-shy)*shfac; w = c.w*shfac; h = c.h*shfac; 
				compplot(cl,k,x,y,w,h,COMPBUT,0);
			}
		}
	}
			
	var dataURL = outcan.toDataURL();	
	
	if(filename == 0){ graphdx = xsize; graphdy = ysize; printimage(dataURL);}
	else{
		dataURL = dataURL.replace(/^data:image\/\w+;base64,/, "");
			
		var buf = new Buffer(dataURL, 'base64');
		var fs = require('fs');
		fs.writeFile(filename, buf, function (err) {
		  if (err) {
			   console.info("There was an error attempting to save your data.");
			   console.warn(err.message);
			   return;
			} else if (callback) { callback();}});
	}
	cv = maincv;
	loading = 0;
}

function saveFileAsText()                                 // Saves an exported file
{
	var st, filename, nsam, n, c, nev;
			
	filename = ById("fileToSave").value;

	switch(exporttype){
	case 0:   // save BICI file
		st = save(0);
		break;

	case 10:  // save BICI file (with results)
		st = save(1);
		break;
		
	case 1:   // export trace
		st = "State\t";
		for(v = 0; v < res.nvar; v++){ st += res.varname[v].replace("→","->"); if(v < res.nvar-1) st += "\t"; else st += "\r\n";}
		numsampmax = 0; for(ch = 0; ch < res.nch; ch++) numsampmax += res.ch[ch].nsamp-res.burnin;
		for(s = 0; s < samplenum; s++){
			ss = Math.floor(s*numsampmax/samplenum);
			ch = 0; while(ss >= (res.ch[ch].nsamp-res.burnin)){ ss -= (res.ch[ch].nsamp-res.burnin); ch++;}
			st += s+"\t"; for(v = 0; v < res.nvar; v++){ st += res.ch[ch].varval[v][res.burnin+ss]; if(v < res.nvar-1) st += "\t"; else st += "\r\n";}
		}
		break;
		
	case 2:   // export state data
		st = "Sample\tID\tBirth time\tState changes\r\n";
		
		if(page == INFERENCEPAGE){ res = infres; indi = res.inddata; smin = res.burninev;}
		else{ res = simres; indi = indsim; smin = 0;}
			
		numsampmax = 0; for(ch = 0; ch < res.nch; ch++) numsampmax += res.ch[ch].nsampev-smin;
		for(s = 0; s < samplenum; s++){
			ss = Math.floor(s*numsampmax/samplenum);
			ch = 0; while(ss >= (res.ch[ch].nsampev-smin)){ ss -= (res.ch[ch].nsampev-smin); ch++;}
			
			ws = res.ch[ch].sampev[ss];
			for(i = 0; i < ws.nind; i++){
				st += s+"\t";
				
				if(indi.ind[i]) st += indi.ind[i].id+"\t"; else st += "New Individual "+i+"\t";  
				wsi = ws.ind[i];
				st += wsi.tbirth +"\t";
				for(e = 0; e < wsi.nev; e++){
					if(e > 0) st += " -> ";
					else st += "Enter ";
					c = wsi.evc[e];
					if(c == NOTALIVE) st += "Leave";
					else{ 
						fl = 0; 
						for(cl = 0; cl < ncla; cl++){ 
							if(cl < ncla-2 || cla[cl].ncomp > 1){
								if(fl == 1) st += ","; fl = 1;
								st += cla[cl].comp[compval[c][cl]].name;
							}
						}
					}
					st += ","+wsi.evt[e];
				}
				st += "\r\n";
			}
		}
		break;
		
	case 3: // diagnostics
		st = "";
		for(c = 0; c < infres.nch; c++){
			st += "***************** RUN "+(c+1)+" *********************\r\n\r\n";
			if(infres.diagnostics[c] == undefined) st += "No output yet given.\n\r";
			else{
				tabs = infres.diagnostics[c].split("|");
				for(j = 1; j < tabs.length; j++) st += tabs[j]+"\r\n";	
				st += "\r\n";
			}
		}
		break;
		
	case 4:
		exportimage(filename);
		return;
		
	case 5:
		var na;
		switch(page){
		case INFERENCEPAGE: na = infpagename[pagesubsub[INFERENCEPAGE][2]]; break;
		case SIMULATEPAGE: na = simpagename[pagesubsub[SIMULATEPAGE][2]]; break;
		}

		switch(na){
		case "Traces":
			st = "State\t";
			for(cu = 0; cu < curvevar.length; cu++){
				st += lablist[cu];
				if(cu < curvevar.length - 1) st += "\t"; else st += "\r\n";
			}
			
			imax = 0; for(cu = 0; cu < curvevar.length; cu++){ ch = curvech[cu]; if(res.ch[ch].nsamp > imax) imax = res.ch[ch].nsamp;}
					
			for(i = 0; i < imax; i++){
				st += i + "\t";
				for(cu = 0; cu < curvevar.length; cu++){
					v = curvevar[cu]; ch = curvech[cu];
					if(i < res.ch[ch].nsamp) st += tpre(res.ch[ch].varval[v][i],5); else st += ".";
					if(cu < curvevar.length - 1) st += "\t"; else st += "\r\n";
				}
			}
			break;
		
		case "Prob. Dist.":
		case "Trans. Dist.":
			max = 0; for(cu = 0; cu < curvevar.length; cu++){ for(j = 0; j < JX; j++){ val = Jbin[cu][j]; if(val > max) max = val;}}
			max *= 1.05;
	
			st = "";
			if(na == "Prob. Dist." && curvevar.length == 1) st += lablist[0]+"\tProbability\r\n";
			else{
				if(na == "Prob. Dist.") st += "Parameter Value\t"; else st += "Duration\t"; 
				for(cu = 0; cu < curvevar.length; cu++){
					st += "Prob. "+lablist[cu];
					if(cu < curvevar.length - 1) st += "\t"; else st += "\r\n";
				}
			}
				
			for(i = 0; i < JX; i++){
				st += tpre((axxmin + (axxmax-axxmin)*(i+0.5)/JX),5) + "\t";
				
				for(cu = 0; cu < curvevar.length; cu++){
					v = curvevar[cu]; ch = curvech[cu];			
					st += tpre((Jbin[cu][i]/max),5);
					if(cu < curvevar.length - 1) st += "\t"; else st += "\r\n";
				}
			}
			break;
			
		case "Correlation":
			if(res.varselx >= 0){  // Scatter
				st = res.varname[res.varselx]+"\t"+res.varname[res.varsely]+"\r\n";	
				
		
				for(ch = 0; ch < res.nch; ch++){
					for(k = 0; k < res.scatnum; k++){
						if(scatsamp[k].ch == ch){
							i = scatsamp[k].i;
							st += tpre(res.ch[ch].varval[res.varselx][i],5)  + "\t" +  tpre(res.ch[ch].varval[res.varsely][i],5)  + "\r\n";
						}
					}
					if(ch < res.nch-1 && res.runsel == res.nch) st += "\r\n\r\n";
				}
			}
			else{             // Correlation
				st = "";
				for(v = 0; v < corlist.length; v++) st += "\t"+res.varname[corlist[v]];
				st += "\r\n";
				for(v = 0; v < corlist.length; v++){
					st += res.varname[corlist[v]];
					for(vv = 0; vv < corlist.length; vv++) st += "\t"+cormat[v][vv].toFixed(2);
					st += "\r\n";
				}
			}
			break;
		
		case "Populations": case "Pos. Pred. Ch.": case "Derived":
			st = "Time\t";
			for(k = 0; k < nquant; k++){
				st += quantname[k] + ":mean\t" + quantname[k]  + ":CImin\t"  + quantname[k]  + ":CImax";
				if(k < nquant-1) st += "\t"; else st += "\r\n";
			}
			
			for(ti = 0; ti < POPX; ti++){
				st += tpre((poptmin + (poptmax-poptmin)*ti/(POPX-1)),5)  +"\t";
				for(k = 0; k < nquant; k++){
					st += tpre(quantmean[ti][k],5)+"\t"+tpre(quantCImin[ti][k],5) +"\t"+tpre(quantCImax[ti][k],5) ;
					if(k < nquant-1) st += "\t"; else st += "\r\n";
				}
			}
			break;
	
		case "Statistics":
			st = "Variable\tMean\t95% Credible Interval\tESS\tR̂\r\n";
			for(vv = vmin; vv < res.filtervar[res.filt].length; vv++){
				v = res.filtervar[res.filt][vv];
				st += res.varname[v]+"\t";
				if(vcalced[vv] == 0) st += "---\t---\t---\r\n";
				else{
					st += tpre(varmean[v],5)+"\t"+tpre(varCImin[v],5)+"  \u2014  "+tpre(varCImax[v],5)+"\t";
					if(varESS[v] == 0) st += "---";
					else st += Math.floor(varESS[v]);
					
					st += "\t"+varGR[v]+"\r\n";
				}
			}
			break;	
		
		case "Dependency":
			st = depopsel+"\t";
			for(d = 0; d < depsel.length; d++){
				h = depsel[d];
				st += histoplot[h].name+" mean\t"+ histoplot[h].name+" 95% Credible Interval";
				if(d < depsel.length-1) st += "\t"; else st += "\r\n";
			}
			
			cl = 0; while(cl < ncla && cla[cl].name != depopsel) cl++;
			imax = cla[cl].ncomp;
		
			for(i = 0; i < imax; i++){
				st += cla[cl].comp[i].name+"\t";
				for(d = 0; d < depsel.length; d++){
					h = depsel[d];
					st += tpre(histoplot[h].stat[i].mean,5)+"\t";
					st += tpre(histoplot[h].stat[i].CImin,5);
					st += "  \u2014  ";
					st += tpre(histoplot[h].stat[i].CImax,5);
					if(d < depsel.length-1) st += "\t"; else st += "\r\n";
				}
			}			
			break;
		
		default: alertp("Export not possible"); return;
		}
		break;
		
	case 6:  // save table
		st = ""; for(i = 0; i < ncol; i++){ st += colname[i]; if(i < ncol-1) st += "\t"; else st +="\r\n";}
		for(j = 0; j < nrow; j++){
			for(i = 0; i < ncol; i++){ st += row[j][i]; if(i < ncol-1) st += "\t"; else st +="\r\n";}
		}
		break;
		
	case 7: exportmodel(filename); return;
	
	default: alert("Export not possible"); return;
	}
	
	var fs = require('fs');
	fs.writeFile(filename, "\ufeff"+st, function (err) {
        if (err) {
            console.info("There was an error attempting to save your data.");
		   return;
        } else if (callback) {
            callback();
        }
    });
	
	loading = 0;
}

function exportoptions(x)                                 // Plots options for export
{
	var na, pngfl=0, txtfl = 0,  parfl = 0, evefl = 0, diafl = 0;
	
	if((page == SIMULATEPAGE && pagesub[page] == 2) || (page == INFERENCEPAGE && pagesub[page] == 2)){ 
		if(page == INFERENCEPAGE) na = infpagename[pagesubsub[INFERENCEPAGE][2]];
		else na = simpagename[pagesubsub[SIMULATEPAGE][2]];
		
		diafl = 1;
		
		if(na != "Start"){
			evefl = 1; 
			if(page == INFERENCEPAGE){ parfl = 1; diafl = 1;}
		}
		
		if(na != "Start" && na != "Export" &&  na != "Individuals" 
			&& !(na == "Populations" && res.popview == "Model") && !(na == "Pos. Pred. Ch." && res.popview == "Model") 
			&& !(na == "Individuals" && res.indview == "Model")){
			if(na != "Statistics") pngfl = 1;
			txtfl = 1;
		}
	}
	
	if(page == MODELPAGE) pngfl = 1;
	
	num = 2*pngfl+txtfl+parfl+evefl+diafl;
	if(num == 0) return;

	if(exporton == 0) addbutton("Export",x,0,75,22,EXPORTAC,LOADBUT,-1,3); 
	else{
		xx = x-37;
		addbutton("Export",xx,0,75+37,25+22*num,EXPORTAC2,LOADBUT,-1,3);
		yy = 24; xx += 8;
		
		pt = "Graph";
		switch(na){
		case "Statistics": case "Correlation": if(res.varselx >= 0) pt = "Plot"; else pt = "Table"; break;
		case "Traces": pt = "Plot"; break;
		}
		
		if(pngfl == 1){
			if(page == MODELPAGE){
				addbutton("Print",xx,yy,95,18,EXPORTMINIBUT,EXPORTMINIBUT,6,-1);
				yy += 22;
			
				addbutton(pt+" (png)",xx,yy,95,18,EXPORTMINIBUT,EXPORTMINIBUT,5,-1);
				yy += 22;
			}
			else{
				addbutton("Print",xx,yy,95,18,EXPORTMINIBUT,EXPORTMINIBUT,4,-1);
				yy += 22;
			
				addbutton(pt+" (png)",xx,yy,95,18,EXPORTMINIBUT,EXPORTMINIBUT,0,-1);
				yy += 22;
			}
		}
		
		if(txtfl == 1){
			addbutton(pt+" (txt)",xx,yy,95,18,EXPORTMINIBUT,EXPORTMINIBUT,1,-1);
			yy += 22;
		}
		
		if(evefl == 1){
			addbutton("Events",xx,yy,95,18,EXPORTSTATEAC,EXPORTMINIBUT,1,-1);
			yy += 22;
		}
		
		if(parfl == 1){
			addbutton("Parameters",xx,yy,95,18,EXPORTPARAMAC,EXPORTMINIBUT,1,-1);
			yy += 22;
		}
		
		if(diafl == 1){
			addbutton("Diagnostics",xx,yy,95,18,EXPORTDIAGNOSTICAC,EXPORTMINIBUT,1,-1);
			yy += 22;
		}
	}
}

function saveoptions(x)                                    // Plots options for saving
{
	if(saveon == 0){ addbutton("Save",x,0,75,22,SAVEAC,LOADBUT,-1,1);}
	else{
		x -= 65; 	
		xx = x+8; yy = 24;
		addbutton("Save",x,0,75+65,25+22*2,-1,LOADBUT,-1,1); 
		addbutton("With results (.bici)",xx,yy,95+30,18,EXPORTMINIBUT,EXPORTMINIBUT,2,0); yy += 22;
		addbutton("W/o results (.bici)",xx,yy,95+30,18,EXPORTMINIBUT,EXPORTMINIBUT,3,0);
	}
	return x-90;
}    

function exportstate()                                     // Exports state samples
{
	if(page == SIMULATEPAGE){
		res = simres;
		samplenum = res.ch[0].nsampev;
		exporttype = 2;
		ById("fileToSave").accept = ".txt";
		savesta();
	}
	else{
		res = infres;
		numsampmax = 0; for(ch = 0; ch < res.nch; ch++) numsampmax += res.ch[ch].nsampev-res.burninev;
	
		if(numsampmax > 1000) numsug = 1000; else numsug = numsampmax;
		var ans = prompt("Please enter the number of samples to be exported (1-"+numsampmax+").\nNote, samples are taken at uniform intervals from burnin\n(which is taken to be 20% of the total sample number).",numsug);
	
		if(ans){
			samplenum = Number(ans);
			if(isNaN(samplenum)){ alertp("Not a number!");}
			else{
				if(samplenum < 1 || samplenum > numsampmax) alertp("Not a valid number!");
				else{
					exporttype = 2;
					ById("fileToSave").accept=".txt";
					savesta();
				}
			}
		} 
	}
}

function tpre(num,p)                                       // Outputs a number to a certain precision
{
	var t1, t2;
	t1 = ""+num;
	t2 = Number(num).toPrecision(p);
	if(t2.length < t1.length) return t2; else return t1;
}

function exportparam()                                     // Exports parameter samples
{
	res = infres;
    numsampmax = 0; for(ch = 0; ch < res.nch; ch++) numsampmax += res.ch[ch].nsamp-res.burnin;

	if(numsampmax > 1000) numsug = 1000; else numsug = numsampmax;
	var ans = prompt("Please enter the number of samples to be exported (1-"+numsampmax+").\nNote, samples are taken at uniform intervals from burnin\n(which is taken to be 20% of the total sample number).",numsug);
	
	if(ans) {
		samplenum = Number(ans);
		if(isNaN(samplenum)){ alertp("Not a number!");}
		else{
			if(samplenum < 1 || samplenum > numsampmax) alertp("Not a valid number!");
			else{
				exporttype = 1;
				ById("fileToSave").accept=".txt";
				savesta();
			}
		}
	} 
}
