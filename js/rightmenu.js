function rightmenu(x,y,type)                              // Draws the right menu
{
	var x1, x2, cl, dy1 = 20, dy2 = 35, yrightbeg;
	x1 = width-140; x2 = width-25;
	
	ybeg = y; ymax = height-80;
	rightheight = height - ybeg - 60;
	
	if(type != "ind" && type != "Statistics") y -= Math.floor(tableyfr*ytotright);
	
	yrightbeg = y;

	switch(type){
	case "seltrans":		
		addbuttonlim("Classification:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,65); y += dy1;
		if(y >= ybeg && y < ymax){
			gdropinfo.push({val:claselect[transcl], x:x, y:y, dx:selbutdx, dy:20, style:4, options:claselect, click:"tranclasel"});		
		}
		y += dy2;
		break;
		
	case "Correlation":
		if(res.varselx == -1){
			y = varchecklistplot(y,RADIOCHOOSE);
			y += 10;
		}
		else{
			y += 10;
			addbutton("Back",x+8,y,84,25,CORBACKAC,SMBACKBUT,-1,-1);
			y += 50;
			
			addbutton("KDE",x,y,110,22,RADIOBUT,RADIOBUT,"KDE",RADIOSCATTY); y += 25;
			addbutton("Scatter",x,y,110,22,RADIOBUT,RADIOBUT,"scatter",RADIOSCATTY); y += 25;
			y += 10;
			
			if(res.scattype == "scatter"){
				addbuttonlim("Number:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,82); y += dy1;
				if(y >= ybeg && y < ymax){
					gdropinfo.push({val:res.scatnum, x:x, y:y, dx:selbutdx, dy:20, style:5, options:scatnumop , click:"scatnum"});
				}
				y +=dy2;
			}
		}
		
		if(res.nch > 1){
			addbuttonlim("Run:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,105); y += dy1;
			if(y >= ybeg && y < ymax){
				gdropinfo.push({val:res.runlab[res.runsel], x:x, y:y, dx:selbutdx, dy:20, style:5, options:res.runlab, click:"runlab"});
			}
			y +=dy2;
		}
		break;
		
	case "Statistics":
		addbuttonlim("Variable type:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,64); y += 20;	
		if(y >= ybeg && y < ymax){
			gdropinfo.push({val:res.filter[res.filt], x:x, y:y, dx:selbutdx, dy:20, style:4, options:res.filter, click:"filter"});
		}
		y += dy2;
		 
		if(res.nch > 1){
			addbuttonlim("Run:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,107); y += dy1;
			if(y >= ybeg && y < ymax){
				gdropinfo.push({val:res.runlab2[res.runsel2], x:x, y:y, dx:selbutdx, dy:20, style:5, options:res.runlab2, click:"runlab2"});
			}
			y += dy2;
		}		
		break;
		
	case "der":
		y = derlistplot(y,RADIODER);
		if(res.nch > 1){
			y += 10;
			addbuttonlim("Run:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,105); y += dy1;
			if(y >= ybeg && y < ymax){
				gdropinfo.push({val:res.runlab2[res.runsel2], x:x, y:y, dx:selbutdx, dy:20, style:5, options:res.runlab2, click:"runlab2"});
			}
			y += dy2;
		}	
		break;

	case "Traces":
	case "Prob. Dist.":
		y = varlistplot(y,RADIOVAR);

		if(type == "Prob. Dist."){
			y += 10;
			addbuttonlim("Graph:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,210); y += 25;	
			addbuttonlim("Mean",x1,y-5,checkdx,15,CHECKBUT6,CHECKBUT6,-1,-1); y += 22;
			addbuttonlim("95% CI",x1,y-5,checkdx,15,CHECKBUT5,CHECKBUT5,-1,-1); y += 22;
		}
		
		if(res.nch > 1){
			y += 10;
			
			if(type == "Prob. Dist.") addbuttonlim("Run:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,105);
			else addbuttonlim("Run:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,106);
			y += dy1;
			if(y >= ybeg && y < ymax){
				if(type == "Prob. Dist."){
					gdropinfo.push({val:res.runlab[res.runsel], x:x, y:y, dx:selbutdx, dy:20, style:5, options:res.runlab, click:"runlab"});
				}
				else{
					gdropinfo.push({val:res.runlab3[res.runsel3], x:x, y:y, dx:selbutdx, dy:20, style:5, options:res.runlab3, click:"runlab3"});
				}
			}
			y += dy2;
		}	
		
		y += 20;
		if(type == "Prob. Dist."){
			if(res.filter[res.filt] != "Trans." && res.filter[res.filt] != "Init. Prob." && res.filter[res.filt] != "Misc."){
				if(curvevar.length == 1){
					addbuttonlim("Bayes factor:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,76); y += dy1;
				
					addbuttonlim("Calculate",x,y,selbutdx,20,BAYESBUT,BAYESBUT,-1,-1);
					y += dy2;
				}
			}
		}
		break;
	
	case "Dependency":
		y = deplistplot(y);     
		
		if(res.nch > 1){
			y += 10;
			addbuttonlim("Run:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,106); y += dy1;
			if(y >= ybeg && y < ymax){
				gdropinfo.push({val:res.runlab2[res.runsel2], x:x, y:y, dx:selbutdx, dy:20, style:5, options:res.runlab2, click:"runlab2"});
			}
			y += dy2;
		}
		break;
		
	case "Trans. Dist.":
		addbuttonlim("Transition:",width-rwid,y,0,0,-1,SMALLTEXTBUT,BLACK,78); y += 20;	
		for(v = 0; v < transpos.length; v++){
			addbuttonlim(transpos[v].name,x1,y,radiodx,15,RADIOBUT,RADIOBUT,v,RADIOTRANS);
			y += 22;
		}
		y += dy2-20;
	
		y = rightotherclass(x,y,x1,dy1,dy2,cla[transpos[transsel[0]].cl].name,"");
			
		if(res.nch > 1){
			addbuttonlim("Run:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,105); y += dy1;
			if(y >= ybeg && y < ymax){
				gdropinfo.push({val:res.runlab[res.runsel], x:x, y:y, dx:selbutdx, dy:20, style:5, options:res.runlab, click:"runlab"});
			}
			y += dy2;
		}	
	
		addbuttonlim("Individual:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,79); y += dy1;
		if(y >= ybeg && y < ymax){
			var pos=[]; pos.push("All"); pos.push("Checkbox");  
			gdropinfo.push({val:transindfilt, x:x, y:y, dx:selbutdx, dy:20, style:4, options:pos, click:"transindfilt"});
		}
				
		if(transindfilt == "Checkbox"){
			y += 26;
			addbuttonlim("Check All",x+3,y,55,15,CHECKALLBUT,CHECKALLBUT,0,-1);
			addbuttonlim("Uncheck All",x+3+55,y,58,15,CHECKALLBUT,CHECKALLBUT,0,-1);
			y += 20;
			
			for(v = 0; v < res.indlist.length; v++){
				addbuttonlim(res.indlist[v],x1,y-5,checkdx,15,CHECKBUT4,CHECKBUT4,v,-1); y += 22;
			}
			y -= 30;
		}
		y += dy2;
		
		break;
		
	case "Manualpop": case "loadinitpop":
		if(type == "Manualpop"){
			addbuttonlim("Defined in:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,67); y += dy1;
			if(y >= ybeg && y < ymax){
				gdropinfo.push({val:cla[simpopinitset].name, x:x, y:y, dx:selbutdx, dy:20, style:7, options:simpopinitlist, click:"simpopinit"}); 
			}
			y += dy2;
		}
		
		addbuttonlim("View:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,65); y += dy1;
				
		if(y >= ybeg && y < ymax){
			gdropinfo.push({val:cla[viewpopinitset].name, x:x, y:y, dx:selbutdx, dy:20, style:4, options:simpopinitlist, click:"viewpopinit"});	
		}
		y += dy2;
		
		addbuttonlim("Style:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,68); y += 30;
		addbuttonlim("Colour",x+15,y-5,radiodx-20,15,RADIOBUT,RADIOBUT,0,RADIOSTYLE); y += 22;
		addbuttonlim("Greyscale",x+15,y-5,radiodx-20,15,RADIOBUT,RADIOBUT,1,RADIOSTYLE); y += 22;
		
		addbuttonlim("Label:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,69); y += 30;
		addbuttonlim("Name",x+15,y-5,radiodx-20,15,RADIOBUT,RADIOBUT,0,RADIONAMLAB); y += 22;
		addbuttonlim("No name",x+15,y-5,radiodx-20,15,RADIOBUT,RADIOBUT,1,RADIONAMLAB); y += 22;
		break;
		
	case "pop": case "ind":
		if(type == "pop"){
			addbuttonlim("View:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,63); y += dy1;
			if(y >= ybeg && y < ymax){
				gdropinfo.push({val:res.popview, x:x, y:y, dx:selbutdx, dy:20, style:7, options:popviewpos, click:"popview"});
			}
			y += dy2;
		}
		
		if(type == "ind"){
			if(page == INFERENCEPAGE && pagesub[page] == 0){ indv = "Timeline"; y += 30;}
			else{
				indv = res.indview;
				addbuttonlim("View:",width-rwid,y,0,0,-1,SMALLTEXTBUT,BLACK,70); y += dy1;
				if(y >= ybeg && y < ymax){
					gdropinfo.push({val:res.indview, x:x, y:y, dx:selbutdx, dy:20, style:7, options:indviewpos, click:"indview"});
				}
				y += dy2;
			}
		
			y += 8;

			switch(indv){
			case "Timeline":
				y -= 10;
				addbuttonlim("Timelines:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,211); y += 26;
				for(cl = 0; cl < ncla; cl++){ addbuttonlim(cla[cl].name,x1,y-5,checkdx,15,CHECKBUT2,CHECKBUT2,cl,-1); y += 22;}
				break;
				
			case "Model":
				addbuttonlim("Individual:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,71); y += dy1;
				
				if(y >= ybeg && y < ymax){
					gdropinfo.push({val:res.indlist[res.indfilt], x:x, y:y, dx:selbutdx, dy:20, style:4, options:res.indlist, click:"indfilt"});
				}
				y += dy2;
				break;
			}
		}
		
		if(!(page == INFERENCEPAGE && pagesub[page] == 0)){
			if(!(type == "ind" && res.indview == "Timeline")){
				addbuttonlim("Classification:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,65); y += dy1;
				if(y >= ybeg && y < ymax){
					gdropinfo.push({val:claselect[res.clasel], x:x, y:y, dx:selbutdx, dy:20, style:4, options:claselect, click:"clasel"});		
				}
			}
		}
		
		if(type == "pop" && res.popview == "Graph" && cla[res.clasel].ncomp > 1){	
			cl = res.clasel;
			y += 26;
			addbuttonlim("Check All",x+3,y,55,15,CHECKALLBUT,CHECKALLBUT,2,cl);
			addbuttonlim("Uncheck All",x+3+55,y,58,15,CHECKALLBUT,CHECKALLBUT,2,cl);
			y += 20;
			
			for(v = 0; v < cla[cl].ncomp; v++){ addbuttonlim(cla[cl].comp[v].name,x1,y-5,checkdx,15,CHECKBUT,CHECKBUT,cl,v); y += 22;}
			y -= 30;
		}
		y += dy2;
	

		if(type == "pop"){
			y = rightotherclass(x,y,x1,dy1,dy2,claselect[res.clasel],"Time");
		}
		
		if(res.nch > 1){
			addbuttonlim("Run:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,105); y += dy1;
			
			if(y >= ybeg && y < ymax){
				gdropinfo.push({val:res.runpoplab[res.runpopsel], x:x, y:y, dx:selbutdx, dy:20, style:5, options:res.runpoplab, click:"runpoplab"});
			}
			y += dy2;
		}	
		
		if(!(page == INFERENCEPAGE && pagesub[page] == 0) && res.ch[0].nsampev > 1){
			var s; 
			var ops=[];
			ops.push("All");
			if(res.runpopsel < res.nch){
				for(s = 0; s < res.ch[res.runpopsel].nsampev; s++) ops.push(s+1);	
			}
			else{
				for(ch = 0; ch < res.nch; ch++){
					for(s = 0; s < res.ch[ch].nsampev; s++) ops.push("R"+(ch+1)+":"+(s+1));	
				}
			}
			
			if(ops.length > 2){
				if(page == SIMULATEPAGE) addbuttonlim("Simulation:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,73);
				else addbuttonlim("Sample:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,74);
				y += dy1;
				if(y >= ybeg && y < ymax){
					gdropinfo.push({val:res.sampfilt, x:x, y:y, dx:selbutdx, dy:20, style:4, options:ops, click:"sampfilt"});
				}
			}
			y += dy2;
		}
		
		if(type == "pop" && res.popview == "Model"){
			addbuttonlim("Style:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,80); y += 30;
			addbuttonlim("Colour",x+15,y-5,radiodx-20,15,RADIOBUT,RADIOBUT,0,RADIOSTYLE2); y += 22;
			addbuttonlim("Greyscale",x+15,y-5,radiodx-20,15,RADIOBUT,RADIOBUT,1,RADIOSTYLE2); y += 22;
			addbuttonlim("Rescaled",x+15,y-5,radiodx-20,15,RADIOBUT,RADIOBUT,2,RADIOSTYLE2); y += 22;
		
			addbuttonlim("Label:",x,y,0,0,-1,SMALLTEXTBUT,BLACK,75); y += 30;
			addbuttonlim("No label",x+15,y-5,radiodx-20,15,RADIOBUT,RADIOBUT,0,RADIOLABEL); y += 22;
			addbuttonlim("Name",x+15,y-5,radiodx-20,15,RADIOBUT,RADIOBUT,1,RADIOLABEL); y += 22;
			addbuttonlim("Value",x+15,y-5,radiodx-20,15,RADIOBUT,RADIOBUT,2,RADIOLABEL); y += 22;
		}
		break;
	}
	
	y -= 10; if(y < yrightbeg) y = yrightbeg;
	if(type != "ind" && type != "Statistics" && type != "Correlation"){
		ytotright = y-yrightbeg;
		y += tableyfr*ytotright;
		
		
		tableyfrac = rightheight/ytotright;
		if(tableyfrac < 1){
			if(tableyfr > 1-tableyfrac) tableyfrac = 1-tableyfrac;
			addbuttonlim("",width-23,ybeg,13,rightheight,SLIDEAC,YSLIDEBUT,-1,-1);
		}
		else talbeyfr = 0;
	}
}

function varlistplot(y,type)                              // Makes a list of variables
{
	var yst, x1, x2;
	
	addbuttonlim("Variable type:",width-rwid,y,0,0,-1,SMALLTEXTBUT,BLACK,64); y += 20;	
	
	if(y >= ybeg && y < ymax){
		gdropinfo.push({val:res.filter[res.filt], x:width-rwid, y:y, dx:selbutdx, dy:20, style:4, options:res.filter, click:"filter"});
	}
	y += 30;
	 
	x1 = width-140; x2 = width-25;
	yst = y;
	
	for(v = 0; v < res.filtervar[res.filt].length; v++){
		addbuttonlim(res.varname[res.filtervar[res.filt][v]],x1,y,radiodx,15,RADIOBUT,RADIOBUT,res.filtervar[res.filt][v],type);
		y += 22;
	}
	
	return y;
}

function varchecklistplot(y,type)                              // Makes a list of variables
{
	var yst, x1, x2;
	
	addbuttonlim("Variable type:",width-rwid,y,0,0,-1,SMALLTEXTBUT,BLACK,64); y += 20;	
	
	if(y >= ybeg && y < ymax){
		gdropinfo.push({val:res.filter[res.filt], x:width-rwid, y:y, dx:selbutdx, dy:20, style:4, options:res.filter, click:"filter"});
	}
	y += 30;
	 
	x1 = width-140; x2 = width-25;
	yst = y;
	if(res.filtervar[res.filt].length > 0){
		addbuttonlim("Check All",x1-2,y,55,15,CHECKALLBUT,CHECKALLBUT,3,res.filt);
		addbuttonlim("Uncheck All",x1-2+55,y,58,15,CHECKALLBUT,CHECKALLBUT,3,res.filt);
		y += 20;
		for(v = 0; v < res.filtervar[res.filt].length; v++){
			addbuttonlim(res.varname[res.filtervar[res.filt][v]],x1,y,radiodx,15,CHECKBUT7,CHECKBUT7,res.filtervar[res.filt][v],-1);
			y += 22;
		}
	}
	
	return y;
}

function deplistplot(y)                                   // Makes a list of dependencies
{
	var yst, x1, x2;
	
	addbuttonlim("Dependency:",width-rwid,y,0,0,-1,SMALLTEXTBUT,BLACK,64); y += 20;	
	
	if(y >= ybeg && y < ymax){
		gdropinfo.push({val:depopsel, x:width-rwid, y:y, dx:selbutdx, dy:20, style:4, options:depop, click:"depop"});
	}
	y += 30;
	 
	x1 = width-140; x2 = width-25;
	yst = y;
	for(n = 0; n < histoplot.length; n++){
		if(histoplot[n].classname == depopsel){
			addbuttonlim(histoplot[n].name,x1,y,radiodx,15,RADIOBUT,RADIOBUT,n,RADIODEP);
			y += 22;
		}
	}
	
	return y;
}

function derlistplot(y,type)                              // Makes a list of variables
{
	var yst, x1, x2;
	
	addbuttonlim("Derived:",width-rwid,y,0,0,-1,SMALLTEXTBUT,BLACK,-1); y += 20;	
	if(page == INFERENCEPAGE) rest = infres; else rest = simres;
	
	x1 = width-140; x2 = width-25;
	yst = y;
	
	for(v = 0; v < rest.nderive; v++){ addbuttonlim(rest.derive[v],x1,y,radiodx,15,RADIOBUT,RADIOBUT,v,type); y += 22;}

	return y;
}

function rightotherclass(x,y,x1,dy1,dy2,claname,claname2) // Gives options for other classifications 
{
	var cl;
	for(cl = 0; cl < ncla; cl++){
		if(cla[cl].name != claname && cla[cl].name != claname2 && cla[cl].ncomp > 1){
			addbuttonlim(cla[cl].name+":",x,y,0,0,-1,SMALLTEXTBUT,BLACK,72); y += dy1;
			var pos=[];
			pos.push("All"); pos.push("Checkbox"); for(c = 0; c < cla[cl].ncomp; c++) pos.push(cla[cl].comp[c].name);
			if(res.popfilt[cl] == -1) val = "All"; 
			else{
				if(res.popfilt[cl] == -2) val = "Checkbox";
				else val = cla[cl].comp[res.popfilt[cl]].name;
			}
			
			if(y >= ybeg && y < ymax){
				gdropinfo.push({val:val, x:x, y:y, dx:selbutdx, dy:20, style:6, options:pos, click:"filtfilt", cl:cl});	
			}
			
			if(val == "Checkbox"){
				y += 26;
				addbuttonlim("Check All",x+3,y,55,15,CHECKALLBUT,CHECKALLBUT,1,cl);
				addbuttonlim("Uncheck All",x+3+55,y,58,15,CHECKALLBUT,CHECKALLBUT,1,cl);
				y += 20;
				for(v = 0; v < cla[cl].ncomp; v++){
					addbuttonlim(cla[cl].comp[v].name,x1,y-5,checkdx,15,CHECKBUT3,CHECKBUT3,cl,v); y += 22;
				}
				y -= 30;
			}
			y += dy2;
		}
	}
	return y;
}
