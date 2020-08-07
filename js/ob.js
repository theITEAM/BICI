function placeob()                                        // Places objects onto a canvas
{
	var x, y, xsh, yst;

	xsh = 0;	
	if(page == INFERENCEPAGE && pagesub[page] == 2 && subna == "Correlation"){
		tablexfrac = corwidth/ytot;
		if(tablexfr < 0) tablexfr = 0; if(tablexfrac < 1 && tablexfr > 1-tablexfrac) tablexfr = 1-tablexfrac;
		xsh = Math.floor(tablexfr*ytot);
	}
	
	tableyfrac = canvasheight/ytot;
	if(tableyfr < 0) tableyfr = 0; if(tableyfrac < 1 && tableyfr > 1-tableyfrac) tableyfr = 1-tableyfrac;
	
	ysh = Math.floor(tableyfr*ytot);
	
	if(page == INFERENCEPAGE && pagesub[page] == 1) ysh = Math.floor(ysh/priordy)*priordy;
	
	for(i = 0; i < nob; i++){
		x = obx[i]-xsh; y = oby[i]-ysh; 
		val = obval[i]; val2 = obval2[i]; val3 = obval3[i]; val4 = obval4[i];
		switch(obty[i]){
		case OBSPEECHOLD:
			addcanbutton("",660,y,25,22,SPEECHBUT,SPEECHBUT,-1,-1);
			if(val2 > 0){
				addcanbutton("",660-2,y+25,28,18,EDITSPEECHBUT,EDITSPEECHBUT,ncanbut+1,-1);
				addcanbutton(datanote,x+5,y,val,val2,PARACANBUT2,PARACANBUT2,-1,-1);
			}
			break;
			
		case OBSPEECH:
			if(val2 > 0) addcanbutton(val3,x+8,y,val,val2,-1,PARACANBUT,-1,-1);	
			break;
			
		case OBSPEECH2:
			addcanbutton("",x,y,val,val2,SPEECHBUT2,SPEECHBUT2,-1,-1);
			break;
		
		case OBSPEECH3:
			if(val2 > 0) addcanbutton(val3,x+8,y,val,val2,-1,PARACANBUT3,-1,-1);	
			break;
			
		case OBEXAMPPIC:
			addcanbutton(val2,x,y-20,25,22,-1,TEXTBUT,-1,-1);
			addcanbutton("",x,y,0,0,-1,EXAMPMODBUT,val,-1);
			addcanbutton("[?]",x+textwidth(val2,"bold 16px arial")+2,y-23,15,20,HELPICONCANBUT,HELPICONCANBUT,val3,-1);
			break;
		
		case OBEXAMP:
			if(val2 == examploaded) addcanbutton(val,x,y,textwidth(val,examplefont2)+20,20,EXAMPBUT,EXAMPBUT,val2,-1);
			else addcanbutton(val,x,y,textwidth(val,examplefont)+20,20,EXAMPBUT,EXAMPBUT,val2,-1);
			break;

		case OBTEXT:
			addcanbutton(val,x,y,150,0,-1,TEXTBUT,-1,-1);
			break;
		
		case OBLINK:
			addcanbutton(val,x,y,150,20,LINKBUT,LINKBUT,val2,val3);
			break;
			
		case OBTEXTEDIT:
			addcanbutton(val,x,y,val2,25,EDITDBUT,EDITDBUT,val3,-1);
			break;
		
		case OBTEXT2:
			addcanbutton(val,x,y,150,0,-1,TEXTBUT2,-1,-1);
			break;
	
		case OBPARAM:
			p = val;	
			addcanbutton(paramsim[p].name,0,y,0,0,-1,PARAMBUT,-1,-1);
			w = textwidth(paramsim[p].name,tablefont);
			x = w+5;
			kmax = paramsim[p].dep.length;
			for(k = 0; k < kmax; k++){
				ww = textwidth(paramsim[p].dep[k],subfont);
				addcanbutton(paramsim[p].dep[k],x,y+7,ww+2,14,SUBBUT2,SUBBUT2,p,k);
				x += ww;
				if(k < kmax-1){ addcanbutton(",",x,y+7,0,14,-1,SUBBUT,-1,-1); x += 4;}
			}
			
			addcanbutton(paramsim[p].sim,250,y-7,100,30,MINMAXBUT,MINMAXBUT,val,2);	
			break;
	
		case OBMINMAX:	
			addcanbutton(val,x,y-7,val3,30,MINMAXBUT,MINMAXBUT,0,val2);
			break;
			
		case OBPD:	
			addcanbutton(val,x,y-7,val3,30,PDBUT,PDBUT,0,val2);
			break;
			
		case OBCOMP:	
			w = textwidth(val3,"20px georgia");
			addcanbutton(val3,x,y,w+20,30,COMPSMALLBUT,COMPSMALLBUT,val,val2);
			break;	

		case OBCOMPPER:	
			w = textwidth(val3,"20px georgia");
			addcanbutton(val3,x,y,w+20,30,COMPSMALLPERBUT,COMPSMALLPERBUT,val,val2);
			break;		

		case OBCOMPPER2:	
			w = textwidth(val3,"20px georgia");
			addcanbutton(val3,x,y,w+20,30,-1,COMPSMALLPERBUT2,val,val2);
			break;						
	
		case OBCOMP2:	
			w = textwidth(val3,"20px georgia");
			addcanbutton(val3,x,y,w+20,30,COMPSMALLBUT2,COMPSMALLBUT2,val,val2);
			break;		

		case OBCOMP3:	
			w = textwidth(val2,"20px georgia");
			addcanbutton(val2,x,y,w+20,30,val3,COMPSMALLBUT3,val,-1);
			break;					
		
		case OBCOMP4:	
			w = textwidth(val3,"16px georgia");
			addcanbutton(val3,x,y,w+20,20,val4,COMPSMALLBUT4,val,val2);
			break;	
		
		case OBCOMP6:	
			w = textwidth(val3,"20px georgia");
			addcanbutton(val3,x,y,w+20,30,COMPSMALLBUT6,COMPSMALLBUT6,val,val2);
			break;	
			
		case OBCOMP7:	
			w = textwidth(val3,"20px georgia");
			addcanbutton(val3,x,y,w+20,30,COMPSMALLBUT7,COMPSMALLBUT7,val,val2);
			break;	
			
		case OBCOMP8:	
			w = textwidth(val3,"20px georgia");
			addcanbutton(val3,x,y,w+20,30,-1,COMPSMALLBUT6,val,val2);
			break;	
			
		case OBDATAHEAD:
			if(data.length == 0) addcanbutton("No data loaded.",x,y,180,0,-1,TABLEHEADBUT,-1,-1); 
			else{
				addcanbutton("Name",x,y,130,0,-1,TABLEHEADBUT,-1,-1); 
				addcanbutton("Type",x+datacol[1],y,130,0,-1,TABLEHEADBUT,-1,-1); 
				addcanbutton("Class.",x+datacol[2],y,130,0,-1,TABLEHEADBUT,-1,-1); 
				addcanbutton("Time range",x+datacol[3],y,130,0,-1,TABLEHEADBUT,-1,-1); 
			}
			break;
		
		case OBDATA:
			addcanbutton(data[val].name,x+1,y,100,24,TABLEBUT3,TABLEBUT3,val,-1); 
			switch(data[val].variety){
			case "state": st = "State"; break;
			case "trans": st = "Transition"; break;
			case "move": st = "Move"; break;
			case "cap": st = "Capture"; break;
			case "capid": st = "Capture ID"; break;
			case "cappd": st = "Capture PD"; break;
			case "pop": st = "Population"; break;
			case "der": st = "Derived"; break;
			case "presence": st = "Presence"; break;
			}
			addcanbutton(st,x+datacol[1],y,130,30,-1,TABLEBUT,-1,-1); 
			if(data[val].variety == "cap" || data[val].variety == "capid" || data[val].variety == "cappd" ||
  			   data[val].variety == "pop" || data[val].variety == "der" || data[val].variety == "presence") st = "---";
			else{
				if(data[val].cl == -1) st = "---"; 
				else st = cla[data[val].cl].name;
			}
			addcanbutton(st,x+datacol[2],y,130,30,-1,TABLEBUT,-1,-1); 
		
			if(data[val].tmin == undefined) te = "---";
			else{
				if(data[val].tmin == -large) te = "-∞"; else te = tpre(data[val].tmin,4);
				te += " — ";
				if(data[val].tmax == large) te += "∞"; else te += tpre(data[val].tmax,4);
			}
			
			addcanbutton(te,x+datacol[3],y,200,30,-1,TABLEBUT,-1,-1); 
			 
			addcanbutton("Data",x+datacol[5]+10,y,70,22,VIEWBUT,VIEWBUT,val,"table"); 
			
			switch(data[val].variety){
			case "state": addcanbutton("Obs. Mod.",x+datacol[6]+13,y,100,22,VIEWBUT,VIEWBUT,val,"obsmodel"); break;
			case "trans": addcanbutton("Filter",x+datacol[6]+13,y,70,22,VIEWBUT,VIEWBUT,val,"transfilt"); break;
			case "pop": addcanbutton("Pop.",x+datacol[6]+13,y,70,22,VIEWBUT,VIEWBUT,val,"pop"); break;
			case "der": addcanbutton("Der.",x+datacol[6]+13,y,70,22,VIEWBUT,VIEWBUT,val,"der"); break;
			case "cap": addcanbutton("Prob. Det.",x+datacol[6]+13,y,100,22,VIEWBUT,VIEWBUT,val,"pdmodel"); break;
			}
			addcanbutton("",x+datacol[7]+10,y,20,20,DELETEBUT,DELETEBUT,val,-1); 
			break;
		
		case OBPRIORTITLE:
			addcanbutton(val,0,y,0,0,-1,PRTITLEBUT,val2,-1);
			break;
		
		case OBAGESMOOTH:
			p = val;
			addcanbutton(paramagesmooth[p].name,0,y,0,0,-1,PARAMBUT,-1,-1);
			if(y >= 0 && y < priorheight-20){
				gdropinfo.push({val:paramagesmooth[p].type, x:cornx+160, y:corny+y, dx:120, dy:20, style:2, 
								options:smoothop, click:"agesmooth", j:p});
			}
			if(paramagesmooth[p].type != "None"){
				addcanbutton("Value: "+paramagesmooth[p].val,350,y-7,dxx,30,MINMAXBUT,MINMAXBUT,p,50);
			}
			break;
		
		case OBTIMESMOOTH:
			p = val;
			addcanbutton(paramtimesmooth[p].name,0,y,0,0,-1,PARAMBUT,-1,-1);
			if(y >= 0 && y < priorheight-20){
				gdropinfo.push({val:paramtimesmooth[p].type, x:cornx+160, y:corny+y, dx:120, dy:20, style:2, 
								options:smoothop, click:"timesmooth", j:p});
			}
			if(paramtimesmooth[p].type != "None"){
				addcanbutton("Value: "+paramtimesmooth[p].val,350,y-7,dxx,30,MINMAXBUT,MINMAXBUT,p,51);
			}
			break;
			
		case OBPRIOR:
			p = val;	
			addcanbutton(param[p].name,0,y,0,0,-1,PARAMBUT,-1,-1);
			w = textwidth(param[p].name,tablefont);
			x = w+5;
			
			kmax = param[p].dep.length;
			for(k = 0; k < kmax; k++){
				ww = textwidth(param[p].dep[k],subfont);
				addcanbutton(param[p].dep[k],x,y+7,ww+2,14,SUBBUT,SUBBUT,p,k);
				x += ww;
				if(k < kmax-1){ addcanbutton(",",x,y+7,0,14,-1,SUBBUT,-1,-1); x += 4;}
			}
			
			if(page == MODELPAGE){
				op = distributiontext; butt = SELPARAMBUT; if(modelsetup == 1) butt = NOTSELPARAMBUT;
			}
			else{ op = priortext; butt = MINMAXBUT; if(param[p].dist == 1) butt = NOTSELPARAMBUT;}
					
			if(param[p].prior == "Dirichlet" || butt == NOTSELPARAMBUT){
				addcanbutton(param[p].prior,140,y,0,0,-1,PARAMBUT,-1,-1);
			}
			else{
				if(y >= 0 && y < priorheight-20){
					gdropinfo.push({val:param[p].prior, x:cornx+140, y:corny+y, dx:120, dy:20, style:2, 
									options:op, click:"prior", j:p});
				}
			}
			
			x = 320; dxx = 160;
			switch(param[p].prior){
			case "Flat":
				addcanbutton("Min.: "+param[p].val[0],x,y-7,dxx,30,butt,butt,p,0);
				addcanbutton("Max.: "+param[p].val[1],x+dxx+20,y-7,dxx,30,butt,butt,p,1);
				break;

			case "Gamma":
				addcanbutton("Mean: "+param[p].val[0],x,y-7,dxx,30,butt,butt,p,30);
				addcanbutton("SD: "+param[p].val[1],x+dxx+20,y-7,dxx,30,butt,butt,p,31);
				break;
				
			case "Normal":
				addcanbutton("Mean: "+param[p].val[0],x,y-7,dxx,30,butt,butt,p,30);
				addcanbutton("SD: "+param[p].val[1],x+dxx+20,y-7,dxx,30,butt,butt,p,31);
				break;
			
			case "Log-Normal":
				addcanbutton("Mean (logscale): "+param[p].val[0],x,y-7,dxx,30,butt,butt,p,32);
				addcanbutton("SD (logscale): "+param[p].val[1],x+dxx+20,y-7,dxx,30,butt,butt,p,33);
				break;
			
			case "Exponential":
				addcanbutton("Rate: "+param[p].val[0],x,y-7,dxx,30,butt,butt,p,34);
				break;
		
			case "Beta":
				addcanbutton("α: "+param[p].val[0],x,y-7,dxx,30,butt,butt,p,36);
				addcanbutton("β: "+param[p].val[1],x+dxx+20,y-7,dxx,30,butt,butt,p,37);
				break;
				
			case "Weibull":
				addcanbutton("λ: "+param[p].val[0],x,y-7,dxx,30,butt,butt,p,38);
				addcanbutton("k: "+param[p].val[1],x+dxx+20,y-7,dxx,30,butt,butt,p,39);
				break;
			
			case "Fix":
				addcanbutton("Value: "+param[p].val[0],x,y-7,dxx,30,butt,butt,p,35);
				break;
				
			case "Dirichlet":
				addcanbutton("α: "+param[p].val[0],x,y-7,dxx,30,butt,butt,p,40);
				break;
			};
			
			if(butt == NOTSELPARAMBUT){
				addcanbutton("",0,y,750,20,INACTIVEBUT,INACTIVEBUT,-1,-1); 
			}
			else{			
				if(page == MODELPAGE) addcanbutton("",x+400,y,20,20,DELDEPAC,DELETEBUT,p,-1); 
			}
			break;

		case OBPRIOR2:
			p = val;	
			b = ncanbut;
			x = 10;
			addcanbutton(param[p].name,1,y,x,23,SELPARBUT,SELPARBUT,p,-1);
		
			addcanbutton(param[p].name,x,y,0,0,-1,PARAMBUT,-1,-1);
			w = textwidth(param[p].name,tablefont);
			x = 10+w+5;
			
			kmax = param[p].dep.length;
			for(k = 0; k < kmax; k++){
				ww = textwidth(param[p].dep[k],subfont);
				addcanbutton(param[p].dep[k],x,y+7,0,0,-1,SUBBUT,p,k);
				x += ww;
				if(k < kmax-1){ addcanbutton(",",x,y+7,0,14,-1,SUBBUT,-1,-1); x += 4;}
			}
			canbutdx[b] = x+20;
		
			break;
		
		case OBSESP:
			addcanbutton("Sensitivity Se:",x,y,150,0,-1,TEXTBUT,-1,-1);
			addcanbutton(Segen,x+110,y-9,80,30,MINMAXBUT,MINMAXBUT,-1,12);
			addcanbutton("Specificity Sp:",x+200,y,150,0,-1,TEXTBUT,-1,-1);
			addcanbutton(Spgen,x+310,y-9,80,30,MINMAXBUT,MINMAXBUT,-1,13);
			break;
			
		case OBDERIVE:
			p = val;	

			addcanbutton(derive[p].name,0,y,0,0,-1,PARAMBUT,-1,-1);
			w = textwidth(derive[p].name,tablefont);
			x = w+5;
			kmax = derive[p].dep.length;
			for(k = 0; k < kmax; k++){
				ww = textwidth(derive[p].dep[k],subfont);
				addcanbutton(derive[p].dep[k],x,y+7,ww+2,14,-1,SUBBUT,p,k);
				x += ww;
				if(k < kmax-1){ addcanbutton(",",x,y+7,0,14,-1,SUBBUT,-1,-1); x += 4;}
			}
		
			addcanbutton(derive[val].eq,150,y-7,500,30,EQBUT,EQBUT,val,0);
			if(modelsetup == 0)	addcanbutton("",725+10,y,20,20,DELDERAC,DELETEBUT,val,-1); 
			else{
				addcanbutton("",0,y-7,750,30,INACTIVEBUT,INACTIVEBUT,-1,-1); 
			}
			break;
			
		case OBNAME:
			st = datatemp.testname;
			addcanbutton(st,x,y,textwidth(st,"18px georgia")+40,20,TESTNAMEBUT,TESTNAMEBUT,-1,-1);
			break;
			
		case OBRADIOWHICH:
			addcanbutton("Entire population",x,y,110,22,CANRADIOBUT,CANRADIOBUT,"all",CANRADIOWHICH);
			addcanbutton("Subpopulation",x+200,y,110,22,CANRADIOBUT,CANRADIOBUT,"sub",CANRADIOWHICH);
			break;
			
		case OBRADIOEVPD:
			addcanbutton("All transitions",x,y,110,22,CANRADIOBUT,CANRADIOBUT,"all",CANRADIOPD);
			addcanbutton("Some transitions",x+200,y,110,22,CANRADIOBUT,CANRADIOBUT,"set",CANRADIOPD);
			break;
			
		case OBRADIOPD:
			addcanbutton("All individuals",x,y,110,22,CANRADIOBUT,CANRADIOBUT,"all",CANRADIOPD);
			addcanbutton("Individuals sampled",x+200,y,110,22,CANRADIOBUT,CANRADIOBUT,"set",CANRADIOPD);
			break;
			
		case OBRADIOPD2:
			addcanbutton("All individuals",x,y,110,22,CANRADIOBUT,CANRADIOBUT,"all",CANRADIOPD2);
			addcanbutton("Individuals sampled",x+200,y,110,22,CANRADIOBUT,CANRADIOBUT,"set",CANRADIOPD2);
			break;
		
		case OBRADIOPD3:
			addcanbutton("The same",x,y,110,22,CANRADIOBUT,CANRADIOBUT,"same",CANRADIOPD3);
			addcanbutton("Different",x+200,y,110,22,CANRADIOBUT,CANRADIOBUT,"dif",CANRADIOPD3);
			break;
			
		case OBRADIO:
			addcanbutton("Single",x,y,110,22,CANRADIOBUT,CANRADIOBUT,"simple",val);
			addcanbutton("Multiple",x+110,y,110,22,CANRADIOBUT,CANRADIOBUT,"binary",val);
			addcanbutton("Diagnostic test",x+220,y,110,22,CANRADIOBUT,CANRADIOBUT,"test",val);
			addcanbutton("General",x+380,y,150,22,CANRADIOBUT,CANRADIOBUT,"expression",val);
			break;
			
		case OBRADIOGEN:
			addcanbutton("Simple",x,y,110,22,CANRADIOBUT,CANRADIOBUT,"simple",val);
			addcanbutton("Edit",x+110,y,110,22,CANRADIOBUT,CANRADIOBUT,"binary",val);
			if(val2 > 1){
				addcanbutton("Diagnostic test",x+220,y,110,22,CANRADIOBUT,CANRADIOBUT,"test",val);
			}
			break;
			
		case OBRADIO2:
			addcanbutton("+ve Test result",x,y,150,22,CANRADIOBUT,CANRADIOBUT,1000*val+1,CANRADIOTESTRES);
			addcanbutton("-ve Test result",x+170,y,150,22,CANRADIOBUT,CANRADIOBUT,1000*val+0,CANRADIOTESTRES);
			addcanbutton("No test",x+340,y,150,22,CANRADIOBUT,CANRADIOBUT,1000*val+2,CANRADIOTESTRES);
			break;
		
		case OBRADIO3:
			addcanbutton(val2,x,y,340,22,CANRADIOBUT,CANRADIOBUT,val,CANRADIODEFINE);
			break;
	
		case OBRADIO5:
			addcanbutton("Manually set initial populations",x,y,350,22,CANRADIOBUT,CANRADIOBUT,"Manualpop",CANRADIOMANUAL);
			addcanbutton("Load initial individuals from file",x+350,y,350,22,CANRADIOBUT,CANRADIOBUT,"Load",CANRADIOMANUAL);
			break;
			
		case OBBRACKET:
			addcanbutton("",x,y,10,val,-1,BRACKETBUT,-1,-1);
			break;
		
		case OBIND:
			var indi;
			if(page == INFERENCEPAGE){
				if(pagesub[page] == 0) indi = inddata.ind[val];
				else indi = infres.inddata.ind[val];
			}
			else indi = indsim.ind[val];
		
			if(!indi){
				if(page == SIMULATEPAGE) addcanbutton("Ind. "+(val+1),x,y,150,20,INDBUT,INDBUT,val,-1);
				else{
					addcanbutton("Unobserved Individual "+(val+1-infres.inddata.nindtotal),x,y,150,20,INDBUT,INDBUT,val,-1);
				}
				yy = y+30;
				for(cl = 0; cl < ncla; cl++){
					if(indshow[cl] == 1){
						addcanbutton(cla[cl].name,tlinexmin,yy,tlinexmax-tlinexmin,1,-1,TIMELINEBUT,val,cl); yy += 25;
					}
				}
			}
			else{
				addcanbutton(indi.id,x,y,150,20,INDBUT,INDBUT,val,-1);
				yy = y+30;
				for(cl = 0; cl < ncla; cl++){
					if(indshow[cl] == 1){
						addcanbutton(cla[cl].name,tlinexmin,yy,tlinexmax-tlinexmin,1,-1,TIMELINEBUT,val,cl); yy += 25;
					}
				}
				
				yy = y+30;

				numon = 0;
				for(cl = 0; cl < ncla; cl++){
					if(indshow[cl] == 1){
						numon++;
						for(e = 0; e < indi.cl[cl].ev.length; e++){
							var ev = indi.cl[cl].ev[e];
							
							r = 6;
							if(ev.t >= axxmin && ev.t < axxmax){
								x = tlinexmin + Math.floor((tlinexmax-tlinexmin)*(ev.t-axxmin)/(axxmax-axxmin));
								switch(ev.variety){
								case "state":
									if(page == INFERENCEPAGE) addcanbutton(e,x-r,yy-r,2*r,2*r,EVBUT,EVBUT,val,cl);
									else addcanbutton(e,x-r,yy-r,2*r,2*r,EVBUT4,EVBUT4,val,cl);
									break;

								case "trans": case "move":					
									if(ev.col[0] >= 0) addcanbutton(e,x-r,yy-r,2*r,2*r,EVBUT2,EVBUT2,val,cl);
									break;
									
								case "transtmin":
									addcanbutton(ev.name,x-4,yy-10,8,20,EVBUT5,EVBUT5,ev.t,-1);
									break;
									
								case "transtmax":
									addcanbutton(ev.name,x-4,yy-10,8,20,EVBUT6,EVBUT6,ev.t,-1);
									break;
								
								case "presence":
									addcanbutton(e,x-r,yy-r,2*r,2*r,EVBUT7,EVBUT7,-1,-1);
									break;
								}
							}
						}
						yy += 25;
					}
				}
				
				cl = ncla-2;      // draws any source sink labels
				if(numon > 0){
					for(e = 0; e < indi.cl[cl].ev.length; e++){
						var ev = indi.cl[cl].ev[e];
								
						r = 6;
						if(ev.t >= axxmin && ev.t < axxmax){
							x = tlinexmin + Math.floor((tlinexmax-tlinexmin)*(ev.t-axxmin)/(axxmax-axxmin));
							if(ev.variety == "trans" || ev.variety == "move"){						
								var col = ev.col;
								if(col[0] < 0){
									if(col[0] == -1) addcanbutton(e,x-r,y+30-r-4,2*r,numon*25+13,EVSOURCEBUT,EVSOURCEBUT,val,cl);
									else addcanbutton(e,x-r,y+30-r-4,2*r,numon*25+13,EVSINKBUT,EVSINKBUT,val,cl);
								}
							}
						}
					}
				}
			}
			break;
		
		case OBINDUO:
			addcanbutton("The system may also contain unobserved individuals.",x+60,y,width-500,0,-1,PARACANBUT2,-1,-1);
			addcanbutton("To observe these select a given sample from the right hand menu.",x+60,y+22,width-500,0,-1,PARACANBUT2,-1,-1);
			break;
			
		case OBUPLOAD:
			addcanbutton(val2,val3,y,80,30,val,CANUPLOADBUT,-1,-1);
			break;
			
		case OBSELIND:
			var ops=[]; for(j = 0; j < cla[val].ncomp; j++) ops.push(cla[val].comp[j].name); ops.push("All");
			
			addcanbutton(cla[val].name+":",x,y,150,20,-1,TEXTBUT2,-1,-1);
			gdropinfo.push({val:clagendata[val], x:cornx+x+val2+3, y:corny+y+2, dx:95, dy:20, style:2, options:ops, click:"gendata", cl:val});
			break;
		
		case OBSELIND2:
			addcanbutton(cla[val].name+":",x,y,150,20,-1,TEXTBUT2,-1,-1);
			addcanbutton(clagendata[val],x+val2+3,y-4,150,20,-1,NOTACTIVEBUT,-1,-1);
			break;
			
		case OBDETECTPROB:
			addcanbutton(detectprob,x,y-7,100,30,MINMAXBUT,MINMAXBUT,val,8);	
			break;
		
		case OBCHOOSECLA:
			if(selclass[val] == 1) st = cla[val].name+" ✓"; else st = cla[val].name+" ✖";
			addcanbutton(st,x,y,val2,30,CANUPLOADBUT,CANUPLOADBUT,val,-1);
			break;
		
		case OBTIMEP:
			addcanbutton("Initial",x,y,140,22,CANRADIOBUT,CANRADIOBUT,1,CANRADIOGENT);
			addcanbutton("Periodic",x+135,y,140,22,CANRADIOBUT,CANRADIOBUT,0,CANRADIOGENT);
			addcanbutton("User defined",x+275,y,140,22,CANRADIOBUT,CANRADIOBUT,2,CANRADIOGENT);
			addcanbutton("At transitions",x+425,y,140,22,CANRADIOBUT,CANRADIOBUT,3,CANRADIOGENT);
			break;
			
		case OBINDUOFLAG:
			addcanbutton("All individuals observed",x+10,y,200,22,CANRADIOBUT,CANRADIOBUT,0,CANRADIOINDUO);
			addcanbutton("Unobserved individuals",x+290,y,200,22,CANRADIOBUT,CANRADIOBUT,1,CANRADIOINDUO);
			break;
			
		case OBHOZLAB: 
			addcanbutton(val,x,y,cordx,cordy,-1,HOZLABBUT,-1,-1);
			break;
			
		case OBVERTLAB:	
			addcanbutton(val,x,y,cordy,cordx,-1,VERTLABBUT,-1,-1);
			break;
			
		case OBCOR:
			if(val2 != val3) addcanbutton(val,x,y,cordy,cordy,CORBUT,CORBUT,val2,val3);
			else addcanbutton(val,x,y,cordy,cordy,-1,CORBUT,val2,val3);
			break;
		}
	}
}

function addob(x,y,ty,val,val2,val3,val4)                // Adds an object to a canvas
{
	obx[nob]=x; oby[nob]=y; obty[nob]=ty; obval[nob]=val; obval2[nob]=val2; obval3[nob]=val3; obval4[nob]=val4; nob++;
}

