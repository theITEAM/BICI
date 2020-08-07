function simulatebuts()                                      // Sets up the buttons on the simulation page
{
	var x, y;
	x = menux+tab; y = 30;
	cv = resultcv;

	res = simres; 
	 
	switch(pagesub[SIMULATEPAGE]){
	case 0:
		addbutton("Set model parameters",x,y,0,0,-1,TITLEBUT,25,-1); y += 50;
		
		if(paramsim.length == 0) addbutton("There are no model parameters.",x+18,y,900,0,-1,PARAGRAPHBUT,1,-1);
		else{
			addbutton("",x+20,y,siminitwidth,simheight,CANVASBUT,CANVASBUT,-1,-1);

			drawparaminit();
				
			tableyfrac = simheight/ytot;
			if(tableyfrac < 1) addbutton("",width-30,y,13,simheight,SLIDEAC,YSLIDEBUT,-1,-1);
		}
		addbutton("Next",width-105,height-45,90,30,NEXTSIMAC,NEXTBUT,0,-1);		
		break;
		
	case 1:
		switch(showsimindinit){
		case 0:
			addbutton("Set initial polulation",x,y,0,0,-1,TITLEBUT,26,-1); y += 50;
			
			cornx = menux+20;
			addbutton("",cornx,y,modeldx,modeldy+50,CANVASBUT,CANVASBUT,-1,-1);
		
			x = 10; y = 0;
			addcanbutton("Manually set",x+40,y,350,22,CANRADIOBUT,CANRADIOBUT,"Manualpop",CANRADIOMANUAL);
			addcanbutton("Load initial individuals",x+350,y,350,22,CANRADIOBUT,CANRADIOBUT,"Load",CANRADIOMANUAL);
			y += 30;
			
			switch(siminit){
			case "":
				addcanbutton("Please select how the initial conditions for the simulation are defined.",x+60,y,150,0,-1,TEXTBUT2,-1,-1);
				break;
				
			case "Manualpop":
				addcanbutton("Specify the initial populations in each of the compartments:",x,y,150,0,-1,TEXTBUT2,-1,-1); y += 30;

				rightmenu(width-rwid,125,"Manualpop");
				
				drawmodel(viewpopinitset,0,y,modeldx,modeldy-y,"frame","popinit");
				
				addbutton("Next",width-105,height-45,90,30,NEXTSIMAC,NEXTBUT,0,-1);		
				break;
		
			case "Load":
				y += 20;
				
				if(simindinit.length == 0){
					addcanbutton("Please upload a file specifying the initial population:",x,y,150,0,-1,TEXTBUT2,-1,-1);
					addbutton("Upload",width-rwid,129,80,30,SIMLOADINITAC,UPLOADBUT,-1,-1);
				}
				else{
					drawmodel(viewpopinitset,0,y,modeldx,modeldy-y,"frame","loadinit");
					addbutton("Reload",width-rwid,125,80,30,SIMLOADINITAC,UPLOADBUT,-1,-1);
					addbutton("Edit",width-rwid,125+40,80,30,SIMVIEWINITAC,GREENBUT,-1,-1);
					
					rightmenu(width-rwid,240,"loadinitpop");
					
					addbutton("Next",width-105,height-45,90,30,NEXTSIMAC,NEXTBUT,0,-1);		
				}
				break;
			}
			break;
			
		case 1:
			addbutton("Initial population",x,y,0,0,-1,TITLEBUT,26,-1); y += 50;
			
			addbutton("The following table defines the initial population:",x+18,y,900,0,-1,PARAGRAPHBUT,1,-1);
			y += 45;
			
			cornx = x+40; corny = y;
			addbutton("",cornx,corny,tablewidth,tableheight,CANVASBUT,CANVASBUT,-1,-1);
			
			drawtable();
			
			tableyfrac = rowmax/nrow;
			if(tableyfrac < 1) addbutton("",x+40+ tablewidth+10,y+22,13,tableheight-22,SLIDEAC,YSLIDEBUT,-1,-1);

			tablexfrac = tablewidth/tottablewidth;
			if(tablexfrac < 1) addbutton("",x+40,y+tableheight+10,tablewidth,13,SLIDEAC,XSLIDEBUT,-1,-1);
			
			addbutton("Cancel",width-205,height-45,90,30,SIMINITCANCELAC,CANCELBUT2,0,-1);
			
			addbutton("Done",width-105,height-45,90,30,SIMINITPOPBACKAC,ADDDATABUT,0,-1);	
			break;
		}
		break;

	case 2:
		if(warning.length > 0){ addwarning(); return;}
	
		na = simpagename[pagesubsub[SIMULATEPAGE][2]];

		switch(na){
		case "Start":
			if(simres.result != 0){
				cornx = menux+tab+20; corny = 80;
				addbutton("",cornx,corny,setupwidth,setupheight,CANVASBUT,CANVASBUT,-1,-1);
		
				if(advop == 0){
					addbutton("Run",x,y,0,0,-1,TITLEBUT,224,-1);
					drawsimsetup();
				}
				else{
					addbutton("Advanced options",x,y,0,0,-1,TITLEBUT,223,-1);
					drawsimadvop();
				}
				loading = 0;
			}
			else{
				if(loading == 0) startloading();
				addbutton("",menux+(width-menux)/2,height/2+60,0,0,-1,PROGRESSBUT,0,-1);
				addbutton("Cancel",width-115,height-45,100,30,STARTAC,ADDDATABUT,0,-1);
			}
			break;
			
		case "Populations": case "Derived":
			if(simres.result == 1) addbutton("Stop!",width-115,height-45,100,30,STARTAC,UPLOADBUT,0,-1);
		
			rightmenu(width-rwid,75,"pop");
		
			showpopulations(2);
			y = height-120;
			
			addbutton("Rerun", width-115,height-45,100,30,STARTAC,GREENBUT,0,-1);
			break;
			
		case "Individuals":
			if(simres.result == 1) addbutton("Stop!",width-115,height-45,100,30,STARTAC,UPLOADBUT,0,-1);
		
			indplot();
			addbutton("Rerun", width-115,height-45,100,30,STARTAC,GREENBUT,0,-1);
			break;
		
		case "Statistics":
			addbutton("Statistics",x,30,0,0,-1,TITLEBUT,48,-1);
			rightmenu(width-rwid,75,"Statistics"); 
			showstatistics();
			break;
			
		case "Gen. Data":
			switch(gentype){
			case -1:
				addbutton("Generate data",x,y,0,0,-1,TITLEBUT,28,-1); y += 40;
	 
				if(simres.nderive == 0){ addbutton("Three basic types of data can be generated:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); y += 95;}
				else{ addbutton("Four basic types of data can be generated:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); y += 45;}
					
				addbutton("State data",x+20,y,800,0,-1,TEXTSUBTITLEBUT,-1,-1); y += 25;
				addbutton("Generate", width-155,y,100,30,GENSTARTAC,UPLOADBUT,0,-1);
				addbutton("Output state measurements of individuals at specified points in time.",x+20,y,600,0,-1,PARAGRAPHBUT,-1,-1);
				y += 90;
				
				addbutton("Transition data",x+20,y,800,0,-1,TEXTSUBTITLEBUT,-1,-1); y += 25;
				addbutton("Generate", width-155,y,100,30,GENSTARTAC,UPLOADBUT,1,-1);
				addbutton("Output the times of a specified transition.",x+20,y,600,0,-1,PARAGRAPHBUT,-1,-1);
				y += 90;
				
				addbutton("Population data",x+20,y,800,0,-1,TEXTSUBTITLEBUT,-1,-1); y += 25;
				addbutton("Generate", width-155,y,100,30,GENSTARTAC,UPLOADBUT,2,-1);
				addbutton("Output population data at specified points in time.",x+20,y,600,0,-1,PARAGRAPHBUT,-1,-1);
				y += 90;
				
				if(simres.nderive > 0){
					addbutton("Derived data",x+20,y,800,0,-1,TEXTSUBTITLEBUT,-1,-1); y += 25;
					addbutton("Generate", width-155,y,100,30,GENSTARTAC,UPLOADBUT,3,-1);
					addbutton("Output derived data at specified points in time.",x+20,y,600,0,-1,PARAGRAPHBUT,-1,-1);
					y += 90;
				}				
				break;
				
			case 0:	
				switch(drawgennum){
				case -2: case -1: addbutton("Individuals observed",x,y,0,0,-1,TITLEBUT,-1,-1); break;
				default: 
					if(drawgennum < sellist.length){
						addbutton("Observation model for '"+cla[sellist[drawgennum]].name+"'",x,y,0,0,-1,TITLEBUT,-1,-1);
					}
					else{
						if(drawgennum == sellist.length) addbutton("State Data",x,y,0,0,-1,TITLEBUT,-1,-1);
						if(drawgennum == sellist.length+1) addbutton("Capture Data",x,y,0,0,-1,TITLEBUT,-1,-1);
					}
				}
				y += 50;
		
				cornx = x+20; corny = y;
				if(drawgennum < sellist.length){
					drawgendata();
					addbutton("",x+20,y,tablewidth,simheight,CANVASBUT,CANVASBUT,-1,-1);
					tableyfrac = simheight/ytot;
					if(tableyfrac < 1) addbutton("",menux+tab+20+tablewidth+10,y,13,simheight,SLIDEAC,YSLIDEBUT,-1,-1);
					
					if(drawgennum >= 0 || drawgennum == -1){
						addbutton("Back",width-305,height-45,90,30,BACKGENAC,BACKBUT,0,-1);
						addbutton("Cancel",width-205,height-45,90,30,GENCANCAC,CANCELBUT2,0,-1);
					}
					else addbutton("Cancel",width-205,height-45,90,30,GENCANCAC,CANCELBUT2,0,-1);
					addbutton("Next",width-105,height-45,90,30,NEXTGENAC,NEXTBUT,0,-1);
				}
				else{
					if(drawgennum == sellist.length){
						addbutton("Here is the state data generated from the simulation and based on the specified observation model:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					}
					else{
						addbutton("This data provides information on when and which individuals are captured:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					}
					y += 40;
					cornx = x+40; corny = y;
					addbutton("",cornx,corny,tablewidth,tableheight,CANVASBUT,CANVASBUT,-1,-1);
					tableyfrac = rowmax/nrow;
					if(tableyfrac < 1) addbutton("",x+40+ tablewidth+10,y+22,13,tableheight-22,SLIDEAC,YSLIDEBUT,-1,-1);

					tablexfrac = tablewidth/tottablewidth;
					if(tablexfrac < 1) addbutton("",x+40,y+tableheight+10,tablewidth,13,SLIDEAC,XSLIDEBUT,-1,-1);
					drawtable();
					
					addbutton("Back",width-305,height-45,90,30,BACKGENAC,BACKBUT,0,-1);
					addbutton("Save",width-205,height-45,90,30,SAVEGENAC,EXPORTBUT,0,-1);
					if(ncoldef == (sellist.length+2) || drawgennum == sellist.length+1){
						addbutton("Done",width-105,height-45,90,30,GENCANCAC,ADDDATABUT,0,-1);
					}
					else addbutton("Next",width-105,height-45,90,30,NEXTGENAC,NEXTBUT,0,-1);
		
				}
				break;
			
			case 1:
				addbutton("Generate event data",x,y,0,0,-1,TITLEBUT,-1,-1); y += 40;
				switch(drawgennum){
				case 0: 
					transop(x,y);
					break;
					
				case 1:
					cornx = x+20; corny = y;

					addbutton("",x+20,y,tablewidth,tableheight,CANVASBUT,CANVASBUT,-1,-1);
					
					addingdata = 0;
					drawfilter(transcl);
					
					tableyfrac = tableheight/ytot;
					if(tableyfrac < 1) addbutton("",menux+tab+20+tablewidth+10,y,13,tableheight,SLIDEAC,YSLIDEBUT,-1,-1);
					addbutton("Back",width-305,height-45,90,30,BACKGENAC,BACKBUT,0,-1);
					addbutton("Cancel",width-205,height-45,90,30,GENCANCAC,CANCELBUT2,0,-1);
					addbutton("Next",width-105,height-45,90,30,NEXTGENEVAC,NEXTBUT,0,-1);
					break;
					
				case 2:
					st = "This data provides information on ";
					switch(transty){
					case "trans": st += "'"+transni+" → "+transnf+"'"; break;
					case "+": st += "source"; break;
					case "-": st += "sink"; break;
					}
					
					st += " transitions";
					var st2=""; for(cl = 0; cl < ncla; cl++) if(clagendata[cl] != "All" && cl != transcl) st2 += clagendata[cl]+",";
					if(st2 != "") st += " (for individuals in "+st2.substr(0,st2.length-1)+")";
					addbutton(st+":",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					y += 40;
					cornx = x+40; corny = y;
					addbutton("",cornx,corny,tablewidth,tableheight,CANVASBUT,CANVASBUT,-1,-1);
					tableyfrac = rowmax/nrow;
					if(tableyfrac < 1) addbutton("",x+40+ tablewidth+10,y+22,13,tableheight-22,SLIDEAC,YSLIDEBUT,-1,-1);

					tablexfrac = tablewidth/tottablewidth;
					if(tablexfrac < 1) addbutton("",x+40,y+tableheight+10,tablewidth,13,SLIDEAC,XSLIDEBUT,-1,-1);
					drawtable();
					
					addbutton("Back",width-305,height-45,90,30,BACKGENAC,BACKBUT,0,-1);
					addbutton("Save",width-205,height-45,90,30,SAVEGENAC,EXPORTBUT,0,-1);
					addbutton("Done",width-105,height-45,90,30,GENCANCAC,ADDDATABUT,0,-1);
					break;
				}
				break;
				
			case 2:  // generate population data
				addbutton("Generate population data",x,y,0,0,-1,TITLEBUT,-1,-1); y += 40;
				switch(drawgennum){
				case 0: 
					drawgenpopdata();
					addbutton("Next",width-105,height-45,90,30,NEXTGENAC3,NEXTBUT,0,-1);
					addbutton("Cancel",width-205,height-45,90,30,GENCANCAC,CANCELBUT2,0,-1);
					break;
					
				case 1:
					addbutton("This table provides population data:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					y += 40;
					
					drawtable();	
					cornx = x+40; corny = y;
					addbutton("",cornx,corny,tablewidth,tableheight,CANVASBUT,CANVASBUT,-1,-1);
					tableyfrac = rowmax/nrow;
					if(tableyfrac < 1) addbutton("",x+40+ tablewidth+10,y+22,13,tableheight-22,SLIDEAC,YSLIDEBUT,-1,-1);

					tablexfrac = tablewidth/tottablewidth;
					if(tablexfrac < 1) addbutton("",x+40,y+tableheight+10,tablewidth,13,SLIDEAC,XSLIDEBUT,-1,-1);
					
					addbutton("Back",width-305,height-45,90,30,BACKGENAC,BACKBUT,0,-1);
					addbutton("Save",width-205,height-45,90,30,SAVEGENAC,EXPORTBUT,0,-1);
					addbutton("Done",width-105,height-45,90,30,GENCANCAC,ADDDATABUT,0,-1);
					break;
				}
				break;
			
			case 3:  // generate derived data
				addbutton("Generate derived data",x,y,0,0,-1,TITLEBUT,-1,-1); y += 40;
				switch(drawgennum){
				case 0: 
					drawgenderdata();
					
					addbutton("Next",width-105,height-45,90,30,NEXTGENAC4,NEXTBUT,0,-1);
					addbutton("Cancel",width-205,height-45,90,30,GENCANCAC,CANCELBUT2,0,-1);
					break;
					
				case 1:
					addbutton("This table provides derived data:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					y += 40;
					
					drawtable();	
					cornx = x+40; corny = y;
					addbutton("",cornx,corny,tablewidth,tableheight,CANVASBUT,CANVASBUT,-1,-1);
					tableyfrac = rowmax/nrow;
					if(tableyfrac < 1) addbutton("",x+40+ tablewidth+10,y+22,13,tableheight-22,SLIDEAC,YSLIDEBUT,-1,-1);

					tablexfrac = tablewidth/tottablewidth;
					if(tablexfrac < 1) addbutton("",x+40,y+tableheight+10,tablewidth,13,SLIDEAC,XSLIDEBUT,-1,-1);
					
					addbutton("Back",width-305,height-45,90,30,BACKGENAC,BACKBUT,0,-1);
					addbutton("Save",width-205,height-45,90,30,SAVEGENAC,EXPORTBUT,0,-1);
					addbutton("Done",width-105,height-45,90,30,GENCANCAC,ADDDATABUT,0,-1);
					break;
				}
				break;
			}		
			break;
		}
		break;
	}
}

function drawsimsetup()
{
	x = 10; y = 10;

	addcanbutton("Time Range",x,y,150,0,-1,TEXTBUT,-1,-1); y += 30;

	addcanbutton("Begin: "+tsimmin,rightval,y-15,135,30,MINMAXBIGBUT,MINMAXBIGBUT,-1,3);
	addcanbutton("End: "+tsimmax,rightval,y+15,135,30,MINMAXBIGBUT,MINMAXBIGBUT,-1,7);

	addcanbutton("This represented the range in time over which simulation is performed.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 25;
	y += 25;
	te = "Advanced options..."; addcanbutton(te,x,y,textwidth(te,"14px arial"),20,ADVOPAC,LINKBUT,-1,-1);
	
	addbutton("Start", width-115,height-45,100,30,STARTAC,GREENBUT,0,-1);
}

function drawsimadvop()
{
	x = 10; y = 10;
	
	addcanbutton("Limit",x,y,150,0,-1,TEXTBUT,-1,-1); y += 30;
	
	addcanbutton("The maximum number of allowable individuals.",x+15,y,150,0,-1,TEXTBUT2,-1,-1);
	addcanbutton(indmaxnumber,rightval,y-8,129,30,MINMAXBIGBUT,MINMAXBIGBUT,0,99); y += 50;

	addcanbutton("Simulation number",x,y,150,0,-1,TEXTBUT,-1,-1); y += 30;
				
	addcanbutton("The total number of simulations to be performed. Since the model is stochastic,",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 25;
	addcanbutton("each simulation will geneate a different result.",x+15,y,150,0,-1,TEXTBUT2,-1,-1); y += 50;

	addcanbutton("Single simulation",x+40,y,180,22,CANRADIOBUT,CANRADIOBUT,0,CANRADIOSIMTY);
	addcanbutton("Multiple simulations",x+240,y,180,22,CANRADIOBUT,CANRADIOBUT,1,CANRADIOSIMTY);
	if(simty == 1) addcanbutton("Number: "+simnumber,rightval,y-8,129,30,MINMAXBIGBUT,MINMAXBIGBUT,0,6); y += 40;

	addbutton("Back",width-105,height-45,90,30,ADVOPAC2,BACKBUT,0,-1);	
}

function transop(x,y)                                      // Page giving options for transitions
{
	var canac;
	if(page == INFERENCEPAGE) canac = CANCELBUT2; else canac = GENCANCAC;
		
	y += 5;
	addbutton("Transition",x+20,y,180,20,RADIOBUT,RADIOBUT,"trans",RADIOTRANSTY);
	addbutton("Source",x+20+180,y,180,20,RADIOBUT,RADIOBUT,"+",RADIOTRANSTY);
	addbutton("Sink",x+20+2*180,y,180,20,RADIOBUT,RADIOBUT,"-",RADIOTRANSTY);
	y += 35;
	
	switch(transty){
	case "trans":
		addbutton("Please select which transition to generate data for:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); y += 50;
		addbutton("",x+20,y,modeldx,modeldytrans,CANVASBUT,CANVASBUT,-1,-1);	

		drawmodel(transcl,0,0,modeldx,modeldytrans,"frame","seltrans");
		addbutton("Cancel",width-105,height-45,90,30,canac,CANCELBUT2,0,-1);
		rightmenu(width-rwid,125,"seltrans");
		break;
		
	case "+":
		addbutton("Source data is generated (i.e. transitions in which individuals enter the system).",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); 
		addbutton("Cancel",width-205,height-45,90,30,canac,CANCELBUT2,0,-1);
		if(page == INFERENCEPAGE && datatype == "move"){
			addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
		}
		else addbutton("Next",width-105,height-45,90,30,NEXTGENAC2,NEXTBUT,0,-1);
		break;
		
	case "-":
		addbutton("Sink data is generated (i.e. transitions in which individuals leave the system).",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); 
		addbutton("Cancel",width-205,height-45,90,30,canac,CANCELBUT2,0,-1);
		if(page == INFERENCEPAGE && datatype == "move"){
			addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
		}
		else addbutton("Next",width-105,height-45,90,30,NEXTGENAC2,NEXTBUT,0,-1);
		break;
	}	
}

function loadsiminit(textFromFile)                       // Loads the initial state for simulation
{
	var li, co;
	var lines = textFromFile.split('\n');
	
	canover = -1; tableyfr = 0;
	ncol = ncla; ncoldef = ncol; ncoldefmax = ncol;
	nrow = 0;

	row = [];
	for(li = 0; li < lines.length; li++){
		temp = lines[li].split('\t');
		for(j = 0; j < temp.length; j++) temp[j] = temp[j].trim();
		if(!(temp.length == 1 && temp[0].length == 0)){
			row[nrow]=[];
			row[nrow][0] = temp[0];
			for(cl = 0; cl < ncla-1; cl++) if(cla[cl].ncomp == 1) row[nrow][1+cl] = cla[cl].comp[0].name;

			for(j = 1; j < temp.length; j++){
				na = temp[j];
				for(cl = 0; cl < ncla-1; cl++){
					for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == na) break;
					if(c < cla[cl].ncomp) break;
				}
				if(cl < ncla-1) row[nrow][cl+1] = na;
				else{ errormsg = "On line "+(li+1)+": '"+na+"' is not a compatment name"; return 1;}
			}
			
			for(cl = 0; cl < ncla-1; cl++){
				if(row[nrow][cl] == undefined){ errormsg = "On line "+(li+1)+": Not all compartments have been specified"; return 1;}
			}
			
			nrow++;
		}
	}
	calcrowwidth();

	colname[0] = "ID";
	for(cl = 0; cl < ncla-1; cl++) colname[1+cl] = cla[cl].name;
	setcolumns();
	showsimindinit = 1;
	return 0;
}

function showsimind()                                   // Shows the individuals in the simulation
{
	canover = -1; tableyfr = 0;
	ncol = ncla; ncoldef = ncol; ncoldefmax = ncol;
	nrow = simindinit.length;
	row=[];
	for(r = 0; r < nrow; r++){
		row[r]=[];
		row[r][0] = simindinit[r].id;
		for(cl = 0; cl < ncla-1; cl++) row[r][cl+1] = simindinit[r].state[cl];
	}
	calcrowwidth();
	
	colname[0] = "ID";
	for(cl = 0; cl < ncla-1; cl++) colname[1+cl] = cla[cl].name;
	setcolumns();
	showsimindinit = 1;
}

function setsimind()                                      // Sets simulated inidividuals from the table
{
	var r, cl, j, i;

	simindinit=[];
	for(r = 0; r < nrow; r++){
		simindinit[r]={id:row[r][0], state:[]};
		for(cl = 0; cl < ncla-1; cl++) simindinit[r].state[cl] = row[r][1+cl];
	}
	
	for(cl = 0; cl < ncla-1; cl++){
		for(j = 0; j < cla[cl].ncomp; j++) cla[cl].comp[j].loadinit = [];
	}
	
	for(i = 0; i < simindinit.length; i++){
		for(cl = 0; cl < ncla-1; cl++){
			val = simindinit[i].state[cl];
			j = 0; while(j < cla[cl].ncomp && cla[cl].comp[j].name != val) j++;
			if(j == cla[cl].ncomp){
				selectelement(i,cl+1,TABLEBUT);
				errmsg = "This is not a compartment"; 
				return 1;
			}
			cla[cl].comp[j].loadinit.push(simindinit[i].id);  
		}
	}
	return 0; 
}

function drawparaminit()                                // Draws the simulation parameter set
{
	var x, y, dyparam = 26;
	
	nob = 0;
	
	x = 0; y = 0;
	
	addob(x,y,OBTEXT2,"Please select the model parameter values used for running the simulation:",0); y += 50; 
	
	if(paramsim.length == 0){ addob(x+20,y,OBTEXT,"The model contains no parameters",0); y += 30;}
	else{
		for(p = 0; p < paramsim.length; p++){ addob(x+20,y,OBPARAM,p,0); y += dyparam;}
	}
	ytot = y;
	placeob();
}

function simsetup()                                      // Sets up the simulation xml file
{
	var j, nsimdata, simfr=[], possum=[];
	
	nchrun = 1;
	simdone = 1;
	 
	nsimdata = 0; simdata=[];
	for(cl = 0; cl < ncla; cl++){
		simdata[nsimdata] = {name:[], variety:"state", id:[], t:[], val:[], cl:cl, type:"binary", pos:[], posbinary:[], posexpression:[], sensitive:[], postestres:[], testname:"T1"};
			
		for(i = 0; i < cla[cl].ncomp; i++){
			simdata[nsimdata].pos[i] = cla[cl].comp[i].name;
			simdata[nsimdata].posbinary[i] = [];
			for(ii = 0; ii < cla[cl].ncomp; ii++){
				if(i == ii) simdata[nsimdata].posbinary[i][ii] = 1;		
				else simdata[nsimdata].posbinary[i][ii] = 0;	
			}
		}
		nsimdata++;
	}
		
	switch(siminit){
	case "Manualpop":			
		cl = simpopinitset;
		nind = 0;
		for(j = cla[cl].ncomp-1; j >= 0; j--){
			for(i = 0; i < cla[cl].comp[j].simpopinit; i++){
				for(cl2 = 0; cl2 < ncla; cl2++){
					if(cl2 == cl){
						simdata[cl2].id.push("Ind. "+(nind+1));
						simdata[cl2].t.push(tsimmin);
						simdata[cl2].val.push(cla[cl].comp[j].name);
					}
				}
				nind++;
			}
		}

		for(cl = 0; cl < ncla; cl++){
			if(cl != simpopinitset){
				if(cla[cl].name == "Time"){
					for(i = 0; i < nind; i++){
						simdata[cl].id.push("Ind. "+(i+1));
						simdata[cl].t.push(tsimmin);
						simdata[cl].val.push(cla[cl].comp[0].name);
					}
				}
				else{
					for(j = 0; j < cla[cl].ncomp; j++) simfr[j] = cla[cl].comp[j].simfracinit;
					
					for(i = 0; i < nind; i++){
						sum = 0; for(j = 0; j < cla[cl].ncomp; j++){ sum += simfr[j]; possum[j] = sum;}
						if(sum == 0){ alertsim("Error code EC13"); return;}
						
						z = Math.random()*sum; j = 0; while(j < cla[cl].ncomp && z > possum[j]) j++;
						if(j == cla[cl].ncomp){ alertsim("Error code EC14"); return;}
						simfr[j] -= 1/nind; if(simfr[j] < 0) simfr[j] = 0;
						
						simdata[cl].id.push("Ind. "+(i+1));
						simdata[cl].t.push(tsimmin);
						simdata[cl].val.push(cla[cl].comp[j].name);
					}
				}
			}
		}
		break;
		
	case "Load":
		nind = simindinit.length;	
		for(i = 0; i < nind; i++){
			for(cl = 0; cl < ncla; cl++){
				if(simindinit[i].id == ""){ alertsim("For initial individual ID blank"); return;}
					
				simdata[cl].id.push(simindinit[i].id);
				simdata[cl].t.push(tsimmin);
				if(cl < ncla-1){
					val = simindinit[i].state[cl];
					if(val == ""){
						changepage(SIMULATEPAGE,1,0);
						alertsim("The value for individual '"+simindinit[i].id+"' in classification '"+cla[cl].name+"' is not set. Please check the population initialisation.");
						return;
					}
					for(c = 0; c < cla[cl].ncomp; c++){ if(cla[cl].comp[c].name == val) break;}
					if(c == cla[cl].ncomp){
						alertsim("The value '"+val+"' for individual '"+simindinit[i].id+"' in classification '"+cla[cl].name+"' is not valid");
						return;
					}
					simdata[cl].val.push(val);
				}
				else simdata[cl].val.push(cla[cl].comp[0].name);	
			}
		}
		break;
	}
	
	converttoobs("sim");
	startinference(1,1);
}

function drawfilter(notcl)                                 // Draws the time filter
{
	var cl;
	
	nob = 0;
	x = 0; y = 20;
	
	for(cl = 0; cl < ncla; cl++) if(cla[cl].ncomp > 1 && cla[cl].name != "Time" && cl != notcl) break;
	if(cl < ncla){
		addob(x,y,OBTEXT,"Population:"); addob(x+110,y,OBRADIOWHICH); y += 30;

		switch(whichind){
		case "all":
			addob(x,y,OBTEXT2,"The selected transition is observed for individuals in all compartments.");
			break;
			
		case "sub":
			addob(x,y,OBTEXT2,"The selected transition is observed for individuals in the following subpopulation:"); y += 25;
			y += 20;
			
			x = 10; 
			for(cl = 0; cl < ncla; cl++){
				if(cla[cl].ncomp > 1 && cla[cl].name != "Time" && cl != notcl){ 
					dx = textwidth(cla[cl].name,"bold 16px arial")+10;
					if(x+dx + 100 > 700){ x = 10; y += 30;}
					//addob(x,y,OBSELIND,cl,dx);
					if(addingdata == 2) addob(x,y,OBSELIND2,cl,dx);
					else addob(x,y,OBSELIND,cl,dx);	
					x += dx+130;
				}
			}
			if(x > 10) y += 30;
			break;
		}
	}
	else{ 
		addob(x,y,OBTEXT,"Population: Entire"); y += 30;
		addob(x,y,OBTEXT2,"The selected transition is observed for individuals in all compartments.");
	}
	
	x = 0; y = 160;
	addob(x,y,OBTEXT,"Detection:"); addob(x+110,y,OBRADIOEVPD); y += 30;
	
	st = ""; if(whichind == "sub"){ for(cl = 0; cl < ncla; cl++) if(clagendata[cl] != "All") st += clagendata[cl]+",";}
	if(st != "") st = st.substring(0,st.length-1);
			
	switch(obspd){
	case "all":
		if(st == "") addob(x,y,OBTEXT2,"All transitions of selected type are observed.");
		else addob(x,y,OBTEXT2,"Transition on all individuals in "+st+" are observed.");
		break;
	
	case "set":
		if(st == "") addob(x,y,OBTEXT2,"Transitions observed with the following probability:");
		else addob(x,y,OBTEXT2,"Transitions in "+st+" are observed with the following probabbility:");
		if(page == SIMULATEPAGE) addob(x+20,y+35,OBMINMAX,"Probability: "+detectprob,14,200);
		else addob(x+20,y+35,OBPD,"Probability: "+datatemp.pd,17,200);
		break;
	}
		
	x = 0; y = 320;
	addob(x,y,OBTEXT,"Times:"); x += 70; y -= 4;
	addob(x,y,OBTEXT2,"From"); x += 54;
	
	dx = textwidth(tgenmin,tablefont);
	addob(x,y+3,OBMINMAX,tgenmin,9,dx+10);
	x += dx+20;
	addob(x,y,OBTEXT2,"to"); x += 30;
	dx = textwidth(tgenmax,tablefont); 
	addob(x,y+3,OBMINMAX,tgenmax,10,dx+10);
	x += dx+20;

	ytot = y;
	
	placeob();
	
	if(addingdata == 2) addcanbutton("",0,0,canvasdx,canvasdy,INACTIVEBUT,INACTIVEBUT,-1,-1);
}

function generatepopdata()                               // Generates population data from simulated data
{
	var st;

	ncol = 3; ncoldef = ncol; ncoldefmax = ncol; datashow = "table";

	nrow = 0; row=[];
	for(ti = 0; ti < tlist.length; ti++){
		row[ti]=[];
		row[ti][0] = tlist[ti]; row[ti][1] = poplist[ti]; row[ti][2] = tpre(errbarlist[ti],4);
		nrow++;
	}
		
	calcrowwidth();
	
	st = getclagenname();
	
	colname[0] = "Time"; colname[1] = st+" population"; colname[2] = "Standard deviation"; 
	setcolumns();
}

function generatederdata()                             // Generates population data from simulated data
{
	var st
	
	d = 0; while(d < derive.length && derive[d].name != dergensel) d++;
	if(d == derive.length) alertp("Error code EC15");
	
	ncol = 3+derive[d].dep.length; ncoldef = ncol; ncoldefmax = ncol; datashow = "table";

	nrow = 0; row=[];
	for(ti = 0; ti < tlist.length; ti++){
		var list = getder(dergensel,tlist[ti]); 
		for(j = 0; j < list.length; j++){
			row[nrow]=[];
			row[nrow][0] = tlist[ti]; row[nrow][1] = tpre(list[j].val,4);

			switch(errbar){
			case 0: eb = errbarscale*Math.sqrt(row[nrow][0]+0.5); break;
			case 1: eb = errbarfix; break;
			}
			row[nrow][2] = tpre(eb,4);
			
			for(k = 0; k < list[j].dep.length; k++) row[nrow][3+k] = list[j].dep[k]; 
			nrow++;
		}
	}
		
	calcrowwidth();
	
	colname[0] = "Time"; colname[1] = dergensel; colname[2] = "Standard deviation";
	for(k = 0; k < derive[d].dep.length; k++) colname[3+k] = derive[d].dep[k];

	setcolumns();
}

function getderval(name,t)
{
	var temp = getder(name,t);
	return temp[0].val;
}

function getder(name,t)                                   // Gets the values for derived quantites
{
	var list=[], d, ch;
	
	ch = 0;
	d = 0;
	i = -1; while(i < simres.DX-1 && simres.ch[ch].derivepl[d][i+1].t < t) i++;
	
	if(i >= 0 && i < res.DX-1){
		f = t - simres.ch[ch].derivepl[d][i].t; ff = simres.ch[ch].derivepl[d][i+1].t - t;			
		for(d = 0; d < simres.nderive; d++){
			var spl = splitname(simres.derive[d]);
			if(spl.name == name){
				list.push({val:(simres.ch[ch].derivepl[d][i].av*ff + simres.ch[ch].derivepl[d][i+1].av*f)/(f+ff), dep:spl.dep});
			}
		}
	}
	return list;
}

function drawpopderdata(type)                           // Draws generated population derived data
{
	var r, x, y, bound=[];
	
	cv = resultcv;
	
	cornx = menux+20; corny = 80;
	
	x = 10; y = 10;
	switch(type){
	case 0:
		addcanbutton("Population:",x,y,150,0,-1,TEXTBUT,-1,-1); x += 100; y -= 4;
		for(cl = 0; cl < ncla; cl++){
			if(cla[cl].ncomp > 1 || cl < ncla-2){
				dx = textwidth(cla[cl].name,"bold 16px arial")+10;
				if(x+dx + 100 > 700){ x = 110; y += 30;}
				
				var ops=[]; for(j = 0; j < cla[cl].ncomp; j++) ops.push(cla[cl].comp[j].name); ops.push("All");
				
				addcanbutton(cla[cl].name,x,y,150,20,-1,TEXTBUT2,-1,-1);
				gdropinfo.push({val:clagendata[cl], x:cornx+x+dx+3, y:corny+y+2, dx:95, dy:20, style:2, options:ops, click:"gendata", cl:cl});
				x += dx+130;
			}
		}
		break;
		
	case 1:
		addcanbutton("Derived quantity:",x,y,150,0,-1,TEXTBUT,-1,-1); x += 150; y -= 4;
		var ops=[]; for(j = 0; j < derive.length; j++) ops.push(derive[j].name);
		gdropinfo.push({val:dergensel, x:cornx+x+3, y:corny+y+2, dx:95, dy:20, style:2, options:ops, click:"dergen"});
		break;
	}
	
	if(x > 10) y += 30;
	
	axxmin = large; axxmax = -large;
	ma = 1;
	for(r = 0; r < nrow; r++){
		t = parseFloat(row[r][0]); if(t < axxmin) axxmin = t; if(t > axxmax) axxmax = t;
		bound[r] = geterrorbar(row[r][1],row[r][2]);
		if(bound[r].max > ma) ma = bound[r].max;
	}
	
	dx = axxmax - axxmin; if(dx == 0) dx = 1;
	axxmin -= 0.05*dx; axxmax += 0.05*dx; setxtics();
	axymin = 0; axymax = ma*1.1; w = setytics();
		
	graphframe(cornx,corny,40+w,50,20,80,w,"Time","Population","genpop");
	
	cv.clearRect(0,0,graphdx,graphdy);
	
	rr = 5;
	for(r = 0; r < nrow; r++){
		x = Math.floor(graphdx*(row[r][0]-axxmin)/(axxmax-axxmin));
		y = Math.floor(graphdy-graphdy*(row[r][1] -axymin)/(axymax-axymin));
		drawline(x-rr,y-rr,x+rr,y+rr,RED,THICKLINE);
		drawline(x-rr,y+rr,x+rr,y-rr,RED,THICKLINE);
	
		ymi = Math.floor(graphdy-graphdy*(bound[r].min-axymin)/(axymax-axymin));
		yma = Math.floor(graphdy-graphdy*(bound[r].max-axymin)/(axymax-axymin));
		drawline(x,ymi,x,yma,RED,THICKLINE);
		drawline(x-rr,ymi,x+rr,ymi,RED,THICKLINE);
		drawline(x-rr,yma,x+rr,yma,RED,THICKLINE);
	}	
}

function drawgenpopdata()                                // Sets up drawing of population data
{
	var poppro=[], timax = 200, ma, bound=[];

	cornx = menux+20; corny = 80;
	
	for(cl = 0; cl < ncla; cl++){ if(cla[cl].ncomp > 1) break;}
	
	x = 10; y = 10;
	if(cl < ncla){
		addcanbutton("Population:",x,y,150,0,-1,TEXTBUT,-1,-1);
		x += 100; y -= 4;
		for(cl = 0; cl < ncla; cl++){
			if(cla[cl].ncomp > 1){
				dx = textwidth(cla[cl].name,"bold 16px arial")+10;
				if(x+dx + 100 > 700){ x = 110; y += 30;}
				
				var ops=[]; for(j = 0; j < cla[cl].ncomp; j++) ops.push(cla[cl].comp[j].name); ops.push("All");
				
				addcanbutton(cla[cl].name,x,y,150,20,-1,TEXTBUT2,-1,-1);
				gdropinfo.push({val:clagendata[cl], x:cornx+x+dx+3, y:corny+y+2, dx:95, dy:20, style:2, 
				                options:ops, click:"gendata", cl:cl});
				x += dx+130;
			}
		}
		if(x > 10) y += 30;
	}

	x = 10; y += 5; 
	addcanbutton("Times:",x,y,150,0,-1,TEXTBUT,-1,-1); x += 90;
	addcanbutton("Periodic",x,y,130,22,CANRADIOBUT,CANRADIOBUT,0,CANRADIOGENT2); x += 130;
	addcanbutton("User defined",x,y,130,22,CANRADIOBUT,CANRADIOBUT,1,CANRADIOGENT2); 
	x += 150; y -= 4;
	
	switch(datagenttype2){
	case 0:
		addcanbutton("From",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 54;
		dx = textwidth(tgenmin,tablefont); addcanbutton(tgenmin,x,y-4,dx+10,30,MINMAXBUT,MINMAXBUT,0,9); x += dx+20;
		addcanbutton("to",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 30;
		dx = textwidth(tgenmax,tablefont); addcanbutton(tgenmax,x,y-4,dx+10,30,MINMAXBUT,MINMAXBUT,0,10); x += dx+20;
		addcanbutton("in steps of",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 90;
		dx = textwidth(dtgen,tablefont); addcanbutton(dtgen,x,y-4,dx+10,30,MINMAXBUT,MINMAXBUT,0,11); x += dx+20;
		break;
	
	case 1:
		addcanbutton("Time points:",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 95;
		dx = textwidth(tgenuserdef,tablefont);
		if(x + dx > canvasdx-20) dx = canvasdx-20-x;
		addcanbutton(tgenuserdef,x,y-5,dx+10,30,GENUSERBUT,GENUSERBUT,-1,-1);	
		break;
	}
	
	x = 10; y += 33;
	addcanbutton("Error bars:",x,y,150,0,-1,TEXTBUT,-1,-1); x += 90;
	addcanbutton("Square root",x,y,130,22,CANRADIOBUT,CANRADIOBUT,0,CANRADIOERRBAR); x += 130;
	addcanbutton("Fixed",x,y,130,22,CANRADIOBUT,CANRADIOBUT,1,CANRADIOERRBAR); 
	x += 150; y -= 4;

	switch(errbar){
	case 0:
		addcanbutton("Scale",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 54;
		dx = textwidth(errbarscale,tablefont); addcanbutton(errbarscale,x,y-4,dx+10,30,MINMAXBUT,MINMAXBUT,0,15); x += dx+20;
		break;
	case 1:
		addcanbutton("Value",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 54;
		dx = textwidth(errbarfix,tablefont); addcanbutton(errbarfix,x,y-4,dx+10,30,MINMAXBUT,MINMAXBUT,0,16); x += dx+20;
		break;
	}
	
	tlist=[];
	switch(datagenttype2){
	case 0: for(tt = parseFloat(tgenmin); tt <= parseFloat(tgenmax); tt += parseFloat(dtgen)) tlist.push(tt); break;
	case 1: tlist = tgenuserdef.split(","); for(j = 0; j < tlist.length; j++) tlist[j] = parseFloat(tlist[j]); break;
	}
	
	ma = 1;
	
	for(ti = 0; ti < tlist.length; ti++){
		t = tlist[ti]; if(t == simres.tmin) t += tiny;
		poplist[ti] = getpop(res.ch[0].sampev[0],t);
		switch(errbar){
		case 0: 
			errbarlist[ti] = errbarscale*Math.sqrt(poplist[ti]+0.5); break;
		case 1: errbarlist[ti] = errbarfix; break;
		}
		
		bound[ti] = geterrorbar(poplist[ti],errbarlist[ti]);
		if(bound[ti].max > ma) ma = bound[ti].max;
	}
	
	for(ti = 0; ti < timax; ti++){
		poppro[ti] = getpop(res.ch[0].sampev[0],simres.tmin+ti*(simres.tmax-simres.tmin)/timax+tiny);
		if(poppro[ti] > ma) ma = poppro[ti];
	}
	
	axxmin = simres.tmin; axxmax = simres.tmax+tiny; setxtics();
	axymin = 0; axymax = ma*1.1; w = setytics();
		
	graphframe(cornx,corny,40+w,50,20,130,w,"Time","Population","genpop");
	
	cv.clearRect(0,0,graphdx,graphdy);
	cv.beginPath(); 
	for(ti = 0; ti < timax; ti++){
		t = simres.tmin+ti*(simres.tmax-simres.tmin)/timax;
		x = Math.floor(graphdx*(t-axxmin)/(axxmax-axxmin));
		y = Math.floor(graphdy-graphdy*(poppro[ti]-axymin)/(axymax-axymin));
		if(ti == 0) cv.moveTo(x,y);
		else cv.lineTo(x,y);
	}
	cv.lineWidth = 3;
	cv.strokeStyle = BLUE;
	cv.stroke();
	
	r = 5;
	for(ti = 0; ti < tlist.length; ti++){
		x = Math.floor(graphdx*(tlist[ti]-axxmin)/(axxmax-axxmin));
		y = Math.floor(graphdy-graphdy*(poplist[ti] -axymin)/(axymax-axymin));
		
		drawline(x-r,y-r,x+r,y+r,RED,THICKLINE);
		drawline(x-r,y+r,x+r,y-r,RED,THICKLINE);
		
		ymi = Math.floor(graphdy-graphdy*(bound[ti].min -axymin)/(axymax-axymin));
		yma = Math.floor(graphdy-graphdy*(bound[ti].max -axymin)/(axymax-axymin));
		drawline(x,ymi,x,yma,RED,THICKLINE);
		drawline(x-r,ymi,x+r,ymi,RED,THICKLINE);
		drawline(x-r,yma,x+r,yma,RED,THICKLINE);
	}	
}

function drawgenderdata()                                // Sets up drawing of derived data
{
	var poppro=[], timax = 200, ma, bound=[];

	cornx = menux+20; corny = 80;
	
	for(cl = 0; cl < ncla; cl++){ if(cla[cl].ncomp > 1) break;}
	
	x = 10; y = 10;
	addcanbutton("Derived quantity:",x,y,150,0,-1,TEXTBUT,-1,-1);
	x += 150; y -= 4;

	var ops=[]; for(j = 0; j < simres.nderive; j++) ops.push(simres.derive[j]);
	gdropinfo.push({val:dergensel, x:cornx+x+3, y:corny+y+2, dx:95, dy:20, style:2, options:ops, click:"dergen"});
	
	y += 30;
	x = 10; y += 5;
	addcanbutton("Times:",x,y,150,0,-1,TEXTBUT,-1,-1); x += 90;
	addcanbutton("Periodic",x,y,130,22,CANRADIOBUT,CANRADIOBUT,0,CANRADIOGENT2); x += 130;
	addcanbutton("User defined",x,y,130,22,CANRADIOBUT,CANRADIOBUT,1,CANRADIOGENT2);
	x += 150; y -= 4;
	
	switch(datagenttype2){
	case 0:
		addcanbutton("From",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 54;
		dx = textwidth(tgenmin,tablefont); addcanbutton(tgenmin,x,y-4,dx+10,30,MINMAXBUT,MINMAXBUT,0,9); x += dx+20;
		addcanbutton("to",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 30;
		dx = textwidth(tgenmax,tablefont); addcanbutton(tgenmax,x,y-4,dx+10,30,MINMAXBUT,MINMAXBUT,0,10); x += dx+20;
		addcanbutton("in steps of",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 90;
		dx = textwidth(dtgen,tablefont); addcanbutton(dtgen,x,y-4,dx+10,30,MINMAXBUT,MINMAXBUT,0,11); x += dx+20;
		break;
	
	case 1:
		addcanbutton("Time points:",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 100;
		dx = textwidth(tgenuserdef,tablefont);
		if(x + dx > canvasdx-20) dx = canvasdx-20-x;
		addcanbutton(tgenuserdef,x,y-5,dx+10,30,GENUSERBUT,GENUSERBUT,-1,-1);	
		break;
	}
	
	x = 10; y += 33;
	addcanbutton("Error bars:",x,y,150,0,-1,TEXTBUT,-1,-1); x += 90;
	addcanbutton("Square root",x,y,130,22,CANRADIOBUT,CANRADIOBUT,0,CANRADIOERRBAR); x += 130
	addcanbutton("Fixed",x,y,130,22,CANRADIOBUT,CANRADIOBUT,1,CANRADIOERRBAR);
	x += 150; y -= 4;

	switch(errbar){
	case 0:
		addcanbutton("Scale",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 54;
		dx = textwidth(errbarscale,tablefont); addcanbutton(errbarscale,x,y-4,dx+10,30,MINMAXBUT,MINMAXBUT,0,15); x += dx+20;
		break;
	case 1:
		addcanbutton("Value",x,y,150,0,-1,TEXTBUT2,-1,-1); x += 54;
		dx = textwidth(errbarfix,tablefont); addcanbutton(errbarfix,x,y-4,dx+10,30,MINMAXBUT,MINMAXBUT,0,16); x += dx+20;
		break;
	}
	
	tlist=[];
	switch(datagenttype2){
	case 0: for(tt = parseFloat(tgenmin); tt <= parseFloat(tgenmax); tt += parseFloat(dtgen)) tlist.push(tt); break;
	case 1: tlist = tgenuserdef.split(","); for(j = 0; j < tlist.length; j++) tlist[j] = parseFloat(tlist[j]); break;
	}
	
	ma = -large;
	for(ti = 0; ti < tlist.length; ti++){
		t = tlist[ti]; if(t == simres.tmin) t += tiny;
		
		poplist[ti] = getderval(dergensel,t);
		switch(errbar){
		case 0: errbarlist[ti] = errbarscale*Math.sqrt(poplist[ti]+0.00001); break;
		case 1: errbarlist[ti] = errbarfix; break;
		}
		
		bound[ti] = geterrorbar(poplist[ti],errbarlist[ti]);
		if(bound[ti].max > ma) ma = bound[ti].max;
	}
	
	for(ti = 0; ti < timax; ti++){
		//poppro[ti] = getpop(res.ch[0].sampev[0],simres.tmin+ti*(simres.tmax-simres.tmin)/timax+tiny);
		poppro[ti] = getderval(dergensel,simres.tmin+ti*(simres.tmax-simres.tmin)/timax+tiny);
		
		if(poppro[ti] > ma) ma = poppro[ti];
	}
	
	axxmin = simres.tmin; axxmax = simres.tmax+tiny; setxtics();
	axymin = 0; axymax = ma*1.1; w = setytics();
		
	graphframe(cornx,corny,40+w,50,20,130,w,"Time",dergensel,"gender");
	
	cv.clearRect(0,0,graphdx,graphdy);
	cv.beginPath(); 
	for(ti = 0; ti < timax; ti++){
		t = simres.tmin+ti*(simres.tmax-simres.tmin)/timax;
		x = Math.floor(graphdx*(t-axxmin)/(axxmax-axxmin));
		y = Math.floor(graphdy-graphdy*(poppro[ti]-axymin)/(axymax-axymin));
		if(ti == 0) cv.moveTo(x,y);
		else cv.lineTo(x,y);
	}
	cv.lineWidth = 3;
	cv.strokeStyle = BLUE;
	cv.stroke();
	
	r = 5;
	for(ti = 0; ti < tlist.length; ti++){
		x = Math.floor(graphdx*(tlist[ti]-axxmin)/(axxmax-axxmin));
		y = Math.floor(graphdy-graphdy*(poplist[ti] -axymin)/(axymax-axymin));
		
		drawline(x-r,y-r,x+r,y+r,RED,THICKLINE);
		drawline(x-r,y+r,x+r,y-r,RED,THICKLINE);
		
		ymi = Math.floor(graphdy-graphdy*(bound[ti].min -axymin)/(axymax-axymin));
		yma = Math.floor(graphdy-graphdy*(bound[ti].max -axymin)/(axymax-axymin));
		drawline(x,ymi,x,yma,RED,THICKLINE);
		drawline(x-r,ymi,x+r,ymi,RED,THICKLINE);
		drawline(x-r,yma,x+r,yma,RED,THICKLINE);
	}	
	
	for(d = 0; d < derive.length; d++) if(derive[d].name == dergensel) break;
	
	x = width-rwid; y = 275;
	for(j = 0; j < derive[d].dep.length; j++){
		addbutton(derive[d].dep[j],x,y,0,0,-1,SMALLTEXTBUT,BLACK,-1); y += 20;
		y += 40;
	}
}

function inserttimes(y)                                 // Inputs for how time sampling is performed
{
	x = 10; y = 330;
	addob(x,y,OBTEXT,"Times:"); addob(x+110,y,OBTIMEP); y += 30;	
	switch(datagenttype){
	case 0:
		addob(x,y,OBTEXT2,"From"); x += 54;
		dx = textwidth(tgenmin,tablefont); addob(x,y+3,OBMINMAX,tgenmin,9,dx+10); x += dx+20;
		addob(x,y,OBTEXT2,"to"); x += 30;
		dx = textwidth(tgenmax,tablefont); addob(x,y+3,OBMINMAX,tgenmax,10,dx+10); x += dx+20;
		addob(x,y,OBTEXT2,"in steps of"); x += 90;
		dx = textwidth(dtgen,tablefont); addob(x,y+3,OBMINMAX,dtgen,11,dx+10);	
		break;
		
	case 1:
		addob(x,y,OBTEXT2,"Initial time: "+simres.tmin);
		break;
		
	case 2:
		addob(x,y,OBTEXT2,"Time points:"); x += 100;
		dx = textwidth(tgenuserdef,tablefont);
		addcanbutton(tgenuserdef,x,y-5,dx+10,30,GENUSERBUT,GENUSERBUT,-1,-1);	
		break;
		
	case 3: 
		addob(x,y,OBTEXT2,"Transition:"); x += 80;
		te = "Select";
		var postrans=[], clop=[], trop=[], cl, tr, ci, cf, transni, transnf;
		for(cl = 0; cl < ncla; cl++){
			for(tr = 0; tr < cla[cl].ntra; tr++){
				ci = cla[cl].tra[tr].i; if(ci < 0) transni = "+"; else transni = cla[cl].comp[ci].name;
				cf = cla[cl].tra[tr].f; if(cf < 0) transnf = "-"; else transnf = cla[cl].comp[cf].name;
				postrans.push(transni+" → "+transnf);
				clop.push(cl); trop.push(tr);
				if(cl == transselcl && tr == transseltr) te = transni+" → "+transnf;
			}	
		}
		gdropinfo.push({val:te, x:cornx+x+10, y:corny+y+3, dx:selbutdx, dy:20, style:4, options:postrans, clop:clop, trop:trop, click:"transsel"});

		break;
	}
	return y;
}

function drawgendata()                                  // Pages related to generating data from simulation
{
	var x, y, j;
	
	nob = 0;
	x = 0; y = 10;

	switch(drawgennum){
	case -2:
		y = 30;
		for(cl = 0; cl < ncla; cl++) if(cla[cl].ncomp > 1) break;
		if(cl < ncla){
			addob(x,y,OBTEXT,"Population:"); addob(x+110,y,OBRADIOWHICH); y += 30;

			switch(whichind){
			case "all": addob(x,y,OBTEXT2,"Individuals in all compartments are observed."); break;
			case "sub":
				addob(x,y,OBTEXT2,"Only individuals from the following subpopulation are observed:"); y += 25;
				y += 20;
				
				x = 10; 
				for(cl = 0; cl < ncla; cl++){
					if(cla[cl].ncomp > 1){
						dx = textwidth(cla[cl].name,"bold 16px arial")+10;
						if(x+dx + 100 > 700){ x = 10; y += 30;}
						addob(x,y,OBSELIND,cl,dx);
						x += dx+130;
					}
				}
				if(x > 10) y += 30;
				break;
			}
		}
		else addob(x,y,OBTEXT,"Population: Entire");
		
		x = 0; y = 180;
		addob(x,y,OBTEXT,"Detection:"); addob(x+110,y,OBRADIOPD); y += 30;
		
		st = ""; if(whichind == "sub"){ for(cl = 0; cl < ncla; cl++) if(clagendata[cl] != "All") st += clagendata[cl]+",";}
		if(st != "") st = st.substring(0,st.length-1);
				
		switch(obspd){
		case "all":
			if(st == "") addob(x,y,OBTEXT2,"All individuals in the entire population are observed.");
			else addob(x,y,OBTEXT2,"All individuals in "+st+" are observed.");
			break;
		
		case "set":
			if(st == "") addob(x,y,OBTEXT2,"Individuals are observed with the following probability:");
			else addob(x,y,OBTEXT2,"Individuals in "+st+" are observed with the following probabbility:");
			
			addob(x+20,y+35,OBMINMAX,"Probability: "+detectprob,14,150);
			break;
		}
	
		y = inserttimes();
		break;
	
	case -1:
		addob(x,y,OBTEXT2,"For observed individuals, which classifications are actually measured? "); y += 25;
		
		x = 20; y += 25;
		
		for(cl = 0; cl < ncla-1; cl++){
			if(cl == 0 || cla[cl].ncomp > 1){
				w = textwidth(cla[cl].name,INPUTFONT)+50;
				if(x+w > width-menux-30){ x = 20; y += 50;} 
				addob(x,y,OBCHOOSECLA,cl,w); 
				x += w + 30;
			}
		}
		if(x > 20){ x = 0; y += 50;}
		y += 10;
		break;
	
	default:
		cl = sellist[drawgennum]; clgl = cl;
		dg = datagen[cl];
	
		addob(x,y,OBTEXT,"Type of observation model:"); addob(x+250,y,OBRADIOGEN,CANRADIOCHECK2,cla[cl].ncomp); y += 30;

		switch(dg.type){
		case "simple":
			addob(x+15,y,OBTEXT2,"Here the data D is the compartment name in which the individual resides."); y += 25;
			y += 10;
			
			dg.pos.length = cla[cl].ncomp;
			for(j = 0; j < cla[cl].ncomp; j++){
				dg.pos[j] = cla[cl].comp[j].name;
				dg.posref[j] = j;
			}
			break;
			
		case "binary":
			addob(x+15,y,OBTEXT2,"Define the data output D for each possible state:"); y += 25;
			y += 10;
			
			dg.pos.length = cla[cl].ncomp;
			for(j = 0; j < cla[cl].ncomp; j++){
				if(dg.posedit[j] == undefined) dg.posedit[j] = cla[cl].comp[j].name;
				dg.posref[j] = j;
			}
			break;
			
		case "test":
			dg.pos.length = 2; dg.pos[0] = "1"; dg.pos[1] = "0";
			
			addob(x+15,y,OBTEXT2,"Generate diagnotic test data assuming the test is sensitive to a defined set of"); y += 25;
			addob(x+15,y,OBTEXT2,"compartments. A sensitivity Se and specificity Sp account for inaccuracies in the test."); y += 25;
			y += 40;
			
			addob(x,y,OBTEXT,"Test sensitive to:");
			x = 160; y -= 6;
			for(k = 0; k < cla[cl].ncomp; k++){
				if(dg.sensitive[k] == undefined){ if(k == 0) dg.sensitive[k] = 1; else dg.sensitive[k] = 0;}
			
				if(dg.sensitive[k] == 1) st = cla[cl].comp[k].name+" ✓";
				else st = cla[cl].comp[k].name+" ✖"; 
				w = textwidth(st,"20px georgia")
			
				if(x+w > tablewidth-40){ x = 160; y += 50;}
				addob(x,y,OBCOMP3,k,st,SENSAC);
				x += w+50;		
			}
			y += 80;
			
			x = 0; addob(x,y,OBSESP); y += 20;
			break;
		}
		
		y += 35;
		
		if(dg.type == "test") istest = 1; else istest = 0;
		
		na = dg.testname;
		
		var wmax = 0;
		for(j = 0; j < dg.pos.length; j++){
			w = textwidth(dg.pos[j],"bold 16px arial");
			if(w > wmax) wmax = w;
		}

		for(j = 0; j < dg.pos.length; j++){	
			ysta = y;
			x = wmax+105;
				
			for(k = 0; k < cla[cl].ncomp; k++){
				st = "Pr(D|"+cla[cl].comp[k].name+") = ";
				
				switch(dg.type){
				case "test":
					switch(j){
					case 0:
						if(dg.sensitive[k] == 1) st += "Se";
						else st += "1-Sp";
						break;
					
					case 1:
						if(dg.sensitive[k] == 1) st += "1-Se";
						else  st += "Sp";
						break;
					}	
					w = textwidth(st,"16px georgia")
					if(x+w > tablewidth-40){ x = wmax+105; y += 30;}
					addob(x,y,OBCOMP4,j,k,st,-1);
					x += w+30;
					break;
				
				case "simple": case "binary":
					if(k == dg.posref[j]){
						st = cla[cl].comp[k].name;
					
						w = textwidth(st,"20px georgia")
						if(x+w > tablewidth-40){ x = wmax+105; y += 50;}
						addob(x,y,OBCOMP8,j,k,st);
						x += w+50;
					}
					break;
				}
			}
			
			if(dg.type == "binary") addob(20,Math.floor((ysta+y)/2),OBTEXTEDIT,'D="'+dg.posedit[j]+'"',wmax+48,j);
			else addob(25,Math.floor((ysta+y)/2+5),OBTEXT,'D="'+dg.pos[j]+'"',0);

			ma = 5;
			addob(wmax+85,ysta-ma,OBBRACKET,y-ysta+30+2*ma);
			
			y += 60;
		}
		break;
	}
	
	ytot = y;
	
	placeob();
}

function generatestatedata()                            // Generate state data from simulated data
{
	var j, ch = 0, s = 0, stat=[], cl, tr;
	
	ncol = 2+sellist.length; ncol++;
	if(obspd == "all"){  // If all individuals detected and no births and deaths then capture not needed.
		for(cl =0; cl < ncla; cl++){
			for(tr = 0; tr < cla[cl].ntra; tr++) if(cla[cl].tra[tr].i < 0 || cla[cl].tra[tr].f < 0) break;
			if(tr < cla[cl].ntra) break;
		}
		if(cl == ncla) ncol--;
	}	
	ncoldef = ncol;
	
	res = simres;

	nrow = 0; row=[];
		
	if(datagenttype == 3){  // Data at transitions
		cl = transselcl; tr = transseltr;
		w = res.ch[ch].sampev[s];
		for(i = 0; i < w.nind; i++){
			w = res.ch[ch].sampev[s];
			wi = w.ind[i];
			c = NOTALIVE;
			for(e = 0; e < wi.nev; e++){
				cf = wi.evc[e];

				fl = 0;
				if(c == NOTALIVE){ if(cla[cl].tra[tr].i < 0) fl++;}
				else{ if(cla[cl].tra[tr].i == compval[c][cl]) fl++;}
					
				if(cf == NOTALIVE){ cc = c; if(cla[cl].tra[tr].f < 0) fl++;}
				else{ cc = cf; if(cla[cl].tra[tr].f == compval[cf][cl]) fl++;}
				
				if(fl == 2 && wi.evt[e] != res.tmin) addstatemeas(i,cc,wi.evt[e]);		
				c = cf;
			}
		}
	}
	else{
		tlist=[];
		switch(datagenttype){
		case 0: for(tt = parseFloat(tgenmin); tt <= parseFloat(tgenmax); tt += parseFloat(dtgen)) tlist.push(tt); break;
		case 1: tlist.push(simres.tmin); break;
		case 2: tlist = tgenuserdef.split(","); for(j = 0; j < tlist.length; j++) tlist[j] = parseFloat(tlist[j]); break;
		}

		for(j = 0; j < tlist.length; j++){
			t = tlist[j];
			w = res.ch[ch].sampev[s];
			for(i = 0; i < w.nind; i++){
				wi = w.ind[i];
				c = NOTALIVE; e = 0; while(e < wi.nev && (wi.evt[e] < t || (wi.evt[e] == t && t == simres.tmin))){ c = wi.evc[e]; e++;} 
				addstatemeas(i,c,t);			
			}
		}
	}
	calcrowwidth()
	colname[0] = "ID"; colname[1] = "Time"; 
	for(k = 0; k < sellist.length; k++){
		var st = cla[sellist[k]].name; if(datagen[sellist[k]].type == "test") st += " test";
		colname[2+k] = st;	
	}
	if(2+sellist.length < ncoldef) colname[2+sellist.length] = "Capture name";
	setcolumns();
}

function addstatemeas(i,c,t)                             // Adds a state measurement onto the generated data
{
	var cl, fl, indi;
	
	fl = 0;
	if(c != NOTALIVE){ // works out if the individual is captured
		for(cl = 0; cl < ncla-1; cl++){
			if(clagendata[cl] != "All" && clagendata[cl] != cla[cl].comp[compval[c][cl]].name) break;
		}
		if(cl == ncla-1){ if(Math.random() < parseFloat(detectprob)) fl = 1;}
	}
	
	if(fl == 1){
		row[nrow]=[];
		indi = indsim.ind[i]; 
		if(!indi) row[nrow][0] = "Ind. "+(i+1); else row[nrow][0] = indi.id; // New individuals
		row[nrow][1] = t;
		for(k = 0; k < sellist.length; k++){
			cl = sellist[k];
			dg = datagen[cl];
			switch(dg.type){
			case "simple": row[nrow][2+k] = cla[cl].comp[compval[c][cl]].name; break;
			case "binary": row[nrow][2+k] = dg.posedit[compval[c][cl]]; break;
			case "test":
				if(dg.sensitive[compval[c][cl]] == 1){
					if(Math.random() < Segen) row[nrow][2+k] = 1;
					else row[nrow][2+k] = 0;
				}
				else{
					if(Math.random() < Spgen) row[nrow][2+k] = 0;
					else row[nrow][2+k] = 1;
				}
				break;
			}
		}
		if(2+sellist.length < ncoldef) row[nrow][2+sellist.length] = getcapname(t);
		nrow++;
	}
}

function getcapname(t)                                   // Makes up the capture name
{
	var cl, st;
	
	st = "Cap:t="+t;
	if(whichind == "sub"){
		for(cl = 0; cl < ncla; cl++){
			if(clagendata[cl] != "All") st += ","+clagendata[cl];
		}
	}
	return st;
}

function getclagenname()                                 // Gets filter name based on clagendata
{
	var st = ""; for(cl = 0; cl < ncla; cl++){ if(clagendata[cl] != "All"){ if(st != "") st += ","; st += clagendata[cl];}}
	if(st == "") st = "All";
	return st;
}

function generatecapturedata()                           // Generates capture data from simulated data
{
	var j, tlist=[]; 
	
	switch(datagenttype){
	case 0: for(tt = parseFloat(tgenmin); tt <= parseFloat(tgenmax); tt += parseFloat(dtgen)) tlist.push(tt); break;
	case 1: tlist.push(simres.tmin); break;
	case 2: tlist = tgenuserdef.split(","); for(j = 0; j < tlist.length; j++) tlist[j] = parseFloat(tlist[j]); break;
	}
	
	st = getclagenname();
	
	ncol = 4; ncoldef = ncol; datashow = "table";
	nrow = 0; row=[];
	for(j = 0; j < tlist.length; j++){
		t = tlist[j];
		row[nrow] = [];
		row[nrow][0] = getcapname(t); row[nrow][1] = st; row[nrow][2] = t; row[nrow][3] = detectprob;  
		nrow++;
	}
	calcrowwidth();
	colname[0] = "Capture name"; colname[1] = "Compartments"; colname[2] = "Time"; colname[3] = "Detection probability";
	setcolumns();
}

function generateeventdata()                             // Generates event data from simulated data
{
	var i, ch = 0, s = 0, c, cf, cc, cl, tr, wi, fl, pd;
	
	cl = transcl;
	ressa = simres.ch[0].sampev[0];
	 
	ncol = 2; ncoldef = ncol; ncoldefmax = ncol; datashow = "table";

	nrow = 0; row=[];

	for(i = 0; i < ressa.nind; i++){
		wi = ressa.ind[i];				
		c = NOTALIVE; 
		for(e = 0; e < wi.nev; e++){
			fl = 0;	
			cf = wi.evc[e];
			if(wi.evt[e] != simres.tmin && wi.evt[e] != simres.tmax){
				switch(transty){
				case "+":
					if(c == NOTALIVE) fl = 2;
					break;
					
				case "-":
					if(cf == NOTALIVE) fl = 2;
					break;
					
				case "trans":
					if(cla[cl].comp[compval[c][cl]].name == transni) fl++;
					if(cla[cl].comp[compval[cf][cl]].name == transnf) fl++;
					break;
				}
			}
			
			if(fl == 2){
				cc = c; if(cc == NOTALIVE) cc = cf;
				if(checkfil(cc,wi.evt[e]) == 1){
					if(obspd == "all") pd = 1;
					else pd = detectprob;
	
					if(Math.random() < pd){ 
						row[nrow] = [];
						indi = indsim.ind[i]; 
						if(!indi) row[nrow][0] = "Ind. "+(i+1); else row[nrow][0] = indi.id; // New individual
						row[nrow][1] = wi.evt[e];
						nrow++;
					}
				}
			}
			c = cf;
		}
	}
		
	calcrowwidth();
	colname[0] = "ID"; colname[1] = "Time";
	setcolumns();	
}

function checkfil(c,t)                                 // Checks compartment c at time t agrees with filter
{
	var cl;
	if(t < tgenmin || t > tgenmax) return 0;
	for(cl = 0; cl < ncla; cl++){ 
		if(clagendata[cl] != "All" && clagendata[cl] != cla[cl].comp[compval[c][cl]].name) return 0;
	}
	return 1;
}

function gendatastart(val)                             // Initialises generation of data
{
	var cl, i;

	for(cl = 0; cl < ncla; cl++) clagendata[cl] = "All";
 
	whichind = "all"; obspd = "all"; detectprob = 1;
 
	tgenmin = simres.tmin; tgenmax = simres.tmax;
 
	switch(val){
	case 0: case 2: case 3:
		tgenuserdef = ""; for(t = tgenmin; t < tgenmax; t += parseFloat(dtgen)) tgenuserdef += t+", ";
		tgenuserdef = tgenuserdef.substr(0,tgenuserdef.length-2);
	}
	
	switch(val){
	case 0:           // state data
		gentype = 0;
		drawgennum = -2;
		for(cl = 0; cl < ncla; cl++){
			datagen[cl] = {type:"simple", pos:[], posedit:[],  posref:[], sensitive:[]};
			if(cl == 0 || cla[cl].ncomp > 1) selclass[cl] = 1; else selclass[cl] = 0;
		}
		dtgen = ((simres.tmax-simres.tmin)/10).toPrecision(3);
		detectprob = 1;
		transselcl = -1;
		break;

	case 1:            // event data
		gentype = 1; drawgennum = 0; transcl = 0; transty = "trans";
		break;
		
	case 2:            // population data	
		gentype = 2; drawgennum = 0;
		break;
	
	case 3:            // Derived data
		gentype = 3; drawgennum = 0; dergensel = simres.derive[0]; 
		break;
	}
}

function changegennum(val)                               // Changes the page under data generation
{
	canover = -1; tableyfr = 0; drawgennum = val;	
}
