function addingdatabuts()                                 // Adds the buttons on the data page
{
	x = menux+tab; y = 30;

	res = datares;

	switch(addingdata){
	case 1: // Adding data
	case 2: // Viewing data
	case 3: // Editing data
		switch(datashow){
		case "pop":
			if(addingdata == 2) addingdata = 3;
			addbutton("Define population",x,y,0,0,-1,TITLEBUT,-1,-1); y += 40;	
			drawpopderdata(0);
			
			if(addingdata == 1) addbutton("Back",width-205,height-45,90,30,BACKAC3,BACKBUT,0,-1);
			else addbutton("Back",width-205,height-45,90,30,BACKBUT,BACKBUT,0,-1); 
			addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
			break;
			
		case "der":
			if(addingdata == 2) addingdata = 3;
			addbutton("Define derived",x,y,0,0,-1,TITLEBUT,-1,-1); y += 40;	
			
			drawpopderdata(1);
			
			if(addingdata == 1) addbutton("Back",width-205,height-45,90,30,BACKAC3,BACKBUT,0,-1);
			else addbutton("Back",width-205,height-45,90,30,BACKBUT,BACKBUT,0,-1); 
			addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
			break;
			
		case "pdmodel":
			addbutton("Detection probability",x,y,0,0,-1,TITLEBUT,-1,-1);

			drawpd();
			cornx = menux+40; corny = 80;
			addbutton("",cornx,corny,tablewidth,tableheight,CANVASBUT,CANVASBUT,-1,-1);
			
			addbutton("Back",width-205,height-45,90,30,BACKBUT,BACKBUT,0,-1); 
			
			if(addingdata == 2) addbutton("Edit",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
			else addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
			break;
			
		case "table":
			switch(addingdata){
			case 1:
				switch(datatype){
				case "cappd": addbutton("Add dection probability data",x,y,0,0,-1,TITLEBUT,-1,-1); break;
				case "capid": addbutton("Add data associating individuals to captures",x,y,0,0,-1,TITLEBUT,-1,-1); break;
				case "state": addbutton("Add state data",x,y,0,0,-1,TITLEBUT,-1,-1); break;
				case "cap": addbutton("Add capture data",x,y,0,0,-1,TITLEBUT,-1,-1); break;
				case "trans": addbutton("Add transition data",x,y,0,0,-1,TITLEBUT,-1,-1); break;
				case "move": addbutton("Add move data",x,y,0,0,-1,TITLEBUT,-1,-1); break;
				case "presence": addbutton("Add presence data",x,y,0,0,-1,TITLEBUT,-1,-1); break;
				case "pop": addbutton("Add population data",x,y,0,0,-1,TITLEBUT,-1,-1); break;
				case "der": addbutton("Add derived data",x,y,0,0,-1,TITLEBUT,-1,-1); break;
				}

				y += 50;
				switch(datatype){
				case "pop":
					if(ncoldef == ncoldefmax){
						addbutton("Please edit entries if needed and press 'next' to continue.",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					
						addbutton("Back",width-305,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
						
						addbutton("Next",width-105,height-45,90,30,NEXTBUT,NEXTBUT,0,-1);
					}
					else{
						switch(ncoldef){
						case 0:
							addbutton("Please select the column representing measurement times:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						case 1:
							addbutton("Please select the column representing population estimates:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						case 2:
							addbutton("Please select the column representing standard deviations in these estimates:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						}
						if(ncoldef > 0) addbutton("Back",width-205,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-105,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
					}
					break;
				
				case "der":
					if(ncoldef == ncoldefmax){
						addbutton("Please edit entries if needed and press 'next' to continue. To add a filter press here.",x+20,y,800,20,ADDFILTAC,TEXTWITHLINKBUT,-1,-1);
									
						addbutton("Back",width-305,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
						
						addbutton("Next",width-105,height-45,90,30,NEXTBUT,NEXTBUT,0,-1);
					}
					else{
						switch(ncoldef){
						case 0:
							addbutton("Please select the column representing measurement times:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						case 1:
							addbutton("Please select the column representing derived quantity estimates:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						case 2:
							addbutton("Please select the column representing standard deviations in these estimates:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						case 3:
							addbutton("Please select the column representing a filter (to cancel click here).",x+20,y,800,20,REMADDFILTAC,TEXTWITHLINKBUT,-1,-1);
							break;
						}
						if(ncoldef > 0) addbutton("Back",width-205,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-105,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
					}
					break;
					
				case "presence":
					switch(ncoldef){
					case 0:
						addbutton("Please select the column representing individual IDs:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
						addbutton("Cancel",width-105,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
						break;
					
					case 1:
						addbutton("Please select the column representing observation times (or choose",x+20,y,550,0,-1,PARAGRAPHBUT,-1,-1);
						addbutton("'All'",x+505,y,28,20,CHOOSEALLAC,CHOOSETIMEBUT,-1,-1);
						addbutton("or a ",x+530,y,450,0,-1,PARAGRAPHBUT,-1,-1);
						addbutton("specific time",x+565,y,90,20,CHOOSETIMEBUT,CHOOSETIMEBUT,-1,-1);
						addbutton(")",x+653,y,450,0,-1,PARAGRAPHBUT,-1,-1);
						
						addbutton("Back",width-205,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-105,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
						break;

					case 2: 
						addbutton("Please edit entries if needed and press 'done' to finish.",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					
						addbutton("Back",width-305,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
						addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
						break;
					}
					break;
					
				case "cap":
					if(ncoldef == ncoldefmax){
						addbutton("Please edit entries if needed and press 'next' to continue.",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					
						addbutton("Back",width-305,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
						addbutton("Next",width-105,height-45,90,30,NEXTBUT,NEXTBUT,0,-1);
					}
					else{
						switch(ncoldef){
						case 0:
							addbutton("Please select the column defining capture names:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;

						case 1:
							addbutton("Please select the column defining compartments (click here for 'All')",x+20,y,800,20,ALLAC,TEXTWITHLINKBUT,-1,-1);
							break;
							
						case 2:
							addbutton("Please select the column defining times:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						}
						if(ncoldef > 0) addbutton("Back",width-205,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-105,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
					}
					break;
				
				case "cappd":
					if(ncoldef == ncoldefmax){
						addbutton("Please edit entries if needed and press 'done' to finish.",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					
						addbutton("Back",width-305,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
						addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
					}
					else{
						switch(ncoldef){
						case 0:
							addbutton("Please select the column defining capture names:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						
						case 1:
							addbutton("Please select the column defining detection probability expressions.",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						}
						if(ncoldef > 0) addbutton("Back",width-205,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-105,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
					}
					break;
					
				case "capid":
					if(ncoldef == ncoldefmax){
						addbutton("Please edit entries if needed and press 'done' to finish.",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					
						addbutton("Back",width-305,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
						addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
					}
					else{
						switch(ncoldef){
						case 0:
							addbutton("Please select the column defining individuals IDs:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						
						case 1:
							addbutton("Please select the column defining capture names:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						}
						if(ncoldef > 0) addbutton("Back",width-205,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-105,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
					}
					break;
					
				default:
					if(ncoldef == ncoldefmax){
						addbutton("Please edit entries if needed and press 'next' to continue.",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
					
						addbutton("Back",width-305,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
					    addbutton("Next",width-105,height-45,90,30,NEXTBUT,NEXTBUT,0,-1);
					}
					else{
						switch(ncoldef){
						case 0:
							addbutton("Please select the column representing individual IDs:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						case 1:
							switch(datatype){
							case "state":	
								addbutton("Please select the column representing observation times (or choose",x+20,y,550,0,-1,PARAGRAPHBUT,-1,-1);
								addbutton("'All'",x+505,y,28,20,CHOOSEALLAC,CHOOSETIMEBUT,-1,-1);
								addbutton("or a ",x+530,y,450,0,-1,PARAGRAPHBUT,-1,-1);
								addbutton("specific time",x+565,y,90,20,CHOOSETIMEBUT,CHOOSETIMEBUT,-1,-1);
								addbutton(")",x+653,y,450,0,-1,PARAGRAPHBUT,-1,-1);
								break;
								
							case "trans": 
								addbutton("Please select the column representing transition times:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
								break;
							
							case "move":
								addbutton("Please select the column representing the times corresponding to moves:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
								break;
							}
							break;

						case 2:
							addbutton("Please select the column giving information on classification:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
							break;
						}
						if(ncoldef > 0) addbutton("Back",width-205,height-45,90,30,TABBACKAC,BACKBUT,0,-1);
						addbutton("Cancel",width-105,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
					}
					break;
				}
				y += 40;
			
				drawtable();
				break;
			
			case 2:
				addbutton("Data Table",x,y,0,0,-1,TITLEBUT,-1,-1); y += 50;
				addbutton("Loaded data:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); y += 40;
					
				drawtable();
				addbutton("Back",width-205,height-45,90,30,BACKBUT,BACKBUT,0,-1);
				addbutton("Edit",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
				break;
			
			case 3:
				addbutton("Edit Data Table",x,y,0,0,-1,TITLEBUT,-1,-1); y += 50;
				addbutton("Edit data and click 'Done' when complete.",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); y += 40;
				
				drawtable();
				addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
				addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
				break
			}
			
			cornx = x+40; corny = y;
			addbutton("",cornx,corny,tablewidth,tableheight,CANVASBUT,CANVASBUT,-1,-1);
			tableyfrac = rowmax/nrow;
			if(tableyfrac < 1) addbutton("",x+40+ tablewidth+10,y+22,13,tableheight-22,SLIDEAC,YSLIDEBUT,-1,-1);

			tablexfrac = tablewidth/tottablewidth;
			if(tablexfrac < 1) addbutton("",x+40,y+tableheight+10,tablewidth,13,SLIDEAC,XSLIDEBUT,-1,-1);
			break;
			
		case "obsmodel": 
			if(datatemp.cl == -1) pickclass();
			else{
				addbutton("Observation Model",x,y,0,0,-1,TITLEBUT,-1,-1); y += 50;
		
				addbutton("Here we define the probability of observing the classifier given a particular compartmental state.",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
				y += 40;
		
				drawmeaning();
			
				x = menux+tab;
				addbutton("",x+20,y,obswidth,obsheight,CANVASBUT,CANVASBUT,-1,-1);
		
				tableyfrac = obsheight/ytot;
				if(tableyfrac < 1) addbutton("",menux+tab+20+obswidth+10,y,13,obsheight,SLIDEAC,YSLIDEBUT,-1,-1);
				
				switch(addingdata){
				case 1:
					addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
					addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
					addbutton("Back",width-305,height-45,90,30,BACKBUT,BACKBUT,0,-1);
					break;
					
				case 2:
					addbutton("Edit",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
					addbutton("Back",width-205,height-45,90,30,BACKBUT,BACKBUT,0,-1);
					break;
					
				case 3:
					addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
					addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
					break;
					
				}
			}
			break;
			
		case "seltrans":
			addbutton("Select transition",x,y,0,0,-1,TITLEBUT,-1,-1); y += 40;
	
			transop(x,y);
			break;
		
		case "transfilt":
			addbutton("Select filter",x,y,0,0,-1,TITLEBUT,-1,-1); y += 50;

			switch(transty){
			case "trans":
				addbutton("This data provides information about the transition: ",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
				addbutton(transni+" → "+transnf,x+390,y,800,0,-1,BOLDMEDBUT,-1,-1);
				break;
			
			case "+":
				addbutton("This data provides information about individual entering the system:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
				break;
				
			case "-":
				addbutton("This data provides information about individual leaving the system:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1);
				break;
			}
			
			y += 30;
			
			/*
			switch(addingdata){
			case 1: case 3: addbutton("Please select which individuals are observed (filtered by compartment) and over which time range:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); break;
			case 2: addbutton("Individuals are observed in the following compartments and time range:",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); break;
			}	
			y += 30;
				*/
				
			cornx = x+20; corny = y;
			addbutton("",x+20,y,tablewidth,tableheight,CANVASBUT,CANVASBUT,-1,-1);

			drawfilter(transcl);
				
			tableyfrac = tableheight/ytot;
			if(tableyfrac < 1) addbutton("",menux+tab+20+tablewidth+10,y,13,tableheight,SLIDEAC,YSLIDEBUT,-1,-1);
			
			switch(addingdata){
			case 1:
				addbutton("Back",width-305,height-45,90,30,BACKEVDATAAC3,BACKBUT,0,-1);
				addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
				addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
				break;
				
			case 2:
				addbutton("Edit",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
				addbutton("Back",width-205,height-45,90,30,BACKBUT,BACKBUT,0,-1);
				break;
					
			case 3:
				addbutton("Done",width-105,height-45,90,30,DONEAC,ADDDATABUT,0,-1);
				addbutton("Cancel",width-205,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
				break;
			}
		}
		break;
	}
}

function addop(te,ty,he)                                  // Places an add data button
{
	var dx = textwidth(te,addfont)+19;
	addbutton(te,x+20,y,dx,20,ADDDATAAC,ADDBUT2,OBSFILE,ty);
	addbutton("[?]",x+20+dx,y,15,20,HELPICONBUT,HELPICONBUT,he,-1);
}

function setpdvar(st)                                     // Sets probability of detection variables in data
{
	datatemp.npdvar = nvarlist;
	for(n = 0; n < nvarlist; n++){ datatemp.pdvar[n] = varlist[n]; datatemp.pdvardep[n] = vardeplist[n];}	
}

function sugname(na)                                      // Sets a suggested name for a data item
{
	var naa;
	if(datatemp.name == ""){
		loop = 1;
		do{
			if(loop == 1) naa = na; else naa = na+""+loop;
			for(d = 0; d < data.length; d++) if(d != dataselected && data[d].name == naa) break;
			if(d == data.length) break;
			loop++;
		}while(1 == 1);
		datatemp.name = naa;
	}		
}

function databuts()                                       // Draws button for the data page
{			
	x = menux+tab; y = 30;

	res = datares;
	
	if(data.length == 0) pagesubsub[INFERENCEPAGE][0] = 0;
	
	switch(pagesubsub[INFERENCEPAGE][0]){
	case 0:
		addbutton("Data sources",x,y,0,0,-1,TITLEBUT,22,-1); y += 50;
	
		if(data.length == 0) addbutton("There is currently no data loaded.",x+18,y,900,0,-1,PARAGRAPHBUT,1,-1);
		else{
			addbutton("",x+20,y,tablewidthbig,tableheightdata,CANVASBUT,CANVASBUT,-1,-1);
			drawdata();

			tableyfrac = tableheightdata/ytot;
			if(tableyfrac < 1) addbutton("",menux+tab+20+tablewidthbig+10,y,13,tableheightdata,SLIDEAC,YSLIDEBUT,-1,-1);
			addbutton("Next",width-105,height-45,90,30,NEXTINFAC,NEXTBUT,0,-1);	
		}
		
		x = menux+tab; y = height-55; dx = 130;
		addop("State","state",90); x += dx;
		addop("Presence","presence",95); x += dx;
		addop("Transition","trans",92); x += dx;
		addop("Move","move",97); x += dx;
		x = menux+tab; y += 20;
		addop("Capture","cap",91); x += dx;
		addop("Capture ID","capid",214); x += dx;
		addop("Capture PD","cappd",215); x += dx;
		addop("Population","pop",93); x += dx;
		addop("Derived","der",94); x += dx;
		break;
	
	case 1:
		indplot();
		addbutton("Next",width-105,height-45,90,30,NEXTINFAC,NEXTBUT,0,-1);	
		break;
	}
}

function indplot()                                        // Plots data on individuals
{
	var x, y;
	
	switch(page){
	case SIMULATEPAGE: addbutton("Individuals",menux+tab,30,0,0,-1,TITLEBUT,50,-1); break;
	case INFERENCEPAGE:
		if(pagesub[page] == 0) addbutton("Individuals",menux+tab,30,0,0,-1,TITLEBUT,51,-1);
		else addbutton("Individuals",menux+tab,30,0,0,-1,TITLEBUT,52,-1);
		break;
	}
	
	if(getnind() == 0 && uoflag == 0){
		addbutton("There is no individual-based data.",menux+tab+10,80,500,0,-1,PARAGRAPHBUT,-1,-1);
	}
	else{
		rightmenu(width-rwid,75,"ind");
		
		if(page == INFERENCEPAGE && pagesub[page] == 0) indv = "Timeline";
		else indv = res.indview;
				
		y = 50;
		switch(indv){
		case "Timeline":
			cornx = menux+20; corny = 80;
			addbutton("",cornx,corny,indtablewidth,indtableheight+indtablemar,CANVASBUT,CANVASBUT,-1,-1);
			canvasheight -= indtablemar;

			drawinddata();
			
			tableyfrac = indtableheight/ytot;
			if(tableyfrac < 1) addbutton("",menux+22+indtablewidth,corny,13,indtableheight,SLIDEAC,YSLIDEBUT,-1,-1);
				
			setxtics();
			break;
			
		case "Model":
			cornx = menux+20;
			addbutton("",cornx,80,modeldx,modeldy,CANVASBUT,CANVASBUT,-1,-1);
		
			drawmodel(res.clasel,0,0,modeldx,modeldy,"frame","indoutput");		
			break;
		}
	}
}	

function pickclass()                                      // Page that allows the user to pick a classification
{
	addbutton("Classifier",menux+tab,30,0,0,-1,TITLEBUT,-1,-1); y += 50;
			
	for(cl = 0; cl < ncla; cl++){ if(cla[cl].ncomp > 1) break;}
	if(cl == ncla){
		addbutton("Unfortunatly, there are no classifiers the data can provide information about!",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); 
	
	}
	else{
		addbutton("Which classifier does the data provide information about?",x+20,y,800,0,-1,PARAGRAPHBUT,-1,-1); y += 60;
		x = menux+tab+60;
		for(cl = 0; cl < ncla-1; cl++){
			if(cla[cl].ncomp > 1){
				w = textwidth(cla[cl].name,INPUTFONT)+50;
				if(x+w > width-30){ x = menux+tab+60; y += 50;} 
				addbutton(cla[cl].name,x,y,w,30,CHOOSECLASSAC,UPLOADBUT,cl,-1);
				x += w + 30;
			}
		}
	}
	addbutton("Cancel",width-105,height-45,90,30,CANCELBUT2,CANCELBUT2,0,-1);
	addbutton("Back",width-205,height-45,90,30,BACKBUT,BACKBUT,0,-1);	
}

function loadobsfile()                                    // Starts when loading up a observation file
{
	canover = -1;
	if(fiformat == ".csv"){ if(loadobsfile2(",") == 0) return 0;}
	else{
		if(loadobsfile2("\t") == 0) return 0;
		if(loadobsfile2(",") == 0) return 0;
		if(loadobsfile2(" ") == 0) return 0;
	}
	return 1;
}

function loadobsfile2(sep)                                // Loads a text file
{
	var lines = textFromFile.split('\n');

	j = 0; 
	do{
		trr = lines[j].trim();
		j++;
		if(trr != "" && trr.substr(0,1) != "#"){
			spli = spl(trr,sep);
			ncol = spli.length;
			for(i = 0; i < ncol; i++) colname[i] = spli[i];
			if(ncol < 2) return 1;
			break;
		}
	}while(j < lines.length);
	if(j == lines.length) return 1;

	nrow = 0;
	do{
		trr = lines[j].trim();
		if(trr != ""){
			spli = spl(trr,sep);
			if(spli.length != ncol){ alertp("Column sizes do not match!"); return 1;}
			row[nrow]=[];
			for(i = 0; i < ncol; i++) row[nrow][i] = spli[i];
			nrow++;
		}
		j++;
	}while(j < lines.length);
	
	calcrowwidth();
	
	obsloaded = 1; 
	setcolumns();
	
	return 0;
}

function spl(inp,sep)                                      // Splits up text line into columns              
{
	if(sep == " "){
		var myRegexp = /[^\s"]+|"([^"]*)"/gi;
		var myArray = [];
		do{
			var match = myRegexp.exec(inp);
			if (match != null){	myArray.push(match[1] ? match[1] : match[0]);}
		}while(match != null);	
		return myArray;
	}
	return inp.split(sep);
}

function calcrowwidth()                                  // Calculates the widths for elements in a table
{
	var r;

	cv.font = tableheadfont;
	rowwidth=[];	
	for(r = 0; r < nrow; r++){
		rowwidth[r]=[];
		for(i = 0; i < ncol; i++) rowwidth[r][i] = Math.floor(cv.measureText(row[r][i]).width);	
	}
}

function dosort(i,type,op)                               // Sorts a column in a table
{
	var sor=[];

	for(j = 0; j < nrow; j++) sor[j] = {si:row[j][i],srow:row[j],swidth:rowwidth};
	
	if(type == "num") sor.sort(function(a, b){return a.si - b.si}); 
	else sor.sort(function(a, b){ return orderstring(a.si,b.si);});

	flag = 0;
	for(j = 0; j < nrow; j++){ if(row[j] != sor[j].srow){ flag = 1; row[j] = sor[j].srow; rowwidth[j] = sor[j].srowwidth;}}
		
	if(op == "swap" && flag == 0){ // if no change sorts in the reverse order 
		for(j = 0; j < nrow; j++){
			row[j] = sor[nrow-j-1].srow;
			rowwidth[j] = sor[nrow-j-1].srowwidth;
		}
	}
}

function orderstring(x,y)                                 // Describes how ordering is performed
{	
	x = (""+x).toLowerCase(), y = (""+y).toLowerCase();
	
	var j = 0; while(j < x.length && j < y.length && x.substr(j,1) == y.substr(j,1) && isNaN(x.substr(j,1))) j++;
	if(j < x.length && j < y.length){
		xi = Number(x.substring(j));
		yi = Number(y.substring(j));
		if(!isNaN(xi) && !isNaN(yi)){
			if(xi < yi) return -1; if(xi > yi) return 1; return 0;
		}
	}
	if(x < y) return -1; if (x > y) return 1; return 0;	
}

function movecol(fr,to)                                  // Moves a column in the table
{
	if(to == fr) return;
	if(to < fr){
		temp = colname[fr];
		for(i = fr; i > to; i--) colname[i] = colname[i-1];
		colname[to] = temp;
		
		for(j = 0; j < nrow; j++){
			temp = row[j][fr];
			tempw = rowwidth[j][fr];
			for(i = fr; i > to; i--){ row[j][i] = row[j][i-1]; rowwidth[j][i] = rowwidth[j][i-1];}
			row[j][to] = temp;
			rowwidth[j][to] = tempw;
		}		
	}
	else{
		temp = colname[fr];
		for(i = fr; i < to; i++) colname[i] = colname[i+1];
		colname[to] = temp;
		
		for(j = 0; j < nrow; j++){
			temp = row[j][fr];
			tempw = rowwidth[j][fr];
			for(i = fr; i < to; i++){
				row[j][i] = row[j][i+1];
				rowwidth[j][i] = rowwidth[j][i+1];
			}
			row[j][to] = temp;
			rowwidth[j][to] = tempw;
		}		
	}
	
	setcolumns();
}

function setcolumns()                                    // Sets the spacing for the columns in the table
{ 
	xx = 0;

	for(i = 0; i < ncol; i++){
		w = 0;
		
		cv.font = tableheadfont; 
		na = colname[i];
		ww = Math.floor(cv.measureText(na).width); if(ww > w) w = ww;
		
		cv.font = tablefont; 
		for(j = 0; j < nrow; j++){
			ww = rowwidth[j][i]; if(ww > w) w = ww;
		}
		w += 25;
		if(w < 100) w = 100;
		
		colx[i] = xx; colw[i] = w;
		xx += w;	
		if(i == ncoldef-1) xx += 50;
	}
	tottablewidth = xx;
}

function drawtable()                                      // Draws a table of data
{
	jmin = Math.floor(tableyfr*nrow); jmax = jmin+rowmax; if(jmax > nrow) jmax = nrow;
	off = Math.floor(tablexfr*tottablewidth);
	
	if(ncoldef == ncoldefmax) onlycol = 0; else onlycol = 1;
	
	if(onlycol == 1 || addingdata == 2) fac = 0; else fac = 1;
	fac2 = fac;

	if(page == SIMULATEPAGE && pagesub[page] == 1){ onlycol = 0; fac2 = 1; fac = 1;}
		
	for(i = 0; i < ncol; i++){
		if(colx[i]-off < tablewidth && colx[i]+colw[i]-off > 0){
			if(i >= ncoldef) fac = 0;
			
			if(onlycol == 1){
				ddy = nrow*tabledy; if(ddy > tableheight) ddy = tableheight;
				addcanbutton("",colx[i]-off,tabledy+3,colw[i],ddy,TABLECOLBUT,TABLECOLBUT,i,-1);
				
				if(row.length == 0){
					addcanbutton("",colx[i]-off,tabledy+3,colw[i],tabledy,TABLECOLBUT,TABLECOLBUT,i,-1);
				}
			}
			addcanbutton(colname[i],colx[i]-off,0,colw[i],fac2*tabledy,TABLEHEADBUT,TABLEHEADBUT,i,-1);
			
			for(j = jmin; j < jmax; j++){		
				if(addingdata == 3 && i == ncol-1 && datatemp.variety == "state"){
					addcanbutton(row[j][i],colx[i]-off,22+(j-jmin)*tabledy,colw[i],fac*tabledy,TABLEDROPBUT,TABLEDROPBUT,i,j);
				}
				else addcanbutton(row[j][i],colx[i]-off,22+(j-jmin)*tabledy,colw[i],fac*tabledy,TABLEBUT,TABLEBUT,i,j);		
			}
			
			if(row.length == 0){
				addcanbutton("No Data",colx[i]-off,22,colw[i],fac*tabledy,-1,TABLEBUT,i,j);		
			}
		}
	}

	if(addingdata == 2) addcanbutton("",0,0,canvasdx,canvasdy,INACTIVEBUT,INACTIVEBUT,0,-1);

	addbutton("# Rows:"+(row.length),210,height-30,0,0,-1,STATBUT2,-1,-1);
	if(row.length > 0) addbutton("# Cols:"+(row[0].length),360,height-30,0,0,-1,STATBUT2,-1,-1);
}

function selectelement(r,i,bub)                           // Selects an element in the table
{
	if(selectbub != -1) buboff(1);
		
	off = Math.floor(tablexfr*tottablewidth);
	jmin = Math.floor(tableyfr*nrow);

	if(r-jmin > rowmax-1 || r-jmin < 0){
		if(r-jmin > rowmax-1) tableyfr = (r - (rowmax-5))/nrow;
		if(r-jmin < 0) tableyfr = (r - 5)/nrow;
		if(tableyfr < 0) tableyfr = 0;
		if(tableyfr > 1-tableyfrac) tableyfr = 1-tableyfrac;
		jmin = Math.floor(tableyfr*nrow);
	}
	
	selectbubx = colx[i]-off; selectbuby = 22+(r-jmin)*tabledy; selectbubdx = colw[i]; selectbubdy = tabledy;
	selectbub = bub;
	selectbubval = i;
	selectbubval2 = r;	
	buttoninit();
}

function selecthead(i,bub)                               // Select the head of a column
{
	off = Math.floor(tablexfr*tottablewidth);
	jmin = Math.floor(tableyfr*nrow);
	
	selectbubx = colx[i]-off; selectbuby = 1; selectbubdx = colw[i]; selectbubdy = tabledy;
	selectbub = bub;
	selectbubval = i;
	selectbubval2 = -1;	
}

function loadsta()                                       // Starts loading a file
{
	startloading();
	ById("fileToLoad").click();	
}

function savesta()                                       // Starts saving a file
{
	startloading();
	ById("fileToSave").click();	
}

function loadfile()                                      // Loads a file
{
	var na = ById("fileToLoad").value, fna, res;	
	fileToLoad = ById("fileToLoad").files[0];
	if(!fileToLoad) return;

	fna = fileToLoad.name;
	if(fna.length > 4) fiformat = fna.substring(fna.length-4).toLowerCase(); else fiformat = "";
		
	var fileReader = new FileReader();
	fileReader.onload = function(fileLoadedEvent){	
		newfile = 1; IDcol = -1;
		textFromFile = fileLoadedEvent.target.result;
		
		res = loadedfile();
	
		loading = 0; 
		if(res == 1){ 
			if(errormsg != "") helptype = 100; 
			else{ alertp("There was a problem loading this file");}
		}
		
		buttoninit();
	};
	fileReader.readAsText(fileToLoad, "UTF-8");
	ById("fileToLoad").value="";

	over = -1; canover = -1;
}

function loadedfile()                                      // Processes a loaded file
{
	errormsg = "";
	
	switch(fitype){
	case OBSFILE: 
		res = loadobsfile(); 
		if(res == 0){
			datashow = "table";
			dataselected = data.length;
			ncoldef = 0;
			initdatatemp();
			if(ncol < ncoldefmax){ helptype = 112; return 2;}
			addingdata = 1; 
		}
		break;
		
	case IMPORTFILE: res = importfile(textFromFile); break;
		
	case SIMINITFILE: res = loadsiminit(textFromFile); break;

	case BICIFILE: examploaded = "Load"; res = load(); if(res != 1) changepage(DESCPAGE,0,0); break;
	}
	return res;
}

function convertdates()                                  // Converts dates into continuous time
{
	var format, nmonthdays =  new Array(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
	
	for(r = 0; r < nrow; r++){  // tries to work out format
		d = row[r][1].split(/[\/\-\.]/, 3);
		if(parseInt(d[0]) > 12){ format = 0; break;}
		if(parseInt(d[1]) > 12){ format = 1; break;}
	}
	if(r == nrow){ alertp("Could not establish format"); return;}
	
	for(r = 0; r < nrow; r++){
		d = getdate(row[r][1],format);
		
		if(d){
			beginyear = new Date(d.getFullYear(),0,1);
			endyear = new Date(d.getFullYear()+1,0,1);
			row[r][1] = d.getFullYear()+(d.getTime()-beginyear.getTime())/(endyear.getTime()-beginyear.getTime());
			row[r][1] = row[r][1].toFixed(3);
		}
	}
	setcolumns();
}

function getdate(mdy,format)                             // Gets the date
{
	d = mdy.split(/[\/\-\.]/, 3);

	if (d.length != 3) return null;

	// Check if date is valid
	if(format == 0){ mon = parseInt(d[1]), day = parseInt(d[0]), year= parseInt(d[2]);}
	if(format == 1){ mon = parseInt(d[0]), day = parseInt(d[1]), year= parseInt(d[2]);}
	
	if (d[2].length == 2) year += 2000;
	if (day <= 31 && mon <= 12) return new Date(year, mon - 1, day);
	return null; 
}
	 
function datapossetup()                                     // Sets up the possibilities a data table has
{
	var cl = datatemp.cl;
	
	switch(datatype){
	case "state":		
		for(k = 0; k < cla[cl].ncomp; k++){
			if(cla[cl].comp[k].name == "I") datatemp.sensitive[k] = 1;
			else datatemp.sensitive[k] = 0;
		}
		
		tnum = 1;
		do{	
			tnam = "T"+tnum;
			for(d = 0; d < data.length; d++){
				if(data[d].type == "test" && data[d].testname == tnam) break;
			}
			if(d == data.length) break;
			tnum++;
		}while(1 == 1);

		datatemp.testname = tnam;
		
		var poslist=[];
		for(i = 0; i < nrow; i++){
			val = row[i][2];
			j = 0; jmax = poslist.length; while(j < jmax && poslist[j] != val) j++;
			if(j == jmax) poslist.push(val);
		}
		poslist.sort( function(a, b){ return orderstring(a,b);});
		
		for(j = 0; j < poslist.length; j++){
			val = poslist[j];
			datatemp.pos[j] = val;
			datatemp.posbinary[j] = [];
			datatemp.posexpression[j] = [];
			for(k = 0; k < cla[cl].ncomp; k++){
				datatemp.posbinary[j][k]=1;
				datatemp.posexpression[j][k]="";
			}
			for(k = 0; k < cla[cl].ncomp; k++) if(cla[cl].comp[k].name == datatemp.pos[j]) break;
			if(k < cla[cl].ncomp){ datatemp.posref[j] = k; datatemp.posrefmore[j] = 0;}
			else{ datatemp.posref[j] = 0; datatemp.posrefmore[j] = 1;}

			if(val == "1" || val == "+" || val == "+ve" || val == "positive") datatemp.postestres[j] = 1;
			else{
				if(val == "0" || val == "-" ||  val == "-ve" || val == "negative") datatemp.postestres[j] = 0;
				else datatemp.postestres[j] = 2;
			}
		}
		break;
	}
}

function drawpd()                                       // Draw the detection probability model
{
	nob = 0;
	
	x = 0; y = 10;
	addob(x,y,OBTEXT,"Detection:"); addob(x+110,y,OBRADIOPD2); y += 30;
	
	st = ""; if(whichind == "sub"){ for(cl = 0; cl < ncla; cl++) if(clagendata[cl] != "All") st += clagendata[cl]+",";}
	if(st != "") st = st.substring(0,st.length-1);
			
	switch(datatemp.obspd){
	case "all":
		if(st == "") addob(x,y,OBTEXT2,"All individuals in the entire population are observed.");
		else addob(x,y,OBTEXT2,"All individuals in "+st+" are observed.");
		break;
	
	case "set":
		addob(x,y,OBTEXT2,"Not all individuals are observed during a capture.");y += 50; 
		addob(x,y,OBTEXT2,"Please select whether the detection probability is the same or different across captures:"); 
		y += 40;
		x += 20;
		addob(x+10,y,OBRADIOPD3); y += 30;
	
		switch(datatemp.pdsame){
		case "same":
			if(st == "") addob(x,y,OBTEXT2,"Individuals are captured with the following probability: ");
			else addob(x,y,OBTEXT2,"Individuals in "+st+" are captured with the following probability:");
			addob(x+404,y+2,OBPD,datatemp.pd,17,150);
			break;
			
		case "dif":
			addob(x,y,OBTEXT2,"Load 'capture PD' data to inform detection probabilities for each capture.")
			break;
		}
		break;
	}

	ytot = y;
	
	placeob();	
	
	if(addingdata == 2) addcanbutton("",0,0,canvasdx,canvasdy,INACTIVEBUT,INACTIVEBUT,1,-1);
}

function drawmeaning()                                  // Associates data to observation model
{
	var x, y, ysta;
	
	cl = datatemp.cl; clgl = cl;
	
	nob = 0;
	x = 10; y = 10;
	addob(x,y,OBTEXT,"Type:"); addob(x+120,y,OBRADIO,CANRADIOCHECK,-1); y += 30;
		
	switch(datatemp.type){
	case "simple":
		addob(x+15,y,OBTEXT2,"Associates a given data value D with a given state."); y += 25;
		y += 30;
		
		addob(x,y,OBTEXT,"Select compartment:");
		break;
		
	case "binary":
		addob(x+15,y,OBTEXT2,"Define multiple compartmental states consistent with a given obsevered data D."); y += 25;
		//addob(x+15,y,OBTEXT2,""); y += 25;
		y += 30;
		
		addob(x,y,OBTEXT,"Select compartments:");
		break;
		
	case "test":
		addob(x+15,y,OBTEXT2,"Diagnotic test data assuming the test is sensitive to a defined set of compartments."); y += 25;
		addob(x+15,y,OBTEXT2,"A sensitivity Se and specificity Sp account for inaccuracies in the test."); y += 25;
		y += 20;
		
		addob(x,y,OBTEXT,"Test sensitive to:");
		
		x = 160; y -= 6;
		for(k = 0; k < cla[cl].ncomp; k++){
			if(datatemp.sensitive[k] == 1) st = cla[cl].comp[k].name+" ✓";
			else st = cla[cl].comp[k].name+" ✖"; 
			w = textwidth(st,"20px georgia")
		
			if(x+w > tablewidth-45){ x = 160; y += 50;}
			addob(x,y,OBCOMP3,k,st,COMPSMALLBUT3);
			x += w+50;		
		}
		y += 50;
		x = 10;
		addob(x,y,OBTEXT,"Test name:"); addob(105,y-3,OBNAME); y += 15;
		break;
		
	case "expression":
		addob(x+15,y,OBTEXT2,"Set abitrary expressions for the observation model. Note,"); y += 25;
		addob(x+15,y,OBTEXT2,"these can involve new model parameters."); y += 25;
		y += 30;
		
		addob(x,y,OBTEXT,"Edit observation probabilities:");
		break;
	}
	
	y += 35;
	
	if(datatemp.type == "test") istest = 1; else istest = 0;
	
	na = datatemp.testname;
	
	var wmax = 0;
	for(j = 0; j < datatemp.pos.length; j++){
		w = textwidth(datatemp.pos[j],"bold 16px arial");
		if(w > wmax) wmax = w;
	}
	
	for(j = 0; j < datatemp.pos.length; j++){	
		ysta = y;
		x = wmax+105;
		
		if(istest == 1){ addob(x,y,OBRADIO2,j,-1); if(datatemp.postestres[j] != 3) y += 30; else y -= 10;}
		
		if(datatemp.type == "simple" && datatemp.posrefmore[j] == 0){
			k = datatemp.posref[j];
			st = cla[cl].comp[k].name;
			w = textwidth(st,"20px georgia")
			addob(x,y,OBCOMP6,j,k,st);
			if(cla[cl].ncomp > 1) addob(x+w+35,y+8,OBLINK,"More...",datatemp,j);
		}
		else{
			kmax = cla[cl].ncomp; if(kmax > 5 && datatemp.posrefmore[j] == 0) kmax = 5;
			for(k = 0; k < kmax; k++){
				st = "Pr(D|"+cla[cl].comp[k].name+") = ";
				
				switch(datatemp.type){
				case "test":
					if(datatemp.postestres[j] != 3){
						switch(datatemp.postestres[j]){
						case 1:
							if(datatemp.sensitive[k] == 1) st += "[Se("+na+")]";
							else  st += "1-[Sp("+na+")]";
							break;
						
						case 0:
							if(datatemp.sensitive[k] == 1) st += "1-[Se("+na+")]";
							else  st += "[Sp("+na+")]";
							break;
						
						case 2:
							st += "1";
							break;
						}	
						w = textwidth(st,"16px georgia")
						if(x+w > tablewidth-40){ x = wmax+105; y += 30;}
						addob(x,y,OBCOMP4,j,k,st,COMPSMALLBUT4);
						x += w+30;
					}
					break;
				
				case "simple":
					if(k == datatemp.posref[j]) st = cla[cl].comp[k].name+" ✓";
					else st = cla[cl].comp[k].name+" ✖";
					
					w = textwidth(st,"20px georgia")
					if(x+w > tablewidth-40){ x = wmax+105; y += 50;}
					addob(x,y,OBCOMP6,j,k,st);
					x += w+50;
					break;
					
				case "binary":
					if(datatemp.posbinary[j][k] == 1) st = cla[cl].comp[k].name+" ✓";
					else st = cla[cl].comp[k].name+" ✖";
					
					w = textwidth(st,"20px georgia")
					if(x+w > tablewidth-40){ x = wmax+105; y += 50;}
					addob(x,y,OBCOMP6,j,k,st);
					x += w+50;
					break;
					
				case "expression":
					switch(datatemp.type){
					case "binary":
						if(datatemp.posbinary[j][k] == 0) st += "0";
						else st += "1";	
						break;
			
					case "expression":
						st += datatemp.posexpression[j][k];
						break;
					}
					w = textwidth(st,"20px georgia")
					if(x+w > tablewidth-40){ x = wmax+105; y += 50;}
					addob(x,y,OBCOMP2,j,k,st);
					x += w+50;
				}
			}
			if(kmax < cla[cl].ncomp){
				w = 50;
				if(x+w > tablewidth-40){ x = wmax+105; if(datatemp.type == "test") y += 30; else y += 50;}

				addob(x,y+8,OBLINK,"More...",datatemp,j);
			}
			else{
				if(cla[cl].ncomp > 5 || (cla[cl].ncomp > 1 && datatemp.type == "simple")){
					w = textwidth("Less...","20px georgia")
					if(x+w > tablewidth-40){ x = wmax+105; y += 50;}
					addob(x,y+8,OBLINK,"Less...",datatemp,j);
				}
			}
		}
		
		addob(25,Math.floor((ysta+y)/2+5),OBTEXT,'D="'+datatemp.pos[j]+'"',0);
		if(istest == 1) ma = 3; else ma = 5;
		addob(wmax+85,ysta-ma,OBBRACKET,y-ysta+30+2*ma);
		
		if(istest == 1) y += 60; else y += 80;
	}
	y -= 10;
	ytot = y;		

	placeob();
	
	if(addingdata == 2) addcanbutton("",0,0,canvasdx,canvasdy,INACTIVEBUT,INACTIVEBUT,-1,-1);
}

function drawdata()                                     // Draws the data sources table
{
	var x, y;
	x = 0; y = 1;
	
	nob = 0;
	
	dx = 650;
	if(datanoteon == 1){
		if(datanote == "") datanote = "Place comment here...";
		alignparagraph(datanote,dx-10);
		addob(x,y,OBSPEECHOLD,dx,nlines*20+6);
		y += nlines*20+20;
	}
	else addob(x,y,OBSPEECHOLD,dx,0);
	
	y += 5;
	addob(x,y,OBDATAHEAD); y += 40;
	for(d = 0; d < data.length; d++){ addob(x,y,OBDATA,d); y += 35;}

	if(data.length > 0){
		y += 10;
		addob(x,y,OBINDUOFLAG);
		y += 20;
	}
	
	ytot = y;
	
	placeob();	
}

function indzoom(x,fac)                                 // Zooms in on a particular individual
{
	if(xaxisauto == 1){ xaxisfixmin = axxmin; xaxisfixmax = axxmax; xaxisauto = 0;}
				
	xaxisfixmin = parseFloat(xaxisfixmin);
	xaxisfixmax = parseFloat(xaxisfixmax);
	d = xaxisfixmax - xaxisfixmin;
	mid = xaxisfixmin + d*x/(tlinexmax-tlinexmin);
	xaxisfixmin = mid - d*fac/2;
	xaxisfixmax = mid + d*fac/2;
	indplotstclear();
}

function selectdataev(isel,nind,clsel,e)                // Select a particualr data event
{
	var dy;
	changepage(INFERENCEPAGE,0,1);
	
	indshow[clsel] = 1;
	dy = 30; for(cl = 0; cl < ncla; cl++){ if(indshow[cl] == 1) dy += 25;}
	tableyfr = isel/nind;
	buttoninit();
	for(i = 0; i < ncanbut; i++){
		if(canbuttype[i] == EVBUT && canbutval[i] == isel && canbutval2[i] == clsel && canbuttext[i] == e){ selbut(i); return;}
	}
}

function drawinddata()                                  // Draws individual data
{
	var x, y, dy, imin, imax;
	
	dy = 30; for(cl = 0; cl < ncla; cl++){ if(indshow[cl] == 1) dy += 25;}

	if(xaxisauto == 1){
		if(page == INFERENCEPAGE){
			if(pagesub[INFERENCEPAGE] == 0){
				tlinetmin = datares.tmin; tlinetmax = datares.tmax;
				if(tlinetmin == tlinetmax){ tlinetmin -= 1; tlinetmax += 1;}
			}
			else{ tlinetmin = infres.tmin; tlinetmax = infres.tmax;}
		}
		else{ tlinetmin = simres.tmin; tlinetmax = simres.tmax;}
		
		axxmin = tlinetmin-0.000001; axxmax = tlinetmax+0.0000001;	
	}
	else{ axxmin = xaxisfixmin; axxmax = xaxisfixmax;}
	
	nind = getnind();
	nob = 0;
	
	ytot = dy*nind;
	ysh = Math.floor(tableyfr*ytot);
	imin = Math.floor(ysh/dy);
	imax = Math.floor((ysh+indtableheight)/dy)+1; if(imax > nind) imax = nind;
	
	x = 0; y = imin*dy;
	for(i = imin; i < imax; i++){ addob(x,y,OBIND,i); y += dy;}

	if(imax == nind && uoflag == 1 && res.sampfilt == "All"){ addob(x,y,OBINDUO); y += dy; ytot += dy;}
	
	placeob();			
	
	addcanbutton("",0,indtableheight,tablewidth,indtablemar,-1,WHITERECTBUT,-1,-1);
	addcanbutton("",tlinexmin,indtableheight+20,tlinexmax-tlinexmin,30,-1,ARROWBUT,-1,-1);
	addcanbutton("X ticks",tlinexmin,indtableheight+20,tlinexmax-tlinexmin,30,XTICKBUT,XTICKBUT,-1,-1);
	addcanbutton("Time",tlinexmin,indtableheight+20+30,tlinexmax-tlinexmin,30,-1,XLABELBUT,-1,-1);

	xx = indtablewidth-30; yy = indtableheight+60; dx = 22; dy = 26; dd = 3;
	
	if(page == INFERENCEPAGE && pagesub[page] == 2){
		addcanbutton("Reload",8,indtableheight+20+22,20,20,RELOADBUT2,RELOADBUT2,-1,-1); 
	}
}

function getnind()                                      // Gets the number of individuals
{
	var nind;
	uoflag = 0;
	switch(page){
	case INFERENCEPAGE: 
		if(pagesub[page] == 0) nind = inddata.nindtotal;
		else{
			if(res.sampfilt == "All") nind = infres.inddata.nindtotal;
			else{ 
				var r = getsampfilt();
				if(r.ch >= infres.ch.length) alertp("chain out of range");
				if(r.s >= infres.ch[r.ch].sampev.length) alertp("sample out of range");
				nind = infres.ch[r.ch].sampev[r.s].nind;
			}
			if(infres.nindmax > nind) uoflag = 1;
		}
		break;
	case SIMULATEPAGE: nind = simres.nindmax; break;
	default: alertp("Error code EC30"); break;
	}
	return nind;
}

function getsampfilt()                                   // Converts the sample filter into chain and sample
{
	var ch, s;
	
	if((""+res.sampfilt).substr(0,1) == "R"){
		var sp = res.sampfilt.split(":");
		ch = parseInt(sp[0].substring(1))-1;
		s = parseInt(sp[1])-1;
	}
	else{
		ch = parseInt(res.runpopsel);
		s = parseInt(res.sampfilt)-1;
	}

	return {ch:ch, s:s};
}

function dosearch()                                     // Does a search in a table
{
	var frag=[], fragsum=[], st;
	
	col = selectbubval;

	frag = searchterm.split("*");
	searchres=[];
	if(frag.length == 1){  // simple search
		for(r = 0; r < nrow; r++){
			st = (""+row[r][col]).trim(); 
			if(st.toLowerCase() == searchterm.toLowerCase()) searchres.push(r);
		}
	}
	else{
		sum = 0; for(j = 0; j < frag.length; j++){ fragsum[j] = sum; sum += frag[j].length;}
		for(j = 0; j < frag.length; j++) fragsum[j] -= sum;
		
		for(r = 0; r < nrow; r++){
			st = (""+row[r][col]).trim();; 
	
			len = st.length;
			if(len >= sum){
				k = 0;
				for(j = 0; j < frag.length; j++){
					len2 = frag[j].length;
					if(j == 0){ k = 0; kmax = 1;}
					else{
						if(j == frag.length-1){ k = len-len2; kmax = k+1;}
						else kmax = len+fragsum[j];
					}
				
					while(k < kmax && st.substr(k,len2) != frag[j]) k++;
					if(k == kmax) break;
					k += len2;
				}
				if(j == frag.length) searchres.push(r);
			}
		}	
	}
	
	if(searchres.length == 0){
		for(i = 0; i < ncanbut; i++) if(canbuttype[i] == TABLEHEADBUT && canbutval[i] == selectbubval) break;
		if(i == ncanbut) alertp("Error code EC31");
		selbut(i);
		selectbub = SEARCHBUB; errmsg = "Sorry no results";
		buttoninit();
	}
	else{
		jmin = Math.floor(tableyfr*nrow);
		searchresnum = 0; while(searchresnum < searchres.length && searchres[searchresnum] < jmin) searchresnum++;
		if(searchresnum == searchres.length) searchresnum = 0;
		
		selectelement(searchres[searchresnum],selectbubval,SEARCHRESBUB);
	}
}

function doreplace()                                    // Replaces elements in a table
{
	var frag=[], fragsum=[], fragrep=[], st, kst=[];
	
	rowcopy = copy(row);
	
	col = selectbubval;
	
	nreplace = 0;
	frag = searchterm.split("*");

	if(frag.length == 1){  // simple replace
		for(r = 0; r < nrow; r++){
			st = (""+row[r][col]).trim(); 
			if(st.toLowerCase() == searchterm.toLowerCase()){ row[r][col] = replaceterm; nreplace++;}
		}
	}
	else{
		fragrep = replaceterm.split("*");
	
		if(fragrep.length != frag.length){ errmsg = "Must have same number of wildcards."; return;}
	
		sum = 0; for(j = 0; j < frag.length; j++){ fragsum[j] = sum; sum += frag[j].length;}
		for(j = 0; j < frag.length; j++) fragsum[j] -= sum;
		
		for(r = 0; r < nrow; r++){
			st = (""+row[r][col]).trim();
	
			len = st.length;
			if(len >= sum){
				k = 0;
				for(j = 0; j < frag.length; j++){
					len2 = frag[j].length;
					if(j == 0){ k = 0; kmax = 1;}
					else{
						if(j == frag.length-1){ k = len-len2; kmax = k+1;}
						else kmax = len+fragsum[j];
					}
				
					while(k < kmax && st.substr(k,len2) != frag[j]) k++;
					if(k == kmax) break;
					kst.push(k);
					k += len2;
				}
				if(j == frag.length){
					for(j = frag.length-1; j >= 0; j--){
						st = st.substr(0,kst[j])+fragrep[j]+st.substring(kst[j]+frag[j].length);
					}
					row[r][col] =  st; nreplace++;
				}
			}
		}	
	}
	
	if(nreplace == 0) errmsg = "Could not find";
	else{
		calcrowwidth(); setcolumns();
		for(i = 0; i < ncanbut; i++) if(canbuttype[i] == TABLEHEADBUT && canbutval[i] == selectbubval) break;
		if(i == ncanbut) alertp("Error code EC");
		selbut(i);
		selectbub = REPLACEDONEBUB;
	}
	buttoninit();
}

function joincols(col1,col2)                              // Joins two columns together
{
	if(col1 == col2){ alertp("Cannot select the same column!"); return;}
	for(r = 0; r < nrow; r++){ row[r][col1] = row[r][col1]+"+"+row[r][col2];}
	colname[col1] = colname[col1]+"+"+colname[col2];
	
	calcrowwidth(); setcolumns();
	selectbubst = -1;
	buttoninit();
}

function dodelrows()                                      // Deletes rows on a table
{
	var frag=[], fragsum=[], st, kst=[];
	
	searchterm = ById("inp").value;
	
	rowcopy = copy(row);
	
	col = selectbubval;
	
	ndel = 0;
	frag = searchterm.split("*");

	if(frag.length == 1){  // simple replace
		for(r = nrow-1; r >= 0; r--){
			st = (""+row[r][col]).trim(); 
			if(st.toLowerCase() == searchterm.toLowerCase()){ row.splice(r,1); rowwidth.splice(r,1); nrow--; ndel++;}
		}
	}
	else{
		sum = 0; for(j = 0; j < frag.length; j++){ fragsum[j] = sum; sum += frag[j].length;}
		for(j = 0; j < frag.length; j++) fragsum[j] -= sum;
		
		for(r = nrow-1; r >= 0; r--){
			st = (""+row[r][col]).trim();
	
			len = st.length;
			if(len >= sum){
				k = 0;
				for(j = 0; j < frag.length; j++){
					len2 = frag[j].length;
					if(j == 0){ k = 0; kmax = 1;}
					else{
						if(j == frag.length-1){ k = len-len2; kmax = k+1;}
						else kmax = len+fragsum[j];
					}
				
					while(k < kmax && st.substr(k,len2) != frag[j]) k++;
					if(k == kmax) break;
					kst.push(k);
					k += len2;
				}
				if(j == frag.length){ row.splice(r,1); rowwidth.splice(r,1); nrow--; ndel++;}
			}
		}	
	}
	
	if(ndel == 0) errmsg = "Could not find";
	else{
		calcrowwidth(); setcolumns();
		for(i = 0; i < ncanbut; i++) if(canbuttype[i] == TABLEHEADBUT && canbutval[i] == selectbubval) break;
		if(i == ncanbut) alertp("Error code EC33");
		selbut(i);
		selectbub = DELROWSDONEBUB;
	}
	buttoninit();
}

function selecttablecol(val)                              // Selects a table column
{
	tablehist.push({colname:copy(colname), ncol:ncol,  ncoldef:ncoldef, row:copy(row), rowwidth:copy(rowwidth)});
	
	movecol(val,ncoldef);
	if(ncoldef == 0){
		switch(datatype){
		case "cap": case "cappd": colname[0] = "Capture name"; break;
		case "pop": colname[0] = "Time"; break;
		case "der": colname[0] = "Time"; break;
		default: colname[0] = "ID"; break;
		}
	}
	if(ncoldef == 1){
		switch(datatype){
		case "capid": colname[1] = "Capture name"; break;
		case "cappd": colname[1] = "Detection prob."; break;
		case "pop": colname[1] = "Population"; break;
		case "der": colname[1] = "Value"; break;
		case "cap": colname[1] = "Compartments"; break;
		default: colname[1] = "Time"; break;
		}
	}
	if(ncoldef == 2){
		switch(datatype){
		case "cap": colname[2] = "Time"; break;
		case "pop": case "der": colname[2] = "Standard deviation"; break;
		}
	}
	
	ncoldef++; setcolumns(); tablexfr = 0;	
	over = -1; canover = -1;
}

function adjustdata(fr,to,type)                          // Changes data to reflect change in compartment or classification
{
	var d, r, dat, cl, i, j;
	
	if(type == "comp"){
		for(i = 0; i < simindinit.length; i++){  // data table which defined initial state for simulation
			for(cl = 0; cl < ncla-1; cl++){
				if(simindinit[i].state[cl] == fr) simindinit[i].state[cl] = to;
			}
		}
	}
	 
	for(d = 0; d < data.length; d++){  // adjusts data
		dat = data[d];
		if(dat.name == fr) dat.name = to;
		
		if(type == "class"){
			for(j = 0; j < dat.ncol; j++){ if(dat.colname[j] == fr) dat.colname[j] = to;}
		}
		
		switch(dat.variety){
		case "state":
			if(type == "comp"){
				for(r = 0; r < dat.nrow; r++){ if(dat.val[r] == fr) dat.val[r] = to;}
				for(j = 0; j < dat.pos.length; j++){ if(dat.pos[j] == fr) dat.pos[j] = to;}
			}
			break;
			
		case "cap":
			if(type == "comp"){
				for(r = 0; r < dat.nrow; r++){
					var sp = dat.comps[r].split(',');
					var ans ="";
					for(j = 0; j < sp.length; j++){
						if(sp[j] == fr) sp[j] = to;
						if(j > 0) ans += ","; ans += sp[j];
					}
					dat.comps[r] = ans;
				}	
			}
			break;
			
		case "der":
			for(r = 0; r < dat.nrow; r++){
				for(j = 0; j < dat.depval[r].length; j++){
					if(dat.depval[r][j] == fr) dat.depval[r][j] = to;
				}
			}
			break;
			
		case "move":
			if(dat.transi == fr) dat.transi = to;
			if(dat.transf == fr) dat.transf = to;
			break;
			
		case "trans":
			if(dat.transi == fr) dat.transi = to;
			if(dat.transf == fr) dat.transf = to;
			for(cl = 0; cl < ncla; cl++){ if(dat.filt[cl] == fr) dat.filt[cl] = to;}
			break;
		
		case "pop":
			for(cl = 0; cl < ncla; cl++){ if(dat.popcl[cl] == fr) dat.popcl[cl] = to;}
			break;
		}
	}
}

// Storing data function:

function initdatatemp()                                 // Initialises data object
{
	switch(datatype){
	case "state":
		datatemp = {name:"", nrow, ncol, colname:[], variety:"state", id:[], t:[], val:[], cl:-1, type:"simple", tmin:0, tmax:0, pos:[], posref:[], posrefmore:[], posbinary:[], posexpression:[], sensitive:[], postestres:[], testname:""};
		ncoldefmax = 3;
		break;
		
	case "presence":
		datatemp = {name:"", nrow, ncol, colname:[], variety:"presence", id:[], t:[]};
		ncoldefmax = 2;
		break;
		
	case "trans":
		datatemp = {name:"", nrow, ncol, colname:[], variety:"trans", id:[], t:[], cl:-1, tmin:0, tmax:0, whichind:"all", obspd:"all", pd:"[pev]", npdvar:0, pdvar:[], pdvardep:[], transi:"", transf:"", transty:"", filt:[]};
		ncoldefmax = 2;
		transty = "trans"; transcl = 0;
		obspd = "all"; whichind = "all";
		break;
		
	case "move":
		datatemp = {name:"", nrow, ncol, colname:[], variety:"move", id:[], t:[], cl:-1, tmin:0, tmax:0, transi:"", transf:"", transty:""};
		ncoldefmax = 2;
		transty = "trans"; transcl = 0;
		break;
		
	case "cap":
		datatemp = {name:"", nrow, ncol, colname:[], variety:"cap", tmin:0, tmax:0, obspd:"all", pdsame:"same", pd:"[p]", npdvar:0, pdvar:[], pdvardep:[], capname:[], comps:[], t:[]};
		ncoldefmax = 3;
		break;
		
	case "capid":
		datatemp = {name:"", nrow, ncol, colname:[], variety:"capid", capname:[], id:[]};
		ncoldefmax = 2;
		break;
		
	case "cappd":
		datatemp = {name:"", nrow, ncol, colname:[], variety:"cappd", npdvar:0, pdvar:[], pdvardep:[], capname:[], pd:[]};
		ncoldefmax = 2;
		break;
		
	case "pop":
		datatemp = {name:"", nrow, ncol, colname:[], variety:"pop", t:[], val:[], sd:[], errbar:[], popcl:[], tmin:0, tmax:0};
		ncoldefmax = 3;
		break;
		
	case "der":
		datatemp = {name:"", nrow, ncol, colname:[], variety:"der", t:[], val:[], sd:[], errbar:[], der:"", dep:[], depval:[], tmin:0, tmax:0};
		ncoldefmax = 3;
		dergensel = "Select";
		break;
	}
	dataselected = data.length;
}

function adddata()                                        // Adds a new data object
{
	var cl, r, t;

	if(checkdata() == 1) return 1;

	switch(datatemp.variety){
	case "pop":
		sugname(getclagenname()+ " pop.");
		
		datatemp.popcl=[];
		for(cl = 0; cl < ncla; cl++) datatemp.popcl[cl] = clagendata[cl];
	 
		datatemp.t=[]; datatemp.val=[]; datatemp.sd=[]; datatemp.errbar=[];
		tmi = large; tma = -large;
		for(r = 0; r < nrow; r++){
			t = parseFloat(row[r][0]); if(t > tma) tma = t; if(t < tmi) tmi = t;
			datatemp.t[r] = t;
			datatemp.val[r] = parseInt(row[r][1]);
			datatemp.sd[r] = parseFloat(row[r][2]);
			datatemp.errbar[r] = geterrorbar(row[r][1],row[r][2]);
		}
		datatemp.tmin = tmi; datatemp.tmax = tma; 
		break;
		
	case "der":
		sugname(dergensel);
	
		d = 0; while(d < derive.length && derive[d].name != dergensel) d++; 
		datatemp.der = derive[d].name;
		datatemp.dep = derive[d].dep;
		
		datatemp.t=[]; datatemp.val=[]; datatemp.sd=[]; datatemp.errbar=[]; datatemp.depval=[];
		tmi = large; tma = -large;
		for(r = 0; r < nrow; r++){
			t = parseFloat(row[r][0]); if(t > tma) tma = t; if(t < tmi) tmi = t;
			datatemp.t[r] = t;
			datatemp.val[r] = parseFloat(row[r][1]);
			datatemp.sd[r] = parseFloat(row[r][2]);
			mean = parseFloat(row[r][1]); sd = parseFloat(row[r][2]); 
			datatemp.errbar[r] = { min:mean-sd, max:mean+sd};
			datatemp.depval[r]=[]; for(j = 0; j < datatemp.dep.length; j++) datatemp.depval[r][j] = row[r][3+j];
		}
		datatemp.tmin = tmi; datatemp.tmax = tma; 
		break;
		
	case "cap":
		sugname("Capture");
		
		datatemp.capname=[]; datatemp.comps=[]; datatemp.t=[];
		tmi = large; tma = -large;
		for(r = 0; r < nrow; r++){
			datatemp.capname[r] = row[r][0];
			datatemp.comps[r] = row[r][1];
			t = parseFloat(row[r][2]); if(t < tmi) tmi = t; if(t > tma) tma = t; 
			datatemp.t[r] = t;
		}
		datatemp.tmin = tmi; datatemp.tmax = tma;
		break;

	case "capid":
		sugname("Capture ID");
		
		datatemp.id=[]; datatemp.capname=[];
		for(r = 0; r < nrow; r++){
			datatemp.id[r] = row[r][0];
			datatemp.capname[r] = row[r][1];
		}
		ncoldefmax = 2;
		break;
		
	case "cappd":
		sugname("Capture PD");
		
		datatemp.capname=[]; datatemp.pd=[]; 
		for(r = 0; r < nrow; r++){
			datatemp.capname[r] = row[r][0];
			datatemp.pd[r] = row[r][1];
		}	
		ncoldefmax = 2;
		break;
		
	case "presence":
		sugname("Presence");
		
		datatemp.id=[]; datatemp.t=[];
		tma = -large; tmi = large;
		for(r = 0; r < nrow; r++){
			if(row[r][1] == "All"){ t = "All"; tmi = -large; tma = large;}
			else{ t = parseFloat(row[r][1]); if(t < tmi) tmi = t; if(t > tma) tma = t;}
			datatemp.id[r] = row[r][0];
			datatemp.t[r] = t;
		}
		datatemp.tmin = tmi; datatemp.tmax = tma; 
		break;
		
	case "state":
		na = cla[datatemp.cl].name; if(datatemp.type == "test") na += " test";
		sugname(na);
		
		datatemp.id=[]; datatemp.t=[]; datatemp.val=[];
		tma = -large; tmi = large;
		for(r = 0; r < nrow; r++){
			datatemp.id[r] = row[r][0];
			if(row[r][1] == "All"){ t = "All"; tmi = -large; tma = large;}
			else{ t = parseFloat(row[r][1]); if(t < tmi) tmi = t; if(t > tma) tma = t;} 
			
			datatemp.t[r] = t;
			datatemp.val[r] = row[r][2];
		}
		datatemp.tmin = tmi; datatemp.tmax = tma; 
		break;
	
	case "trans": case "move":
		switch(transty){
		case "+": na = "Source"; break;
		case "-": na = "Sink"; break;
		case "trans": na = transni+" → "+transnf; break;
		}
		for(cl = 0; cl < ncla; cl++) if(clagendata[cl] != "All" && cl != transcl) na += ","+clagendata[cl];
		sugname(na);
		
		datatemp.id=[]; datatemp.t=[]; 
		for(r = 0; r < nrow; r++){
			datatemp.id[r] = row[r][0];
			datatemp.t[r] = parseFloat(row[r][1]);
		}
		
		datatemp.cl = transcl;
		datatemp.transi = transni; 
		datatemp.transf = transnf;
		datatemp.transty = transty;
		
		if(datatemp.variety == "trans"){
			datatemp.tmin = tgenmin; datatemp.tmax = tgenmax;
			for(cl = 0; cl < ncla; cl++) datatemp.filt[cl] = clagendata[cl];
			datatemp.obspd = obspd; datatemp.whichind = whichind;
		}
		else{
			tma = -large; tmi = large;
			for(r = 0; r < nrow; r++){ t = parseFloat(row[r][1]); if(t < tmi) tmi = t; if(t > tma) tma = t;} 
			datatemp.tmin = tmi; datatemp.tmax = tma; 
		}
		break;
	}	
	var colnamenew=[]; for(i = 0; i < ncoldefmax; i++) colnamenew[i] = colname[i];
	datatemp.nrow = nrow; datatemp.ncol = ncoldefmax; datatemp.colname = colnamenew;
	
	data[dataselected] = datatemp;

	addingdata = 0;
	return 0;
}

function checkdata(fl)                                    // Checks data
{
	var r, c, cl, v, i, p, j, rownew=[], rowwidthnew=[], id, fl, d, nrownew, ty, sortt;
	
	for(r = 0; r < nrow; r++){                  // Checks for any empty elements
		for(i = 0; i < ncol; i++){
			if(row[r][i] == "" && row[r][i].length == 0){ errormsg = "An element is empty"; selectelement(r,i,EMPTYBUB); return 1;}
		}
	}
	
	sortt = -1;
	for(i = 0; i < ncol; i++){                  // Checks the table has the correct values
		ty = "";
		
		switch(datatemp.variety){ 
		case "cap": if(i == 1) ty = "comps"; if(i == 2) ty = "numt"; break;
		case "capid": case "cappd": break;	
		case "pop": ty = "num"; break;
		case "der": if(i <= 2) ty = "num"; break;
		case "state": case "trans": case "move": case "presence": if(i == 1) ty = "numt"; break;
		}
	
		switch(ty){
		case "num":		
			for(r = 0; r < nrow; r++){  
				if(!(row[r][i] == "All" && datatemp.variety == "state")){
					v = parseFloat(row[r][i]);
					if(isNaN(v)){ errormsg = "An element is not a number"; selectelement(r,i,NANBUB); return 1;}
				}
			}
			break;
		
		case "numt":	
			sortt = i;		
			for(r = 0; r < nrow; r++){  
				if(!(row[r][i] == "All" && datatemp.variety == "state")){
					v = parseFloat(row[r][i]);
					if(isNaN(v)){ 
						parts = st.split('/');
						if(parts.length == 3) selectelement(r,1,CONVERTDATEBUB);
						else selectelement(r,1,NANBUB);
						errormsg = "An element is not a number";
						return 1;
					}
				}
			}
			break;

		case "comps":
			for(r = 0; r < nrow; r++){  
				var sp = row[r][i].split(',');
				for(j = 0; j < sp.length; j++){
					for(cl = 0; cl < ncla; cl++){
						for(c = 0; c < cla[cl].ncomp; c++){
							if(cla[cl].comp[c].name == sp[j]) break;
						}
						if(c < cla[cl].ncomp) break;
					}
					if(cl == ncla){ errormsg = "The compartments do not match"+row[r][i]; selectelement(r,i,PROBBUB); return 1;}
				}
			}
			break;
		}
	}

	switch(datatemp.variety){             // Checks the filter clagendata
	case "pop": case "trans":
		if(clagendata.length != ncla){ alertp("There is a problem with the filter"); return 1;}
		for(cl = 0; cl < ncla; cl++){
			if(clagendata[cl] != "All"){
				for(c = 0; c < cla[cl].ncomp; c++){
					if(cla[cl].comp[c].name == clagendata[cl]) break;
				}
				if(c == cla[cl].ncomp){ alertp("There is a problem with the filter"); return 1;}
			}
		}
		break;
	}
	
	switch(datatemp.variety){
	case "state":
		cl = datatemp.cl;
		switch(datatemp.type){
		case "binary":
			for(p = 0; p < datatemp.pos.length; p++){
				for(j = 0; j < cla[cl].ncomp; j++) if(datatemp.posbinary[p][j] != 0) break;
				if(j == cla[cl].ncomp){
					alertp('D="'+datatemp.pos[p]+'" must have at least one non-zero observation probability.');
					return 1;
				}
			}
			break;
		}
		break;
	
	case "cap":
		if(datatemp.pdsame == "same"){
			nvarlist = 0; varlist=[]; vardeplist=[];   
			checkeqn(datatemp.pd);
			if(errmsg != ""){ alertp("Problem with probability of detection"); return 1;}
			setpdvar();
		}
		break;
	
	case "cappd":
		nvarlist = 0; varlist=[]; vardeplist=[];   
		for(r = 0; r < nrow; r++){
			checkeqn(row[r][1]);
			if(errmsg != ""){ errormsg = "Expression is not valid"; selectelement(r,1,PDBUT); return 1;}
		}
		setpdvar();	
		break;
		
	case "der":
		if(fl != 1){
			d = 0; while(d < derive.length && derive[d].name != dergensel) d++; 
			if(d == derive.length){ alertp("The derived quantity must be set."); return 1;}
		
			for(j = 0; j < derive[d].dep.length; j++){
				cl = 0; while(cl < ncla && cla[cl].name != derive[d].dep[j]) cl++;
				if(cl == ncla){ alertp("The derived dependancy cannot be found."); return 1;}
				
				for(r = 0; r < nrow; r++){
					for(c = 0; c < cla[cl].ncomp; c++) if(cla[cl].comp[c].name == row[r][3+j]) break;
					if(c == cla[cl].ncomp){ errormsg = "Compartment is not valid"; selectelement(r,2+j,PROBBUB); return 1;}
				}
			}
		}
		break;
		
	case "trans": 
		if(fl != 1){
			for(r = 0; r < nrow; r++){
				if(row[r][1] < tgenmin || row[r][1] > tgenmax){ alertp("Data time is outside of the range!"); return 1;}
			}
		}
		
		if(obspd == "set"){
			nvarlist = 0; varlist=[]; vardeplist=[];   
			checkeqn(datatemp.pd);
			if(errmsg != ""){ alertp("Problem with probability of detection"); return 1;}
			setpdvar();
		}
		break;
	}
	
	// removes repeated records
	
	if(sortt >= 0) dosort(sortt,"num");  
	if(datatemp.variety != "pop" && datatemp.variety != "der") dosort(0,"alph"); else dosort(0,"num"); 
	
	nrownew = 0;
	for(r = 0; r < nrow; r++){
		id = row[r][0]; 
		fl = 0;
		if(nrownew == 0) fl = 1;
		else{ for(i = 0; i < ncoldef; i++) if(row[r][i] != rownew[nrownew-1][i]) fl = 1;}
		
		if(fl == 1){
			rownew[nrownew] = row[r]; rowwidthnew[nrownew] = rowwidth[r]; nrownew++;
		}
		else{
			switch(datatemp.variety){
				case "state":
					if(row[r][0] == rownew[nrownew-1][0] && row[r][1] == rownew[nrownew-1][1] && row[r][2] != rownew[nrownew-1][2]){ 
						errormsg = "Individual '"+row[r][0]+"' has multiple states associated with it";
						selectelement(r,0,MULTISTATEBUB); return 1;
					}
					break;
					
				case "trans": case "move": break;
			}
		}
	}
	
	row = rownew; rowwidth = rowwidthnew; nrow = nrownew;
	return 0;
}

function reloaddata()                                    // Reloads data from a data object
{
	var r;
	
	nrow = datatemp.nrow; ncol = datatemp.ncol; colname = datatemp.colname;
	
	row = []; for(r = 0; r < nrow; r++) row[r]=[];
	switch(datatemp.variety){
	case "cap":
		for(r = 0; r < nrow; r++){ 
			row[r][0] = datatemp.capname[r]; 
			row[r][1] = datatemp.comps[r]; 
			row[r][2] = datatemp.t[r];
		}
		break;
		
	case "capid":
		for(r = 0; r < nrow; r++){ 
			row[r]=[];
			row[r][0] = datatemp.id[r]; 
			row[r][1] = datatemp.capname[r];
		}
		break;
		
	case "cappd":
		for(r = 0; r < nrow; r++){
			row[r][0] = datatemp.capname[r]; 
			row[r][1] = datatemp.pd[r];
		}
		break;
		
	case "pop":
		for(r = 0; r < nrow; r++){
			row[r][0] = datatemp.t[r]; 
			row[r][1] = datatemp.val[r]; 
			row[r][2] = datatemp.sd[r];
		}
		clagendata=[]; for(cl = 0; cl < ncla; cl++) clagendata[cl] = datatemp.popcl[cl];
		break;
	
	case "der":
		dergensel = datatemp.der;	
		for(r = 0; r < nrow; r++){
			row[r][0] = datatemp.t[r];
			row[r][1] = datatemp.val[r];
			row[r][2] = datatemp.sd[r];
			for(j = 0; j < datatemp.dep.length; j++) row[r][3+j] = datatemp.depval[r][j];
		}
		break;

	case "state":
		for(r = 0; r < nrow; r++){
			row[r][0] = datatemp.id[r]; 
			row[r][1] = datatemp.t[r]; 
			row[r][2] = datatemp.val[r];
		}
		break;
		
	case "trans": case "move":
		for(r = 0; r < nrow; r++){ 
			row[r][0] = datatemp.id[r]; 
			row[r][1] = datatemp.t[r];
		}
		
		transcl = datatemp.cl;
		transni = datatemp.transi; 
		transnf = datatemp.transf;
		transty = datatemp.transty;
		if(datatemp.variety == "trans"){
			tgenmin = datatemp.tmin;
			tgenmax = datatemp.tmax;
			clagendata=[]; for(cl = 0; cl < ncla; cl++) clagendata[cl] = datatemp.filt[cl];
			obspd = datatemp.obspd; whichind = datatemp.whichind;
		}
		break;
		
	case "presence":
		for(r = 0; r < nrow; r++){
			row[r][0] = datatemp.id[r]; 
			row[r][1] = datatemp.t[r];
		}
		break;
	}	
}
