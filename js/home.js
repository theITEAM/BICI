function homebuts()                                       // The home page
{
	x = menux+tab; y = 30;

	addbutton("About",x,y,0,0,-1,TITLEBUT,-1,-1); y += 35;

	dx = (width-menux-300)/2;
	addbutton("BICI stands for 'Bayesian individual-based compartmental inference'. To start BICI either begin a new analysis (right), load an existing analysis (top right) or try one of the examples below. Click on the [?] icons for help.",x+20,y+5,dx+100,0,-1,PARAGRAPHBUT,1,-1);
	
	addbutton("",x+dx+120,y+20,80,80,NEWMODBUT,NEWMODBUT,-1,-1);
	
	addbutton("A description of how to use this software is provided in the attached manual. Further details are given in a paper hosted on bioRxiv.",x+dx+200,y+5,dx-20+5,0,-1,PARAGRAPHBUT,1,-1);

	addbutton("",x+width-260,y,60,60,PDFBUT,PDFBUT,-1,-1);
	addbutton("",x+width-260,y+65,60,60,PDFBUT2,PDFBUT2,-1,-1);

	y = 200;	
	addbutton("Examples",x,y,0,0,-1,TITLEBUT,1,-1); y += 50;

	cornx = x+20, corny = y;
	exampwidth = width-230;
	addbutton("",cornx,corny,exampwidth,exampheight,CANVASBUT,CANVASBUT,-1,-1);

	drawexample();

	tableyfrac = exampheight/ytot;
	if(tableyfrac < 1) addbutton("",menux+tab+20+exampwidth+10,y,13,exampheight,SLIDEAC,YSLIDEBUT,-1,-1);
}

function drawexample()                                     // Draws the example on the openning page
{
	var x, y, gap = 15, dyy, num, dy = 23;
	
	x = 0; y = 6; nob = 0;
	
	xx = x + Math.floor(exampwidth*0.45);
	xxx = Math.floor(exampwidth*0.25)-180; if(xxx < 0) xxx = 0;

	dyy = 44; num = 7; dyyy = 0; if(dyy+20 > num*dy) dyyy = Math.floor(dyy+20 - num*dy)/2;
	y += dyyy;
	yst = y;

	//addob(xx,y,OBEXAMP,"atemp","atemp"); y += dy;
	addob(xx,y,OBEXAMP,"EX 1: Complete knowledge of events","EX1"); y += dy;
	addob(xx,y,OBEXAMP,"EX 2: Known initial state and recoveries","EX2"); y += dy;
	addob(xx,y,OBEXAMP,"EX 3: Known recoveries only","EX3"); y += dy;
	addob(xx,y,OBEXAMP,"EX 4: Periodic disease status data","EX4"); y += dy;
	addob(xx,y,OBEXAMP,"EX 5: Diagnostic test results","EX5"); y += dy;
	addob(xx,y,OBEXAMP,"EX 6: Future and past predition","EX6"); y += dy;
	addob(xx,y,OBEXAMP,"EX 7: Time classification","EX7"); y += dy;
	
	addob(xx-20,yst,OBBRACKET,y-yst);
	addob(xxx,Math.floor((y+yst)/2 - dyy/2),OBEXAMPPIC,0,"SIR model",2);
	y += dyyy+gap;
	
	dyy = 44; num = 2; dyyy = 0; if(dyy+20 > num*dy) dyyy = Math.floor(dyy+20 - num*dy)/2; y += dyyy;
	yst = y;
	addob(xx,y,OBEXAMP,"EX 8: Non-Markovian latent period","EX8"); y += dy;
	addob(xx,y,OBEXAMP,"EX 9: Disease diagnostic tests","EX9"); y += dy;
	addob(xx,y,OBEXAMP,"EX 10: Environmental and diagnostic tests","EX10"); y += dy;
	addob(xx,y,OBEXAMP,"EX 11: Estimating test Se and Sp","EX11"); y += dy;

	addob(xx-20,yst,OBBRACKET,y-yst);
	addob(xxx,Math.floor((y+yst)/2 - dyy/2),OBEXAMPPIC,1,"SEIR model",3);
	y += dyyy+gap;
	
	dyy = 110; num = 7; dyyy = 0; if(dyy+20 > num*dy) dyyy = Math.floor(dyy+20 - num*dy)/2; y += dyyy;
	yst = y;
	addob(xx,y,OBEXAMP,"EX 12: Known infection times and initial state","EX12"); y += dy;
	addob(xx,y,OBEXAMP,"EX 13: Diagnostic test results","EX13"); y += dy;
	addob(xx,y,OBEXAMP,"EX 14: Initial / final disease status","EX14"); y += dy;
	addob(xx,y,OBEXAMP,"EX 15: Incorporating vaccination status","EX15"); y += dy;
	addob(xx,y,OBEXAMP,"EX 16: Incorporating group effect","EX16"); y += dy;
	addob(xx,y,OBEXAMP,"EX 17: SNP effects on susceptiblity / infectivity","EX17"); y += dy;
	addob(xx,y,OBEXAMP,"EX 18: SIR with non-Markovian recovery","EX18"); y += dy;
	
	addob(xx-20,yst,OBBRACKET,y-yst);
	addob(xxx,Math.floor((y+yst)/2 - dyy/2),OBEXAMPPIC,3,"Disease transmission experiment",4);
	y += dyyy+gap;
	
	dyy = 60; num = 5; dyyy = 0; if(dyy+20 > num*dy) dyyy = Math.floor(dyy+20 - num*dy)/2; y += dyyy;
	yst = y;
	addob(xx,y,OBEXAMP,"EX 19: Periodic population estimates","EX19"); y += dy;
	addob(xx,y,OBEXAMP,"EX 20: Death times","EX20"); y += dy;
	addob(xx,y,OBEXAMP,"EX 21: Captures","EX21"); y += dy;
	addob(xx,y,OBEXAMP,"EX 22: Age dependent mortality","EX22"); y += dy;
	addob(xx,y,OBEXAMP,"EX 23: Captures and a fraction of death times","EX23"); y += dy;
	addob(xx-20,yst,OBBRACKET,y-yst);
	addob(xxx,Math.floor((y+yst)/2 - dyy/2),OBEXAMPPIC,4,"Population with births and deaths",5);	
	y += dyyy+gap;
	
	dyy = 148; num = 4; dyyy = 0; if(dyy+20 > num*dy) dyyy = Math.floor(dyy+20 - num*dy)/2; y += dyyy;
	yst = y;
	addob(xx,y,OBEXAMP,"EX 24: Periodic disease status measurements","EX24"); y += dy;
	addob(xx,y,OBEXAMP,"EX 25: Disease status from captures","EX25"); y += dy;
	addob(xx,y,OBEXAMP,"EX 26: Disease diagnostic tests from captures","EX26"); y += dy;
	addob(xx,y,OBEXAMP,"EX 27: Disease induced mortality","EX27"); y += dy;
	addob(xx-20,yst,OBBRACKET,y-yst);
	addob(xxx,Math.floor((y+yst)/2 - dyy/2),OBEXAMPPIC,5,"SIS model with births and deaths",6);
	y += dyyy+gap;
		
	dyy = 90; num = 2; dyyy = 0; if(dyy+20 > num*dy) dyyy = Math.floor(dyy+20 - num*dy)/2; y += dyyy;
	yst = y;
	addob(xx,y,OBEXAMP,"EX 28: Population measurements","EX28"); y += dy;
	addob(xx,y,OBEXAMP,"EX 29: Captures and death times","EX29"); y += dy;
	
	addob(xx-20,yst,OBBRACKET,y-yst);
	addob(xxx,Math.floor((y+yst)/2 - dyy/2),OBEXAMPPIC,6,"Predator-prey model",8);
	y += dyyy+gap;
	
	dyy = 127; num = 2; dyyy = 0; if(dyy+20 > num*dy) dyyy = Math.floor(dyy+20 - num*dy)/2; y += dyyy;
	yst = y;
	addob(xx,y,OBEXAMP,"EX 30: State measurements at all locations","EX30"); y += dy;
	addob(xx,y,OBEXAMP,"EX 31: Captures at some locations","EX31"); y += dy;

	addob(xx-20,yst,OBBRACKET,y-yst);
	addob(xxx,Math.floor((y+yst)/2 - dyy/2),OBEXAMPPIC,2,"Spatial diffusion model",7);
	y += dyyy+gap;
	
	dyy = 129; num = 2; dyyy = 0; if(dyy+20 > num*dy) dyyy = Math.floor(dyy+20 - num*dy)/2; y += dyyy;
	yst = y;
	addob(xx,y,OBEXAMP,"EX 32: Location and disease status","EX32"); y += dy;
	addob(xx,y,OBEXAMP,"EX 33: Captures with disease diagnostic tests","EX33"); y += dy;
	
	addob(xx-20,yst,OBBRACKET,y-yst);
	addob(xxx,Math.floor((y+yst)/2 - dyy/2),OBEXAMPPIC,8,"Spatial disease model",9);
	y += dyyy+gap;
	
	dyy = 144; num = 2; dyyy = 0; if(dyy+20 > num*dy) dyyy = Math.floor(dyy+20 - num*dy)/2; y += dyyy;
	yst = y;
	addob(xx,y,OBEXAMP,"EX 34: Power-law spatial kernel","EX34"); y += dy;
	addob(xx,y,OBEXAMP,"EX 35: Exponential spatial kernel","EX35"); y += dy;
	
	addob(xx-20,yst,OBBRACKET,y-yst);
	addob(xxx,Math.floor((y+yst)/2 - dyy/2),OBEXAMPPIC,7,"Location-based spatial disease model",13);
	y += dyyy+gap;
	
	dyy = 137; num = 1; dyyy = 0; if(dyy+20 > num*dy) dyyy = Math.floor(dyy+20 - num*dy)/2; y += dyyy;
	yst = y;
	addob(xx,y,OBEXAMP,"EX 36: A model of badger social groups","EX36"); y += dy;
	
	addob(xx-20,yst,OBBRACKET,y-yst);
	addob(xxx,Math.floor((y+yst)/2 - dyy/2),OBEXAMPPIC,9,"Spatial disease model with births and deaths",14);
	y += dyyy+gap;
	
	ytot = y;
	
	placeob();		
}
