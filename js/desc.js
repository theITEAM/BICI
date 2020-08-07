function descbuts()                                        // Add buttons for the description page
{
	x = menux+tab; y = 30;
	addbutton("Description",x,y,0,0,-1,TITLEBUT,20,-1); y += 50;
	
	addbutton("",x+20,y,tablewidthdesc,tableheightdesc,CANVASBUT,CANVASBUT,-1,-1);
	
	drawdesc();
	
	tableyfrac = tableheightdesc/ytot;
	if(tableyfrac < 1) addbutton("",menux+tab+20+tablewidthdesc+10,y,13,tableheightdesc,SLIDEAC,YSLIDEBUT,-1,-1);
	
	addbutton("Edit",width-205,height-45,90,30,EDITDESCAC,ADDDATABUT,0,-1);
	addbutton("Next",width-105,height-45,90,30,NEXTDESCAC,NEXTBUT,0,-1);
}
	 
function drawdesc()                                         // Draws the description canvas
{
	var x, y, j, dxx, dx;
	x = 0; y = 1;
	
	nob = 0;
	dx = width-menux-80;
	
	var splw = descnote.split('\n');
	j = 0;
	while(j < splw.length){
		splw[j] = splw[j].trim();
		if(splw[j].length == 0) splw.splice(j,1);
		else j++;
	}		
	
	for(j = 0; j < splw.length; j++){
		if(y < 150) dxx = dx-120;
		else  dxx = dx;
		
		if(j == 0){
			alignparagraph(splw[j],dxx-10,"bold 18px arial");
			addob(x+5,y,OBSPEECH3,dxx,nlines*23+5,splw[j]);
			y += 5;
		}
		else{
			alignparagraph(splw[j],dxx-10);
			addob(x+5,y,OBSPEECH,dxx,nlines*23+5,splw[j]);
		}
		y += nlines*23+10;
	}
	ddy = y; if(ddy > tableheightdesc) ddy = tableheightdesc-2;
	if(ddy < 100) ddy = 100;
	addob(x,1,OBSPEECH2,dx,ddy);
	y += 10;
	
	ytot = y;
	
	placeob();
	addcanbutton("",tablewidthdesc-100,0,0,0,-1,DESCBUT,20,-1);
}
