function init()                                           // Intialises the software
{
	var sum;

	document.onselectstart = function() { return false; };

	if(fin == 0 && nwon == 1) require('nw.gui').Window.get().showDevTools();

	if(nwon == 1){
		var os = require('os');
		ncpu = os.cpus().length; 
	}

	initpages();
		
	logopic = new Image(); 
	logopic.src = "pics/bicilogo.png";
	logopic.onload = function(){ buttonplot();};

	descpic = new Image(); 
	descpic.src = "pics/desc.png";
	descpic.onload = function(){ buttonplot();};

	examppic[0] = new Image(); 
	examppic[0].src = "pics/SIR.png";
	examppic[0].onload = function(){ buttonplot();};
	
	examppic[1] = new Image(); 
	examppic[1].src = "pics/SEIR.png";
	examppic[1].onload = function(){ buttonplot();};
	
	examppic[2] = new Image(); 
	examppic[2].src = "pics/grid2.png";
	examppic[2].onload = function(){ buttonplot();};
	
	examppic[3] = new Image(); 
	examppic[3].src = "pics/transexp.png";
	examppic[3].onload = function(){ buttonplot();};

	examppic[4] = new Image(); 
	examppic[4].src = "pics/Simple population.png";
	examppic[4].onload = function(){ buttonplot();};
	
	examppic[5] = new Image(); 
	examppic[5].src = "pics/SIbd.png";
	examppic[5].onload = function(){ buttonplot();};
	
	examppic[6] = new Image(); 
	examppic[6].src = "pics/LV.png";
	examppic[6].onload = function(){ buttonplot();};
	
	examppic[7] = new Image(); 
	examppic[7].src = "pics/spatialSIR.png";
	examppic[7].onload = function(){ buttonplot();};
	
	examppic[8] = new Image(); 
	examppic[8].src = "pics/gridSI.png";
	examppic[8].onload = function(){ buttonplot();};
	
	examppic[9] = new Image(); 
	examppic[9].src = "pics/badger.png";
	examppic[9].onload = function(){ buttonplot();};
	
	
	speechpic = new Image(); 
	speechpic.src = "pics/speech.png";
	speechpic.onload = function(){ buttonplot();};
	
	newmodpic = new Image(); 
	newmodpic.src = "pics/newmod.png";
	newmodpic.onload = function(){ buttonplot();};
	
	manualpic = new Image(); 
	manualpic.src = "pics/manual.png";
	manualpic.onload = function(){ buttonplot();};
	
	paperpic = new Image(); 
	paperpic.src = "pics/paper.png";
	paperpic.onload = function(){ buttonplot();};
	
	warnpic = new Image(); 
	warnpic.src = "pics/warn.png";
	
	myc = ById("myCanvas");
	maincv = myc.getContext("2d");
	cv = maincv;

	mycover = ById("overlay");
	cvover = mycover.getContext("2d");

	cornx = menux + 80;
	corny = 120;
  
    graphcan = document.createElement('canvas');
    graphcv = graphcan.getContext('2d');

	resultcan = document.createElement('canvas');
    resultcv = resultcan.getContext('2d');

	a = ById("main");
	
	a.addEventListener('mousemove', function(evt) {
		var mousePos = getMousePos(myc, evt);
		mousemove(mousePos.x,mousePos.y);
	}, false);

	a.addEventListener('mousedown', function(evt) {
		var mousePos = getMousePos(myc, evt);
		mdown(mousePos.x,mousePos.y); 
	}, false);

	a.addEventListener ("mouseout", function(evt) {
		drag = 0; 
		if(addtransmode == 1) deletetransmode();
		addcompmode = 0; addsourcemode = 0;
		buttoninit();
	}, false);
	
	a.addEventListener('mouseup', function(evt) {
		var mousePos = getMousePos(myc, evt); 
		ctrlkey = evt.ctrlKey;
		if(evt.altKey) location.reload(true);
		mup(mousePos.x,mousePos.y);
	}, false);

	nticksize = 0;
	for(sh = -10; sh <= 10; sh++){
		for(i = 0; i < 3; i++){   
			f = tickpo[i];
			if(sh < 0){ for(j = 0; j < -sh; j++) f /= 10;}
			else{ for(j = 0; j < sh; j++) f *= 10;}
			ticksize[nticksize] = f; nticksize++;
		}
	}
	
	initvar();
	 
	modelsetup = 0;
	setageon(); settimeon();
	
	setsize();

	ById("bod").style.visibility = "visible";	
	
	xmlhttp = new XMLHttpRequest();    // Used for testing
	xmlhttp.onreadystatechange = function(){
		if(xmlhttp.status == 200 && xmlhttp.readyState == 4){
			textFromFile = xmlhttp.responseText;
			load();
			buttoninit();
			loading = 0;
			//changepage(MODELPAGE,0,0);
			//changepage(MODELPAGE,-3,0);
			//changepage(SIMULATEPAGE,2,0); //changepage(SIMULATEPAGE,2,1); 
			//changepage(INFERENCEPAGE,0,0);
			//changepage(INFERENCEPAGE,1,0);
			changepage(DESCPAGE,0,0);
			//changepage(INFERENCEPAGE,2,0);
			
			//modelsetup = 0; changepage(MODELPAGE,-3,0);
			buttoninit();
		}
	};
	//xmlhttp.open("GET","Examples/EX1.bici",true); xmlhttp.send(); 
	//xmlhttp.open("GET","Examples/atemp.bici",true); xmlhttp.send(); 
}

function initpages()                                       // Initialises all the pages and scrollbars
{
	var i;
	
	for(i = 0; i < 8; i++){ pagesub[i] = 0; pagesubsub[i]=[]; for(j = 0; j < 10; j++) pagesubsub[i][j] = 0;}
	for(i = 0; i < 8; i++){ pageyst[i]=[]; for(j = 0; j < 10; j++){ pageyst[i][j]=[]; for(k = 0; k < 10; k++) pageyst[i][j][k]=0;}}	
}

function setsize()                                         // Sets the size of page objects (activates when window resized)
{
	var w = window.innerWidth;
	var h = window.innerHeight;
	width = w-25; height = h-45;
	
	if(height < 538) height = 538;
	if(width < 917) width = 917;
	
	if(height < 600) pagesize = "small"; else pagesize = "big";
	
	ById("overlay").style.width = width+"px";
	ById("overlay").style.height = height+"px";
	ById("pw").style.width = width+"px";
	
	mycover.width = width;
	mycover.height = height;
	
	graphcan.width = width;
    graphcan.height = height;
 
    resultcan.width = width;
    resultcan.height = height;
	
	myc.width = width;
    myc.height = height;

	if(widthold != undefined){
		for(cl = 0; cl < ncla; cl++){
			zoomin(cl,menux+(widthold-menux)/2,heightold/2,(width-menux)/(widthold-menux),(width-widthold)/2,(height-heightold)/2);
		}
	}

	ById("overlay").style.left = ById("pw").offsetLeft+"px";
	ById("overlay").style.top = ById("pw").offsetTop+"px";

    canw = width; canh = height;

    modeldx = width-330; modeldy = height-135; modelplaydx = modeldx-modelplaybar;
	modeldytrans = height-195;
	 
	tablewidth = width-246; tableheight = height-203; 
	obswidth = width-226; obsheight = height-180; 

	tablewidthbig = width - 232; tableheightbig = height-163; 
	derwidth = width - 232; derheight = height-180; 
	
	tableheightdata = height-143; 
	tablewidthbigger = width-212; 
	
	corwidth = width-360; corheight = height-193; 	
	
	simheight = height-153;
	rowmax = Math.floor(tableheight/tabledy)-1;
  
    indtablewidth = width-345; indtableheight = height-165;
	tablewidthdesc = width-220; tableheightdesc = height-145;
    siminitwidth = width-356;
	initmodeldx = width-356;
	graphdy = height - 255;
    rightval = width-340;
	tlinexmax = indtablewidth-20;
	exampheight = height - 280;
		
	setupwidth = tablewidthbig+30; setupheight = height-133; 
	priorwidth = width-226;
	priorheight = Math.floor((height-153)/priordy)*priordy+2;

	indplotst=[];
	
	settimeon(); setageon();
	
	buttoninit();
	widthold = width; heightold = height;
	
	if(tableyfr+tableyfrac > 1){
		tableyfr = 1-tableyfrac; if(tableyfr < 0) tableyfr = 0;
		buttoninit();
	}
	
	if(tablexfr+tablexfrac > 1){
		tablexfr = 1-tablexfrac; if(tablexfr < 0) tablexfr = 0;
		buttoninit();
	}
 }

function getMousePos(canvas, evt)                         // Gets the mouse position
{
	var rect = canvas.getBoundingClientRect();
	return {
		x: evt.clientX - rect.left,
		y: evt.clientY - rect.top
	};
}

function ById(a){ return document.getElementById(a);}     // Gets an element in DOM

function mdown(xx,yy)                                     // Fires if mouse button is clicked down
{
    var d = new Date(); timedown = d.getTime();
	
	dragged = 0;
	if(timedown-timeup < 300 && timeclick < 300){ mousedblclick(); timedown = 0; return;}
	
	drag = 0; 
	if(selectbub != -1 && (over == -1 || (buttype[over] != CANVASBUT && buttype[over] != GDROPBUT && buttype[over] != GDROPSELBUT) || canover == -1)) buboff();

	if(over >= 0){
		switch(buttype[over]){
		case GDROPSELBUT:
			if(gdropslider == 1){
				gdropmyst = my; drag = 101; gdropfracst = gdropfrac;
			}
			break;
			
	    case XSLIDEBUT:
			i = over;
			x = butx[i]; dx = butdx[i];
			x1 = x+tablexfr*dx; x2 = x+dx*(tablexfr+tablexfrac);
			if(mx >= x1 && mx <= x2){ mxst = mx; myst = my; tablexfrst = tablexfr; drag = 1;}
			else{
				if(mx > x2){ tablexfr += 0.9*tablexfrac; if(tablexfr > 1-tablexfrac) tablexfr = 1-tablexfrac; buttoninit();}
				else{
					if(mx < x1){ tablexfr -= 0.9*tablexfrac; if(tablexfr < 0) tablexfr = 0; buttoninit();}
				}
			}
			break;
		
        case YSLIDEBUT:
			i = over;
			if(my >= sliy1 && my <= sliy2){ mxst = mx; myst = my; tableyfrst = tableyfr; coheight = butdy[over]; drag = 2;}
			else{
				if(my > sliy2){ tableyfr += 0.9*tableyfrac; if(tableyfr > 1-tableyfrac) tableyfr = 1-tableyfrac; buttoninit();}
				else{
					if(my < sliy1){ tableyfr -= 0.9*tableyfrac; if(tableyfr < 0) tableyfr = 0; buttoninit();}
				}
			}
			break;	
			
		case CANVASBUT:
			if(canover >= 0){
				switch(canbuttype[canover]){
				case COMPBUT:
					dragval = canbutval[canover];
					dragval2 = canbutval2[canover];
					mxst = mx; myst = my; mxlast = mx; mylast = my;
					drag = 3;
					break;
				
				case TRANSPBUT:
					dragval = canbutval[canover];
					dragval2 = canbutval2[canover];
					mxst = mx; myst = my; mxlast = mx; mylast = my;
					drag = 5;
					break;
					
				case TRANSPLBUT:
					dragval = canbutval[canover];
					dragval2 = canbutval2[canover];
					mxst = mx; myst = my;
					drag = 6;
					break;
					
				case CLASSVALBUT:
					mxst = mx; myst = my; drag = 8;
					dragcol = cla[canbutval[canover]].comp[canbutval2[canover]].col;
					dragval = canbutval[canover];
					dragtext = canbuttext[canover]; dragx = canbutx[canover]; dragy = canbuty[canover]; dragdx = canbutdx[canover]; dragdy = canbutdy[canover];
					break;
					
				case TABLECOLBUT: case TABLEHEADBUT:
					mxst = butx[over]; myst = buty[over]; 
					dragval = canbutval[canover]; drag = 9;
					break;
					
					
				case SLIDERBUT:
					mxst = mx; myst = my; kdest = kde; drag = 10;
					break;
				}
			}
			else{
				if(arrow == 2){ 
					if(xaxisauto == 1){ xaxisfixmin = axxmin; xaxisfixmax = axxmax; xaxisauto = 0;}
					mxst = mx; mxlast = mx; myst = my; mylast = my;
					drag = 11;
				}
				if(arrow == 1){ mxst = mx; myst = my; mxlast = mx; mylast = my; drag = 4;}
			}
			break;
		}
	}
	
	if(drag > 0 && selectbub >= 0) buboff(1);
}

function mup(xx,yy)                                        // Fires when mouse button is released
{
    var i, j, time;
	
	timeup = (new Date()).getTime();
	timeclick  = timeup-timedown; 
	if(drag != 0){ drag = 0; if(timeclick > 200) buttoninit();}
	if(timeclick < 500){ mouseclick(xx,yy);}
}

function pr(te){ console.log(te);}                         // Prints to the console
