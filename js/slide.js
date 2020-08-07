var nslide = 10;
var pic=[], picload=[];
var GX = 20, GY = 12, TX = 35, TY = 25;
var s = 1, cv, width=700, height=300, stat = -2, loop, loopmax = 70, dir = 0;

function slideshow()
{
	myc = ById("can");
	cv = myc.getContext("2d");

	for(i = 0; i < nslide; i++) picload[i] = 0;		
	load(1); load(2);
	
	setInterval(function(){ plot();},50);
}

function load(i)
{
	if(picload[i] == 0){		
		pic[i] = new Image(); 
		pic[i].src = "slides/slide"+i+"b.png";
		pic[i].onload = function(){
			num = parseInt(this.src.substr(this.src.substr.length-8,1));
			picload[num] = 1;
		}
	}
}

function plot()
{
	switch(stat){
	case -2:
		if(picload[s] == 1 && picload[(s+1)%nslide] == 1){ stat = -1; frac = 0;}		
		break;
	
	case -1:
		x1 = (pic[s].width-width)/2; y1 = (pic[s].height - height)/2;
		cv.globalAlpha = frac;
		cv.drawImage(pic[s],x1,y1,width,height,0,0,width,height);
		cv.globalAlpha = 1;
		frac += 0.1;
		if(frac > 1) stat = 0;
		break;
		
	case 0:
		if(picload[s] == 1 && picload[(s+1)%nslide] == 1){
			stat = 1;
			x1 = (pic[s].width-width)/2; y1 = (pic[s].height - height)/2; x2 = pic[s].width-width; y2 = pic[s].height - height;
			x3 = 0; y3 = 0; x4 = (pic[(s+1)%nslide].width-width)/2; y4 = (pic[(s+1)%nslide].height - height)/2;
		
			x1 = (pic[s].width-width)/2; y1 = (pic[s].height - height)/2; x2 = x1; y2 = pic[s].height - height;
			x3 = (pic[(s+1)%nslide].width-width)/2; y3 = 0; x4 = x3; y4 = (pic[(s+1)%nslide].height - height)/2;
		
			load((s+2)%nslide);
			loop = 1;
		}
		break;
	
	case 1:
		xxa = Math.floor(x1+(x2-x1)/(loopmax/loop)); yya = Math.floor(y1+(y2-y1))/(loopmax/loop);
		xxb = Math.floor(x3+(x4-x3)/(loopmax/loop)); yyb = Math.floor(y3+(y4-y3))/(loopmax/loop);
		//xxa = (x1+(x2-x1)*loop/loopmax); yya = (y1+(y2-y1)*loop/loopmax);
		//xxb = (x3+(x4-x3)*loop/loopmax); yyb =(y3+(y4-y3)*loop/loopmax);
		
		xxx = xxa;
		xxxx = 0;
		for(gx = 0; gx < GX; gx++){
			 yyyy = 0; yyy = yya;
			for(gy = 0; gy < GY; gy++){
				if((gx+gy)%2 == 0){
					frac = 2 + (gx+gy)/32 - 6/(loopmax/loop);
				}
				else{
					frac = 2.5 + (gx+y/2)/32 - 6/(loopmax/loop);
				}
				if(frac > 1) frac = 1; if(frac < 0) frac = 0;
		
				if(frac == 1) cv.drawImage(pic[s],xxx,yyy,TX,TY,xxxx,yyyy,TX,TY);
				else{
					if(frac == 0) cv.drawImage(pic[(s+1)%nslide],xxx,yyy,TX,TY,xxxx,yyyy,TX,TY);
					else{
						cv.beginPath();
						cv.rect(gx*TX,gy*TY,TX,TY);
						cv.fillStyle = "#aaaaaa";
						cv.fill();
						
						cv.globalAlpha = frac;
						cv.drawImage(pic[s],xxx,yyy,TX,TY,xxxx,yyyy,TX,TY);
						cv.globalAlpha = 1-frac;
						cv.drawImage(pic[(s+1)%nslide],xxx,yyy,TX,TY,xxxx,yyyy,TX,TY);
						cv.globalAlpha = 1;
					}
				}
				yyy += TY; yyyy += TY;
			}
			xxx += TX; xxxx += TX;
		}
		loop++; if(loop == loopmax){ s++; if(s == nslide) s = 0; stat = 0;}

		
		break;
	}
	drawcorners(0,0,width,height,30,"#ddddff","#9999ff");
}

function drawcorners(x,y,dx,dy,r,col,col2)
{
	var th, i, nth = Math.floor(r/3);
	if(nth < 1) nth = 1;

	cv.lineWidth = 2;
	cv.beginPath();
	cv.fillStyle = col; 

	cv.moveTo(x+dx,y);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	} 
	cv.closePath();
	cv.fill();

	cv.moveTo(x+dx,y+dy);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();

	cv.moveTo(x,y+dy);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}
	cv.closePath();
	cv.fill();

	cv.moveTo(x,y);
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();
	
	cv.beginPath();
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	} 

	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}
	
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}
	
	for(i = 0; i <= nth; i++){
		th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}

	cv.closePath();
	cv.strokeStyle = col2;
	cv.stroke();
}

function ById(a){ return document.getElementById(a);}

function pr(te){ console.log(te);}
