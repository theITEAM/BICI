"use strict";


/// Changes colour format
function hex_to_RGB(hex) 
{                 
	let result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
	return result ? {
		r: parseInt(result[1], 16),
		g: parseInt(result[2], 16),
		b: parseInt(result[3], 16)
	} : null;
}


/// Converts decimal to hexidecimal
function hex(c) 
{
	let hex = (Math.floor(c)).toString(16);
	return hex.length == 1 ? "0" + hex : hex;
}


/// Returns a darker version of a colour
function dark_colour(col)
{ 
	if(col == BLACK) return BLACK;

	const frac = 0.7;
	let bigint = parseInt(col.substring(1), 16);	
	let r = (bigint >> 16) & 255;
	let g = (bigint >> 8) & 255;
	let b = bigint & 255;
	
	return "#" + hex(frac*r) + hex(frac*g) + hex(frac*b);
}


/// Returns a lighter version of a colour
function light_colour(col) 
{ 
	const frac = 0.13;
	let bigint = parseInt(col.substring(1), 16);
	let r = 255-(255-((bigint >> 16) & 255))*frac;
	let g = 255-(255-((bigint >> 8) & 255))*frac;
	let b = 255-(255-(bigint & 255))*frac;
	return "#" + hex(r) + hex(g) + hex(b);
}


/// Print left-aligned text
function plot_text(te,x,y,font,col,width)               
{
	x = ro_int(x); y = ro_int(y);

	if(te == undefined){ te = "undefined"; error("undefined text");}
	if(te.length == 0) return;

	te = restrict_width(te,font,width);

	cv.font = font;		
	cv.textAlign = 'left';
	cv.fillStyle = col;
	cv.fillText(te,x,y);
}


/// Centers some text in an area
function center_text(te,x,y,font,col,width)               
{
	let w = text_width(te,font);
	if(w >= width) plot_text(te,x-width/2,y,font,col,width);
	else plot_text(te,x-w/2,y,font,col);
}


/// Right-aligns some text in an area
function right_text(te,x,y,font,col,width)               
{
	te = restrict_width(te,font,width);
	
	let w = text_width(te,font);
	plot_text(te,x-w,y,font,col,width);
}


// Center-aligns verticle text
function center_vert_text(te,x,y,font,col,width)   
{
	te = restrict_width(te,font,width);
	x = ro_int(x); y = ro_int(y);

	let th = Math.PI/2;
	
	cv.font = font;
	cv.save();
	cv.translate(x, y);
	cv.rotate(-th);
	cv.textAlign = 'center';
	cv.fillStyle = col;
	cv.fillText(te, 0, 0);
	cv.restore();
}


/// Vertical text
function vert_text(te,x,y,font,col,width)   
{
	te = restrict_width(te,font,width);
	
	x = ro_int(x); y = ro_int(y);

	let th = Math.PI/2;
	
	cv.font = font;
	cv.save();
	cv.translate(x, y);
	cv.rotate(-th);
	cv.textAlign = 'left';
	cv.fillStyle = col;
	cv.fillText(te, 0, 0);
	cv.restore();
}


/// Ensures that text is not beyond a given width
function restrict_width(te,font,width)
{
	if(width){		
		cv.font = font;	
	
		width *= inter.sca;

		let rule = cv.measureText("l").width;
		let max = Math.floor(width/rule);

		if(te.length > max) te = te.substr(0,max);

		if(cv.measureText(te).width > width){
			let rule2 = cv.measureText("Z").width;
			
			while(te.length > 0){
				let ww = cv.measureText(te+"...").width;
				if(ww <= width) break;
				
				let num = Math.floor(0.9*(ww-width)/rule2);
				if(num == 0) num = 1;
				te = te.substr(0,te.length-num);
			}
				
			if(te == "") error("Space is too small");
			te += "...";
		}
	}	
	return te;
}


/// Vertical text
function vert_text_right(te,x,y,font,col,width)   
{
	te = restrict_width(te,font,width);
	
	x = ro_int(x); y = ro_int(y);

	let th = Math.PI/2;
	
	cv.font = font;
	cv.save();
	cv.translate(x, y);
	cv.rotate(-th);
	cv.textAlign = 'right';
	cv.fillStyle = col;
	cv.fillText(te, 0, 0);
	cv.restore();
}

/// Find the width of text 
function text_width(te,font)                             
{ 
	cv.font = font;
	return (cv.measureText(te).width)/inter.sca;
}


/// Draw an image into a rectangular region
function draw_image(pic,x,y,dx,dy)
{
	let x2 = ro(x+dx), y2 = ro(y+dy);
	x = ro(x); y = ro(y); 
	cv.drawImage(pic,x,y,x2-x,y2-y);
}


/// Draws image directly to the screen
function draw_image_simple(pic,x,y)
{
	x = ro(x); y = ro(y);

	cv.drawImage(pic,x,y);
}


/// Draws a plus shape
function draw_plus(x, y, r, col, col2)
{
	x = ro(x); y = ro(y); r = r*inter.sca;
	
	cv.beginPath();
	cv.arc(x,y,r,0,2*Math.PI);
	cv.fillStyle = col;
	cv.fill();
	
	let frac = 0.7;
	cv.lineWidth = 2;
	cv.strokeStyle = col2;
	cv.beginPath();
	cv.moveTo(x-0.7*r,y);
	cv.lineTo(x+0.7*r,y);
	cv.moveTo(x,y-0.7*r);
	cv.lineTo(x,y+0.7*r);
	cv.stroke();
}	

/// Draws a rectangle
function draw_rect(x, y, dx, dy, col, col2, style)            
{
	let x1 = ro(x), y1 = ro(y), x2 = ro(x+dx), y2 = ro(y+dy);
	//let x1 = ro_down(x), y1 = ro_down(y), x2 = ro_up(x+dx), y2 = ro_up(y+dy);

	cv.beginPath();
	cv.lineWidth = style;
	cv.fillStyle = col;
	cv.rect(x1,y1,x2-x1,y2-y1);
	cv.fill();
	cv.strokeStyle = col2;
	cv.stroke();
}

/// Draws a rectangle
function draw_rectangle(x, y, dx, dy, col, style)            
{
	//let x1 = ro(x), y1 = ro(y), x2 = ro(x+dx), y2 = ro(y+dy);
	let x1 = ro_down(x), y1 = ro_down(y), x2 = ro_up(x+dx), y2 = ro_up(y+dy);

	cv.beginPath();
	cv.lineWidth = style;
	
	cv.rect(x1,y1,x2-x1,y2-y1);
	cv.strokeStyle = col;
	cv.stroke();
}


/// Draws a filled rectangle
function fill_rectangle(x, y, dx, dy, col)                   
{
	let x1 = ro_down(x), y1 = ro_down(y), x2 = ro_up(x+dx), y2 = ro_up(y+dy);

	cv.beginPath();
	cv.rect(x1,y1,x2-x1,y2-y1);
	cv.fillStyle = col;
	cv.fill();
}


/// Draw rectangle checkbox
function check_rectangle(x, y, dx, dy, col, style)            
{
	let x1 = ro(x), y1 = ro(y), x2 = ro(x+dx), y2 = ro(y+dy);

	cv.beginPath();
	cv.lineWidth = style;
	
	cv.rect(x1,y1,x2-x1,y2-y1);
	cv.strokeStyle = col;
	cv.stroke();
	
	let mar = 2;
	cv.beginPath();
	cv.rect(x1+mar,y1+mar,x2-x1-2*mar,y2-y1-2*mar);
	cv.fillStyle = col;
	cv.fill();
}


/// Draws reload symbol
function reload_sign(x,y,col)                              
{
	x = ro(x); y = ro(y);

	let dx = 20, dy = 20;
	cv.lineWidth = 3; 
	cv.beginPath();
	cv.arc(x+dx/2, y+dy/2+1, 8, -0.1, 1.5*Math.PI);
	cv.strokeStyle = col;
	cv.stroke(); 
	draw_arrow(x+dx/2+7,y+3,x-100,y+6,8,col);
	center_text("Reload",x+dx/2,y+dy+9,"bold 10px georgia",col);
}


/// Rounds a position to a pixel
function ro(val)
{
	return Math.round(val*inter.sca)+0.5;
}

/// Rounds a position to a pixel
function ro_int(val)
{
	return Math.round(val*inter.sca);
}

function ro_down(val)
{
	return Math.floor(val*inter.sca);
}

function ro_up(val)
{
	return Math.floor(val*inter.sca+0.9999999);
}

function nearest_pixel(val)
{
	return Math.round(val*inter.sca)/inter.sca;
}


/*
function fill_rectangle_large(x1, y1, x2, y2, col)                    // Draws a filled rectangle
{
	x1 = ro_down(x1); y1 = ro_down(y1); x2 = ro_up(x2); y2 = ro_up(y2);

	cv.beginPath();
	cv.rect(x1,y1,x2,y2);
	cv.fillStyle = col;
	cv.fill();
}
*/


/// Draws a dimming rectangle (used as background for help)
function dim_rectangle(x,y,dx,dy,col)
{
	if(col == undefined) col = BLACK;
	
	clear_rectangle(x,y,dx,dy);
	cv.globalAlpha = 0.4;
	fill_rectangle(x,y,dx,dy,col);
	cv.globalAlpha = 0.7;
	fill_rectangle(x,y,dx,dy,WHITE);
	cv.globalAlpha = 1;
}
	
	
/// Clears rectangle
function clear_rectangle(x, y, dx, dy)                    
{
	let x1 = ro(x), y1 = ro(y), x2 = ro(x+dx), y2 = ro(y+dy);
	cv.clearRect(x1,y1,x2-x1,y2-y1);
}


/// Draws a line with a colour gradient
function draw_gradient_line(x1, y1, x2, y2, f1, f2)        
{
	x1 = ro(x1); y1 = ro(y1); x2 = ro(x2); y2 = ro(y2);

	let grad= cv.createLinearGradient(x1,y1,x2,y2);
	let v = Math.floor(f1*225); grad.addColorStop(0, "rgb("+v+","+v+","+v+")");
	v = Math.floor(f2*225); grad.addColorStop(1, "rgb("+v+","+v+","+v+")");

	cv.lineWidth = 4; 
	cv.beginPath();
	cv.moveTo(x1,y1);
	cv.lineTo(x2,y2);
	cv.strokeStyle = grad;
	cv.stroke();
}


/// Sets the dash type
function set_dash_type(dash)                                    
{	
	let scale = inter.sca/10;
	let sm = 3*scale, md = 9*scale, bi = 14*scale;
	
	switch(dash%10){
	case 0: cv.setLineDash([]); break;
	case 1: cv.setLineDash([md, sm]); break;
	case 2: cv.setLineDash([sm, sm]); break;
	case 3: cv.setLineDash([md, sm, sm, sm]); break;
	case 4: cv.setLineDash([bi, md]); break;
	case 5: cv.setLineDash([bi, sm, sm, sm]); break;
	case 6: cv.setLineDash([bi, sm, md, sm]); break;
	case 7: cv.setLineDash([bi, sm, sm, sm, sm, sm]); break;
	case 8: cv.setLineDash([bi, sm, md, sm, sm, sm]); break;
	case 9: cv.setLineDash([bi, sm, md, sm, md, sm]); break;
	default: error("Option not recognised 101"); break;
	}
}


/// Draws the cross used for source
function draw_cross(x,y,r,col)
{
	x = ro(x); y = ro(y); 
	let rx = ro(r), ry = ro(0.3*r);
	
	cv.beginPath();
	cv.rect(x-rx,y-ry,rx*2,ry*2);
	cv.rect(x-ry,y-rx,ry*2,rx*2);
	cv.fillStyle = col;
	cv.fill()
}


/// Draws the minus used for sink
function draw_minus(x,y,r,col)
{
	x = ro(x); y = ro(y); 
	let rx = ro(r), ry = ro(0.3*r);
	
	cv.beginPath();
	cv.rect(x-rx,y-ry,rx*2,ry*2);
	cv.fillStyle = col;
	cv.fill()
}


/// Draws a line
function draw_line(x1, y1, x2, y2, col, style, dash)       
{// zzq
	x1 = ro(x1); y1 = ro(y1); x2 = ro(x2); y2 = ro(y2);

	cv.lineWidth = style;

	if(dash) set_dash_type(dash);

	cv.beginPath();
	cv.moveTo(x1,y1);
	cv.lineTo(x2,y2);
	cv.strokeStyle = col;
	cv.stroke();
	
	if(dash) cv.setLineDash([]);
}


/// Draws a line with no rounding 
function draw_line_noround(x1, y1, x2, y2, col, style, dash)       
{
	x1 = x1*inter.sca; y1 = y1*inter.sca; x2 = x2*inter.sca; y2 = y2*inter.sca;

	cv.lineWidth = style;

	if(dash) set_dash_type(dash);

	cv.beginPath();
	cv.moveTo(x1,y1);
	cv.lineTo(x2,y2);
	cv.strokeStyle = col;
	cv.stroke();
	
	if(dash) cv.setLineDash([]);
}


/// Draws a line with an corner arrow (used for expand symbol
function corner_arrow(x1,y1,x2,y2,col,style)
{
	x1 = ro(x1); y1 = ro(y1); x2 = ro(x2); y2 = ro(y2);
	let d = ro(0.3);

	cv.beginPath();
	cv.moveTo(x1,y1);
	cv.lineTo(x2,y2);
	cv.lineStyle = col;
	cv.strokeStyle = style;
	cv.stroke();

	let dx = d; if(x2 < x1) dx = -d;
	let dy = d; if(y2 < y1) dy = -d;

	cv.moveTo(x2,y2);
	cv.lineTo(x2-dx,y2);
	cv.lineTo(x2,y2-dy);
	cv.fillStyle = col;
	cv.fill();
}


/// Draws a polygon
function draw_polygon(polypoint,col,col2,style)             
{
	cv.lineWidth = style;

	cv.beginPath();
	cv.moveTo(ro(polypoint[0].x),ro(polypoint[0].y));
	for(let i = 1; i < polypoint.length; i++){
		cv.lineTo(ro(polypoint[i].x),ro(polypoint[i].y));
	}
	
	cv.closePath();
	cv.fillStyle = col;
	cv.fill();
	cv.strokeStyle = col2;
	cv.stroke();
}


// Draws a polygon
function draw_feature(x,y,dx,dy,poly,col,col2,style)             
{
	cv.lineWidth = style;

	for(let p = 0; p < poly.length; p++){
		let cor = poly[p].cor;
		
		cv.beginPath();
		for(let i = 0; i < cor.length; i++){
		//for(let i = 0; i < cor.length; i+=100){
			let xx = ro(x+dx*cor[i][0]);
			let yy = ro(y+dy*cor[i][1]);
			if(i == 0) cv.moveTo(xx,yy);
			else cv.lineTo(xx,yy);
		}
		cv.closePath();
		cv.fillStyle = col;
		cv.fill();
		cv.strokeStyle = col2;
		cv.stroke();
	}
}


// Draws a rectangle with rounded corners
function draw_round_rectangle(x,y,dx,dy,r,col,col2,thick)            
{
	if(thick == undefined) thick = 1;
	
	cv.lineWidth = thick;	
	
	//let x2 = ro_up(x+dx), y2 = ro_up(y+dy);
	//x = ro_down(x); y = ro_down(y); r = ro(r);
	let x2 = ro(x+dx), y2 = ro(y+dy);
	x = ro(x); y = ro(y); r = ro(r);
	dx = x2-x; dy = y2-y;

	let nth = Math.floor(r/2);
	if(nth < 1){
		cv.beginPath();
		cv.beginPath();
		cv.moveTo(x+r,y);
		cv.lineTo(x+dx-r,y);
		cv.lineTo(x+dx,y+r);
		cv.lineTo(x+dx,y+dy-r);
		cv.lineTo(x+dx-r,y+dy);
		cv.lineTo(x+r,y+dy);
		cv.lineTo(x,y+dy-r);
		cv.lineTo(x,y+r);
		cv.closePath();
		if(col != "no_fill"){
			cv.fillStyle = col; 
			cv.fill();
		}		
		cv.strokeStyle = col2;
		cv.stroke();
		return;
	}
	
	cv.beginPath(); 
	cv.moveTo(x+r,y);

	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	}

	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}

	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}

	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}

	cv.closePath();
	if(col != "no_fill"){
		cv.fillStyle = col; 
		cv.fill();
	}
	cv.strokeStyle = col2;
	cv.stroke();
}


// Draws the rounded corners around a preexisting rectangle
function draw_inverse_round_rectangle(x,y,dx,dy,r,col)            
{
	let x2 = ro(x+dx), y2 = ro(y+dy);
	x = ro(x); y = ro(y); r = ro(r);
	x--; x2++; y--; y2++; 
	dx = x2-x; dy = y2-y;

	let nth = Math.floor(r/2);

	cv.strokeStyle = col; 
	cv.fillStyle = col; 
	cv.lineWidth = 1;

	cv.beginPath(); 
	cv.moveTo(x+dx,y);
	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	}
	cv.closePath();
	cv.fill();
	cv.stroke();

	cv.beginPath(); 
	cv.moveTo(x+dx,y+dy);
	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();
	cv.stroke();
	
	cv.beginPath(); 
	cv.moveTo(x,y+dy);
	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}
	cv.closePath();
	cv.fill();
	cv.stroke();
	
	cv.beginPath(); 
	cv.moveTo(x,y);
	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();
	cv.stroke();
}


// Draws a rectangle with rounded corners at the bottom
function draw_lower_round_rectangle(x,y,dx,dy,r,col,col2,thick)            
{
	if(thick == undefined) thick = 0;
	cv.lineWidth = thick;	

	//let x2 = ro_up(x+dx), y2 = ro_up(y+dy);
	//x = ro_down(x); y = ro_down(y); r = ro(r);
	let x2 = ro(x+dx), y2 = ro(y+dy);
	x = ro(x); y = ro(y); r = ro(r);
	dx = x2-x; dy = y2-y;

	var nth = Math.floor(r/2);
	
	cv.beginPath(); 
	cv.moveTo(x,y);
	cv.lineTo(x+dx,y);
	
	if(col){
		cv.strokeStyle = col;
		cv.stroke();
	}
	
	cv.beginPath(); 
	cv.moveTo(x+dx,y);

	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}

	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}
	cv.lineTo(x,y);
	//cv.closePath();

	if(col != "no_fill"){
		cv.fillStyle = col; 
		cv.fill();
	}
	
	if(col2){
		cv.strokeStyle = col2;
		cv.stroke();
	}
	
	
	/*
	cv.beginPath(); 
	cv.moveTo(x,y);
	cv.lineTo(x+dx,y);
	
	if(col){
		cv.strokeStyle = col;
		cv.stroke();
	}
	*/
	//cv.moveTo(x,y);
	
	
}


/// Draws rounded rectangle used in dropdown buttons
function draw_round_menu_top(x,y,dx,dy,r,col,col2)          
{
	//let x2 = ro_up(x+dx), y2 = ro_up(y+dy);
	//x = ro_down(x); y = ro_down(y); r = ro(r);
	
	let x2 = ro(x+dx), y2 = ro(y+dy);
	x = ro(x); y = ro(y); r = ro(r);
	
	dx = x2-x; dy = y2-y;

	let nth = Math.floor(r/2);
	if(nth < 1) nth = 1;

	cv.lineWidth = 1;
	cv.beginPath();
 
	cv.moveTo(x+r,y);

	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	}
	cv.lineTo(x+dx,y+dy);
	cv.lineTo(x,y+dy);
	
	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}

	cv.closePath();
	cv.fillStyle = col; 
	cv.fill();
	cv.strokeStyle = col2;
	cv.stroke();
}


/// Draws rounded rectangle used in dropdown buttons going up
function draw_round_menu_bottom(x,y,dx,dy,r,col,col2)          
{
	let x2 = ro(x+dx), y2 = ro(y+dy);
	x = ro(x); y = ro(y); r = ro(r);
	
	dx = x2-x; dy = y2-y;

	let nth = Math.floor(r/2);
	if(nth < 1) nth = 1;

	cv.lineWidth = 1;
	cv.beginPath();
 
	cv.moveTo(x+dx,y);

	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}
	
	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}

	cv.lineTo(x,y);
		cv.lineTo(x+dx,y);
		
	cv.closePath();
	cv.fillStyle = col; 
	cv.fill();
	cv.strokeStyle = col2;
	cv.stroke();
}


/// Draws circle
function circle(x,y,r,col,style)                        
{
	x = ro(x); y = ro(y); r = ro(r);

	cv.lineWidth = style;
	cv.beginPath();
	cv.arc(x,y,r,0,2*Math.PI);
	cv.strokeStyle = col;
	cv.stroke();
}


/// Draws filled circle
function fill_circle(x,y,r,col,col2,style)                
{
	x = ro(x); y = ro(y); r = ro(r);
	if(r < 0.5) r = 0.5; 
	
	cv.lineWidth = style;
	cv.beginPath();
	cv.arc(x,y,r,0,2*Math.PI);
	cv.fillStyle = col;
	cv.fill();
  
	cv.strokeStyle = col2;
	cv.stroke();
}


/// Draws the corners of the main page
function draw_corners(x,y,dx,dy,r,col,col2)               
{
	let x2 = ro_up(x+dx), y2 = ro_up(y+dy); 
	x = ro_down(x); y = ro_down(y); r = ro(r);
	dx = x2-x; dy = y2-y;
	
	let nth = Math.floor(r/3);
	if(nth < 1) nth = 1;

	cv.lineWidth = 2;
	cv.beginPath();
	cv.fillStyle = col; 

	cv.moveTo(x+dx,y);
	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.sin(th),y+r-r*Math.cos(th));
	} 
	cv.closePath();
	cv.fill();

	cv.moveTo(x+dx,y+dy);
	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+dx-r+r*Math.cos(th),y+dy-r+r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();

	cv.moveTo(x,y+dy);
	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.sin(th),y+dy-r+r*Math.cos(th));
	}
	cv.closePath();
	cv.fill();

	cv.moveTo(x,y);
	for(let i = 0; i <= nth; i++){
		let th = i*Math.PI/(2*nth);
		cv.lineTo(x+r-r*Math.cos(th),y+r-r*Math.sin(th));
	}
	cv.closePath();
	cv.fill();
}


/// Converts a text string into annotated text
function text_convert_annotation(te,si,lh,dx,align,col)      
{
	let subsi = 0.7*si;
	
	let col_basic = col;
	
	/// This allows for links to be added (these have format ['Text','Action'] and must be separated from other words)
	/// Additionally information can passed through the link using ['Text','Action|Info']
	let link=[];
	for(let i = 0; i < te.length-4; i++){
		if(te.substr(i,2) == "['"){
			let j = i+2;
			while(j < te.length && te.substr(j,1) != "'") j++;
			if(j < te.length-3){
				let tex = te.substr(i+2,j-(i+2));
				
				if(te.substr(j,3) != "','") error("Link not right");
				else{
					j += 3;
					let k = j;
					while(k < te.length-1 && te.substr(k,2) != "']") k++;
					
					link.push({te:tex,ac:te.substr(j,k-j)});
					
					te = te.substr(0,i)+"\r"+te.substr(k+2);
				}
			}					
		}
	}
	
	let para = te.split("\n");
	
	let word=[];
	
	let x = 0, y = 0.7*lh;
	
	let width_first;
	
	let li = 0;
	let li_width=[];
	
	let link_num = 0;
	for(let p = 0; p < para.length; p++){
		let wo_basic = para[p].split(" ");
		let wo = [];
		
		let it_flag = false, bo_flag = false, eq_flag = false;
		let col = col_basic;
		for(let i = 0; i < wo_basic.length; i++){
			let te = wo_basic[i];
			
			if(eq_flag == true && wo.length > 0){
				let j = 0; while(j < te.length && te.substr(j,4) != "</e>") j++;
				wo[wo.length-1].te += " "+te.substr(0,j);
				te = te.substr(j);
			}
			
			let k = 0; 
			let te2 = "";
			
			let space = true;
			while(k < te.length){
				let flag = false;
				
				switch(te.substr(k,3)){
				case "<b>": case "<i>": case "<c>": case "<t>": case "<g>": case "<e>": 
				case "</b": case "</i": case "</c": case "</t": case "</g": case "</e": 
					flag = true;
					break;				
				}
				
				if(flag == false){
					te2 += te.substr(k,1);
					k++;
				}
				else{
					if(te2 != ""){
						wo.push({ te:te2, space:space, col:col, bo_flag:bo_flag, it_flag:it_flag, eq_flag:eq_flag});
						te2 = "";
						space = false;
					}
		
					switch(te.substr(k,3)){
					case "<b>": bo_flag = true; k += 3; break;
					case "<i>": it_flag = true; k += 3; break;
					case "<c>": col = BLUE; k += 3; break;
					case "<t>": col = DGREEN; k += 3; break;
					case "<g>": col = DGREY; k += 3; break;
					case "<e>": eq_flag = true; k += 3; break;
					case "</b": bo_flag = false; k += 4; break;
					case "</i": it_flag = false; k += 4; break;
					case "</c": col = col_basic; k += 4; break;
					case "</t": col = col_basic; k += 4; break;
					case "</g": col = col_basic; k += 4; break;
					case "</e": eq_flag = false; k += 4; break;
					}
				}
			}
			
			if(te2 != ""){
				wo.push({ te:te2, space:space, col:col, bo_flag:bo_flag, it_flag:it_flag, eq_flag:eq_flag});	
			}
		}	
		
		let contain_word = false;
		for(let j = 0; j < wo.length; j++){
			let wor = wo[j];
			
			let bold = ""; 
			if(wor.bo_flag == true) bold = "bold";
			if(wor.it_flag == true) bold += " italic";
		
			let fo = get_font(si,bold);
				
			if(j != 0 && x != 0 && wor.space == true){
				let te = " ";
						
				let w = text_width(te,fo);
				if(x+w < dx) x += w;
				else{ li_width[li] = x; li++; x = 0; y += lh;}
			}
			
			let te = wor.te;
			
			if(te == "\r"){
				let te3 = link[link_num].te;
				let w = text_width(te3,fo);
				if(x+w >= dx){ li_width[li] = x; li++; x = 0; y += lh;}
			
				let lin = link[link_num];
				lin.x = x; lin.y = y; lin.w = w; lin.font = fo; lin.si = si;
				x += w;		
				link_num++;	
			}
			else{
				if(wor.eq_flag == true){
					let end;
					if(te.length > 3 && (te.substr(te.length-3,3) == "(t)" || te.substr(te.length-3,3) == "(a)")){
						end = te.substr(te.length-3,3);
						te = te.substr(0,te.length-3);
					}
					
					let main = "", sub = "", sup = "";
					
					let spl = te.split("^");
					if(spl.length == 1){
						let spl2 = spl[0].split("_");
						if(spl2.length == 1) main = spl2[0];
						if(spl2.length == 2){ main = spl2[0]; sub = spl2[1];}
						if(spl2.length > 2) error("Problem with split");
					}
					
					if(spl.length == 2){
						main = spl[0];
						let spl2 = spl[1].split("_");
						sup = spl2[0];
						if(spl2.length == 2) sub = spl2[1];
						if(spl2.length > 2) error("Problem with split");
					}
					
					if(spl.length > 2) error("Problem with split");
					
					bold = ""; if(wor.bo_flag == true) bold = "bold";
					fo = get_font(si,bold,"times");
					
					let fo_sub = get_font(subsi,"italic "+bold,"times");
					let fo_sup = get_font(subsi,bold,"times");
					
					sub = remove_bracket(sub);
					sup = remove_bracket(sup);
					
					let end_bracket = false;
					
					if(sub.substr(sub.length-1,1) == ")"){
						sub = sub.substr(0,sub.length-1);
						end_bracket = true;
					}
					
					if(sup.substr(sup.length-1,1) == ")"){
						sup = sup.substr(0,sup.length-1);
						end_bracket = true;
					}
					
					let w = text_width(main,fo);
					
					let w2 = 0;
					if(sub != "") w2 = text_width(sub,fo_sub);
					if(sup != ""){ let ww = text_width(sup,fo_sup); if(ww > w2) w2 = ww;}
					
					if(x+w+w2+0.1 >= dx){ li_width[li] = x; li++; x = 0; y += lh;}
					
					word.push({te:main, x:x, y:y, font:fo, col:col});		
					if(sub != "") word.push({te:sub, x:x+w+0.1, y:y+0.2*si, font:fo_sub, col:wor.col});		
					if(sup != "") word.push({te:sup, x:x+w+0.1, y:y-0.4*si, font:fo_sup, col:wor.col});

					contain_word = true;
					x += w+w2+0.1;	

					if(end_bracket == true){
						word.push({te:")", x:x, y:y, font:fo, col:col});		
						x += text_width(")",fo);
					}
					
					if(end != undefined){
						word.push({te:end, x:x, y:y, font:fo, col:col});		
						x += text_width(end,fo);
					}
				}		
				else{		
					let w = text_width(wor.te,fo);
						
					if(x+w >= dx){ li_width[li] = x; li++; x = 0; y += lh;}
					
					word.push({te:te, x:x, y:y, font:fo, col:wor.col, li:li});		
					
					contain_word = true;
					x += w;
				}
			}
		}
		
		width_first = x;
		
		if(x != 0){ li_width[li] = x; li++; x = 0; y += lh;}
		
		if(contain_word == true && p < para.length-1){
			y += 0.3;
		}
	}
	
	if(x == 0 && y == 0.7*lh){ li++; y += lh;}
	
	//if(y == 0.7*lh && x < dx && 
	if(align == "center"){
		for(let i = 0; i < word.length; i++){
			word[i].x += (dx-li_width[word[i].li])/2;
		}
	}
	
	return { width_first:width_first,  height:y-lh+si*0.4, word:word, link:link};
}


/// Selects white or black depending on the colour of the background
function white_black(col)
{
	let bigint = parseInt(col.substring(1), 16);	
	let r = (bigint >> 16) & 255;
	let g = (bigint >> 8) & 255;
	let b = bigint & 255;
	
	if(r+g+b > 200*3) return BLACK;
	return WHITE;
}


/// Draws an arrow
function draw_arrow(x,y,x2,y2,size,col)                   
{
	let nx = x2-x, ny = y2-y;
	let r = Math.sqrt(nx*nx+ny*ny);
	//if(size > r/5) size = r/5;
	nx *= size/r; ny *= size/r; 
	let px = 0.5*ny, py = -0.5*nx;

	let polypoint=[];
	polypoint.push({x:x, y:y});
	polypoint.push({x:x+0.8*nx, y:y+0.8*ny});
	polypoint.push({x:x+nx+px, y:y+ny+py});
	draw_polygon(polypoint,col,col,NOLINE);
	
	polypoint=[];
	polypoint.push({x:x, y:y});
	polypoint.push({x:x+0.8*nx, y:y+0.8*ny});
	polypoint.push({x:x+nx-px, y:y+ny-py});
	
	draw_polygon(polypoint,col,col,NOLINE);
}


/// Plots the loading symbol
function plot_loading_symbol(x,y,r)         
{
	x = ro(x), y = ro(y), r = ro(r);
	let	rr = ro(0.2);

	cv.fillStyle = BLACK; 
	for(let i = 0; i < 12; i++){
		let th = (i-inter.loading_symbol.offset)*2*Math.PI/12.0;
		cv.globalAlpha = Math.exp(-i*0.4);
		cv.beginPath();
		cv.arc(Math.floor(x+r*Math.sin(th)),Math.floor(y+r*Math.cos(th)),rr,0,2*Math.PI);
		cv.fill();
	}
	cv.globalAlpha = 1;
}


function draw_menu_arrow_down(x,y,si,col)
{
	x = ro(x); y = ro(y); si *= inter.sca;
	
	cv.beginPath();
	cv.fillStyle = col;

	cv.moveTo(x-si,y-si/2);
	cv.lineTo(x+si,y-si/2);
	cv.lineTo(x,y+si/2);
	cv.fill();
}


/// Converts a label to individual text elements
function label_convert(te,si,wmax)
{
	let fo_big = get_font(1.4*si,"","Times");
	let fo = get_font(si,"","Times");
	let sub_fo = get_font(0.6*si,"italic","Times");
	let sup_fo = get_font(0.6*si,"","Times");
	
	let frag = [];
	
	let prob_flag = false;
	 
	let gap = 0.05;

	let x = 0;
	let i = 0;
	while(i < te.length){
		let i_store = i;
	
		let type;
		do{
			let ch = te.substr(i,1);
			if(ch == "{"){ type = "pop"; break;}
			if(ch == "["){ type = "ie"; break;}
			if(ch == "〈"){ type = "fe"; break;}
			if(!notparam_list.includes(ch)){ type = "param"; break;}				
			i++;
		}while(i < te.length);

		if(i_store != i){
			let tex = te.substr(i_store,i-i_store).replace(/ /g,"");

			frag.push({ te:tex, col:BLACK, x:x, y:0*si, fo:fo});
			x += text_width(tex,fo)+gap;  
		}
		
		if(i < te.length){
			let pair, col;
			
			let i_store = i;
			let tex;
			
			if(type == "param"){
				i = param_end(te,i);
				if(typeof i == 'string'){ prob_flag = true; break;}
				
				tex = te.substr(i_store,i-i_store);	
				if(i < te.length && te.substr(i,1) == "$") i++;
			}
			else{
				switch(te.substr(i,1)){
				case "{": pair = "}"; col = BLUE; break;
				case "[": pair = "]"; col = DGREEN; break;
				case "〈": pair = "〉"; col = DRED; break;
				case "Σ": pair = "("; col = BLACK; break;
				default: alertp("No pair"); break;
				}	
				
				while(i < te.length && te.substr(i,1) != pair) i++;
			
				if(i < te.length) i++;
	
				if(pair == "("){
					let num = 1;
					for(let j = i; j < te.length; j++){
						if(te.substr(j,1) == "(") num++;
						if(te.substr(j,1) == ")"){
							num--;
							if(num == 0){
								te = te.substr(0,j)+te.substr(j+1);
								break;
							}
						}
					}
				}
				
				tex = te.substr(i_store+1,i-i_store-2);
			}
			
			if(type == "param" || pair == "}" || pair == "("){
				let sub = "", sup = "", time = "";
			
				if(tex.length > 3){
					let end = tex.substr(tex.length-3,3);
					if(end == "(t)" || end == "(a)"){
						time = end;
						tex = tex.substr(0,tex.length-3);
					}
				}
				
				let k = 0; while(k < tex.length && tex.substr(k,1) != "^" && tex.substr(k,1) != "_") k++;
				let main = tex.substr(0,k);
				if(pair == "}"){
					if(main == "") main = "All";
					main = "{"+main+"}";
				}
			
				if(pair == "("){
					main = "Σ"+main;
				}
				
				for(let loop = 0; loop < 2; loop++){
					if(k < tex.length){
						let k_store = k;
						k++; 
						while(k < tex.length && tex.substr(k,1) != "^" && tex.substr(k,1) != "_" ) k++;
						let st = tex.substr(k_store+1,k-(k_store+1));
						switch(tex.substr(k_store,1)){
							case "^": sup = st; break;
							case "_": sub = st; break;
							default: error("Option not recognised 103"); break;
						}
					}
				}			
				
				sup = remove_bracket(sup);
				sub = remove_bracket(sub);

				let ysh = 0;
				let fo_main = fo; 
				if(main == "Σ" && pair == "("){ fo_main = fo_big; ysh = 0.2;}

				frag.push({ te:main, col:col, x:x, y:ysh*si, fo:fo_main});
				x += text_width(main,fo_main);  
				if(sub != "" || sup != "") x += gap/2;
				
				let w_sub = 0, w_sup = 0;
				if(sub != ""){
					frag.push({ te:sub, col:col, x:x, y:0.15*si, fo:sub_fo});
					w_sub = text_width(sub,sub_fo)
				}
				
				if(sup != ""){
					frag.push({ te:sup, col:col, x:x, y:-0.3*si, fo:sup_fo});
					w_sup = text_width(sup,sup_fo)
				}
				
				if(w_sub > w_sup) x += w_sub; 
				else x += w_sup;

				if(time != ""){
					frag.push({ te:time, col:col, x:x, y:0*si, fo:fo});
					x += text_width(time,fo);
				}

				x += gap;
			}
			else{
				//if(pair == "}") tex = "{"+tex+"}";
			
				if(pair == "]"){
					tex = "["+tex+"]";
				}
				
				if(pair == "〉"){
					tex = "〈"+tex+"〉";
				}
					
				frag.push({ te:tex, col:col, x:x, y:0*si, fo:fo});
				x += text_width(tex,fo)+gap;  
			}
		}
	}
	
	if(prob_flag == true){
		frag = [];
		frag.push({ te:te, col:BLACK, x:0, y:0, fo:fo})
	}
	
	let w = 0;
	for(let i = 0; i < frag.length; i++){
		let x = frag[i].x + text_width(frag[i].te,frag[i].fo);
		if(x > w) w = x; 
	}
	
	if(w > wmax){
		let te = "...";
		let tew = text_width(te,fo);
		
		let i;
		for(i = frag.length-1; i >= 0; i--){
			let fo =  frag[i].fo;
			if(frag[i].x < wmax-tew && frag[i].fo != sub_fo && frag[i].fo != sup_fo ) break;
		}
		
		if(i < 1) frag.length = 0;
		else{
			let x = frag[i].x;
			frag.length = i;
			frag.push({te:te, col:BLACK, x:x, y:0, fo:fo});
		}
		
		w = 0;
		for(let i = 0; i < frag.length; i++){
			let x = frag[i].x + text_width(frag[i].te,frag[i].fo);
			if(x > w) w = x; 
		}
	}

	return {frag:frag, w:w};
}


/// Plots an equation
function plot_equation(te,x,y,si,w,col)
{
	plot_angle_equation(te,x,y,1,0,si,w,"left",col);
}


/// Plots equation at an angle
function plot_angle_equation(te,x,y,nx,ny,si,w,orient,col)
{
	if(col == undefined) col = BLACK;
	let label_conv = label_convert(te,si,w);

	if(w == undefined) w = LARGE;

	if(label_conv.w > w) return;
	
	let px = -ny, py = nx;
	
	for(let i = 0; i < label_conv.frag.length; i++){
		let frag = label_conv.frag[i];
		
		let xf = frag.x; if(orient != "left") xf -= label_conv.w/2;
		let yf = frag.y;
		
		let xx = x + xf*nx + yf*px;
		let yy = y + xf*ny + yf*py;
		
		xx = ro_int(xx); yy = ro_int(yy); 

		let th = Math.acos(nx);
		if(ny > 0) th *= -1;
				
		cv.font = frag.fo;
		cv.save();	
		cv.translate(xx,yy);
		cv.rotate(-th);
		cv.textAlign = 'left';
		//cv.fillStyle = frag.col;
		cv.fillStyle = col;
		cv.fillText(frag.te, 0, 0);
		cv.restore();
	}
}


/// Draws a left arrow 
function draw_left_arrow(x,y,dx,dy,col,col2,style)
{
	let x1 = ro(x), y1 = ro(y), x2 = ro(x+dx), y2 = ro(y+dy);

	let ym = ro(y+0.5*dy);
	let xm = ro(x+0.4*dx);
	let d = ro(0.3*dy);
	
	cv.lineWidth = style;
	cv.beginPath();
	cv.moveTo(x1,ym);
	cv.lineTo(xm,y1);
	cv.lineTo(xm,y1+d);
	cv.lineTo(x2,y1+d);
	cv.lineTo(x2,y2-d);
	cv.lineTo(xm,y2-d);
	cv.lineTo(xm,y2);
	cv.lineTo(x1,ym);
	cv.strokeStyle = col;
	cv.stroke();
	
	if(col2 != "no_fill"){
		cv.fillStyle = col2;
		cv.fill();
	}
}



/// Plots a parameter label
function plot_paramlabel(x,y,dx,dy,tcol,info)
{		
	let xx = x+dx/2-info.dx/2;
	let yy = y+dy/2+0.3*si_big;
	plot_text(info.name,xx,yy,info.fo_big,tcol,dx);

	let xx2 = xx+text_width(info.name,info.fo_big)+0.02;
	plot_text(info.sup,xx2,yy-0.6,info.fo_sup,tcol,dx);
	plot_text(info.sub,xx2,yy+0.1,info.fo_sub,tcol,dx);

	if(info.age_time != ""){
		plot_text(info.age_time,xx+info.age_time_pos,yy,info.fo_big,tcol,dx);
	}
}
		

/// Gets parameter info
function plot_param_info(te,si)
{
	let frag = [];
	
	let name, sub, sup, end;
	
	if(te.length >= 3){
		let pos = te.substr(te.length-3,3);
		if(pos == "(t)" || pos == "(a)"){ te = te.substr(0,te.length-3); end = pos;}
	}		
			
	let spl = te.split("_");
	if(spl.length == 2) sub = spl[1];
		
	let spl3 = spl[0].split("^");
	if(spl3.length == 2) sup = spl3[1];
	
	name = spl3[0];
	
	let font = get_font(si);
	let x = 0;
	let w = text_width(name,font);
	frag.push({te:name, x:x, y:0, font:font});
		
	if(sub || sup || end){
		x += w;
	
		w = 0;
		if(sub){
			let fo_sub = get_font(si*0.7,"italic","Times");
			let ww = text_width(sub,fo_sub);
			frag.push({te:sub, x:x, y:0.3*si, font:fo_sub});
			if(ww > w) w = ww;
		}
		
		if(sup){
			let fo_sup = get_font(si*0.7,undefined,"Times");
			let ww = text_width(sup,fo_sup);
			frag.push({te:sup, x:x, y:-0.3*si, font:fo_sup});
			if(ww > w) w = ww;
		}
		
		x += w+0.05;
		if(end){
			frag.push({te:end, x:x, y:0, font:font});
			x += text_width(end,font);
		}
	}
	
	return { w:x, frag:frag};
}


/// Plots parameter text
function plot_param_text(te,x,y,si,col,width)
{
	let info = plot_param_info(te,si);

	for(let i = 0; i < info.frag.length; i++){
		let fr = info.frag[i];
		if(fr.x < width) plot_text(fr.te,x+fr.x,y+fr.y,fr.font,col,width-fr.x);
	}
}


/// Plots parameter text
function center_param_text(te,x,y,si,col,width)
{
	let info = plot_param_info(te,si);
	let sh = info.w/2; if(sh > width/2) sh = width/2;
	let xx = x-sh;
	for(let i = 0; i < info.frag.length; i++){
		let fr = info.frag[i];
		if(fr.x < width) plot_text(fr.te,xx+fr.x,y+fr.y,fr.font,col,width-fr.x);
	}
}


/// Plots verticle parameter text
function center_vert_param_text(te,x,y,si,col,width) 
{
	let info = plot_param_info(te,si);

	let sh = info.w/2; if(sh > width/2) sh = width/2;
	let yy = y+sh;
	for(let i = 0; i < info.frag.length; i++){
		let fr = info.frag[i];
		if(fr.x < width){
			vert_text(fr.te,x+fr.y,yy-fr.x,fr.font,col,width-fr.x);
		}			
	}
}


/// Plots right-aligned parameter text
function right_param_text(te,x,y,si,col,width)
{
	let info = plot_param_info(te,si);
	
	let sh = info.w; if(sh > width) sh = width;
	let yy = y+sh;
	for(let i = 0; i < info.frag.length; i++){
		let fr = info.frag[i];
		if(fr.x < width) vert_text(fr.te,x+fr.y,yy-fr.x,fr.font,col,width-fr.x); 
	}
}
