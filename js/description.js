"use strict";
// Functions related to description of model

/// Adds buttons on the description page
function add_description_buts(lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	let yy = cy;
	
	if(model.description.edit == true){
		yy += 2;
		let nline = Math.round((lay.dy-yy-5.5)/get_font_info().lh);
		
		lay.add_textbox(cx,yy+0.5,lay.dx-2*cx,nline,{type:"description", wrap:"word"});
		
		lay.add_corner_button([["Cancel","Grey","CancelEditDescription"],["Done","Grey","DoneEditDescription"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
	else{
		lay.add_button({x:lay.dx-6, y:cy+2.5, dx:5, dy:5, type:"DescPic"}); 
	
		add_layer("DescriptionContent",lay.x+cx,lay.y+yy,lay.dx-2*cx-5,lay.dy-yy-3.5,{});
		
		lay.add_corner_button([["Edit","Grey","EditDescription"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
}


/// Replaces characters in markdown
function desc_replace(te,st,bef,aft)
{
	let len = st.length;
	let imax = te.length+1-len;
	
	let i = 0;
	while(i < imax){
		while(i < imax && te.substr(i,len) != st) i++;
		if(i < imax){
			let j = i+len;
			while(j < imax && te.substr(j,len) != st) j++;
			if(j < imax){
				let tebef = te.substr(0,i)+bef+te.substr(i+len,j-(i+len))+aft; 
				te = tebef+te.substr(j+len)
				imax = te.length+1-len;
				i = tebef.length;
			}
			else i = j;
		}
	}
	return te;
}


/// Adds content in the scrollable box
function add_description_content(lay)
{
	let mar = 0;
	let bullet = 1;
	let cy = 0.;
		
	let te = model.description.te;
	if(te){
		let para = te.split("\n");
	
		for(let p = 0; p < para.length; p++){
			let gapfrac = 1;
			
			let te = para[p].trim();
			if(te != ""){
				let si = para_si;
				let type = "";
				let col = BLACK;

				if(te.length > 0 && te.substr(0,1) == "#"){
					let j = 0; while(j < te.length && te.substr(j,1) == "#") j++;
					if(j > 0){
						switch(j){
						case 1: si *= 1.2; col = DBLUE; type = "bold"; break; 
						default: gapfrac = 0; si *= 1; col = DDBLUE; type = "bold"; break; 
						}
						te = te.substr(j).trim();
						if(type == "bold") te = "<b>"+te+"</b>";
					}
				}
				
				let x = mar;
				if(te.length > 2 && te.substr(0,2) == "- "){
					te = te.substr(2);
					x += bullet;
					lay.add_button({x:mar, y:cy+0.06, dx:bullet, dy:1.3, type:"BulletPoint", col:BLACK}); 
					gapfrac = 0.5;
				}
				
				// Converts equations
				te = desc_replace(te,"$","<e>","</e>");
				te = desc_replace(te,"**","<b>","</b>");
				te = desc_replace(te,"*","<i>","</i>");
				
				cy = lay.add_paragraph(te,lay.inner_dx-mar-x,x,cy,col,si,1.4*si);
				
				cy += 0.5*gapfrac;
			}
		}
	}
}
