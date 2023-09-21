"use strict";

/// Adds buttons on the description page
function add_description_buts(lay)
{
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("Description",cx,cy,{te:desc_text});
	
	let yy = cy;
	
	lay.add_button({x:lay.dx-7, y:cy+0.5, dx:5, dy:5, type:"DescPic"}); 
	
	if(model.description.edit == true){
		let nline = Math.round((lay.dy-yy-5.5)/textbox_linesi);
		
		lay.add_textbox(cx,yy+0.5,lay.dx-2*cx,nline,{type:"description", wrap:"word"});
		
		lay.add_corner_button([["Cancel","Grey","CancelEditDescription"],["Done","Grey","DoneEditDescription"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
	else{
		add_layer("DescriptionContent",lay.x+cx,lay.y+yy,lay.dx-2*cx-6,lay.dy-yy-3.5,{});
		
		lay.add_corner_button([["Edit","Grey","EditDescription"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
}


/// Adds content in the scrollable box
function add_description_content(lay)
{
	let mar = 0.;
	let bullet = 1;
	let cy = 0.;
		
	let te = model.description.te;
	if(te){
		let para = te.split("\n");
		for(let p = 0; p < para.length; p++){
			let te= para[p];
			if(te != ""){
				lay.add_button({x:mar, y:cy+0.06, dx:bullet, dy:1.3, type:"BulletPoint", col:BLACK}); 
		
				cy = lay.add_paragraph(te,lay.inner_dx-2*mar-bullet-7,mar+bullet,cy,BLACK,para_si,para_lh);
				cy += 0.5;
			}
		}
	}
}
