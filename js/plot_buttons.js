"use strict";
// Functions for plotting buttons

/// Plots individual button
Layer.prototype.plot_button = function (bu,ov)
{ 	
	if(!bu){ error("No button"); return;}
	if(!bu.x && bu.x != 0) error("no button");
	
	let x = bu.x-this.x_shift;
	let y = bu.y-this.y_shift;
	
	if(inter.export_image == true){ x = bu.x; y = bu.y;}
		
	let dx = bu.dx, dy = bu.dy, te = bu.te;

	cv = this.cv;

	switch(bu.type){  // Extra fast plotting functions for results
	case "CompLatLng":
		{
			let col = bu.col; 
			if(ov) col = shift_colour(col);	
			
			if(dx*inter.sca < 1) plot_pixel(x,y,col);
			else{
				fill_circle(x+dx/2,y+dy/2,0.9*dx/2,col,bu.col_dark,NORMLINE);
			}
		}
		return;
		
		case "CompGraph2": 
		{
			let col = bu.col; if(ov) col = shift_colour(col);
			let te_col = white_black(col);
			
			if(dx*inter.sca < 1) plot_pixel(x,y,col);
			else{
				clear_rectangle(x,y,dx,dy);
				draw_round_rectangle(x,y,dx,dy,COMPARTMENT_CURVE*dy,col,BLUE);
				let si = dy*si_cogr_text_frac;

				let tsa = text_sup_anno(te,si,dx);
				center_text_tsa(tsa,x+dx/2,y+dy/2+0.1*si,te_col);
				
				let val = bu.value;
				if(val != undefined){
					let si_num = si/2; 
					let wmax = dx*0.9;
					do{
						let fo = get_font(si_num,"","times")
						let w = text_width(val,fo);
						if(w > wmax){
							si_num *= 0.95*wmax/w;
						}
						else{
							center_text(val,x+dx/2,y+dy/2+0.6*si,fo,te_col);
							break;
						}
					}while(si_num > 0.2);
				}
			}
		}
		return;
		
	case "CompLatLngGraph": case "CompLatLngGraph2":
		{
			let col = bu.col;	if(ov) col = shift_colour(col);
			if(dx*inter.sca < 1) plot_pixel(x,y,col);
			else{
				fill_circle(x+dx/2,y+dy/2,dx/2,col,bu.col_dark,NORMLINE);
		
				let val = bu.value;
				if(val != undefined){
					let si_num = dx/2; 
					if(si_num > 0.2){
						while(si_num > 0.2 && text_width(val,get_font(si_num,"","times")) > dx){
							si_num *= 0.7;
						}
						center_text(val,x+dx/2,y+dy/2+0.3*si_num,get_font(si_num,"","times"),white_black(col));
					}
				}
			}
		}
		return;
		
	case "Compartment": case "CompGraph": case "CompPop":
		{
			let col = bu.col;
			if(ov && bu.ac != "NotEditCompPopPercent") col = shift_colour(col);
		
			if(dx*inter.sca < 1) plot_pixel(x,y,col);
			else{
				let sh = 0.3;
				if(bu.type == "CompPop"){ sh = 0.35; clear_rectangle(x,y,dx,dy);}
			
				draw_round_rectangle(x,y,dx,dy,COMPARTMENT_CURVE*dy,col,bu.col_dark,MEDIUMLINE);
				let si = dy*si_comp_text_frac;
			
				if(si > 0.3){
					let tsa = text_sup_anno(te,si,dx);
					center_text_tsa(tsa,x+dx/2,y+dy/2+sh*si,white_black(col));
				}
				
				if(inter.mode.type == undefined && bu.show != true){
					if(bu.ac == undefined || bu.ac == "NotEditCompPopPercent"){
						cv.globalAlpha = 0.5;
						draw_round_rectangle(x,y,dx,dy,COMPARTMENT_CURVE*dy,WHITE,WHITE);
						cv.globalAlpha = 1;
					}
				}
			}
		}
		return;
	}
	
	x = nearest_pixel(x); y = nearest_pixel(y); 
	dx = nearest_pixel(dx); dy = nearest_pixel(dy); 
		
	let val = bu.val, val2 = bu.val2, ac = bu.ac;
	
	switch(bu.type){ 
	case "Rect":
		fill_rectangle(x,y,dx,dy,val); 
		break;
		
	case "WarningStore":
		{	
			let mar = 0.05;		
			let col = DDGREY; if(ov) col = DGREY;
			let sh = 0.7;
			draw_round_rectangle(x+mar,y-sh,dx-2*mar,dy+sh-mar,sh,col,col);  
			draw_image(find_pic("warn_black"),x+mar+0.4,y+0.01,1.3,1.1);
		
			let si = 0.8;
			let font = get_font(si);
			center_text("View model errors",x+dx/2+0.6,y+si+0.1,font,WHITE,dx);
		}
		break;

	case "Timebar":
		{
			clear_rectangle(x,y,dx,dy);
			let r = 0.3, r2 = 0.4;
	
			let anim = inter.graph.animation;
			let frac = anim.playframe/anim.playframe_max;
			let xx = x+dy/2+frac*(dx-dy);
		
			fill_rectangle(x+dy/2,y+dy/2-0.15,dx-dy,0.3,LLGREY);
			
			fill_rectangle(x+dy/2,y+dy/2-0.15,xx-(x+dy/2),0.3,RED);
			
			if(ov) fill_circle(xx,y+dy/2,r2,LRED,LRED);
		
			fill_circle(xx,y+dy/2,r,RED,RED);
		}
		break;
	
	case "DensityBut":
		inter.graph.density_plot(x,y,dx,dy);
		break;
		
	case "Circle":
		fill_circle(x+dx/2,y+dy/2,dx/2,bu.col,bu.col2);
		break;
		
	case "PlayButton":
		{
			clear_rectangle(x,y,dx,dy);
			let col = BLACK; if(ov) col = GREY;	
			if(inter.graph.animation.playing == false){
				let ddx = 0.5, ddy = 0.7;
				let polypoint=[];    
				polypoint.push({x:x+dx/2-ddx, y:y+dy/2+ddy});
				polypoint.push({x:x+dx/2-ddx, y:y+dy/2-ddy});
				polypoint.push({x:x+dx/2+ddx, y:y+dy/2});
				draw_polygon(polypoint,col,col,0);
			}
			else{
				let ddx = 0.15, ddx2 = 0.4, ddy = 0.7;
				fill_rectangle(x+dx/2-ddx-ddx2,y+dy/2-ddy,ddx2,2*ddy,col);
				fill_rectangle(x+dx/2+ddx,y+dy/2-ddy,ddx2,2*ddy,col);
			}
		}
		break;
		
	case "PlayBackward": case "PlayForward":
		{
			clear_rectangle(x,y,dx,dy);
			let col = BLACK; if(ov) col = GREY;
			
			let ddx = 0.2, ddy = 0.35;
			if(bu.type == "PlayBackward") ddx *= -1;
			for(let lo = 0; lo < 2; lo++){
				let shift = 0.6*(lo-0.5);
				let polypoint=[];    
				polypoint.push({x:x+dx/2-ddx+shift, y:y+dy/2+ddy});
				polypoint.push({x:x+dx/2-ddx+shift, y:y+dy/2-ddy});
				polypoint.push({x:x+dx/2+ddx+shift, y:y+dy/2});
				draw_polygon(polypoint,col,col,0);
			}
		}
		break;
	
	case "Slider":
		{
			let d = slider_dx/2;
			let info = bu.info;
			let fr = (info.value-info.min)/(info.max-info.min);
			draw_line(x+d,y+dy/2,x+dx-d,y+dy/2,BLACK,2);  
			let xx = x+d+fr*(dx-2*d);
			let col = BLACK; if(ov) col = GREY;
			fill_rectangle(xx-d,y,2*d,dy,col);  
		}
		break;
		
	case "EmptyRect":
		draw_rectangle(x,y,dx,dy,val,1);  
		break;
		
	case "ExpandIcon":
		{
			fill_rectangle(x,y,dx,dy,WHITE);
			let col = BLACK; if(ov) col = GREY;
			let xmid = x+dx/2,ymid = y+dy/2;
			let d = 0.4*dx;
			corner_arrow(xmid,ymid,xmid+d,ymid+d,col,THICKLINE);
			corner_arrow(xmid,ymid,xmid-d,ymid+d,col,THICKLINE);
			corner_arrow(xmid,ymid,xmid+d,ymid-d,col,THICKLINE);
			corner_arrow(xmid,ymid,xmid-d,ymid-d,col,THICKLINE);
		}
		break;

	case "Box":
		{
			clear_rectangle(x,y,dx,dy);
			let ruler = model.get_ruler();
			let scale = model.get_scale();
		
			let mar = 0.2*bu.si;
			let li = NORMLINE; 
			
			cv.globalAlpha = 0.1;
			draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,scale*ruler,bu.col,dark_colour(bu.col),li);  
			cv.globalAlpha = 1;
			let fo = get_font(bu.si);
			center_text(te,x+dx/2,y+bu.si*1.2,fo,bu.col,dx-1);             
		}
		break;
		
	case "Tab":
		if(bu.ac != undefined || val == inter.pa){
			let col = WHITE, col_back = BLACK;
			
			let si = 0.8;
			let font = get_font(si);
			let arrowdown = false;

			if(bu.ac == undefined){
				if(inter.page[inter.pa].sub.length > 0) arrowdown = true;
				col_back = WHITE;
				col = BLACK;
			}
			else{
				if(ov) col_back = "#888888";
			}
			
			if(bu.ac == undefined && val != inter.pa) cv.globalAlpha = 0.3;
			fill_rectangle(x,y,dx,dy,col_back);
			
			plot_text(te,x+0.4,y+dy/2+0.4*si,font,col);
			
			if(arrowdown == true){
				draw_menu_arrow_down(x+dx-0.7,y+dy/2,0.25,BLACK);
			}			
			cv.globalAlpha = 1;
		}
		break;
		
	case "PageSub":
		{	
			let col = WHITE, col_back = "#333333";
			
			let si = 0.8;
			let font = get_font(si);
			let arrowdown = false;

			let pag = inter.page[inter.pa];
			
			if(bu.ac == undefined){
				if(pag.sub[pag.index].sub.length > 0) arrowdown = true;
				else fill_rectangle(x+dx-0.3,y+dy/2-0.3,0.3,0.6,WHITE); 
	
				col_back = "#eeeeee";
				col = BLACK;
			}
			else{
				if(ov) col_back = "#777777";
			}
			
			fill_rectangle(x,y,dx,dy,col_back);
			
			plot_text(te,x+0.4,y+dy/2+0.4*si,font,col);
			
			if(arrowdown == true){
				draw_menu_arrow_down(x+dx-0.7,y+dy/2,0.25,BLACK);
			}			
			
		}
		break;
	
	
	case "AddSpecies":
		{
			let si = 0.7;
			let font = get_font(si,"bold");
			
			let col = WHITE; if(ov == 1) col = GREY;
			 
			fill_rectangle(x, y, dx, dy, BACKGROUND);
			
			draw_plus(x+1.0,y+dy/2,0.32,col,BLACK);
			
			plot_text(te,x+1.5, y+dy/2+0.42*si,font,col);
		}
		break;
		
		case "PageSubSub":
		{	
			let col = WHITE, col_back = BLACK;
			
			let si = 0.8;
			let font = get_font(si);
		
			let pag = inter.page[inter.pa];
			let sub = pag.sub[pag.index];
	
			fill_rectangle(x,y,dx,dy,col_back);
			
			if(val == sub.index){
				fill_rectangle(x+dx-0.3,y+dy/2-0.3,0.3,0.6,WHITE); 
			}
			else{
				if(ov) col = GREY;
			}
			plot_text(te,x+0.4,y+dy/2+0.4*si,font,col,dx-1);	
		}
		break;
		
	case "PageSubSubSub":
		{	
			let col = WHITE, col_back = BLACK;
			
			let si = 0.8;
			let font = get_font(si);
		
			let pag = inter.page[inter.pa];
			let sub = pag.sub[pag.index];
			let subsub = sub.sub[sub.index];
	
			fill_rectangle(x,y,dx,dy,col_back);
			
			if(val == subsub.index){
				let dyy = 0.5;
				rem_triangle(x+0.9,y+dy/2-dyy/2,-0.3,dyy,col) 
			}
			else{
				if(ov) col = GREY;
			}
			
			plot_text(te,x+1.1,y+dy/2+0.4*si,font,col,dx-1.4);	
		}
		break;
	
	case "Selected":
		{
			let mar = 0.1;
			clear_rectangle(x,y,dx,dy);
			x += mar; y += mar; dx -= 2*mar; dy -= 2*mar;
			cv.globalAlpha = 0.05;
			fill_rectangle(x,y,dx,dy,BLACK); 
			cv.globalAlpha = 0.1;
			draw_rectangle(x,y,dx,dy,BLACK,NORMLINE); 
			cv.globalAlpha = 1;
		}
		break;
		
	case "Corners":
		if(inter.mode && inter.mode.type == "Drag_Selector"){	
			let box = select_box();
			
			if(box.x1 < x || box.x2 > x+dx || box.y1 < y || box.y2 > y+dy){
				mode_off();
			}
			else{	
				if(box.on){			
					cv.globalAlpha = 0.1;
					fill_rectangle(box.x1,box.y1,box.x2-box.x1,box.y2-box.y1,BLACK); 
					cv.globalAlpha = 0.3;
					draw_rectangle(box.x1,box.y1,box.x2-box.x1,box.y2-box.y1,BLACK,NORMLINE); 
					cv.globalAlpha = 1;
				}
			}
		}
		
		draw_corners(x,y,dx,dy,1.5,BACKGROUND,RED)
		break;
		
	case "Logo":
		{
			draw_image(find_pic("logo_big"),x,y,dx,dy);
			let mar = 0.2;
			if(ov) draw_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,WHITE,2); 
		}
		break;
		
	case "Settings":
		if(ov){ cv.globalAlpha = 0.4; clear_rectangle(x,y,dx,dy);}
		draw_image(find_pic("settings"),x,y,dx,dy);
		if(ov) cv.globalAlpha = 1;
		break;
		
	case "BayesFactor":
		{
			let mar = 0.05;
			if(ov){ cv.globalAlpha = 0.4; clear_rectangle(x,y,dx,dy);}
			draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,0.2,WHITE,BLACK);  
			center_text("BF",x+dx/2+0.05,y+1.05,get_font(0.8,"bold"),BLACK,dx); 
			if(ov) cv.globalAlpha = 1;
		}
		break;
		
	case "Title":
		{	
			let word = bu.word;
			for(let i = 0; i < word.length; i++){
				let wo = word[i];
				plot_text(wo.te,x+wo.x,y+wo.y+0.12,wo.font,wo.col);
			}
			draw_line(x,y+1.5,x+dx,y+1.5,bu.col_line,NORMLINE);
		}
		break;
		
	case "MiniMenu":
		plot_text(te,x+0.2,y+1.1,get_font(0.8),DBLUE);
		draw_line(x,y+1.5,x+dx,y+1.5,DDBLUE,NORMLINE);
		break;
		
	case "SubTitle":	
		plot_text(te,x+0.2,y+1.1,get_font(1,"bold"),bu.col);
		break;
		
	case "Paragraph":
		{
			let word = bu.word;
			
			for(let i = 0; i < word.length; i++){
				let wo = word[i];
				plot_text(wo.te,x+wo.x,y+wo.y,wo.font,wo.col);
			}
		}
		break;

	case "BubbleLink":
		{
			fill_rectangle(x,y,dx,dy,BUBBLE_COL);
			let col = BLUE; if(ov) col = LBLUE;
			plot_text(te,x,y+bu.si,bu.font,col,dx);
		}
		break;
		
	case "NameLink":
		{
			fill_rectangle(x,y,dx,dy,WHITE);
			let col = BLUE; if(ov) col = LBLUE;
			plot_text(te,x,y+bu.si,bu.font,col,dx);
		}
		break;
		
	case "Link":
		{
			fill_rectangle(x,y,dx,dy,WHITE);
			let col = BLUE; if(ov) col = LBLUE;
			plot_text(te,x,y+bu.si,bu.font,col,dx);
		}
		break;
		
	case "LinkBlue":
		{
			fill_rectangle(x,y,dx,dy,BLUE_BACK);
			let col = LRED; if(ov) col = RED;
			plot_text(te,x,y+bu.si,bu.font,col,dx);
		}
		break;
		
	case "LinkPara":
		{
			fill_rectangle(x,y,dx,dy,WHITE);
			let col = BLUE; if(ov) col = LBLUE;
			let word = bu.word;
		
			for(let i = 0; i < word.length; i++){
				let wo = word[i];
				plot_text(wo.te,x+wo.x+0.3,y+wo.y+0.15,wo.font,col);
			}
		}
		break;

	case "ParamText":
		if(bu.back_col) fill_rectangle(x,y,dx,dy,bu.back_col);
		else clear_rectangle(x,y,dx,dy);
	
		plot_param_text(te,x,y+dy/2+0.35*bu.si,bu.si,bu.col,dx-0.2);
		break;
	
	case "LoadMessage":
		{
			if(inter.loading_symbol) te = inter.loading_symbol.message;
			if(te){
				fill_rectangle(x,y,dx,dy,WHITE);
				let col = BLACK; if(ov) col = GREY;
				let si = 0.7;
				center_text(te,x+dx/2,y+dy/2+0.4*si,get_font(si),col);     
			}
		}
		break;
		
	case "Stop":
		{
			fill_rectangle(x,y,dx,dy,WHITE);
			if(!ov) cv.globalAlpha = 0.3;
			let mar = 0.1;
			draw_image(find_pic("stop"),x+mar,y+mar,dx-2*mar,dy-2*mar);
			cv.globalAlpha = 1;     
		}
		break;
		
	case "Message":
		{
			fill_rectangle(x,y,dx,dy,WHITE);
			let col = BLACK; if(ov) col = GREY;
			let si = 0.7;
			center_text(bu.te,x+dx/2+0.5,y+dy/2+0.4*si,get_font(si),col);     
		}
		break;
			
	case "Script Line":
		{
			if(bu.back_col) fill_rectangle(x,y,dx,dy,bu.back_col);
			else clear_rectangle(x,y,dx,dy);
			
			let text_anno = bu.text_anno;
			
			let word = text_anno.word;
			for(let i = 0; i < word.length; i++){
				let wo = word[i];
				plot_text(wo.te,x+wo.x+0.3,y+wo.y,wo.font,wo.col);
			}
		}
		break;
		
	case "Letter":
		{
			fill_circle(x+dx/2,y+dy/2,dx/2-0.1,"#FFFF00",BLACK,MEDIUMLINE);
			//fill_circle(x+dx/2,y+dy/2,dx/2-0.1,WHITE,BLACK);
			
			let si = 1.3;
			let sh = 0;
			switch(bu.te){
			case "A": sh = 5; break;
			case "B": sh = 8; break;
			case "C": sh = -1; break;
			case "D": sh = 7; break;
			case "E": sh = 5; break;
			case "F": sh = 1; break;
			case "G": sh = 4; break;
			case "H": sh = 4; break;
			case "I": sh = 3; break;
			}
			
			center_text(bu.te,x+dx/2+sh*0.01,y+dy/2+0.34*si,get_font(si),BLACK,dx);
			//plot_text(te,x,y+dy/2+0.35*bu.si,bu.font,bu.col,ddx);
		}
		break;
	
	case "Text": case "LabelText":
		{
			if(bu.back_col) fill_rectangle(x,y,dx,dy,bu.back_col);
			else clear_rectangle(x,y,dx,dy);
		
			let ddx = dx; if(ddx != undefined) ddx -= 0.2;
			
			if(bu.comp_col){
				let tsa = text_sup_anno(te,bu.si,ddx,"arial");
				plot_text_tsa(tsa,x,y+dy/2+0.35*bu.si,bu.col);
			}
			else{
				plot_text(te,x,y+dy/2+0.35*bu.si,bu.font,bu.col,ddx);
			}
		}
		break;
		
	case "SubTitle":
		if(bu.back_col) fill_rectangle(x,y,dx,dy,bu.back_col);
		else clear_rectangle(x,y,dx,dy);
		
		plot_text(te,x,y+dy/2+0.35*si,fo,DBLUE,dx);
		break;
		
	case "RightText":
		if(bu.back_col) fill_rectangle(x,y,dx,dy,bu.back_col);
		else clear_rectangle(x,y,dx,dy);
		
		plot_text(te,x+dx-text_width(te,bu.font),y+dy/2+0.35*bu.si,bu.font,bu.col,dx);
		break;
		
	case "LatLngLabel":
		{
			if(bu.si > 0.3){
				let tsa = text_sup_anno(te,bu.si);
				center_text_tsa(tsa,x,y+0.35*bu.si,bu.col);
			}
		}
		break;
		
	case "RowNumber":
		let col = WHITE; if(ov) col = LLRED;
		fill_rectangle(x,y,dx,dy,col);
	
		let www = text_width(te,bu.font);
	
		plot_text(te,x+dx-www-0.4,y+dy/2+0.35*bu.si,bu.font,bu.col);
		break;

	case "Column":
		{
			let col = WHITE; if(ov) col = LLBLUE;
			fill_rectangle(x,y,dx,dy,col);
				
			let space = bu.space;
			plot_text(te[0],x+0.3,y+0.7*space, get_font(si_table,"bold"),BLACK,dx);
			
			let rmi = Math.floor((0-0.7*space-y)/space); if(rmi < 1) rmi = 1; 
			let rma = Math.floor(2+(this.dy-0.7*space-y)/space); if(rma > te.length) rma = te.length;	
			for(let r = rmi; r < rma; r++){
				plot_text(te[r],x+0.3,y+r*space+0.7*space,get_font(si_table),BLACK,dx);
			}
		}
		break;
		
	case "Element":
	case "ParamElement": case "ParamElementConst": case "ParamElementEq": 
	case "ParamFactorConst": case "ParamWeightConst":
	case "PriorSplitElement": case "DistSplitElement":
	case "ParamSimElement": case "DistSimElement": 
	case "ReparamElement": case "ReparamTableElement":  
	case "PriorElement": case "DistElement": case "SplineKnots": 
	case "SmoothValue": case "Amatrix": case "Xvector": case "AmatrixElement": case "XvectorElement":
		if(bu.ac == undefined && (bu.type == "ParamElement" || bu.type == "ReparamElement" || bu.type == "ParamElementConst")){
			let d = 0.4, d2 = 0.2;
			fill_rectangle(x+d2,y+d2,dx-2*d2,dy-2*d2,LLLGREY);
			draw_line(x+d,y+d,x+dx-d,y+dy-d,GREY,NORMLINE);
			draw_line(x+d,y+dy-d,x+dx-d,y+d,GREY,NORMLINE);
		}
		else{
			let col = WHITE; if(ov) col = LLBLUE;	
			
			if(selected(bu)) col = LLRED;
			else{
				if(te == inter.bubble.find && inter.bubble.mode == "SearchResult"){
					if(selected_column(bu)) col = LLLGREY;
				}
			}
			
			fill_rectangle(x,y,dx,dy,col);
			let colte = BLACK; 
			if(bu.type != "Element") colte = BLUE;
			
			if(ac == undefined){
				colte = GREY;			
				if(bu.type == "ParamSimElement" && bu.source != model) colte = BLACK;
			}
		
			let ddy = dy; if(bu.type == "Element") ddy = dy_table;
		
			if(bu.head){
				plot_text(te,x+0.3,y+0.7*ddy,bu.font,colte); 
			}
			else{
				switch(bu.type){
				case "ReparamElement": case "ReparamTableElement":  
					{
						let tsa = text_sup_anno(te,si_table,dx-0.5);
						plot_text_tsa(tsa,x+0.3,y+0.7*ddy,colte);
					}
					break;
					
				default:
					plot_text(te,x+0.3,y+0.7*ddy,bu.font,colte,dx);
					break;
				}
			}
		}
		break;
		
	case "TooBigElement": case "No element":	
		plot_text(String(te),x+0.3,y+0.7*dy,bu.font,BLACK,dx-0.5);
		break;
	
	case "ReparamEqn":
		{	
			let col = WHITE; if(ov) col = LLBLUE;	
			if(selected(bu)) col = LLRED;
			fill_rectangle(x,y,dx,dy,col);
			let si = 1.3;		
			plot_equation(bu.te,x+0.3,y+dy/2+0.3*si,si,x+dx-0.6,BLUE);
		}
		break;
		
	case "ProbEqn":
		{
			let col = BUBBLE_COL; if(ov) col = LLBLUE;
			fill_rectangle(x,y,dx,dy,col);
			plot_text(te,x+0.3,y+0.7*dy,bu.font,BLUE,dx-0.5);
		}
		break;
		
		
	case "ElementTooBig":
		fill_rectangle(x,y,dx,dy,WHITE);
		plot_text(te,x+0.3,y+0.7*dy,bu.font,BLACK,dx-0.5);
		break;
		
	case "CurvedOutline":
		{
			let mar = 0.2;
			draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,1,"no_fill",bu.col,6.5); 
			let fo = get_font(0.7);
			plot_text(te,x+2*mar+0.5,y+mar+1,fo,bu.col2);				
		}
		break;
		
	case "Outline":
		{
			let mar = 0.1;
			let back = bu.back_col; if(back == undefined) back = "no_fill";
			draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,0.5,back,bu.col);  
		}
		break;
	
	case "Char":
		{
			let cur = inter.cursor;
			let curi = cur.i;
			if(curi != undefined && cur.l == this.index){
				let sel = cur.select_pos;
				let i = bu.i;
			
				let back = false;
				if(sel != undefined){
					if(sel < curi){
						if(i >= sel && i < curi) back = true; 
					} 
					else{
						if(i >= curi && i < sel) back = true; 
					}
				}
				
				if(back == true) fill_rectangle(x,y,dx,dy,LLBLUE,LLBLUE); 
				
				if(sel == undefined){
					let ch_cu_back = bu.ch_cu_back;
					if(ch_cu_back != undefined && find_in(ch_cu_back,curi) != undefined){
						fill_rectangle(x,y,dx,dy,LGREY,LGREY); 
					}
				}
				
				if(bu.underline == true) fill_rectangle(x,y,dx,dy,LLRED,LLRED);
			
				if(curi == i){
					if(inter.cursor.on == 1){
						draw_line(x,y,x,y+dy,BLACK,NORMLINE);
					}
				}
			}
			
			let ymid = y+dy/2;

			if(ymid > 0 && ymid < this.dy){
				let yy = ymid+bu.si*0.35; if(bu.sh) yy += bu.sh;
				plot_text(te,x+textbox_space/2,yy,bu.font,bu.col);
			}
		}
		break;
	
	case "Cross":
		fill_rectangle(x,y,dx,dy,WHITE);
		draw_line(x,y,x+dx,y+dy,BLACK,NORMLINE);
		draw_line(x+dx,y,x,y+dy,BLACK,NORMLINE);
		break;
	
	case "Tick":
		{
			fill_rectangle(x,y,dx,dy,WHITE);
			let d = 0.3*dx;
			draw_line(x,y+dy-d,x+d,y+dy,BLACK,NORMLINE);
			draw_line(x+d,y+dy,x+dx,y,BLACK,NORMLINE);
		}
		break;
		
	case "NewModel":
		fill_rectangle(x,y,dx,dy,WHITE);
		if(ov) cv.globalAlpha = 0.3;
		draw_image(find_pic("newmod"),x,y,dx,dy);
		if(ov) cv.globalAlpha = 1;
		break;
		
	case "Pdf":
		fill_rectangle(x,y,dx,dy,WHITE);
		if(ov) cv.globalAlpha = 0.3;
		draw_image(find_pic("manual"),x,y,dx,dy);
		if(ov) cv.globalAlpha = 1;
		break;
	
	case "Pdf2":
		fill_rectangle(x,y,dx,dy,WHITE); 
		if(ov) cv.globalAlpha = 0.3;
		draw_image(find_pic("paper"),x,y,dx,dy);
		if(ov) cv.globalAlpha = 1;
		break;
	
	case "Example":
		{
			fill_rectangle(x,y,dx,dy,WHITE); 
				
			let si = 0.8;

			let col = DBLUE; if(bu.sel == true) col = DRED;
			if(ov) col = LBLUE;
			
			plot_text(te,x+dy*0.45, y+dy/2+0.35*si,bu.font,col);
		
			let ddx = dy*0.04;
			let ddx2 = dy*0.16, ddy = ddx2;
				
			draw_line(x+ddx,y+dy/2-ddy,x+ddx+ddx2,y+dy/2,col,THICKLINE)
			draw_line(x+ddx,y+dy/2+ddy,x+ddx+ddx2,y+dy/2,col,THICKLINE)
		}
		break;
	
	case "ShowModel":
		{
			fill_rectangle(x,y,dx,dy,WHITE); 
			let si = 0.8;
			let col = DBLUE; if(ov) col = LBLUE;		
			plot_text(te,x+0.1, y+dy/2+0.35*si,bu.font,col);
		
		}
		break;
		
	case "Bracket":
		let line = MEDIUMLINE;
		draw_line(x,y,x,y+dy,BLACK,line);
		draw_line(x,y,x+dx,y,BLACK,line);
		draw_line(x,y+dy,x+dx,y+dy,BLACK,line);
		draw_line(x-0.7*dx,y+dy/2,x,y+dy/2,BLACK,line);	
		break;
		
	case "ExampleModel":
		{
			fill_rectangle(x,y,dx,dy,WHITE);
			
			let W = dx;
			let H = 0.75*dy;
			
			let w = bu.pic.width;
			let h = bu.pic.height;
			if(w/h > W/H){
				let HH = W*(h/w);
				draw_image(bu.pic,x,y+H/2-HH/2,W,HH);
			}
			else{
				let WW = H*(w/h);
				draw_image(bu.pic,x+W/2-WW/2,y,WW,H);
			}
			
			let mar = 0.2;
				
			{
				let ma = 0.4;
				let TX = 1.4;
				let TY = TX*bu.mod_ty.height/bu.mod_ty.width;
				draw_image(bu.mod_ty,x+W-TX-ma,y+ma,TX,TY);
			}
			
			if(bu.mod_ty2){
				let ma = 0.4;
				let TX = 1.4;
				let TY = TX*bu.mod_ty.height/bu.mod_ty.width;
				draw_image(bu.mod_ty2,x+W-TX-ma,y+ma+1.7,TX,TY);
			}
			
			let word = bu.word;
			
			for(let i = 0; i < word.length; i++){
				let wo = word[i];
				let col = wo.col; if(bu.sel) col = DRED;
				plot_text(wo.te,x+wo.x,y+wo.y,wo.font,col);
			}
			
			if(ov){
				cv.globalAlpha = 0.2;
			
				let r = 1;
				let col = "#cccccc";
				draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,r,col,col,0);		
				cv.globalAlpha = 1;				
			}
		}
		break;
		
	case "PlotWord":
		{
			let word = bu.word;
			
			for(let i = 0; i < word.length; i++){
				let wo = word[i];
				plot_text(wo.te,x+wo.x,y+wo.y,wo.font,wo.col);
			}
		}
		break;
			
			
	case "CenterText":
		fill_rectangle(x,y,dx,dy,WHITE);
		center_text(te,x+dx/2,y+0.8*dy,bu.font,BLACK,dx-1);
		break;
		
	case "LeftText":
		fill_rectangle(x,y,dx,dy,WHITE);
		center_text(te,x+dx/2,y+0.8*dy,bu.font,BLACK);
		break;
		
	case "HelpIcon":
		{
			let back_col = WHITE; if(bu.help.back_col) back_col = bu.help.back_col;
			
			fill_rectangle(x,y,dx,dy,back_col);
			let col;
			if(back_col == WHITE || back_col == BUBBLE_COL){
				col = BLUE; if(ov) col = LBLUE;
			}
			else{
				col = WHITE; if(ov) col = LRED;
			}
			plot_text("[?]",x-0.05,y+0.4,get_font(0.6,"bold"),col);
		}
		break;
		
	case "HelpBackground":	case "EquationBackground":
		{
			cv.globalAlpha = 0.5;
			fill_rectangle(x,y,dx,dy,WHITE);
			fill_rectangle(x,y,dx,dy,BLACK);
			cv.globalAlpha = 1;		
		}
		break;
		
	case "HelpBack":
		draw_round_rectangle(x,y,dx,dy,1.5,BLUE_BACK,"#000000");  
		break;
	
	case "EquationBack":
		{
			let mar = 0.1;
			draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,1.5,GREY_BACK,"#000000");  
		}
		break;
	
	case "EquationAddBack":
		{
			draw_round_rectangle(x,y,dx,dy,1.4,BLACK,"#000000");  
			let xx = x + 1., yy = y + 2.1;
	
			plot_text(te,xx+0.2,yy-0.5,get_font(1.0),WHITE);
			draw_line(xx,yy,xx+dx-2.,yy,WHITE,NORMLINE);
		}
		break;

	case "WhiteCross":
		{
			let col = WHITE; if(ov) col = RED;
			let mar = 0.2;
			let xx = x+mar, yy = y+mar, dxx = dx-2*mar, dyy = dy-2*mar;
			draw_line(xx,yy,xx+dxx,yy+dyy,col,THICKLINE);
			draw_line(xx+dxx,yy,xx,yy+dyy,col,THICKLINE);
		}
		break;
		
	case "Blue":
		{
			let col = DDBLUE; col2 = DDBLUE; col3 = WHITE;  if(ov){ col = DBLUE;}
			let font = get_font(1,"bold");
			
			draw_round_rectangle(x,y,dx,dy,0.7,col,col2);
			center_text(te,x+dx/2,y+dy/2+0.35,font,col3);
		}
		break;
		
	case "White":
		{
			let col = WHITE, col2 = BLACK, col3 = BLACK; if(ov){ col = LLBLUE;}
			let font = get_font(1,"bold");
			
			draw_round_rectangle(x,y,dx,dy,0.7,col,col2);
			center_text(te,x+dx/2,y+dy/2+0.4,font,col3);
		}
		break;

	case "SmallButton":
		{
			let col = WHITE, col2 = BLACK, col3 = BLACK; if(ov){ col = LLBLUE;}
		
			let si = 0.9;
			let font = get_font(si);
			
			draw_round_rectangle(x,y,dx,dy,0.4,col,col2);
			center_text(te,x+dx/2,y+dy/2+si*0.4,font,col3);
		}
		break;
		
	case "Grey":
		{ 
			let col = "#333333"; if(ov) col = "#777777"; if(bu.ac == undefined) col = "#dddddd";
			let si = 1;
			let font = get_font(si,"bold");
			draw_round_rectangle(x,y,dx,dy,0.5*dy,col,col);
			center_text(te,x+dx/2,y+dy/2+0.37*si,font,WHITE);
		}
		break;
		
	case "GreyEdit":
		{ 
			let col = DDGREEN; if(ov) col = DGREEN; if(bu.ac == undefined) col = LLGREY;
			let si = 0.8;
			let font = get_font(si,"bold");
			draw_round_rectangle(x,y,dx,dy,0.5*dy,col,col);
			center_text(te,x+dx/2,y+dy/2+0.37*si,font,WHITE);
		}
		break;
		
	case "GreyView":
		{ 
			let col = DBLUE; if(ov) col = BLUE; if(bu.ac == undefined) col = LLGREY;
			let si = 0.8;
			let font = get_font(si,"bold");
			
			draw_round_rectangle(x,y,dx,dy,0.5*dy,col,col);
			center_text(te,x+dx/2,y+dy/2+0.37*si,font,WHITE);
		}
		break;
		
	
	case "ViewSmall": case "CombineIE":
		{ 
			let col = OK_butcol, col2 = OK_butcol2; 
			if(ov){ col = OK_butcol3; col2 = OK_butcol4;} 
			if(ac == undefined){ col = "#aaaaff"; col2 = "#aaaaff";}
			
			let si = 0.8;
			let font = get_font(si);
			draw_round_rectangle(x,y,dx,dy,0.4*dy,col,col2);
			center_text(te,x+dx/2,y+dy/2+0.4*si,font,WHITE);
		}
		break;
		
	case "RedInvalid":
		{ 
			let col = RED; if(ov) col = LRED;
			let si = 0.9;
			let font = get_font(si,"bold");
			draw_round_rectangle(x,y,dx,dy,0.5*dy,col,col);
			center_text(te,x+dx/2,y+dy/2+0.4*si,font,WHITE);
		}
		break;

	case "GreyBack":
		{ 
			let col = "#333333"; if(ov) col = "#777777";
			let si = 1;
			let font = get_font(si,"bold");
			draw_round_rectangle(x,y,dx,dy,0.5*dy,col,col);
			plot_text(te,x+1.35,y+dy/2+0.37*si,font,WHITE);
			let xx = x+0.55, mid = y+dy/2, d = 0.4;
			draw_line(xx,mid,xx+d,mid+d,WHITE,THICKLINE);	
			draw_line(xx,mid,xx+d,mid-d,WHITE,THICKLINE);	
		}
		break;
		
	case "GreyNext":
		{ 
			let col = "#333333"; if(ov) col = "#777777"; if(bu.ac == undefined) col = "#dddddd";
			let si = 1;
			let font = get_font(si,"bold");
			draw_round_rectangle(x,y,dx,dy,0.5*dy,col,col);
			plot_text(te,x+0.85,y+dy/2+0.37*si,font,WHITE);
			let xx = x+dx-0.55, mid = y+dy/2, d = -0.4;
			draw_line(xx,mid,xx+d,mid+d,WHITE,THICKLINE);	
			draw_line(xx,mid,xx+d,mid-d,WHITE,THICKLINE);	
		}
		break;
		
	case "Delete":
		{ 
			if(bu.back_col) fill_rectangle(x,y,dx,dy,bu.back_col);
			else clear_rectangle(x,y,dx,dy);
			
			let col = RED; if(ov) col = LRED;
			let mar = 0.2;			
			draw_line(x+mar,y+mar,x+dx-mar,y+dy-mar,col,MEDIUMLINE);	
			draw_line(x+mar,y+dy-mar,x+dx-mar,y+mar,col,MEDIUMLINE);	
		}
		break;
		
	case "MenuBut":
		{
			let col = LLLBLUE; if(ov) col = LLBLUE;
			let col2 = BLACK; if(ac == undefined) col2 = GREY
			
			fill_rectangle(x,y,dx,dy,col);
			let si = 0.8;
			let font = get_font(si);
			plot_text(te,x+0.5,y+dy/2+0.34*si,font,col2);
		}
		break;
		
	case "BubbleEndBut":
		{ 
			clear_rectangle(x,y,dx,dy);
			let col = OK_butcol, col2 = OK_butcol2; 
			if(ov){ col = OK_butcol3; col2 = OK_butcol4;} 
		
			if(!ac) cv.globalAlpha = 0.5;
			let si = 0.8;
			let font = get_font(si,"bold");
			draw_round_rectangle(x,y,dx,dy,0.5,col,col2);
			
			switch(te){
			case "Back":
				{
					let xx = x+0.40, mid = y+dy/2, d = 0.30;
					draw_line(xx,mid,xx+d,mid+d,WHITE,THICKLINE);	
					draw_line(xx,mid,xx+d,mid-d,WHITE,THICKLINE);	
					center_text(te,x+dx/2+0.2,y+dy/2+0.4*si,font,WHITE);
				}
				break;
			
			case "Next":
				{
					let xx = x+dx-0.40, mid = y+dy/2, d = -0.30;
					draw_line(xx,mid,xx+d,mid+d,WHITE,THICKLINE);	
					draw_line(xx,mid,xx+d,mid-d,WHITE,THICKLINE);	
					center_text(te,x+dx/2-0.2,y+dy/2+0.4*si,font,WHITE);
				}
				break;
				
			default:			
				center_text(te,x+dx/2,y+dy/2+0.4*si,font,WHITE);
				break;
			}
			if(!ac) cv.globalAlpha = 1;
		}
		break;
	
	case "Green":
		{
			let col = DGREEN, col2 = DDGREEN, col3 = WHITE; if(ov){ col = DDGREEN; col2 = DDDGREEN;}
			let font = get_font(1,"bold");
			
			draw_round_rectangle(x,y,dx,dy,0.7,col,col2);
			center_text(te,x+dx/2,y+dy/2+0.37,font,col3);
		}
		break;
	
	case "Red":
		{
			let col = RED, col2 = DRED, col3 = WHITE; if(ov){ col = DRED; col2 = DDRED;}
			let font = get_font(1,"bold");
			
			draw_round_rectangle(x,y,dx,dy,0.7,col,col2);
			center_text(te,x+dx/2,y+dy/2+0.37,font,col3);
		}
		break;
		
	case "YscrollUp":
		{
			let scol = get_scroll_col(bu.style,ov);
			let c1 = scol.c1, c2 = scol.c2, c3 = scol.c3;
		
			let col = c1; if(ov) col = c2;
			let col_line = c3;
			
			if(col == "none") clear_rectangle(x,y,dx,dy);
			else fill_rectangle(x,y,dx,dy,col); 
		
			let ddx = 0.3*dx;
			let ddy = -0.7*ddx;
			draw_line_noround(x+dx/2-ddx,y+dy/2-ddy,x+dx/2,y+dy/2+ddy,col_line,MEDIUMLINE);
			draw_line_noround(x+dx/2+ddx,y+dy/2-ddy,x+dx/2,y+dy/2+ddy,col_line,MEDIUMLINE);
		}	
		break;
		
	case "YscrollDown":
		{
			let scol = get_scroll_col(bu.style,ov);
			let c1 = scol.c1, c2 = scol.c2, c3 = scol.c3;
		
			let col = c1; if(ov) col = c2;
			let col_line = c3;
			
			if(col == "none") clear_rectangle(x,y,dx,dy);
			else fill_rectangle(x,y,dx,dy,col); 
		
			let ddx = 0.3*dx;
			let ddy = 0.7*ddx;
			draw_line_noround(x+dx/2-ddx,y+dy/2-ddy,x+dx/2,y+dy/2+ddy,col_line,MEDIUMLINE);
			draw_line_noround(x+dx/2+ddx,y+dy/2-ddy,x+dx/2,y+dy/2+ddy,col_line,MEDIUMLINE);
		}	
		break;
		
	case "XscrollLeft":
		{
			let scol = get_scroll_col(bu.style,ov);
			let c1 = scol.c1, c2 = scol.c2, c3 = scol.c3;
			
			let col = c1; if(ov) col = c2;
			let col_line = c3;
			
			if(col == "none") clear_rectangle(x,y,dx,dy);
			else fill_rectangle(x,y,dx,dy,col); 
		
			let ddy = 0.3*dy;
			let ddx = -0.7*ddy;
			draw_line(x+dx/2-ddx,y+dy/2-ddy,x+dx/2+ddx,y+dy/2,col_line,MEDIUMLINE);
			draw_line(x+dx/2-ddx,y+dy/2+ddy,x+dx/2+ddx,y+dy/2,col_line,MEDIUMLINE);
		}	
		break;
		
	case "XscrollRight":
		{
			let scol = get_scroll_col(bu.style,ov);
			let c1 = scol.c1, c2 = scol.c2, c3 = scol.c3;
			
			let col = c1; if(ov) col = c2;
			let col_line = c3;
			
			if(col == "none") clear_rectangle(x,y,dx,dy);
			else fill_rectangle(x,y,dx,dy,col); 
		
			let ddy = 0.3*dy;
			let ddx = 0.7*ddy;
			draw_line(x+dx/2-ddx,y+dy/2-ddy,x+dx/2+ddx,y+dy/2,col_line,MEDIUMLINE);
			draw_line(x+dx/2-ddx,y+dy/2+ddy,x+dx/2+ddx,y+dy/2,col_line,MEDIUMLINE);
		}	
		break;
		
	case "VerticalSlider": case "HorizontalSlider":
		{
			let c1 = LGREY, c2 = GREY, c3 = DGREY;
			switch(bu.style){
				case "White": c1 = WHITE; c2 = LLBLUE; c3 = LBLUE; break;
				case "Blue": c1 = HELP_BLUE; c2 = DBLUE; c3 = BLUE; break;
				case "Bubble": c1 = DBLUE; c2 = BLUE; c3 = BLUE; break;
				case "TabSlider": c1 = HELP_BLUE; c2 = DBLUE; c3 = BLUE; break;
				case "TextBox": 
					{
						c1 = HELP_BLUE; c2 = DBLUE; c3 = BLUE;
						
						if(bu.type == "VerticalSlider"){ let ddx = 0.4*dx; x += ddx*0.5; dx -= ddx;}
						else{ let ddy = 0.4*dy; y += ddy*0.5; dy -= ddy;}
					}
					break;
			}

			let col = c1; if(ov) col = c2;
			if(inter.mode.type == "Drag_VerticalSlider" || inter.mode.type == "Drag_HorizontalSlider"){
				if(inter.mode.pn == bu.pn) col = c3;
			}

			fill_rectangle(x,y,dx,dy,col); 
		}
		break;
		
	case "BackSlider":
		{
			let col = LLGREY;
			switch(bu.style){
			case "White": col = HELP_BLUE; break;
			case "Blue": col = BACKGROUND; break;
			case "Bubble": col = WHITE; break;
			case "TextBox": col = WHITE; break;
			case "TabSlider": col = WHITE; break;
			}

			fill_rectangle(x,y,dx,dy,col);
		}
		break;
			
	case "TextboxTab":
		{
			let col = GREY;
			let xx = x+0.9*dx, yy = y+dy/2;
			draw_line(x+0.5*dx,y+dy/2,xx,yy,col,NORMLINE);
			let r = 0.2*dx;
			draw_line(xx,yy,xx-r,yy-r,col,NORMLINE);
			draw_line(xx,yy,xx-r,yy+r,col,NORMLINE);
		}
		break;

	case "BulletPoint":
		{ 
			let ddx = 0.3, ddy = 0.3;
			fill_rectangle(x+dx/2-ddx/2,y+dy/2-ddy/2,ddx,ddy,bu.col);
		}
		break;
		
	case "CompAlpha": case "CompPrior":
		{
			let sh = 0.35; clear_rectangle(x,y,dx,dy);
			let col = bu.col;	if(ov) col = shift_colour(col);
			draw_round_rectangle(x,y,dx,dy,COMPARTMENT_CURVE*dy,col,dark_colour(col));
		
			let si = dy*si_comp_text_frac;
			let xx = x+1;
			let coll = white_black(col);
			
			plot_text(te,xx,y+dy/2+sh*si,bu.fo,coll);
			xx += text_width(te,bu.fo)+0.1;
			plot_text_tsa(bu.sup_tsa,xx,y+dy/2+0.05*si,coll);
			xx += bu.sup_tsa.w;
			plot_text(bu.te2,xx,y+dy/2+sh*si,bu.fo,coll);
		}
		break;
		
	case "MultiCompAlpha":
		{
			let sh = 0.35; clear_rectangle(x,y,dx,dy);
			
			let n = bu.col_list.length;
			for(let i = 0; i < n; i++){
				let col = bu.col_list[i];
				if(ov) col = shift_colour(col);
				fill_rectangle(x+(i*dx/n),y,dx/n,dy,col);
			}
			
			draw_inverse_round_rectangle(x,y,dx,dy,COMPARTMENT_CURVE*dy,WHITE);
			
			draw_round_rectangle(x,y,dx,dy,COMPARTMENT_CURVE*dy,"no_fill",BLACK);
			
			let si = dy*si_comp_text_frac;
			let xx = x+1;
			
			let coll = WHITE;
			
			plot_text(te,xx,y+dy/2+sh*si,bu.fo,coll);
			xx += text_width(te,bu.fo)+0.1;
			plot_text_tsa(bu.sup_tsa,xx,y+dy/2+0.05*si,coll);
			xx += bu.sup_tsa.w;
			plot_text(bu.te2,xx,y+dy/2+sh*si,bu.fo,coll);
		}
		break;
		
	case "CompLink":
		inter.graph.draw_complink(bu.lay);
		break;
		
	
	case "MultiComp":
		{
			clear_rectangle(x,y,dx,dy);
			
			let n = bu.col_list.length;
			for(let i = 0; i < n; i++){
				let col = bu.col_list[i]; if(ov) col = shift_colour(col);
				fill_rectangle(x+(i*dx/n),y,dx/n,dy,col);
			}
			
			draw_inverse_round_rectangle(x,y,dx,dy,COMPARTMENT_CURVE*dy,WHITE);
			
			draw_round_rectangle(x,y,dx,dy,COMPARTMENT_CURVE*dy,"no_fill",BLACK);
			
			let si = dy*si_comp_text_frac;
			center_text_tsa(bu.tsa,x+dx/2,y+dy/2+0.3*si,WHITE);
		}
		break;
		
	case "AddCompartmentMode":
		{
			clear_rectangle(x,y,dx,dy);
		
			cv.globalAlpha = 0.1;
			if(dx == dy) fill_circle(x+dx/2,y+dy/2,(dx/3),BLACK,BLACK);
			else draw_round_rectangle(x,y,dx,dy,dy/4,BLACK,BLACK);
			cv.globalAlpha = 1;
		}
		break;
			
	case "Transition":
		{
			let col = BLACK; 
			if(ov){
				if(inter.mode.type != "Drag_TransitionPoint") col = DGREY;
			}
			
			let cam = bu.cam;
			
			let scale = cam.scale;
			
			let si = 0.5*scale*cam.ruler; 
			
			let len = 0;
			
			let liw = MEDIUMLINE;
			if(big_eqn){ liw = THICKLINE; si *= 1.5;}
			
			for(let i = 0; i < bu.points.length - 1; i++){
				let colp = col;
				if(bu.tr.type == "" && i == bu.points.length - 2) colp = GREY;
				
				let po = bu.points[i];
				let po2 = bu.points[i+1];
				len += Math.sqrt((po2.x-po.x)*(po2.x-po.x) + (po2.y-po.y)*(po2.y-po.y));
				
				draw_line(po.x,po.y,po2.x,po2.y,colp,liw);
				if(bu.tr.f != undefined && bu.tr.f != SINK){
					if(i == bu.points.length-2){
						if(si > len/5) si = len/5;
						draw_arrow(po2.x,po2.y,po.x,po.y,si,colp);
					}
				}
			}
		
		
			let si_lab = 2.5*si;
			if(si_lab >  si_limit_label){
				let label = bu.tr.label;
		
				if(label != ""){
					let cen = bu.center;
				
					let x = cen.x, y = cen.y, nx = cen.nx, ny = cen.ny;
					let px = -ny, py = nx;
				
					let shift = -0.4;
					if(py < 0){
						shift = 0.4;
						if(bu.tr.pair == true) shift = -0.90;
						
						nx *= -1; ny *= -1;
					}
			
					plot_angle_equation(label,x+px*shift*si_lab,y+py*shift*si_lab,nx,ny,si_lab,cen.w*0.95*cam.scale);
				}
			}
		}
		break;
	
	case "Derive":
		{
			let tcol = bu.col;
			let col = WHITE; if(ov) col = LLBLUE;
			if(selected(bu)) col = LLRED;
			fill_rectangle(x,y,dx,dy,col); 
	
			if(bu.par != undefined){
				let info = bu.par.label_info;
				plot_paramlabel(bu.x_eq-bu.par.label_info.dx-0.7,y,info.dx,dy,BLACK,info);
			}
			else{
				let fo = get_font(bu.si,"","times");
				plot_text(bu.te_eqn1,bu.x_eq-text_width(bu.te_eqn1,fo)-0.7,y+dy/2+0.35*bu.si,fo,BLACK,dx);
			}
			
			let si = 1.5;
			plot_text("=",bu.x_eq,y+dy/2+0.35*si,get_font(si),BLACK);
			
			si = 1.4;
			plot_equation(bu.te_eqn2,bu.x_eq+1.6,y+dy/2+0.3*si,si,x+dx-(bu.x_eq+1.6),tcol);
		}
		break;
		
	case "Table":
		{
			let tab = bu.tab;
			let col = BLACK;//WHITE;
			
			//fill_rectangle(x,y,dx,dy,DBLUE); 
			fill_rectangle(x,y,dx,dy,LLBLUE); 
			
			for(let i = 0; i <= bu.ncol; i++){
				let xx = x+i*dx/bu.ncol;
				draw_line(xx,y,xx,y+dy,col,NORMLINE);
			}
			
			for(let i = 0; i <= bu.nrow; i++){
				let yy = y+i*dy/bu.nrow;
				draw_line(x,yy,x+dx,yy,col,NORMLINE);
			}
			
			let fo_bold = get_font(1,"bold");
			let fo_norm = get_font(1);
			
			let w = dx/bu.ncol;
			for(let j = 0; j < bu.nrow; j++){
				let fo = fo_norm; if(j == 0) fo = fo_bold;
				for(let i = 0; i < bu.ncol; i++){
					let xx = x+i*w;
					let yy = y+j*dy/bu.nrow;
					plot_text(tab[j][i],xx+0.5,yy+1.1,fo,col,w-1); 
				}
			}
		}
		break;
		
	case "PlotEquation":
		{
			let col = WHITE; if(ov) col = LLBLUE;
			if(selected(bu)) col = LLRED;
			fill_rectangle(x,y,dx,dy,col); 
			plot_equation(bu.te,x+0.3,y+dy/2+0.3*bu.si,bu.si);
		}
		break;
		
	case "TransitionPoint":
		{
			switch(bu.variety){
			case MIDPOINT:
				if(ov) fill_circle(x+dx/2,y+dy/2,0.3,DGREY,DGREY,NORMLINE); 
				break;
				
			case SOURCE: case SINK:
				{
					let col = BLACK; if(ov) col =  DGREY; 
					let r = dx/2-0.1;
					if(thick_line) r *= 2;
					let xx = x+dx/2, yy = y+dy/2;
					fill_circle(xx,yy,r,col,col,NORMLINE);
					let r2 =0.6*r;
					
					if(bu.variety == SOURCE) draw_cross(xx,yy,0.6*r,WHITE);
					else draw_minus(xx,yy,0.6*r,WHITE);
				}
				break;
			}
		}
		break;
		 
	case "AddButton":
		{
			if(bu.op.back == undefined)	clear_rectangle(x,y,dx,dy);
			else draw_round_rectangle(x+0.1,y+0.1,dx-0.2,dy-0.2,0.3,bu.op.back,bu.op.back);
		
			let si = 0.8;
			
			let col = BLACK;
			if(bu.ac == undefined) col = LLBLUE;
			else{
				if(ov == 1) col = BLUE;
			}
				
			plot_text(te,x+1,y+dy/2+0.44*bu.si,bu.font,col);		
			draw_plus(x+0.5,y+dy/2,0.35,col,WHITE);
		}
		break;
		
		case "Grid":
			plot_grid(bu.p,bu.cl,dx,dy)
			break;
			
		case "GridIcon":
		{
			clear_rectangle(x,y,dx,dy);
			
			let thick = NORMLINE;
			if(!bu.ac){
				cv.globalAlpha = 0.3;
			}
			else{
				let grid = model.species[bu.p].cla[bu.cl].camera.grid;
			
				if(grid=="on"){				
					if(ov) thick = MEDIUMLINE;
				}
				else{
					cv.globalAlpha = 0.3;
					if(ov) cv.globalAlpha = 0.6;
				}
			}
			
			let d = 0.3;
			let kmax = 3;
			for(let k = 0; k <= kmax; k++){
				let xx = x+d+k*(dx-2*d)/kmax;
				draw_line(xx,y+d,xx,y+dy-d,BLACK,thick);
				
				let yy = y+d+k*(dy-2*d)/kmax;
				draw_line(x+d,yy,x+dx-d,yy,BLACK,thick);
			}
			cv.globalAlpha = 1;	
		}
		break;
		
	case "ZoomIn":
		clear_rectangle(x,y,dx,dy);
		if(ov || ac == undefined) cv.globalAlpha = 0.3;
		draw_image(find_pic("zoomin"),x,y,dx,dy);
		cv.globalAlpha = 1;	
		break;
		
	case "ZoomOut":
		clear_rectangle(x,y,dx,dy);
		if(ov || ac == undefined) cv.globalAlpha = 0.3;
		draw_image(find_pic("zoomout"),x,y,dx,dy);
		cv.globalAlpha = 1;	
		break;
		
	case "ClassTab":
		{
			cv.globalAlpha = 1;
			clear_rectangle(x,y,dx,dy);
			
			let col = WHITE, col2 = BLACK;
			if(bu.cl == bu.cl_sel){ col = BLACK; col2 = WHITE;}
			else{
				if(ov) col = LLGREY;
			}
			draw_lower_round_rectangle(x+0.1,y,dx-0.2,dy-0.1,0.5,col,BLACK,NORMLINE);
			plot_text(te,x+1,y+dy/2+0.25*bu.si,bu.font,col2);
			let w = text_width(te,bu.font);
			
			plot_text(bu.index,x+1.5+w,y+dy/2+0.25*bu.si,bu.font_it,col2);
		}
		break;	
		
	case "File":
		{
			let si = 1;
			clear_rectangle(x,y,dx,dy);
			let col = DDGREY; if(ov) col = GREY;
			draw_lower_round_rectangle(x+0.1,y,dx-0.2,dy-0.1,0.5,col,col,NORMLINE);
			center_text(te,x+dx/2,y+dy/2+0.25*0.8,get_font(si),WHITE);
		}
		break;	
		
	case "Menu":
		{
			clear_rectangle(x,y,dx,dy);
			let col = BLACK; if(ov) col = DGREY;
			let padx = 0.2, pady = 0.5;
			let thick = 3; 
		
			let y1 = nearest_pixel(y+pady);
			let y2 = nearest_pixel(y+pady+((dy-2*pady)/2));
			
			let sca = inter.sca;
			if(inter.printing) sca /= print_scale_factor;
			
			let w = Math.floor(sca/5);
		
			if(w < 1) w = 1;
			for(let i = 0; i < 3; i++){
				let yy = y1+(y2-y1)*i;
				draw_line(x+padx,yy,x+dx-padx,yy,col,w);
			}				
		}
		break;
	
	case "InputBoxName":
		{
			if(bu.back_col != undefined) fill_rectangle(x,y,dx,dy,bu.back_col);
			else clear_rectangle(x,y,dx,dy);
			
			let si = 0.7;
			plot_text(te,x,y+dy/2+0.35*si,get_font(si,"bold"),BUBBLE_TEXT_COL);
		}
		break;
	
	case "RadioButton":
		{
			if(bu.back_col != undefined) fill_rectangle(x,y,dx,dy,bu.back_col);
			else clear_rectangle(x,y,dx,dy);
			
			let col = DGREY, col2 = WHITE;
			if(bu.ac == undefined) col = LGREY;
		
			if(bu.source.value == bu.value){
				col = DBLUE; col2 = DBLUE;
				if(bu.ac == undefined){ col = GREY; col2 = GREY;}
			}
			else{
				if(ov) col = GREY;
			}
			
			fill_circle(x+dx/2,y+dy/2,0.4,WHITE,col,THICKLINE);	 
			fill_circle(x+dx/2,y+dy/2,0.2,col2,col2,0);
		}
		break;
		
	case "RadioWhiteButton":
		{
			if(bu.back_col != undefined) fill_rectangle(x,y,dx,dy,bu.back_col);
			else clear_rectangle(x,y,dx,dy);
			
			let col = WHITE, col2 = BLACK;
			
			if(bu.source.value == bu.value){
				col = WHITE; col2 = WHITE;
			}
			else{
				if(ov) col = GREY;
			}
			
			fill_circle(x+dx/2,y+dy/2,0.4,BLACK,col,THICKLINE);	 
			fill_circle(x+dx/2,y+dy/2,0.2,col2,col2,0);
		}
		break;

	case "CheckboxButton":
		{
			if(bu.back_col != undefined) fill_rectangle(x,y,dx,dy,bu.back_col);
			else clear_rectangle(x,y,dx,dy);

			let col = DDGREY, col2 = WHITE;
		
			if(ov) col = GREY;
			
			let mar = 0.1;
			
			fill_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,WHITE);	

			check_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,col,bu.source.check,NORMLINE);	
		}
		break;

	case "CheckboxWhiteButton":
		{
			if(bu.back_col != undefined) fill_rectangle(x,y,dx,dy,bu.back_col);
			else clear_rectangle(x,y,dx,dy);
			
			let col = WHITE, col2 = BLACK;
		
			if(ov) col = GREY;
			
			let mar = 0.1;
			
			check_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,col,bu.source.check,NORMLINE);	
		}
		break;
		
	case "RadioButtonText":	case "CheckboxButtonText":
		{
			if(bu.back_col != undefined) fill_rectangle(x,y,dx,dy,bu.back_col);
			else clear_rectangle(x,y,dx,dy);
			 
			let fo = get_font(si_radio); if(bu.fo != undefined) fo = bu.fo;
		
			if(bu.source.param == true){
				plot_param_text(te,x+0.2,y+dy/2+0.4*si_radio,si_radio,bu.col,dx+5);
			}
			else{
				let tsa = text_sup_anno(te,si_radio,dx+5,"arial");
				plot_text_tsa(tsa,x+0.2,y+dy/2+0.4*si_radio,bu.col);
			}
		}
		break;
		
	case "BubbleTitle":
		{
			clear_rectangle(x,y,dx,dy);
			let si = si_bubble_title;
			plot_text(te,x,y+dy/2+0.37*si,get_font(si,"bold"),BUBBLE_TEXT_COL);
		}
		break;
		
	case "WarnPic":
		draw_image(find_pic("warn"),x,y,dx,dy);
		break;
		
	case "WarnPicBlack":
		draw_image(find_pic("warn_black"),x,y,dx,dy);
		break;
		
	case "CopyPic":
		clear_rectangle(x,y,dx,dy);
		if(ov) cv.globalAlpha = 0.5;
		draw_image(find_pic("copy_icon"),x,y,dx,dy);
		if(ov) cv.globalAlpha = 1;
		break;
		
	case "DescPic":
		draw_image(find_pic("desc"),x,y,dx,dy);
		break;
		
	case "Warning":
		{
			if(bu.back_col != undefined) fill_rectangle(x,y,dx,dy,bu.back_col);
			else clear_rectangle(x,y,dx,dy);
			let si = 0.8;
			center_text(te,x+dx/2,y+dy/2+0.37*si,get_font(si),RED);
		}
		break;
		
		
	case "ColourSelect":
		{
			fill_rectangle(x,y,dx,dy,BUBBLE_COL);
	
			let col = bu.col; 
	
			let sel = false;
			let col2 = bu.sel_bu.col;
			if(bu.sel_bu.type == "AddButton")	col2 = inter.bubble.radio.col;
			if(col == col2) sel = true;
				
			if(ov && sel == false) col = shift_colour(col);
			
			let gap = 0.2;
			
			fill_rectangle(x+gap,y+gap,dx-2*gap,dy-2*gap,col);
			if(sel == true) draw_rectangle(x+gap,y+gap,dx-2*gap,dy-2*gap,BLACK,2);
		}
		break;
	
	case "BubbleBack":
		{
			let col_back = BUBBLE_COL;
			let col = BUBBLE_SURROUND;

			draw_round_rectangle(x,y,dx,dy,1,col_back,col);
		
			let ta = bu.tag;
		
			draw_line(ta.x1,ta.y1,ta.x2,ta.y2,col_back,NORMLINE);	
		
			let polypoint=[];    
			polypoint.push({x:ta.x0, y:ta.y0});
			polypoint.push({x:ta.x1, y:ta.y1});
			polypoint.push({x:ta.x2, y:ta.y2});
			draw_polygon(polypoint,col_back,col_back,0);
			
			draw_line(ta.x0,ta.y0,ta.x1,ta.y1,col,NORMLINE);	
			draw_line(ta.x0,ta.y0,ta.x2,ta.y2,col,NORMLINE);
		}
		break;
		
	case "Feature":
		draw_feature(x,y,dx,dy,bu.polygon,LGREY,WHITE,NORMLINE);
		break;
		
	case "CompMap":
		{
			let col = bu.col; if(ov) col = shift_colour(col);
			if(bu.ac == undefined){ 
				if(col == MAP_DEFAULT) col = LRED; else col = dark_colour(col);
			}
			draw_feature(x,y,dx,dy,bu.polygon,col,WHITE,NORMLINE);
		}
		break;
		
	case "CompMapGraph":
		{
			let col = bu.col; if(ov) col = shift_colour(col);
			draw_feature(x,y,dx,dy,bu.polygon,col,WHITE,NORMLINE);
		}
		break;
		
	case "LeftArrow":
		clear_rectangle(x,y,dx,dy);
		draw_left_arrow(x,y,dx,dy,BLACK,DGREY,THICKLINE);
		break;
		
	case "CloseBubble":
		{
			cv.globalAlpha=1;
			fill_rectangle(x,y,dx,dy,BUBBLE_COL);
		
			let col = BLACK; if(ov) col = RED;//DGREY;
			let mar = 0.2;			
			draw_line(x+mar,y+mar,x+dx-mar,y+dy-mar,col,THICKLINE);	
			draw_line(x+mar,y+dy-mar,x+dx-mar,y+mar,col,THICKLINE);	
		}
		break;
	
	case "MessageMenu":
		dim_rectangle(x,y,dx,dy,BLUE);
		let si = 1;
		center_text(te,x+dx/2,y+dy/2+0.4*si,get_font(si),BLACK);
		break;
	
	case "MessageClose":
		{
			let col = BLACK; if(ov) col = RED;
			let mar = 0.2;			
			draw_line(x+mar,y+mar,x+dx-mar,y+dy-mar,col,THICKLINE);	
			draw_line(x+mar,y+dy-mar,x+dx-mar,y+mar,col,THICKLINE);	
		}
		break;
	
	case "PopFilt": case "SliceTime":
		{
			clear_rectangle(x,y,dx,dy);

			let mar = 0.05;
			let col = WHITE; if(ov) col = LLLBLUE;
		
			let col2 = DDBLUE, col3 = DDBLUE;
			draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,0.4,col,col2);
			
			let si = si_drop;
			let font = get_font(si);
			plot_text(te,x+0.5,y+dy/2+0.4*si,font,col3,dx-2);
		}
		break;
		
	case "Dropdown":
		{
			let mar = 0.05;
			
			clear_rectangle(x,y,dx,dy);
			
			let col = DDBLUE, col2 = DDBLUE, col3 = WHITE, col4 = DBLUE; 
			if(bu.style == "white_style"){
				col = WHITE; col2 = WHITE; col3 = BLACK; col4 = LLGREY;
				if(dy >= 1.1*dropdown_height){ col = LGREY; col2 = LGREY;}
			}
			
			let si = si_drop;
			let font = get_font(si);
		
			let yy = y;
			let ddy = 0.15;
			
			if(ac == undefined) cv.globalAlpha = 0.5;
			if(dy < 1.1*dropdown_height){
				if(ov){ col = col4;}
				draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,0.4,col,col2);
			}
			else{
				draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,0.4,WHITE,col2);
				if(bu.direction == "Up"){
					draw_round_menu_bottom(x+mar,y+dy-dropdown_height+mar,dx-2*mar,dropdown_height-2*mar,0.4,col,col2); 
					yy = y+dy-dropdown_height;
					ddy *= -1;
				}
				else{
					draw_round_menu_top(x+mar,y+mar,dx-2*mar,dropdown_height-2*mar,0.4,col,col2); 
				}
			}
			if(bu.source.param == true){
				plot_param_text(te,x+0.5,yy+dropdown_height/2+0.4*si,si,col3,dx-2);
			}
			else{
				let tsa = text_sup_anno(te,si,dx-2,"arial");
				plot_text_tsa(tsa,x+0.5,yy+dropdown_height/2+0.4*si,col3);
			}
		
			let x1 = x+dx-1.2, ddx = 0.25, y1 = yy+dropdown_height/2-ddy, y2 = yy+dropdown_height/2+ddy;
			draw_line(x1,y1,x1+ddx,y2,col3,THICKLINE);	
			draw_line(x1+2*ddx,y1,x1+ddx,y2,col3,THICKLINE);	

			if(ac == undefined) cv.globalAlpha = 1;
		}
		break;
		
	case "DropdownOption":
		{
			let col = WHITE, col2 = BLACK;
			if(ov) col = LLBLUE;
			fill_rectangle(x,y,dx,dy,col);
			let si = si_drop;
			let font = get_font(si);
	
			if(bu.source.param == true){
				plot_param_text(bu.pos.te,x+0.2,y+dy/2+0.4*si,si,col2,dx-0.4);
			}
			else{
				let tsa = text_sup_anno(bu.pos.te,si,dx-0.4,"arial");
				plot_text_tsa(tsa,x+0.2,y+dy/2+0.4*si,col2);
			}
		}			
		break;
		
	case "Calculator_But":
		{
			fill_rectangle(x,y,dx,dy,GREY_BACK);
			
			let col = bu.col; if(ov == 1) col = bu.colov;
			draw_round_rectangle(x,y,dx,dy,0.3,col,dark_colour(col));
			let si = 0.9;
			let font = get_font(si);
			center_text(te,x+dx/2,y+dy/2+0.35*si,font,WHITE);
		}
		break;
	
	case "Toolbar_But":
		{
			fill_rectangle(x,y,dx,dy,GREY_BACK);
			
			let col = "#333333"; 
			if(bu.te2 == inter.equation.toolbar) col = BLACK;
			else{
				if(ov == 1) col = DGREY;
			}
			
			draw_round_rectangle(x,y,dx,dy,0.3,col,dark_colour(col));
			let font = get_font(si_toolbut);
			center_text(te,x+dx/2,y+dy/2+0.35*si_toolbut,font,WHITE);
			
			if(ac == undefined){
				cv.globalAlpha = 0.6;
				fill_rectangle(x,y,dx,dy,GREY_BACK);
				cv.globalAlpha = 1;
			}
		}
		break;
	
	case "Letter_But":
		{
			fill_rectangle(x,y,dx,dy,BLACK);
			let col = bu.col, col2 = dark_colour(col); 
	
			if(inter.equation.letter == te){
				col2 = WHITE;
			}
			else{
				if(ov == 1) col = bu.colov;
			}
			
			draw_round_rectangle(x+0.05,y+0.05,dx-0.1,dy-0.1,0.1,col,col2,THICKLINE);
			let si = 1.1;
			let font = get_font(si);
			center_text(te,x+dx/2,y+dy/2+0.25*si,font,WHITE);
		}
		break;
	
	case "LoadingSymbol":
		fill_rectangle(x,y,dx,dy,WHITE);
		
		let per = inter.loading_symbol.percent;
	
		if(per != undefined){
			center_text(per,x+dx/2+0.02,y+dy/2+0.26,get_font(0.8,"bold"),BLACK);
		}
		plot_loading_symbol(x+dx/2,y+dy/2,1);
		break;
		
	case "TableSlider":
		break;
		
	case "TableSlider2":
		{
			let ddx = 0.2, ddx2 = 0.2, ddy = 1.3;
			let col = LGREY; if(ov) col = GREY;
			fill_rectangle(x+dx/2-ddx/2,y+dy/2-ddy/2,ddx,ddy,col);
		}
		break;
	
	case "Polygon":
		break;
		
	case "Select":
		let mo = inter.mode;
		if(mo.mx != undefined){
			mo.x1 = inter.mode.mx - bu.lay_x;
			mo.y1 = inter.mode.my - bu.lay_y;
			mo.x2 = inter.mx - bu.lay_x;
			mo.y2 = inter.my - bu.lay_y;
			cv.globalAlpha = 0.3;
			fill_rectangle(mo.x1,mo.y1,mo.x2-mo.x1,mo.y2-mo.y1,LBLUE);
			draw_rectangle(mo.x1,mo.y1,mo.x2-mo.x1,mo.y2-mo.y1,BLUE,NORMLINE);
			cv.globalAlpha = 1;
		}
		break;
		
	case "IEcombine":
		{
			let colb = WHITE; if(bu.back_col) colb = bu.back_col
			fill_rectangle(x,y,dx,dy,colb);
			let col = BLACK, bord = BLACK; if(ov){ col = DGREY; bord = GREY;} 
			
			let mar = 0.1;
			draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,dy*0.2,col,bord,NORMLINE);  
	
			center_text(bu.te,x+dx/2,y+0.7*dy,bu.font,WHITE);
		}
		break;
	
	case "ParamLabel":
		{
			let colb = WHITE; if(bu.back_col) colb = bu.back_col
			fill_rectangle(x,y,dx,dy,colb);
		
			let mar = 0.1;

			let col = LLBLUE;
			if(ov) col = LBLUE; 
			let bord = dark_colour(col);
			
			if(bu.col == WHITE){
				col = BLACK; bord = BLACK; if(ov){ col = DGREY; bord = GREY;} 
			}
		
			if(ac != undefined){
				draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,dy*0.2,col,bord,NORMLINE);  
			}

			plot_paramlabel(x,y,dx,dy,bu.col,bu.info);
		}
		break;
	
	case "ParamLabel2":
		{
			let mar = 0.1;
			let co = BLACK; if(ov) co = DGREY;
			let par = model.param[bu.i];
	
			switch(par.variety){
			case "normal":
				if(par.type == "derive_param"){ co = RED; if(ov) co = LRED;} 
				if(par.type == "fixed effect"){ co = DGREY; if(ov) co = GREY;} 
				
				let i = find(eqn_types,"name",par.type);
				if(i != undefined){
					if(eqn_types[i].obs_model == true){ co = DBROWN; if(ov) co = BROWN;}
				}					
				if(par.type == "variance" || par.type == "correlation"){ co = DCYAN; if(ov) co = CYAN;} 
				break;
			case "const": co = DBLUE; if(ov) co = BLUE; break;
			case "dist": co = DPURPLE; if(ov) co = PURPLE; break;
			case "reparam": co = DGREEN; if(ov) co = "#66bb66"; break;
			}
			
			co = BLACK; if(ov) co = DGREY;
			
			draw_round_rectangle(x+mar,y+mar,dx-2*mar,dy-2*mar,dy*0.2,co,WHITE,NORMLINE); 
			let col = WHITE;
			
			plot_paramlabel(x,y,dx,dy,col,bu.info);
		}
		break;
	
	case "x-axis":
		draw_line(x,y,x+dx,y+dy,BLACK,NORMLINE);    
		draw_arrow(x+dx-0.1,y+dy,x-0.1,y,0.7,BLACK);
		break;
		
	case "y-axis":
		draw_line(x,y,x+dx,y+dy,BLACK,NORMLINE);    
		draw_arrow(x,y,x+dx,y+dy,0.7,BLACK);
		break;
	
	case "x-tick":
		inter.graph.draw_xtick(x,y,dx,dy);
		break;
	
	case "x-tick-label":
		inter.graph.draw_xtick_label(x,y,dx,dy,bu.mar,bu.x_lab,bu.x_vert,bu.x_param,bu.x_si,ov);
		break;
		
	case "x-label":
		if(bu.param) center_param_text(te,x+dx/2,y+0.75*dy,si_graph_param_label,BLACK,20);
		else{
			let fo = get_font(si_graph_label); if(bu.italic) fo = get_font(si_graph_label,"italic","Times");
			center_text(te,x+dx/2,y+0.75*dy,fo,BLACK);
		}
		break;
		
	case "y-tick":
		inter.graph.draw_ytick(x,y,dx,dy);
		break;
		
	case "y-tick-label":
		inter.graph.draw_ytick_label(x,y,dx,dy,bu.mar,bu.y_lab,bu.y_hor,bu.y_param,bu.y_si,ov);
		break;
	
	case "y-label":
		if(bu.param) center_vert_param_text(te,x+1.1,y+dy/2,si_graph_param_label,BLACK,20);  
		else{
			let fo = get_font(si_graph_label); if(bu.italic) fo = get_font(si_graph_label,"italic","Times");
			center_vert_text(te,x+0.9,y+dy/2,fo,BLACK);  
		}
		break;
		
	case "VertText":
		vert_text(te,x+0.7*dx,y+dy-0.3,bu.font,bu.col);  
		break;
	
	case "GraphContentBut":
		inter.graph.draw_content_button(x,y,dx,dy);
		break;
		
	case "MatrixEleBut":
		{
			let col = bu.col;
			fill_rectangle(x,y,dx,dy,col);
			if(dx > 3){
				let si = 1.0; let fo = get_font(si,"bold");
				center_text(te,x+dx/2,y+dy/2+0.3*si,fo,white_black(col),dx-0.0);
			}
			
			let mar = 0.1;	
			if(ov){
				cv.globalAlpha = 0.1;
				fill_rectangle(x,y,dx,dy,WHITE);
				fill_rectangle(x,y,dx,dy,BLACK);
				cv.globalAlpha = 1;
			}
		}
		break;
		
	case "MatrixKey":
		{
			clear_rectangle(x,y,dx,dy);
			let ddx = 2;
			let mary = 0.5;
			let grad = bu.mkey.grad;
			let imax = grad.length;
			
			for(let i = 0; i < imax; i++){
				fill_rectangle(x,y+mary+(dy-2*mary)*(i/imax),ddx,dy/imax,grad[imax-1-i]);
			}
			draw_rectangle(x,y+mary,ddx,dy-2*mary,BLACK,THINLINE);
			
			let tick = bu.mkey.tick;
			let min = bu.mkey.min, max = bu.mkey.max; 
			for(let j = 0; j < tick.length; j++){
				let yy = y+mary+(dy-2*mary)*(1-(tick[j].value-min)/(max-min+TINY));
				draw_line(x+ddx,yy,x+ddx+0.3,yy,BLACK,NORMLINE);
				plot_text(tick[j].te,x+ddx+0.5,yy+0.3,get_font(0.8),BLACK,dx-(ddx+0.5)); 
			}
		}
		break

	case "IndTimelineBut":
		inter.graph.draw_ind_timeline_button(x,y,dx,dy,bu.name,bu.info,bu.col_timeline,bu.obs);
		break;
		
	case "InfBarBut":
		inter.graph.draw_ind_infbar_button(x,y,dx,dy,bu.name,bu.inf_line,bu.tmin,bu.tmax,bu.left_shift);
		break;
	
	case "TransArrow":
		inter.graph.draw_transarrow_button(x,y,dx,dy,bu.barh,bu.transmission);
		break;
		
	case "PhyloButton":
		inter.graph.draw_phylo_button(x,y,dx,dy,bu.barh,bu.pnode);
		break;
		
	case "Key":
		{
			clear_rectangle(x,y,dx,dy);
			let th = THICKLINE; if(bu.thick) th = bu.thick;
			draw_line(x+0.8,y+dy/2,x+2.5,y+dy/2,bu.col,th,bu.dash); 
			plot_equation(te,x+3,y+dy/2+0.4,1,dx-(x+3));
		}
		break;
		
	case "KeyRect":
		{
			clear_rectangle(x,y,dx,dy);
			let ddy = 0.7; 
			fill_rectangle(x+0.8,y+dy/2-ddy/2,1.7,ddy,bu.col); 
			plot_equation(te,x+3,y+dy/2+0.35,1,dx-(x+3));
		}
		break;
		
	case "KeyErrBar":
		{
			clear_rectangle(x,y,dx,dy);
			let xx = x+1.7, dxx = 0.4, mar = 0.1;
			draw_line(xx,y+mar,xx,y+dy-mar,bu.col,NORMLINE,bu.dash); 
			draw_line(xx-dxx,y+mar,xx+dxx,y+mar,bu.col,NORMLINE,bu.dash); 
			draw_line(xx-0.5*dxx,y+dy/2,xx+0.5*dxx,y+dy/2,bu.col,NORMLINE,bu.dash); 
			draw_line(xx-dxx,y+dy-mar,xx+dxx,y+dy-mar,bu.col,NORMLINE,bu.dash); 
			plot_equation(te,x+3,y+dy/2+0.4,1,dx-(x+3));
		}
		break;
		
	case "KeyCross":
		{
			clear_rectangle(x,y,dx,dy);
			let d = 0.3;
			let xp = x+1.7;
			let yp = y+dy/2+0.05;
			let thick = THICKLINE;
			draw_line(xp-d,yp-d,xp+d,yp+d,bu.col,thick);
			draw_line(xp-d,yp+d,xp+d,yp-d,bu.col,thick);
			plot_text(te,x+3,y+dy/2+0.4,get_font(0.8),BLACK); 
		}
		break;
	
	case "KeySymb":
		{
			let xx = x+1.7;
			let si = 0.5, si2 = 0.5;
			let siw = 0.5, sih = 0.7;
			clear_rectangle(x,y,dx,dy);
			switch(bu.symb){
			case "AddObs":
				{			
					let list=[]; list.push({col:WHITE, frac:1});
					multi_triangle(xx-0.3,y+0.6*dy-sih/2,siw,sih,list,BLACK,NORMLINE);
				}
				break;
			
			case "RemObs": rem_triangle(xx-0.3,y+0.6*dy-sih/2,siw,sih,BLACK); break;
				
			case "MoveObs": fill_semicircle(xx,y+0.55*dy,0.32,WHITE,BLACK,NORMLINE); break;
			
			case "TransObs":
				{
					draw_rectangle(xx-si2/2,y+0.60*dy-si2/2,si2,si2,BLACK,NORMLINE); 
					draw_line(xx,y+0.6*dy-si2/2,xx,y+0.6*dy+si2/2,BLACK,NORMLINE);
				}
				break;
				
			case "CompObs": circle(xx,y+dy/2,0.32,BLACK,NORMLINE); break;
			
			case "DiagObs": draw_rect(xx-si/2,y+0.6*dy-si/2,si,si,WHITE,BLACK,NORMLINE); break;
			
			case "GeneticObs": draw_X(x+1.7,y+dy/2,0.45*dy); break;
				
			}
			plot_text(te,x+3,y+dy/2+0.4,get_font(0.8),BLACK); 
		}
		break;
	
	case "AddFilter":
		{
			fill_rectangle(x,y,dx,dy,WHITE);
			let si = 0.7;
			let font = get_font(si,"bold");
			
			let col = BLACK; if(ov == 1) col = GREY;
			 
			draw_plus(x+0.5,y+dy/2,0.32,col,WHITE);
			
			plot_text(te,x+1, y+dy/2+0.45*si,font,col);
		}
		break;

	case "BannerWarn":
		{
			//fill_rectangle(x,y,dx,dy,BLUE); 
			let fo = get_font(0.9,"bold");
			let w = text_width(te,fo);
			let sh = 1;
			center_text(te,x+dx/2+sh,y+0.9,fo,BLACK,dx); 
			draw_image_sc(find_pic("warn"),x+dx/2-w/2-1.3+sh,y,1.2);
		}
		break;

	case "ModType":
		fill_rectangle(x,y,dx,dy,WHITE);
		if(ov) cv.globalAlpha = 0.5;
		draw_image(bu.mod_ty,x,y,dx,dy);
		if(ov) cv.globalAlpha = 1;
		break;

	case "BubblePic":
		draw_image(bu.pic,x,y,dx,dy);
		break;
		
	case "Nothing":
		break;
		
	case "ClassBack":
		break;
		
	default: error(bu.type+" Button not recognised"); break;
  }
}
