/// Layers describe canvases which are sequentially plotted to generate the final interface

"use strict";

class Layer
{
	x;
	y; 
	dx; 
	dy;
	but = [];
	name;
	can;
	cv;
	
	constructor(na,x,y,dx,dy,op)	
	{
		this.name = na;
		this.x = x;
		this.y = y;
		this.dx = dx;
		this.dy = dy;
		this.op = op;
		this.x_shift = 0;
		this.y_shift = 0;	
		this.index = inter.layer.length;

		this.x_scroll = false;
		this.y_scroll = false;
	}


	/// Adds buttons to the main area (this is to the right of the menu bar)
	add_main_buts()
	{
		this.add_button({x:0, y:0, dx:this.dx, dy:this.dy, type:"Rect", val:WHITE});
	
		let tree = inter.page_name.split("->");
		
		if(model.warn.length > 0){
			model.add_warning_buts(this);
			return;			
		}
		
		switch(tree[0]){
		case "Home": add_home_page_buts(this); break;
		case "Model":
			switch(tree[1]){
			case "Compartments":
				model.add_classback_buts(this);
				break;
			case "Parameters":
				if(model.warn.length > 0) model.add_warning_buts(this); 
				else add_model_param_buts(this); 
				break;
			case "Description": add_description_buts(this); break;
			}
			break;
		
		case "Simulation":
			switch(tree[1]){
			case "Population": add_data_buts(this,"sim"); break;
			
			case "Parameters":
				if(model.warn.length > 0) model.add_warning_buts(this); 
				else add_param_value_buts(this); 
				break;
				
			case "Run":
				if(model.warn.length > 0) model.add_warning_buts(this); 
				else{
					if(inter.running_status == false){
						add_sim_start_buts(this);
					}
				}					
				break;
					
			case "Results":
				switch(tree[2]){
				case "Populations":
					add_pop_buts(sim_result,this);
					break;
				case "Transitions":
					add_trans_buts(sim_result,this);
					break;
				case "Individuals":
					add_individual_buts(sim_result,this);
					break;
				case "Splines":
					add_spline_buts(sim_result,this);
					break;
				case "Parameters":
					add_parameter_buts(sim_result,this);
					break;
				}
				break;
				
			case "Generate Data": add_data_buts(this,"gen"); break;
			}
			break;
			
			
			
		case "Inference":
			switch(tree[1]){
			case "Data": add_data_buts(this,"inf"); break;
			case "Prior":
				if(model.warn.length > 0) model.add_warning_buts(this); 
				else add_param_prior_buts(this); 
				break;
			case "Run":
				if(model.warn.length > 0) model.add_warning_buts(this); 
				else add_inf_start_buts(this);	
				break;
	
			case "Results":
				switch(tree[2]){
				case "Populations":
					add_pop_buts(inf_result,this);
					break;
				case "Transitions":
					add_trans_buts(inf_result,this);
					break;
				case "Individuals":
					add_individual_buts(inf_result,this);
					break;
				case "Splines":
					add_spline_buts(inf_result,this);
					break;
				case "Parameters":
					add_parameter_buts(inf_result,this);
					break;
				}
				break;
			}
			break;
			
		default: error("Page "+inter.page_name+" not set up"); break;
		}
	}
		
		
	/// Adds a button onto a page
	add_button(bu)      
	{
		if(!bu.type) alertp("'type' must be specified");
		if(!bu.dx && bu.type == "Text") bu.dx = text_width(bu.te,bu.font)+0.3;
		if(!bu.dx) bu.dx = 0; 
		if(!bu.dy && bu.dy != 0) bu.dy = 0;
		bu.ind = this.but.length;
		this.but.push(bu);
	}


	/// Plots all the buttons on a layer
	plot_buttons(not_plot) 
	{	
		if(this.but.length == 0){ 
			if(this.can != undefined) this.can = undefined;
			return;
		}
		
		this.can = document.createElement('canvas');
    this.cv = this.can.getContext('2d');
		
		let wid = Math.round(this.inner_dx*inter.sca);
		let hei = Math.round(this.inner_dy*inter.sca);
		
		this.can.width = wid;
		this.can.height = hei;

		if(this.background != undefined){
			this.cv.beginPath();
			this.cv.rect(0,0,wid,hei);
			this.cv.fillStyle = this.background;
			this.cv.fill()
		}
		
		let x1 = this.x_shift, y1 = this.y_shift; 
		if(inter.export_image == true){ x1 = 0; y1 = 0;}
		
		let x2 = x1+this.inner_dx, y2 = y1+this.inner_dy;
		for(let i = 0; i < this.but.length; i++){
			let bu = this.but[i];
			
			if(bu.x < x2 && bu.x+bu.dx > x1 && bu.y < y2 && bu.y+bu.dy > y1){
				if(!(not_plot && find_in(not_plot,bu.type) != undefined)){				
					let dx_can = this.dx, dy_can = this.dy;
	
					let ov = false;
			
					if(this.index == inter.over.layer && i == inter.over.i) ov = true;
					this.plot_button(bu,ov);
				}
			}
		}
	}


	/// Adds a title to the layer
	add_title(te,cx,cy,op)
	{ 
		let col = DBLUE, col_line = LBLUE;
		if(op == undefined) op = {};
		if(op.col != undefined) col = op.col;
		if(op.col_line != undefined) col_line = op.col_line;
		
		let text_anno = text_convert_annotation(te,si_title,1.4*si_title,100,"",col);
		
		this.add_button({word:text_anno.word, x:cx, y:cy, dx:this.dx-cx*2, dy:1.5, col:col, col_line:col_line, type:"Title"}); 
		
		if(op.te){
			let title = te; if(op.title) title = op.title;
			this.add_help_button(cx+text_anno.width_first+0.2,cy+1.5,{title:title, te:op.te, back_col:op.back_col});
		}
		
		return cy + 2;
	}


	/// Adds a subtitle
	add_subtitle(te,cx,cy,back_col,op)
	{	
		let si = 0.8;
		let fo = get_font(si,"Bold");
	
		this.add_button({te:te, x:cx,dx:this.dx-cx*2, y:cy, dy:si, type:"Text", font:fo, si:si, col:DBLUE, back_col:back_col});
		
		if(op && op.te){
			let title = te; if(op.title != undefined) title = op.title;
			let w = text_width(te,fo);
			this.add_help_button(cx+w+0.1,cy+1.2,{title:title, te:op.te, back_col:op.back_col});
		}
		
		return cy+1.2;
	}
	
	
	/// Adds a link to the layer
	add_example_link(te,file,cx,cy)
	{ 
		let font = get_font(0.9);
		
		this.add_button({te:te, x:cx, y:cy, dx:text_width(te,font)+1, dy:1, ac:"Example", type:"Example", font:font, file:file});
		return cy + 1.4;
	}
	
	/// Adds some simple text
	add_text(te,cx,cy,col,si)
	{
		this.add_button({te:te, x:cx, y:cy, dy:si, type:"Text", font:get_font(si), si:si, col:BLACK});
		
		//this.add_paragraph(te,text_width(te,get_font(si)),cx,cy,col,si);
	}
	
	
	/// Adds some blank space
	add_space(cy,dy)
	{
		this.add_button({x:0, y:cy, dx:0, dy:dy, type:"Nothing"});
		return cy+dy;
	}
	
	
	/// Adds a paragraph to the layer
	add_paragraph(te,dx,cx,cy,col,si,lh,back_col,align)
	{
		if(lh == undefined) l = 1.4*si;
		
		let text_anno = text_convert_annotation(te,si,lh,dx,align,col);
	
		let dy = text_anno.height;
		
		if(back_col != undefined) this.add_button({x:cx, y:cy, dx:dx, dy:dy, type:"Rect", val:back_col});
		
		this.add_button({x:cx, y:cy, dx:dx, dy:dy, type:"Paragraph", word:text_anno.word}); 
		
		for(let i = 0; i < text_anno.link.length; i++){
			let li = text_anno.link[i];
			let spl = li.ac.split("|");
			
			this.add_button({te:li.te, x:cx+li.x, y:cy+li.y-1*li.si, dx:li.w, dy:1.4*li.si, si: li.si, ac:spl[0], info:spl[1], type:"Link", font:li.font});
		}
		
		return cy+dy;
	}
	
	
	/// Adds an example picture
	example_picture(te,x,y,dx,dy,pic,help_title,help_text)
	{
		this.add_button({x:x, y:y, dx:dx, dy:dy, type:"ExampleModel", pic:pic}); 
		let si = 0.9;
		let font = get_font(si,"bold");
		let w = text_width(te,font)+0.3;
		this.add_button({te:te, x:x+dx/2-w/2, y:y-si, dx:w, dy:si, type:"CenterText", font:font}); 
		
		this.add_help_button(x+dx/2+w/2,y+0.2,help_title,help_text);
	}
	
	
	/// Adds a help icon
	add_help_button(x,y,help)
	{
		this.add_button({x:x, y:y-1.2, dx:0.75,  dy:0.6, ac:"HelpIcon", type:"HelpIcon", help:help}); 
	}
	
	
	/// Adds a check box on a black background
	add_checkbox_white(x,y,value,te,source)
	{
		this.add_button({value:value, x:x+0.5, y:y, dx:1, dy:1, ac:"CheckboxButton", type:"CheckboxWhiteButton", source:source});
		this.add_button({te:te, x:x+1.6, y:y, dx:text_width(te,get_font(si_radio))+0.4, dy:1, type:"CheckboxButtonText", source:source, col:WHITE});
	}

	/// Adds a check box on a white background
	add_checkbox(x,y,value,te,source,back_col,op)
	{
		this.add_button({value:value, x:x+0.5, y:y, dx:1, dy:1, ac:"CheckboxButton", type:"CheckboxButton", source:source});
		let dx = text_width(te,get_font(si_radio))+0.4;
		this.add_button({te:te, x:x+1.6, y:y, dx:dx, dy:1, type:"CheckboxButtonText", source:source, col:BLACK, back_col:back_col});
		
		if(op != undefined && op.title != undefined){
			let w = text_width(te,get_font(si_title));
			this.add_help_button(x+dx+1.6,y+1.3,{title:op.title, te:op.te, back_col:back_col});
		}
	}


	/// Adds a radio button on a black background
	add_radio_white(x,y,value,te,source,op)
	{
		let back_col; 
		let fo = get_font(si_radio);
		if(op){
			back_col = op.back_col;
			if(op.fo != undefined) fo = op.fo;
		}
			
		this.add_button({value:value, x:x+0.5, y:y, dx:1, dy:1, ac:"RadioButton", type:"RadioWhiteButton", source:source, back_col:back_col});
		this.add_button({te:te, x:x+1.6, y:y, dx:text_width(te,fo)+0.4, dy:1, type:"RadioButtonText", source:source, fo:fo, col:WHITE, back_col:back_col});
	}
	
	
	/// Adds a radio button on white background
	add_radio(x,y,value,te,source,op)
	{
		let back_col; 
		let fo = get_font(si_radio);
		let ac = "RadioButton"; 
		if(op){
			back_col = op.back_col;
			if(op.fo != undefined) fo = op.fo;
			if(op.disable == true) ac = undefined;
		}
	
		this.add_button({value:value, x:x+0.5, y:y, dx:1, dy:1, ac:ac, type:"RadioButton", source:source, back_col:back_col});
		this.add_button({te:te, x:x+1.6, y:y, dx:text_width(te,get_font(si_radio))+0.4, dy:1, type:"RadioButtonText", source:source, col:BLACK, fo:fo, back_col:back_col});
		
		return x + text_width(te,get_font(si_radio))+3;
	}

	/// Adds yscroll buttons
	add_yscroll_buts()
	{
		let lay = this.op.lay;
		let style = this.op.style;

		let box = get_but_box(lay.but);
		let height = box.ymax;
		
		let frac = lay.inner_dy/height;
		let frac_trunc = frac; if(frac_trunc < 0.05) frac_trunc = 0.05;

		let dx = this.dx;
		let si = dx;
		if(lay.name == "DropdownOptions") si = 0;
	
		let pn = inter.page_name+"_"+lay.name+"_Y";
		
		let scroll_ref = 0; while(scroll_ref < inter.scroll_position.length && pn != inter.scroll_position[scroll_ref].page) scroll_ref++;
		
		if(scroll_ref == inter.scroll_position.length){
			inter.scroll_position.push({page:pn, shift:0, max_create:height}); 
		}
		
		this.scroll_ref = scroll_ref;
		
		let dy_tot = lay.inner_dy-2*si;
		
		let pos = inter.scroll_position[scroll_ref];
			
		pos.max = height;
		pos.scale = height/dy_tot;
		
		pos.frac = frac;
		pos.scroll_layer = this.index;
		
		//if(height != pos.max_create) pos.shift = 0;
		if(pos.shift < 0) pos.shift = 0;
		if(pos.shift > pos.max*(1-pos.frac)) pos.shift = pos.max*(1-pos.frac);

		let ysh = pos.shift;
		let ysh_disc = ysh;

		if(lay.name == "TextBox"){
			ysh_disc = Math.round(ysh_disc/textbox_linesi)*textbox_linesi;
			pos.discretise = true;
		}
		
		if(lay.name == "DropdownOptions"){
			ysh_disc = Math.round(ysh_disc/dropdown_opheight)*dropdown_opheight;
			pos.discretise = true;
		}
	
		let fr = ysh_disc/(height*(1-frac));
		lay.y_shift = ysh_disc;
	
		let y0 = si;
		let y1 = si+fr*(1-frac_trunc)*dy_tot;
		let y2 = y1+frac_trunc*dy_tot;
		let y3 = si+dy_tot;

		if(si != 0) this.add_button({x:0, y:0, dx:dx, dy:si, ac:"ScrollUp", type:"YscrollUp", scroll_ref:scroll_ref, frac:frac, style:style});
	
		this.add_button({x:0, y:y0, dx:dx, dy:y1-y0, ac:"ScrollUp", type:"BackSlider", scroll_ref:scroll_ref, frac:frac, style:style});
		
		this.add_button({x:0, y:y1, dx:dx, dy:y2-y1, ac:"VerticalSlider", type:"VerticalSlider", scroll_ref:scroll_ref, scale:dy_tot*(1-frac_trunc), pn:pn, style:style});
	
		this.add_button({x:0, y:y2, dx:dx, dy:y3-y2, ac:"ScrollDown", type:"BackSlider", scroll_ref:scroll_ref, frac:frac, style:style});
		
		if(si != 0) this.add_button({x:0, y:si+dy_tot, dx:dx, dy:si, ac:"ScrollDown", type:"YscrollDown", scroll_ref:scroll_ref, frac:frac, style:style});
	}


	/// Adds xscroll buttons
	add_xscroll_buts()
	{
		let lay = this.op.lay;
		let style = this.op.style;

		let box = get_but_box(lay.but);
		let width = box.xmax;
		
		let frac = lay.inner_dx/width;
		let frac_trunc = frac; if(frac_trunc < 0.05) frac_trunc = 0.05;

		let si = this.dy;
	
		let pn = inter.page_name+"_"+lay.name+"_X";
		
		let scroll_ref = 0; while(scroll_ref < inter.scroll_position.length && pn != inter.scroll_position[scroll_ref].page) scroll_ref++;
		
		if(scroll_ref == inter.scroll_position.length){
			inter.scroll_position.push({page:pn, shift:0, max_create:width}); 
		}
		this.scroll_ref = scroll_ref;
	
		let dx_tot = lay.inner_dx-2*si;
		
		let pos = inter.scroll_position[scroll_ref];
		pos.max = width;
		pos.scale = width/dx_tot;
		pos.frac = frac;
		pos.scroll_layer = this.index;
		
		//if(width != pos.max_create) pos.shift = 0;
		if(pos.shift < 0) pos.shift = 0;
		if(pos.shift > pos.max*(1-pos.frac)) pos.shift = pos.max*(1-pos.frac);

		let xsh = pos.shift;

		let fr = xsh/(width*(1-frac));
		lay.x_shift = xsh;
	 
		let x0 = si;
		let x1 = si+fr*(1-frac_trunc)*dx_tot;
		let x2 = x1+frac_trunc*dx_tot;
		let x3 = si+dx_tot;
	
		this.add_button({x:0, y:0, dx:si, dy:si, ac:"ScrollUp", type:"XscrollLeft", scroll_ref:scroll_ref, frac:frac, style:style});
	
		this.add_button({x:x0, y:0, dx:x1-x0, dy:si, ac:"ScrollUp", type:"BackSlider", scroll_ref:scroll_ref, frac:frac, style:style});
		
		this.add_button({x:x1, y:0, dx:x2-x1, dy:si, ac:"HorizontalSlider", type:"HorizontalSlider", scroll_ref:scroll_ref, scale:dx_tot*(1-frac_trunc), pn:pn, style:style});
	
		this.add_button({x:x2, y:0, dx:x3-x2, dy:si, ac:"ScrollDown", type:"BackSlider", scroll_ref:scroll_ref, frac:frac, style:style});
		
		this.add_button({x:si+dx_tot, y:0, dx:si, dy:si, ac:"ScrollDown", type:"XscrollRight", scroll_ref:scroll_ref, frac:frac, style:style});
	}


	/// Adds corner scroll buttons
	add_corner_scroll_buts()
	{
		let si = this.dx;
		//this.add_button({x:0, y:0, dx:si, dy:si, type:"Rect", val:LLGREY, style:this.op.style});
	}
	

	/// Adds loading symbol buttons
	add_loadingsymbol_buts()
	{
		this.add_button({x:this.dx/2-loading_si/2, y:0, dx:loading_si, dy:loading_si, type:"LoadingSymbol"});
		
		if(inter.loading_symbol.processing == true) this.add_button({te:"Processing", x:0, y:loading_si, dx:this.dx, dy:1, type:"Stop"});
		else this.add_button({te:"Stop", x:0, y:loading_si, dx:this.dx, dy:1, type:"Stop", ac:"Stop"});
	}
	
	
	/// Adds buttons for input box
	add_input_buts()
	{
		//this.add_button({x:0, y:0, dx:this.dx, dy:this.dy, type:"Rect", val:BLUE});				
		
		this.copy_from_source();
		
		let te = this.get_text_from_source();

		if(te == undefined) return;

		let marx = input_margin;
	
		let font;
		if(this.op.source.eqn == true) font = get_font(inputbox_fontsi,"","times");
		else font = get_font(inputbox_fontsi);
		
		let box_dx = this.inner_dx;
		
		let cx = marx;
		let cx_click = 0;
	
		let lh = inputbox_linesi;
		
		let eqn = this.op.source.eqn;
	
		let textcol = BLACK;
		let textcol_param_store = BLACK;
		let textcol_pop_store = BLACK;
		let textcol_ie_store = BLACK;
		let textcol_fe_store = BLACK;

		for(let i = 0; i < te.length; i++){
			let ch = te.substr(i,1);
			let w = text_width(ch,font)+textbox_space;	

			this.add_button({x:cx_click, y:0, dx:(cx+w/2)-cx_click, dy:lh, ac:"PositionCursor", type:"Nothing", i:i, eqn:eqn});
			cx_click = (cx+w/2);
			
			if(ch == "{"){ textcol_pop_store = textcol; textcol = BLUE;}
			if(ch == "["){ textcol_ie_store = textcol;  textcol = DGREEN;}
			if(ch == "〈"){ textcol_fe_store = textcol; textcol = RED;}
			this.add_button({te:ch, x:cx, y:0, dx:w, dy:lh, type:"Char", col:textcol, i:i, font:font, si:inputbox_fontsi});
			if(ch == "}") textcol = textcol_pop_store;
			if(ch == "]") textcol = textcol_ie_store;
			if(ch == "〉") textcol = textcol_fe_store;
					
			cx += w;
		}
		
		let i = te.length;
		this.add_button({te:"", x:cx, y:0, dx:0, dy:lh, type:"Char", col:BLACK, i:i, font:font, si:inputbox_fontsi});				
			
		this.add_button({x:cx_click, y:0, dx:box_dx-cx_click, dy:lh, ac:"PositionCursor", type:"Nothing", i:i, eqn:eqn});
						
		this.outline = {style:"round", col:INPUT_OUTLINE_COL, line:NORMLINE, marx:outline_marx, mary:outline_mary}
	}
	
	add_expand_buts()
	{
		this.add_button({x:0, y:0, dx:expand_dx, dy:inputbox_linesi, type:"ExpandIcon", ac:"ExpandIcon"});				
	}

	/// Adds buttons for textbox
	add_textbox_buts()
	{
		this.copy_from_source();
		
		let te = this.get_text_from_source();
		if(te == undefined) return;
		
		let textcol = BLACK;
		
		let lh = textbox_linesi;
		let marx = 0.1;
		
		let cx = marx;
		let cy = 0;
		let cx_click = 0;
	
		let font;
		if(this.op.source.type == "equation") font = get_font(textbox_fontsi,"","times");
		else font = get_font(textbox_fontsi);
			
		let box_dx = this.inner_dx;
		let box_dy = this.dy;
		
		let wrap = this.op.source.wrap;
	
		let warn = this.op.source.warn;
		
		let textcol_param_store = BLACK;
		let textcol_pop_store = BLACK;
		let textcol_ie_store = BLACK;
		let textcol_fe_store = BLACK;

		let text_nrow = 1;

		let list=[];
		for(let i = 0; i < te.length; i++){
			let ch = te.substr(i,1);
			let check_word_wrap = false;
			if(wrap == "word" && ch != " " && (i == 0 ||  list[i-1].ch == " " || list[i-1].ch == "\n")) check_word_wrap = true;
			list[i] = {ch:ch, w:text_width(ch,font)+textbox_space, check_word_wrap:check_word_wrap};
		}
	
		let eqn = this.op.source.eqn;
		for(let i = 0; i <= te.length; i++){
			let char_begin_cx = cx, char_begin_cy = cy;
			
			if(i < te.length){
				let li = list[i];
				let ch = li.ch;
				
				if(ch == "\n"){
					this.add_button({x:cx_click, y:cy, dx:box_dx-cx_click, dy:lh, ac:"PositionCursor", type:"Nothing", i:i, eqn:eqn});
					cx_click = 0;

					if(cx == marx){
						let w = text_width(" ",font);
						this.add_button({te:" ", x:cx, y:cy, dx:w, dy:lh, type:"Char", col:BLACK, i:i, font:font, si:textbox_fontsi });
					}
					else{
						this.add_button({te:"", x:cx, y:cy, dx:0, dy:lh, type:"Char", col:BLACK, i:i, font:font, si:textbox_fontsi});				
					}
				
					cx = marx; cy += lh; text_nrow++;
				}
				else{
					if(li.check_word_wrap == true && cx != marx){  // This checks to make sure that word does no tneed to wrap
						let xx = cx; let ii = i; 
						while(ii < te.length && xx < box_dx-marx && list[ii].ch != " " && list[ii].ch != "\n"){
							xx += list[ii].w;
							ii++;
						}
						if(xx >= box_dx-marx){
							cx = marx; cy += lh; text_nrow++;
							char_begin_cx = cx; char_begin_cy = cy;

							cx_click = 0;
						}
					}					
					
					let w = li.w;
					if(cx+w >= box_dx-marx){ 
						if(wrap == "word") cx = marx;
						else{
							cx = textbox_tab; 
							this.add_button({x:cx, y:cy, dx:textbox_tab-2*marx, dy:lh, type:"TextboxTab", col:BLACK, i:i, font:font, si:textbox_fontsi });
						}
			
						cy += lh; text_nrow++;
						
						char_begin_cx = cx; 
						char_begin_cy = cy;

						cx_click = 0;
					}

					this.add_button({x:cx_click, y:cy, dx:(cx+w/2)-cx_click, dy:lh, ac:"PositionCursor", type:"Nothing", i:i, eqn:eqn});
					cx_click = (cx+w/2);

					if(ch == "{"){ textcol_pop_store = textcol; textcol = BLUE;}
					if(ch == "["){ textcol_ie_store = textcol; textcol = DGREEN;}
					if(ch == "〈"){ textcol_fe_store = textcol; textcol = RED;}
					
					let underline = false;
					if(warn != undefined){
						let len = warn.len; if(len == 0) len = 1;
						if(i >= warn.cur && i < warn.cur+len) underline = true;
					}
					
					this.add_button({te:ch, x:cx, y:cy, dx:w, dy:lh, type:"Char", col:textcol, i:i, font:font, si:textbox_fontsi, underline:underline});
					
					if(ch == "}") textcol = textcol_pop_store;
					if(ch == "]") textcol = textcol_ie_store;
					if(ch == "〉") textcol = textcol_fe_store;
					cx += w;
				}
			}

			if(i == te.length){
				this.add_button({te:"", x:cx, y:cy, dx:0, dy:lh, type:"Char", col:BLACK, i:i, font:font, si:textbox_fontsi});				
			
				this.add_button({x:cx_click, y:cy, dx:box_dx-cx_click, dy:lh, ac:"PositionCursor", type:"Nothing", i:i, eqn:eqn});
				cx = marx; cy += lh; 
			
				this.add_button({x:0, y:cy, dx:box_dx, dy:box_dy-cy, ac:"PositionCursor", type:"Nothing", i:i, eqn:eqn});
			}
		}
		

		this.text_nrow = text_nrow;
		this.outline = {style:"round", col:BLACK, line:NORMLINE, marx:outline_marx, mary:outline_mary}
	}


	/// Adds textbox
	add_textbox(x,y,dx,nrow,source)	
	{		
		add_ref(source);
		add_layer("TextBox",this.x+x,this.y+y,dx,nrow*textbox_linesi,{nrow:nrow, source:source});
	}


	/// Adds input box
	add_input(x,y,dx,source)	
	{		
		add_ref(source);

		if(source.eqn == true){
			add_layer("Input",this.x+x,this.y+y,dx-expand_dx,inputbox_linesi,{source:source});
			add_layer("Expand",this.x+x+dx-expand_dx,this.y+y,expand_dx,inputbox_linesi);
		}
		else{
			add_layer("Input",this.x+x,this.y+y,dx,inputbox_linesi,{source:source});
		}
	}
	
	
	/// Adds dropdown menu
	add_dropdown(x,y,dx,dymax,source,pos,style)	
	{		
		add_layer("Dropdown",this.x+x,this.y+y,dx,dropdown_height,{source:source, pos:pos, style:style});
	}

	
	/// Adds dropdown buttons
	add_dropdown_buts()
	{
		let op = this.op;

		let ac; if(op.pos != undefined) ac = "Dropdown";
		this.add_button({te:op.source.te, x:0, y:0, dx:this.dx, dy:this.dy, ac:ac, type:"Dropdown", source:op.source, direction:op.direction, style:op.style});
	}

	
	/// Adds dropdown options
	add_dropdownoptions_buts()
	{
		let pos = this.op.pos;
		for(let i = 0; i < pos.length; i++){
			this.add_button({pos:pos[i], x:0, y:i*dropdown_opheight, dx:this.inner_dx, dy:dropdown_opheight, ac:"DropdownOption", type:"DropdownOption", source:this.op.source});
		}
	}


	/// Finds the button with a given action
	search_button_ac(ac)
	{
		let i = 0; while(i < this.but.length && this.but[i].ac != ac) i++;
		if(i == this.but.length) error("Cannot find button"+ac);
		return i;
	}


	/// Gets the text store for this layer 
	get_text_box_store()
	{
		let ref = this.op.source.ref;

		let sto = inter.textbox_store;
		let i = 0; while(i < sto.length && sto[i].ref != ref) i++;
		if(i == sto.length){ error("PROBLEM getting store"); return "";}
		
		return sto[i];
	}


	/// Gets the text associated with this layer
	get_text_from_source()
	{
		return this.get_text_box_store().te;
	}


	/// Puts a given set of text into the store
	put_text_in_source(st,remember)
	{
		let tbs = this.get_text_box_store();
		
		if(tbs.te != st){
			if(remember.on == true){
				let cur = inter.cursor;
				tbs.store_text.push({text:tbs.te, i:remember.i, select_pos:remember.select_pos});			
			}			
			
			tbs.te = st;
		
			if(tbs.source.update == true){ // If update is true then copies back to source
				copy_back_to_source2(tbs);
			}
		}
	}


	/// Displays a parameter
	display_param(x,y,info)
	{
		this.add_button({x:x, y:y, dx:info.dx, dy:1.2*si_big, type:"ParamLabel", info:info, col:BLACK});
	}
	
	/// Add corner buttons
	add_corner_button(arr,pos,op)
	{
		let gap = 0.8;
		let cx = pos.x-arr.length*but_width- (arr.length-1)*gap-0.3, cy = pos.y-but_height-0.3;

		for(let i = 0; i < arr.length; i++){
			let row = arr[i]; if(row.length != 3) alert("Problem with button");
			
			this.add_button({te:row[0], x:cx, y:cy, dx:but_width, dy:but_height, op:op, ac:row[2], type:row[1]});

			cx += but_width+gap;
		}
	}


	/// Adds a gap with nothing in it
	add_gap(cy,dy)
	{
		this.add_button({x:0, y:cy, dx:0, dy:dy, type:"Nothing"});
		return cy+dy;
	}
	
	
	/// Deals with when the up and down keys are pressed within a text box
	arrow_up_down(i,dir)
	{
		let x, y;
		
		for(let j = 0; j < this.but.length; j++){  
			let bu = this.but[j];
			if(bu.type == "Char" && bu.i == i){ x = bu.x; y = bu.y;}
		}

		let cur = inter.cursor;
		if(cur.xstore != undefined) x = cur.xstore;
		else cur.xstore = x;
		
		let inew, dist = LARGE;
		for(let j = 0; j < this.but.length; j++){
			let bu = this.but[j];
			if((dir == "up" && bu.y < y) || (dir == "down" && bu.y > y) ){
				let d = (bu.y-y)*(bu.y-y)*1000+(bu.x-x)*(bu.x-x);
				if(d < dist){ dist = d; inew = bu.i;}
			}
		}
		
		if(inew == undefined) return i;

		return inew;
	}


	/// Ensures that if the cursor goes out of view the screen is repositioned
	ensure_cursor_in_view()
	{
		let cur = inter.cursor;
		
		if(cur.i == undefined) return;
		
		switch(this.name){
		case "TextBox":
			{
				if(this.y_scroll == false) return;
				
				let y;
				
				for(let j = 0; j < this.but.length; j++){  
					let bu = this.but[j];
					if(bu.type == "Char" && bu.i == cur.i) y = bu.y;
				}
				
				if(y == undefined){
					let te = inter.layer[cur.l].get_text_from_source();
				}
					
				let yscroll_lay = inter.layer[this.index+1];
				
				if(y-this.y_shift < 0.5){
					change_scroll(-1,yscroll_lay.but[0]);
					yscroll_lay.initialise_buttons();
					yscroll_lay.plot_buttons();
					this.plot_buttons();
				}

				if(y-this.y_shift > this.inner_dy-1.5){
					change_scroll(1,yscroll_lay.but[0]);
					
					yscroll_lay.initialise_buttons();
					yscroll_lay.plot_buttons();
					this.plot_buttons();
				}
			}
			break;
			
		case "Input":
			{
				if(this.x_scroll == false) return;
				
				let xscroll_lay = inter.layer[this.index+1];
				
				let x;
				
				for(let j = 0; j < this.but.length; j++){  
					let bu = this.but[j];
					if(bu.type == "Char" && bu.i == cur.i) x = bu.x;
				}
				if(x == undefined){
					let te = inter.layer[cur.l].get_text_from_source();
				}				
				
				let sh = 0;
				if(x-this.x_shift > this.dx-input_margin){
					sh = (x-this.x_shift) - (this.dx-input_margin);
				}
				
				if(x-this.x_shift < input_margin){
					sh = (x-this.x_shift) - input_margin;
				}
				
				if(sh != 0){
					let bu = xscroll_lay.but[0];
					
					let pos = inter.scroll_position[bu.scroll_ref];
		
					this.x_shift += sh;
					pos.shift += sh;
					
					xscroll_lay.initialise_buttons();
					xscroll_lay.plot_buttons();
					this.plot_buttons();
					plot_screen();
				}	
			}			
			break;

		default: error("Option not recognised 65"); break;
		}
	}	


	/// Sets up all the buttons associated with a layer
	initialise_buttons() 
	{
		this.but=[];
	
		switch(this.name){
		case "Screen": add_screen_buts(this); break;
		case "Logo": add_logo_buts(this); break;
		case "Menu": add_menu_buts(this); break;
		case "Main": this.add_main_buts(); break;
		case "EditTable": create_edit_table(this); break;	
		case "EditInitPop": create_initial_pop(this); break;	
		case "EditInitPopPrior": create_initial_pop_prior(this); break;	
		case "EditParam": create_edit_param(this); break;	
		case "EditAmatrix": create_edit_Amatrix(this); break;	
		case "EditXvector": create_edit_Xvector(this); break;	
		case "ViewGraph": create_view_graph(this); break;	
		case "Compartment": model.add_compartment_buts(this); break;
		case "Select": add_select_buts(this); break;
		case "Transition": model.add_transition_buts(this); break;
		case "Annotation": add_annotation_buts(this); break;
		case "AnnotationMap": add_annotation_map_buts(this); break;
		case "BubbleBack": add_bubble_back_buts(this); break;
		case "Bubble": add_bubble_buts(this); break;
		case "LowerMenu": model.add_lower_menu_buts(this); break;
		case "UpperMenu": model.add_upper_menu_buts(this); break;
		case "AddCompartment": model.add_addcompartment_buts(this); break;
		case "AddLabel": model.add_addlabel_buts(this); break;
		case "Examples": add_examples_buts(this); break;
		case "Frame": add_frame_buts(this); break;
		case "Xscroll": this.add_xscroll_buts(); break;
		case "Yscroll": this.add_yscroll_buts(); break;
		case "CornerScroll": this.add_corner_scroll_buts(); break;
		case "HelpBackground": add_help_background_buts(this); break;
		case "Help": add_help_buts(this); break;
		case "HelpContent": add_help_content_buts(this); break;
		case "EquationBackground": add_equation_background_buts(this); break;
		case "Equation": add_equation_buts(this); break;
		case "EquationAddQuantity": add_quantity_content_buts(this); break
		case "TextBox": this.add_textbox_buts(); break;
		case "Input": this.add_input_buts(); break;
		case "Expand": this.add_expand_buts(); break;
		case "DescriptionContent": add_description_content(this); break;
		case "InitialPopulationContent": add_initial_population_content(this); break;
		case "InitPopPriorContent": add_init_pop_prior_content(this); break;
		case "TableContent": add_table_content(this,this.op.table,WHITE,1.7); break;
		case "WarningContent": model.add_warning_content(this); break;
		case "ParamValueContent": add_param_value_content(this); break;
		case "ParamPriorContent": add_param_prior_content(this); break;
		case "ModelParamContent": add_model_param_content(this); break;
		case "Dropdown": this.add_dropdown_buts(); break;
		case "DropdownOptions": this.add_dropdownoptions_buts(); break;
		case "LoadingSymbol": this.add_loadingsymbol_buts(); break;
		case "CreateEditTableContent": add_create_edit_table_buts(this); break;
		case "CreateEditXvectorContent": add_Xvector_buts(this); break;
		case "CreateEditParamContent": add_create_edit_param_buts(this); break;
		case "CreateEditAmatrixContent": add_Amatrix_buts(this); break;
		case "BubbleScrollable": add_bubble_scrollable_buts(this); break;
		case "GraphContent": inter.graph.content(this); break;
		case "Axes": inter.graph.axes(this); break;
		case "TimeAxis": inter.graph.time_axis(this); break;
		case "RightMenu": right_menu_buts(this); break;
		case "RightBotMenu": rightbot_menu_buts(this); break;
		case "RightMidMenu": rightmid_menu_buts(this); break;
		case "GraphCompartments": inter.graph.add_compartment_buts(this); break;
		case "GraphTransitions": inter.graph.add_transition_buts(this); break;
		case "AnimControls": inter.graph.add_animcontrol_buts(this); break;
		
		default: error("Layer "+this.name+" not recognised"); break;
		}
	}
	
	
	/// Copied stored text from the source
	copy_from_source()
	{
		let so = this.op.source;
	
		for(let i = 0; i < inter.textbox_store.length; i++){
			if(inter.textbox_store[i].ref == so.ref){
				inter.textbox_store[i].source.hidden = so.hidden;
				inter.textbox_store[i].used = true;
				return;
			}
		}
	
		let te = "";
		let eqn;
		switch(so.type){ 
		case "equation": te = inter.equation.te; break;
		case "description": te = model.description.te; break;
		case "compartment": te = model.species[so.p].cla[so.cl].comp[so.i].name; break;
		
		case "trans_mean":
			eqn = model.species[so.p].cla[so.cl].tra[so.i].value.mean_eqn; 
			break;
		
		case "trans_rate":
			eqn = model.species[so.p].cla[so.cl].tra[so.i].value.rate_eqn; 
			break;
		
		case "trans_shape": 
			eqn = model.species[so.p].cla[so.cl].tra[so.i].value.shape_eqn; 
			break;
		
		case "trans_scale":
			eqn = model.species[so.p].cla[so.cl].tra[so.i].value.scale_eqn; 
			break;
		
		case "trans_cv": 
			eqn = model.species[so.p].cla[so.cl].tra[so.i].value.cv_eqn; 
			break;
		
		case "trans_bp": 
			eqn = model.species[so.p].cla[so.cl].tra[so.i].value.bp_eqn; 
			break;
		
		case "trans_shape_erlang":
			let value = model.species[so.p].cla[so.cl].tra[so.i].value;
			te = value.shape_erlang.te; 				
			break;
			
		case "add_species_name": te = ""; break;
		case "add_classification_name": te = ""; break;
		case "classification_name": te = model.species[so.p].cla[so.cl].name; break;
		case "species_name": te = model.species[so.p].name; break;
		case "init_population": te = String(inter.edit_source.cla[so.cl].comp_init_pop[so.c].pop); break;
		case "init_per": te = String(inter.edit_source.cla[so.cl].comp_init_pop[so.c].pop_per); break;
		case "init_globpopulation": te = String(inter.edit_source.glob_comp[so.c].pop);	break;
		case "element": te = String(inter.edit_source.table.ele[so.r][so.c]);	break;
		case "element_param": te = String(get_element(inter.edit_param.value,so.pindex)); break;
		case "element_param_const": te = String(get_element(inter.edit_param.value,so.pindex)); break;
		case "element_Amatrix": te = String(inter.edit_Amatrix.value[so.j][so.i]); break;
		case "element_Xvector": te = String(inter.edit_Xvector.value[so.i]); break;
		case "element_eqn":
			{
				let st = String(get_element(inter.edit_param.value,so.pindex));
				eqn = create_equation(st,"reparam");
			}
			break;
		case "reparam_eqn": 
			{
				let st = String(model.param[inter.bubble.th].value);
				eqn = create_equation(st,"reparam"); 
			}
			break;
		case "find": te = String(inter.bubble.find); break;
		case "burnin": te = String(inter.bubble.burnin); break;
		case "replace": te = String(inter.bubble.replace); break;
		case "time_start": te = inter.edit_source.spec.time_start; break;
		case "time_end": te = inter.edit_source.spec.time_end; break;
		case "Se": eqn = inter.edit_source.spec.Se_eqn; break;
		case "Sp": eqn = inter.edit_source.spec.Sp_eqn; break;
		case "trap_prob": eqn = inter.edit_source.spec.trap_prob_eqn; break;
		case "pos_result": te = inter.edit_source.spec.pos_result; break;
		case "neg_result": te = inter.edit_source.spec.neg_result; break;
		case "percent": te = inter.edit_source.spec.percent; break;
		case "sd": te = inter.edit_source.spec.sd; break;
		case "param_val": te = model.param[so.val].value; break; 
		case "label": te = inter.bubble.label.te; break;
		case "label_anno": te = model.species[so.p].cla[so.cl].annotation[so.i].te; break;
		case "prior_min": te = inter.bubble.prior.value.min_eqn.te; break;
		case "prior_max": te = inter.bubble.prior.value.max_eqn.te; break;
		case "prior_mean": te = inter.bubble.prior.value.mean_eqn.te; break;
		case "prior_shape": te = inter.bubble.prior.value.shape_eqn.te; break;
		case "prior_sd": te = inter.bubble.prior.value.sd_eqn.te; break;
		case "prior_cv": te = inter.bubble.prior.value.cv_eqn.te; break;
		case "prior_alpha": te = inter.bubble.prior.value.alpha_eqn.te; break;
		case "prior_beta": te = inter.bubble.prior.value.beta_eqn.te; break;
		case "prior_dist_min": eqn = inter.bubble.prior.value.min_eqn; break;
		case "prior_dist_max": eqn = inter.bubble.prior.value.max_eqn; break;
		case "prior_dist_mean": eqn = inter.bubble.prior.value.mean_eqn; break;
		case "prior_dist_shape": eqn = inter.bubble.prior.value.shape_eqn; break;
		case "prior_dist_sd": eqn = inter.bubble.prior.value.sd_eqn; break;
		case "prior_dist_cv": eqn = inter.bubble.prior.value.cv_eqn; break;
		case "prior_dist_alpha": eqn = inter.bubble.prior.value.alpha_eqn; break;
		case "prior_dist_beta": eqn = inter.bubble.prior.value.beta_eqn; break;
		case "derive_eqn1": eqn = inter.bubble.derived.eqn1; break;
		case "derive_eqn2": eqn = inter.bubble.derived.eqn2; break;
		case "derive_eqn": eqn = model.derive[so.val].eqn2; break;
		case "deriveparam_eqn": eqn = model.derive[so.val].eqn1; break;
		case "sim_t_start": te = String(model.sim_details.t_start); break;
		case "sim_t_end": te = String(model.sim_details.t_end); break;
		case "sim_number": te = String(model.sim_details.number); break;
		case "sim_timestep": te = String(model.sim_details.timestep); break;
		
		case "inf_t_start": te = String(model.inf_details.t_start); break;
		case "inf_t_end": te = String(model.inf_details.t_end); break;
		case "inf_timestep": te = String(model.inf_details.timestep); break;
		case "inf_sample": te = String(model.inf_details.sample); break;
		case "inf_chain": te = String(model.inf_details.nchain); break;
		case "inf_abcsample": te = String(model.inf_details.abcsample); break;
		case "inf_thinparam": te = String(model.inf_details.thinparam); break;
		case "inf_thinstate": te = String(model.inf_details.thinstate); break;
		case "inf_accfrac": te = String(model.inf_details.accfrac); break;
		case "inf_accfracsmc": te = String(model.inf_details.accfracsmc); break;
		case "inf_numgen": te = String(model.inf_details.numgen); break;
		case "inf_kernelsize": te = String(model.inf_details.kernelsize); break;
		case "inf_indmax": te = String(model.inf_details.indmax); break;
		case "sim_indmax": te = String(model.sim_details.indmax); break;

		case "knot_times": te = stringify(model.param[so.val].spline.knot); break;
		case "smooth_value": te = String(model.param[so.i].spline.smooth.value); break;
		case "min": te = String(inter.bubble.min); break;
		case "max": te = String(inter.bubble.max); break;
		case "fixed_time": te = inter.bubble.fixed_time; break;
		case "const_column": te = inter.bubble.const_column; break;
		case "time_gen": te = String(inter.edit_source.time_gen); break;
		case "time_step": te = String(inter.edit_source.time_step); break;
		default: error("SOURCE PROBLEM: "+so.type); break;
		}

		if(eqn != undefined){
			te = eqn.te; 
			eqn = copy(eqn);
		}
		
		if(te == set_str) te = "";
	
		if(te == undefined) error("Problem getting from source: "+so.type);
	
		inter.textbox_store.push({ref:so.ref, source:so, te:te, eqn:eqn, store_text:[], used:true});
	}
}
