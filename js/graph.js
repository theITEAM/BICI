"use strict";

class Graph
{
	variety;                         // The thing being plotted (e.g. "Population", "Transition")
	type;                            // The type of graph (e.g. "Graph", "Compartment")
	data = [];                       // Data for graph
	view;                            // Shows how the data is viewed (graph or compartments)
	op;                              // Options for the graph
	range;                           // Range for axes 
	range_default;                   // The default range
	animation = {speed:{te:"×1", value:1}};   // Stores animation information
	matrix_key;                               // Information about a matrix key
	tick = {value:[], wmax:0, x:[], y:[]};    // Potential tick marks
	table_width = [];                // Table width for statistics tables
	warn;                            // Any warning to show why the graph could not be plotted
	
	constructor()
	{
		this.initialise_ticks();
	}
	
	initialise_ticks()            // Intialises tick marks used when drawing graphs
	{
		let tickpo = [1,2,5];
		
		for(let sh = -10; sh <= 10; sh++){
			for(let i = 0; i < 3; i++){   
				let f = tickpo[i];
				if(sh < 0){ for(let j = 0; j < -sh; j++) f /= 10;}
				else{ for(let j = 0; j < sh; j++) f *= 10;}
				this.tick.value.push(f);
			}
		}
	}

	
	/// Defines the type and data
	define(variety,type,data,op)    
	{
		if(data.length == 0){
			let te = "No data";
			if(variety == "Transition") te = "No transitions in classification";	
			variety = "No Graph"; type = "No Graph",op = {label:te, key:[]};
		}
		
		this.variety = variety;
		this.type = type;
		this.data = data;
		this.op = op;
		this.set_range();
		this.expand_boundaries();
		
		this.initialise_axes();
		
		this.range_default = copy(this.range);
	}
	
	initialise_axes()
	{
		switch(plot_variety(this.type)){
		case "Line plot": this.set_axes(); break;
		case "Scatter plot": this.set_axes(); break;
		case "Comp plot": this.comp_init(); break;
		case "Individual plot": this.set_ind_time(); break;
		case "Histogram plot": this.set_axes(); break;
		case "Matrix plot": this.set_axes(); this.set_matrix_colour(); break;
		case "Stat table plot": break;
		case "No graph plot": break;
		default: error("Graph variety not recognised:"+this.type+" "+plot_variety(this.type)); break;
		}
	}


	/// Sets the axes for the plot
	set_range()           
	{
		let xmin = LARGE, xmax = -LARGE, ymin = LARGE, ymax = -LARGE;
		
		for(let li = 0; li < this.data.length; li++){
			let da = this.data[li];
		
			switch(da.type){
			case "Line CI":
				for(let i = 0; i < da.point.length; i++){
					let po = da.point[i];
					if(po.x < xmin) xmin = po.x; if(po.x > xmax) xmax = po.x;
					if(typeof po.CImin == "number"){
						if(po.CImin < ymin) ymin = po.CImin; 
					}
					
					if(typeof po.CImax == "number"){
						if(po.CImax > ymax) ymax = po.CImax;
					}
				}
				break;
				
			case "Individual":
				if(da.tmin < xmin) xmin = da.tmin; if(da.tmax > xmax) xmax = da.tmax;
				break;
			
			case "Bar":
				if(da.x-0.5 < xmin) xmin = da.x-0.5; if(da.x+0.5 > xmax) xmax = da.x+0.5;
				if(0 < ymin) ymin = 0; if(da.y > ymax) ymax = da.y;
				break;
				
			case "ErrorBar":
				if(da.x < xmin) xmin = da.x; if(da.x > xmax) xmax = da.x;
				if(da.ymin < ymin) ymin = da.ymin; if(da.ymax > ymax) ymax = da.ymax;
				break;
				
			case "Matrix": case "Correlation":
				xmin = 0; xmax = da.xlab.length;
				ymin = 0; ymax = da.ylab.length;
				break;
		
			case "Table":
				break;
			
			case "VertLine": case "VertDashLine":
				{
					let x = da.x;
					if(x < xmin) xmin = x; if(x > xmax) xmax = x;
				}
				break;
			
			case "PriorLine":
				break;
				
			case "BurninLine":
				break;
				
			default:
				for(let i = 0; i < da.point.length; i++){
					let po = da.point[i];
					if(po.x < xmin) xmin = po.x; if(po.x > xmax) xmax = po.x;
					if(po.y < ymin) ymin = po.y; if(po.y > ymax) ymax = po.y;
				}
				break;
			}
		}
		
		if(xmax < xmin+TINY) xmax = xmin;
		
		if(xmin == xmax){     // Ensures thay xmin and xmax are not the same
			if(xmax > 0){ xmin = 0; xmax *= 1.1;}
			else{
				if(xmin < 0){ xmax = 0; xmin *= 1.1;}
				else xmax = 1;
			}
		}
		else{
			if(this.variety != "Scatter"){
				if(xmin > 0 && xmax > 0 && xmin < 0.8*xmax) xmin = 0; // Puts lower bound at zero
			}
		}
		
		if(ymax < ymin+TINY) ymax = ymin;
		
		if(ymin == ymax){     // Ensures thay ymin and ymax are not the same
			if(ymax > 0){ ymin = 0; ymax *= 1.1;}
			else{
				if(ymin < 0){ ymax = 0; ymin *= 1.1;}
				else ymax = 1;
			}
		}
		else{
			if(this.variety != "Scatter"){
				if(ymin > 0 && ymax > 0 && ymin < 0.8*ymax) ymin = 0; // Puts lower bound at zero
			}
		}
	
		if(this.op.yaxis == false) ymin = 0;
		
		this.range = {xmin:xmin, ymin:ymin, xmax:xmax, ymax:ymax};
	}
	
	
	/// Sets the axes for the plot
	set_axes()           
	{
		let ra = this.range;
		
		this.tick.x = this.get_tick_marks(ra.xmin,ra.xmax,3);
		
		let wmax = 0;
		if(this.op.yaxis == false) this.tick.y = [];
		else{
			this.tick.y = this.get_tick_marks(ra.ymin,ra.ymax,3);
		
			// Works out the width of the y-ticks 
			let fo = get_font(si_graph_tick);
			for(let i = 0; i < this.tick.y.length; i++){
				let w = text_width(this.tick.y[i].te,fo);
				if(w > wmax) wmax = w;
			}
		}
		this.tick.wmax = wmax;
	}
	
	
	/// Given a range and a minimum number of ticks works out how ticks should be distributed
	get_tick_marks(min,max,nummin)
	{
		let tv = this.tick.value;
		let i = tv.length-1; while(i >= 0 && Math.floor((max-min)/tv[i]) < nummin) i--;
		
		let si = tv[i];

		let val = Math.floor(min/si + ALMOST_ONE)*si;
		
		let tick = [];
		while(val < max+TINY){
			tick.push({value:val, te:this.rn(val)}); 
			val += si;
		}
		
		return tick;
	}
	
	
	/// Looks to expand boundaries do they neatly fit with ticks
	expand_boundaries()
	{
		let ra = this.range;
		if(ra.xmin != LARGE){
			let b = this.expand_bound(ra.xmin,ra.xmax,4);
			ra.xmin = b.minb; ra.xmax = b.maxb;
		}
		if(ra.ymin != LARGE){
			let b = this.expand_bound(ra.ymin,ra.ymax,4);
			ra.ymin = b.minb; ra.ymax = b.maxb;
		}
	}
	
	
	/// Looks to expand a boundary to neatly fit with ticks
	expand_bound(min,max,nummin)
	{
		let tv = this.tick.value;
		let i = tv.length-1; while(i >= 0 && Math.floor((max-min)/tv[i]) < nummin) i--;
		
		let si = tv[i];

		let val = Math.floor(min/si + ALMOST_ONE)*si;
		let minb = min;
		if(val > minb) minb = val-si;
		
		val = Math.floor(max/si)*si;
		let maxb = max;
		if(val < maxb) maxb = val+si;		
		
		return {minb:minb, maxb:maxb};
	}
	
	
	/// Sets the timeline for individual plots
	set_ind_time()
	{
		let ra = this.range;
		
		let tv = this.tick.value;
		let i = tv.length-1; while(i >= 0 && Math.floor((ra.xmax-ra.xmin)/tv[i]) < 4) i--;
		
		let si = tv[i];

		let val = Math.floor(ra.xmin/si + ALMOST_ONE)*si;
		
		this.tick.x = [];
		while(val < ra.xmax+TINY){
			this.tick.x.push({value:val, te:this.rn(val)}); 
			val += si;
		}
	}
	
	
	/// Sets the colours used to represent the matrix
	set_matrix_colour()
	{
		let da = this.data[0];
		let mat = da.mat;
		da.mat_col = copy(mat);
		let mat_col = da.mat_col;
	
		let min = LARGE, max = -LARGE;
		for(let j = 0; j < mat.length; j++){
			for(let i = 0; i < mat[j].length; i++){
				let val = Number(mat[j][i]);
				if(val < min) min = val; if(val > max) max = val;
			}
		}
		
		if(da.type == "Correlation"){ min = -1; max = 1;}
		
		if(min == max){
			if(min == 0){ min = -1; max = 1;}
			else{
				if(min > 0) min = 0;
				else max = 0;
			}
		}
		
		let mag = max; if(-min > mag) mag = -min;
		
		for(let j = 0; j < mat.length; j++){
			for(let i = 0; i < mat[j].length; i++){
				let val = Number(mat[j][i]);
				
				mat_col[j][i] = this.get_matrix_col(val,mag);
			}
		}
		
		let grad = [];
		let imax = 100;
		for(let i = 0; i <= imax; i++){
			grad.push(this.get_matrix_col(min+(max-min)*i/imax,mag));
		}
		
		this.matrix_key = {min:min, max:max, tick:this.get_tick_marks(min,max,3), grad:grad};
	}
	
	
	/// Gets a matrix colour from a value
	get_matrix_col(val,mag)
	{
		if(val > 0){
			let f = val/mag;
			return "rgb("+Math.floor(255*(1-f))+","+Math.floor(255*(1-f))+","+255+")";
		}
		else{
			let f = -val/mag;
			return "rgb("+255+","+Math.floor(255*(1-f))+","+Math.floor(255*(1-f))+")";
		}
	}
	
				
	/// Initialises the view of the comparments
	comp_init()
	{
		let species = this.op.species;
		let anim = this.animation;
		
		anim.playframe = 0;
		anim.playframe_max = this.data[0].point.length-1;
		anim.playing = false;
	}
	
	
	/// Gets the clock time (in milliseconds)
	clock()
	{
		return (new Date()).getTime();
	}
	
	
	/// Activates when the play button is pressed
	press_play()
	{
		let anim = this.animation;
		
		if(anim.playing == false){
			anim.playing = true;
			if(anim.playframe == anim.playframe_max) anim.playframe = 0;
			anim.playframe_start = anim.playframe;
			anim.clock = this.clock();
			this.playanim();
		}
		else{
			anim.playing = false;
		}			
	}
	
	
	/// Activates when the timebar is pressed
	press_timebar(bu)
	{
		let lay = get_lay("AnimControls");
			
		let x = inter.mx-bu.x-lay.x;
		let frac = (x-bu.dy/2)/(bu.dx-bu.dy);
		
		let anim = inter.graph.animation;
		anim.playing = false;
		anim.playframe = Math.round(frac*anim.playframe_max);
		if(anim.playframe < 0) anim.playframe = 0;
		if(anim.playframe > anim.playframe_max) anim.playframe = anim.playframe_max;
	}
	
	
	/// Changes the frame (forward or backward)
	change_frame(d)
	{
		let anim = inter.graph.animation;
		anim.playing = false;	
		anim.playframe += d;
		if(anim.playframe < 0) anim.playframe = 0;
		if(anim.playframe > anim.playframe_max) anim.playframe = anim.playframe_max;
	}
	
	
	/// Generates an animation
	playanim()
	{
		let anim = this.animation;
	
		if(anim.playing == true){
			anim.playframe = Math.round(anim.playframe_start + (anim.speed.value*(this.clock()-anim.clock)/4000)*anim.playframe_max);
			if(anim.playframe >= anim.playframe_max){
				anim.playframe = anim.playframe_max;
				anim.playing = false;
			}
			
			replot_layer("GraphCompartments");
			replot_layer("AnimControls");
			plot_screen();
			
			if(anim.playframe < anim.playframe_max){
				setTimeout(function(){ inter.graph.playanim();}, 10);
			}
		}
	}
	

	/// Rounds a number when put in interface
	rn(n)                                           
	{
		n = parseFloat(n);

		let x = 0;
		do{
			let num = n.toFixed(x);
			if((num-n)*(num-n) < 0.000000000000001) return num
			x++;
		}while(true);
	}


	/// Plots the content of the graph
	content(lay)
	{
		switch(this.variety){
		case "Individual":
			{
				let y = 0, dy = 1.8, gap = 0.2; 
				for(let i = 0; i < this.data.length; i++){
					let da = this.data[i];
					lay.add_button({x:0, y:y, dx:lay.dx-scrollw, dy:dy, name:da.name, info:da.info, col_timeline:da.col_timeline, obs:da.obs, type:"IndTimelineBut"});
					y += dy+gap;
				}
				lay.add_button({x:0, y:0, dx:lay.dx-scrollw, dy:y, ac:"TimelineGrab", type:"Nothing"});
			}
			break;
			
		case "Histogram":
			lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"GraphContentBut"});
			break;
			
		case "Matrix":
			{
				let da = this.data[0];
				for(let j = 0; j < da.ylab.length; j++){
					let fracj = j/ da.ylab.length;
					for(let i = 0; i < da.xlab.length; i++){	
						let fraci = i/da.xlab.length; 
						
						let dx = lay.dx/da.xlab.length;
						let dy = lay.dy/da.ylab.length;
						let marx = 0.03*dx;
						let mary = 0.03*dy;
						lay.add_button({te:da.mat[j][i], x:fraci*lay.dx+marx, y:fracj*lay.dy+mary, dx:dx-2*marx, dy:dy-2*mary, value:da.mat[j][i], col:da.mat_col[j][i], type:"MatrixEleBut", ac:"MatrixEleBut"});
					}
				}
			}
			break;
			
		default:
			lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"GraphContentBut"});
			lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"GraphGrab", type:"Nothing"});
			break;
		}
	}
	
	
	/// Draws the content of the graph
	draw_content_button(x,y,dx,dy)
	{
		fill_rectangle(x,y,dx,dy,WHITE); 
		
		let ra = this.range;
		
		for(let li = 0; li < this.data.length; li++){
			let da = this.data[li];
		
			switch(da.type){
			case "Line": case "BurninLine": this.draw_graph_line(x,y,dx,dy,ra,da); break;
			case "PriorLine": this.draw_graph_line(x,y,dx,dy,ra,da); break;
			case "Line CI": this.draw_graph_CI(x,y,dx,dy,ra,da); this.draw_graph_line(x,y,dx,dy,ra,da); break;
			case "Bar": this.draw_bar(x,y,dx,dy,ra,da); break;
			case "ErrorBar": this.draw_errorbar(x,y,dx,dy,ra,da); break;
			case "Distribution": this.draw_distribution(x,y,dx,dy,ra,da); break;
			case "VertLine": this.draw_vert_line(x,y,dx,dy,ra,da); break;
			case "VertDashLine": this.draw_vert_dash_line(x,y,dx,dy,ra,da); break;
			case "Points": this.draw_point(x,y,dx,dy,ra,da); break;
			default: error("Option not recongnised: "+da.type); break;
			}
		}
	}
	
	
	/// Draws and individual timeline
	draw_ind_timeline_button(x,y,dx,dy,name,info,col_timeline,obs)
	{
		let bar = 1;
		//fill_rectangle(x,y,dx,dy,RED);
		plot_text(name,x+0.2,y+0.65,get_font(0.8),BLACK); 
		right_text(info,x+dx-0.2,y+0.65,get_font(0.8),BLACK); 
		
		let ra = this.range;
			
		let pos = [];
		for(let i = 0; i < col_timeline.length; i++){
			pos.push(ro_down(x+dx*((col_timeline[i].t-ra.xmin)/(ra.xmax-ra.xmin))));
		}
	
		let y1 = ro_down(y+dy/2), y2 = ro_down(y+dy);

		for(let i = 0; i < col_timeline.length-1; i++){
			cv.globalAlpha = col_timeline[i].alpha;
			cv.beginPath();
			cv.rect(pos[i],y1,pos[i+1]-pos[i],y2-y1);
			cv.fillStyle = "rgb("+col_timeline[i].col+")";
			cv.fill();
		}
		cv.globalAlpha = 1;
	
		let r = 0.32, si = 0.34, si2 = 0.4;
		for(let i = 0; i < obs.length; i++){
			let ob = obs[i];
			let xx = x+dx*((ob.t-ra.xmin)/(ra.xmax-ra.xmin));
				
			switch(ob.type){
			case "Transition":
				fill_rectangle(xx-si2/2,y+0.70*dy-si2/2,si2/2,si2,ob.col[0],NORMLINE); 
				fill_rectangle(xx,y+0.70*dy-si2/2,si2/2,si2,ob.col[1],NORMLINE); 
			
				draw_rectangle(xx-si2/2,y+0.70*dy-si2/2,si2,si2,BLACK,NORMLINE); 
				break;
			
			case "Compartment":
				fill_circle(xx,y+0.70*dy,r,ob.col,BLACK,NORMLINE);   
				break;
			
			case "Diag. Test":
				if(ob.res == true) draw_rect(xx-si/2,y+0.70*dy-si/2,si,si,BLACK,BLACK,NORMLINE); 
				else draw_rect(xx-si/2,y+0.70*dy-si/2,si,si,WHITE,BLACK,NORMLINE);  
				break;
			}
		}
	}
	
	
	/// Draws a graph line
	draw_graph_line(x,y,dx,dy,ra,da)
	{
		cv.beginPath();

		if(da.dash) set_dash_type(da.dash);

		for(let i = 0; i < da.point.length; i++){
			let po = da.point[i];
			let xp = ro(x+dx*((po.x-ra.xmin)/(ra.xmax-ra.xmin)));
			let yp = ro(y+dy-dy*((po.y-ra.ymin)/(ra.ymax-ra.ymin)));

			if(i == 0) cv.moveTo(xp,yp);
			else cv.lineTo(xp,yp);
		}
		
		cv.lineWidth = 2;
		if(da.view == "Graph (lines)") cv.lineWidth = 1;
		cv.strokeStyle = da.col;
		cv.stroke();
		
		if(da.dash) set_dash_type(0);
	}
	
	
	/// Draws points
	draw_point(x,y,dx,dy,ra,da)
	{
		cv.fillStyle = da.col;
	
		cv.beginPath();
		let r = ro(0.1);
		
		for(let i = 0; i < da.point.length; i++){
			let po = da.point[i];
			let xp = ro(x+dx*((po.x-ra.xmin)/(ra.xmax-ra.xmin)));
			let yp = ro(y+dy-dy*((po.y-ra.ymin)/(ra.ymax-ra.ymin)));
			
			cv.beginPath();
			cv.arc(xp,yp,r,0,2*Math.PI);
			cv.fill();
		}
	}
	
	
	/// Draws a distribution
	draw_distribution(x,y,dx,dy,ra,da)
	{
		let vec=[];
		for(let i = 0; i < da.point.length; i++){
			let po = da.point[i];
			let xp = ro(x+dx*((po.x-ra.xmin)/(ra.xmax-ra.xmin)));
			let yp = ro(y+dy-dy*((po.y-ra.ymin)/(ra.ymax-ra.ymin)));
			vec.push({ xp:xp, yp:yp});
		}
		
		// Draw background
		cv.beginPath();
		for(let i = 0; i < vec.length; i++){
			if(i == 0) cv.moveTo(vec[i].xp,vec[i].yp);
			else cv.lineTo(vec[i].xp,vec[i].yp);
		}
		cv.fillStyle = da.col;
		cv.globalAlpha = 0.2;
		cv.fill();
		cv.globalAlpha = 1;
		
		// Plots lower CI 
		if(da.CImin){
			let xp_min = ro(x+dx*((da.CImin-ra.xmin)/(ra.xmax-ra.xmin)));
				
			cv.beginPath();
			let xp_old, yp_old, flag = false;
			for(let i = 0; i < vec.length; i++){
				let po = da.point[i];
				
				let xp = vec[i].xp, yp = vec[i].yp;
				
				if(i > 0 && xp > xp_min){
					let f = (xp_min-xp_old)/(xp-xp_old);
					yp = yp_old + (yp-yp_old)*f;
					xp = xp_min;
					flag = true;
				}
			
				if(i == 0) cv.moveTo(xp,yp);
				else cv.lineTo(xp,yp);
				
				if(flag == true){
					yp = ro(y+dy-dy*((0-ra.ymin)/(ra.ymax-ra.ymin)));
					cv.lineTo(xp,yp);
					break;
				}
				xp_old = xp; yp_old = yp;
			}
			cv.fillStyle = da.col;//WHITE;
			cv.globalAlpha = 0.5;
			cv.fill();
			cv.globalAlpha = 1;
		}
		
		// Plots upper CI 
		if(da.CImax){
			let xp_max = ro(x+dx*((da.CImax-ra.xmin)/(ra.xmax-ra.xmin)));
				
			cv.beginPath();
			let xp_old, yp_old, flag = false;
			for(let i = vec.length-1; i >= 0; i--){
				let po = da.point[i];
				
				let xp = vec[i].xp, yp = vec[i].yp;
				
				if(i != vec.length-1 && xp < xp_max){
					let f = (xp_max-xp_old)/(xp-xp_old);
					yp = yp_old + (yp-yp_old)*f;
					xp = xp_max;
					flag = true;
				}
			
				if(i == 0) cv.moveTo(xp,yp);
				else cv.lineTo(xp,yp);
				
				if(flag == true){
					yp = ro(y+dy-dy*((0-ra.ymin)/(ra.ymax-ra.ymin)));
					cv.lineTo(xp,yp);
					break;
				}
				xp_old = xp; yp_old = yp;
			}
			cv.fillStyle = da.col;//WHITE;
			cv.globalAlpha = 0.5;
			cv.fill();
			cv.globalAlpha = 1;
		}
		
		// Draws line
		cv.beginPath();
		for(let i = 0; i < vec.length; i++){
			if(i == 0) cv.moveTo(vec[i].xp,vec[i].yp);
			else cv.lineTo(vec[i].xp,vec[i].yp);
		}
		cv.lineWidth = 1;
		cv.strokeStyle = da.col;
		cv.stroke();	
	}
	
	
	/// Draws a verticle line
	draw_vert_line(x,y,dx,dy,ra,da)
	{
		let xx = x+dx*((da.x-ra.xmin)/(ra.xmax-ra.xmin));
		let xp = ro(xx);
		cv.beginPath();
		
		cv.moveTo(xp,ro(y));
		cv.lineTo(xp,ro(y+dy));
		cv.lineWidth = da.thick;
		cv.strokeStyle = da.col;
		cv.stroke();	
		
		let fo = get_font(0.8);
		let wid = text_width(da.te,fo);
		vert_text(da.te,xx-0.3,y+wid+0.2,fo,da.col,10); 
	}
	
	
	/// Draws a verticle dashed line
	draw_vert_dash_line(x,y,dx,dy,ra,da)
	{
		set_dash_type(2);
		this.draw_vert_line(x,y,dx,dy,ra,da);
		set_dash_type(0);
	}
	
	
	/// Shades a credible interval
	draw_graph_CI(x,y,dx,dy,ra,da)
	{
		cv.beginPath();

		for(let i = 0; i < da.point.length; i++){
			let po = da.point[i];
			let xp = ro(x+dx*((po.x-ra.xmin)/(ra.xmax-ra.xmin)));
			let yp = ro(y+dy-dy*((po.CImax-ra.ymin)/(ra.ymax-ra.ymin)));

			if(i == 0) cv.moveTo(xp,yp);
			else cv.lineTo(xp,yp);
		}
		
		for(let i = da.point.length-1; i >= 0; i--){
			let po = da.point[i];
			let xp = ro(x+dx*((po.x-ra.xmin)/(ra.xmax-ra.xmin)));
			let yp = ro(y+dy-dy*((po.CImin-ra.ymin)/(ra.ymax-ra.ymin)));
			cv.lineTo(xp,yp);
		}
		
		cv.globalAlpha = 0.1;
		cv.fillStyle = da.col;
		cv.fill();
		cv.globalAlpha = 1;
	}
	
	
	/// Draws a bar on a histogram
	draw_bar(x,y,dx,dy,ra,da)
	{
		let x1 = x+dx*((da.x-da.thick/2-ra.xmin)/(ra.xmax-ra.xmin));
		let x2 = x+dx*((da.x+da.thick/2-ra.xmin)/(ra.xmax-ra.xmin));
		
		let y1 = y+dy-dy*((da.y-ra.ymin)/(ra.ymax-ra.ymin));
		let y2 = y+dy-dy*((0-ra.ymin)/(ra.ymax-ra.ymin));

		fill_rectangle(x1,y1,x2-x1,y2-y1,da.col);
		draw_rectangle(x1,y1,x2-x1,y2-y1,dark_colour(da.col),NORMLINE);
	}
	
	
	/// Draws an errorbar 
	draw_errorbar(x,y,dx,dy,ra,da)
	{
		let d = 0.3;
		let xp = x+dx*((da.x-ra.xmin)/(ra.xmax-ra.xmin));
		let yp = y+dy-dy*((da.y-ra.ymin)/(ra.ymax-ra.ymin));
		let ymin = y+dy-dy*((da.ymin-ra.ymin)/(ra.ymax-ra.ymin));
		let ymax = y+dy-dy*((da.ymax-ra.ymin)/(ra.ymax-ra.ymin));

		let thick = NORMLINE;
		draw_line(xp,ymin,xp,ymax,da.col,thick);
		draw_line(xp-d,ymin,xp+d,ymin,da.col,NORMLINE);
		draw_line(xp-d/2,yp,xp+d/2,yp,da.col,NORMLINE);
		draw_line(xp-d,ymax,xp+d,ymax,da.col,NORMLINE);
	}
	
	
	/// Draws the x tick marks
	draw_xtick(x,y,dx,dy)
	{
		let ra = this.range;
		
		let ti = this.tick.x;
		for(let i = 0; i < ti.length; i++){
			let xx = x + dx*(ti[i].value-ra.xmin)/(ra.xmax-ra.xmin);
			
			draw_line(xx,y,xx,y+dy,BLACK,NORMLINE);
		}
	}
	
	
	/// Draws the x tick marks
	draw_xtick_label(x,y,dx,dy,mar,x_lab,x_vert,x_param,ov)
	{
		fill_rectangle(x,y,dx,dy,WHITE); 
	
		let ra = this.range;
		
		let col = BLACK; if(ov) col = DGREY;
		
		if(x_lab){
			let wid = (dx-mar.left-mar.right)/x_lab.length;
			
			for(let i = 0; i < x_lab.length; i++){
				let xx = x+mar.left + (dx-mar.left-mar.right)*(x_lab[i].x-ra.xmin)/(ra.xmax-ra.xmin);
		
				if(x_param == true){
					if(x_vert == true) right_param_text(x_lab[i].name,xx,y+0.6,si_histo_label,col,wid);	
					else center_param_text(x_lab[i].name,xx,y+0.6,si_histo_label,col,wid);
				}	
				else{
					if(x_vert == true) right_text(x_lab[i].name,xx,y+0.6,get_font(si_histo_label),col,wid);	
					else center_text(x_lab[i].name,xx,y+0.6,get_font(si_histo_label),col,wid);
				}
			}
		}
		else{
			let fo = get_font(si_graph_tick);
			let ti = this.tick.x;
			for(let i = 0; i < ti.length; i++){
				let xx = x+mar.left + (dx-mar.left-mar.right)*(ti[i].value-ra.xmin)/(ra.xmax-ra.xmin);
			
				center_text(ti[i].te,xx,y+0.6,fo,col);  
			}
		}
	}
	
	
	/// Draws the y tick marks
	draw_ytick(x,y,dx,dy)
	{
		let ra = this.range;
		
		let ti = this.tick.y;
		for(let i = 0; i < ti.length; i++){
			let yy = y + dy - dy*(ti[i].value-ra.ymin)/(ra.ymax-ra.ymin);
			
			draw_line(x,yy,x+dx,yy,BLACK,NORMLINE);
		}
	}
	
	
	/// Draws the y tick marks
	draw_ytick_label(x,y,dx,dy,mar,y_lab,y_hor,y_param,ov)
	{
		fill_rectangle(x,y,dx,dy,WHITE); 
		
		let ra = this.range;
		
		let col = BLACK; if(ov) col = DGREY;
		
		if(y_lab){
			let fo = get_font(si_histo_label);
			for(let i = 0; i < y_lab.length; i++){
				let yy = y + dy - 0.2 - (dy-yaxisgap-mar.top)*(y_lab[i].y-ra.ymin)/(ra.ymax-ra.ymin);
				
				if(y_param == true){
					if(y_hor == true) right_param_text(y_lab[i].name,x+dx,yy,si_histo_label,col,mar.left-2,20);
					else center_vert_param_text(y_lab[i].name,x+dx,yy,si_histo_label,col,20); 
				}
				else{
					if(y_hor == true) right_text(y_lab[i].name,x+dx,yy,fo,col,mar.left-2);
					else center_vert_text(y_lab[i].name,x+dx,yy,fo,col); 
				}
			}
		}
		else{
			let fo = get_font(si_graph_tick);
			let ti = this.tick.y;
			for(let i = 0; i < ti.length; i++){
				let yy = y + dy - 0.2 - (dy-yaxisgap-mar.top)*(ti[i].value-ra.ymin)/(ra.ymax-ra.ymin);
			
				right_text(ti[i].te,x+dx,yy,fo,col); 
			}
		}
	}
	
	
	/// Plots the axes
	axes(lay)
	{
		let mar = lay.op.mar;
		
		lay.add_button({x:mar.left, y:lay.dy-mar.bottom, dx:lay.dx-mar.left, dy:0, type:"x-axis"});
		
		let x_lab = lay.op.x_label;
		let x_vert = lay.op.x_label_vert;
		let x_param = lay.op.x_param;
		
		let yl = lay.dy-mar.bottom;
		if(x_lab){
			yl += 0.8;
			lay.add_button({x:0, y:yl, dx:lay.dx, dy:0.8, mar:mar, x_lab:x_lab, x_vert:x_vert, x_param:x_param, type:"x-tick-label"});
			yl += 1;
		}
		else{
			lay.add_button({x:mar.left, y:yl, dx:lay.dx-mar.left-mar.right, dy:tick_si, type:"x-tick"});
			yl += tick_si+0.3;
			lay.add_button({x:0, y:yl, dx:lay.dx, dy:0.8, mar:mar, type:"x-tick-label", ac:"X Tick"});
			yl += 1;
			
			let xz = lay.dx-4, yz = lay.dy-1.7, f= 0.8;
			lay.add_button({x:xz, y:yz+0.15, dx:1.8*f, dy:2*f, ac:"ZoomInGraph", type:"ZoomIn"});

			lay.add_button({x:xz+2*f, y:yz+0.15, dx:1.8*f, dy:2*f, ac:"ZoomOutGraph", type:"ZoomOut"});
			
			if(this.variety == "Distribution"){
				lay.add_button({x:mar.left, y:yz+0.2, dx:1.4, dy:1.4, ac:"GraphSettings", type:"Settings"});
			}
			
			if(this.variety == "Trace"){
				lay.add_button({x:mar.left, y:yz, dx:1.4, dy:1.4, ac:"TraceSettings", type:"Settings"});
			}
		}
		
		lay.add_button({te:this.op.x_label, x:mar.left, y:lay.dy-si_graph_label-0.21, dx:lay.dx-mar.left-mar.right, dy:si_graph_label+0.1, mar:mar, param:this.op.x_param, type:"x-label"});
	
	
		lay.add_button({x:mar.left, y:0, dx:0, dy:lay.dy-mar.bottom, param:this.op.param, type:"y-axis"});
		
		let y_lab = lay.op.y_label;
		let y_hor = lay.op.y_label_hor;
		let y_param = lay.op.y_param;
	
		if(y_lab){
			let xl = 1;
			lay.add_button({x:xl, y:0, dx:mar.left-1.5, dy:lay.dy-mar.bottom+yaxisgap, mar:mar, y_lab:y_lab, y_hor:y_hor, y_param:y_param, type:"y-tick-label"});
			xl += 1;
		}
		else{
			let xl = mar.left-0.2-this.tick.wmax;
			xl -= tick_si;
			lay.add_button({x:mar.left-tick_si, y:mar.top, dx:tick_si, dy:lay.dy-mar.top-mar.bottom, type:"y-tick"});
			lay.add_button({x:xl, y:0, dx:this.tick.wmax, dy:lay.dy-mar.bottom+yaxisgap, mar:mar, type:"y-tick-label", ac:"Y Tick"});
			xl -= si_graph_label+0.4;
		}
		
		lay.add_button({te:this.op.y_label, x:0, y:mar.top, dx:si_graph_label+0.2, dy:lay.dy-mar.top-mar.bottom, mar:mar, param:this.op.y_param, type:"y-label"});
	}
	
	
	/// Plots the time axis
	time_axis(lay)
	{
		let mar = lay.op.mar;
		
		lay.add_button({x:mar.left, y:lay.dy-mar.bottom, dx:lay.dx-mar.left, dy:0, type:"x-axis"});
		
		let yl = lay.dy-mar.bottom;
		lay.add_button({x:mar.left, y:yl, dx:lay.dx-mar.left-mar.right, dy:tick_si, type:"x-tick"});
		yl += tick_si+0.3;
		lay.add_button({x:0, y:yl, dx:lay.dx, dy:0.8, mar:mar, type:"x-tick-label", ac:"X Tick"});
		yl += 1;
		lay.add_button({te:this.op.x_label, x:mar.left, y:yl, dx:lay.dx-mar.left-mar.right, dy:si_graph_label+0.1, mar:mar, param:this.op.x_param, type:"x-label"});
		
		let xz = lay.dx-6, yz = lay.dy-1.7, f= 0.8;
		lay.add_button({x:xz, y:yz, dx:1.8*f, dy:2*f, ac:"ZoomInTimeline", type:"ZoomIn"});

		lay.add_button({x:xz+2*f, y:yz, dx:1.8*f, dy:2*f, ac:"ZoomOutTimeline", type:"ZoomOut"});
	}
	
	
	/// Zooms into the timeline
	zoom_timeline_factor(fac,f)
	{
		if(f == undefined) f = 0.5;
		let ra = this.range;
		let mean = ra.xmax*f+ra.xmin*(1-f);
		
		ra.xmin = mean - (mean-ra.xmin)/fac;
		ra.xmax = mean + (ra.xmax-mean)/fac;
		this.set_ind_time();
	}
	
	
	/// Zooms into the timeline
	zoom_graph_factor(fac,fx,fy)
	{
		if(fx == undefined) fx = 0.5;
		if(fy == undefined) fy = 0.5;
		
		let ra = this.range;
		let meanx = ra.xmax*fx+ra.xmin*(1-fx);
		ra.xmin = meanx - (meanx-ra.xmin)/fac;
		ra.xmax = meanx + (ra.xmax-meanx)/fac;
		
		let meany = ra.ymax*fy+ra.ymin*(1-fy);
		ra.ymin = meany - (meany-ra.ymin)/fac;
		ra.ymax = meany + (ra.ymax-meany)/fac;
		if(this.op.yaxis == false){
			ra.ymax -= ra.ymin;
			ra.ymin = 0;
		}
		
		this.set_axes();
	}
	
	
	/// Replots layers associated with the graph
	replot()
	{
		switch(this.variety){
		case "Population": case "Transition":	case "Distribution": case "Trace": case "Samples":
			this.set_axes();
			replot_layer("Axes");
			break;
			
		case "Individual":
			this.set_ind_time();
			replot_layer("TimeAxis");
			break;
		}
		
		replot_layer("GraphContent");
		plot_screen();
	}
	
	
	/// Adds all the compartment buttons
	add_compartment_buts(lay)
	{
		let p = this.op.p;
		let cl = this.op.cl;
		
		lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, p:p, cl:cl, ac:"ClassificationBack", type:"Nothing"});
	
		let claa = this.op.species[p].cla[cl];
		
		let cam = claa.camera;
		
		let ac = "Graph Comp";
	
		if(this.data.length != claa.ncomp) error("number of data not right");

		let fr = this.animation.playframe;
	
		for(let k = 0; k < claa.ncomp; k++){
			let c = claa.comp[k];
		
			let value = this.data[k].point[fr].y;
			let val = value/this.range.ymax;
			if(val > 1) val = 1; if(val < 0) val = 0;
			let valc = 255-Math.floor(255*val);
			let col = "rgb("+valc+","+valc+","+valc+")";
					
			switch(c.type){
			case "box":
				{
					let pt = trans_point(c.x,c.y,cam,lay);
					let w =	c.w*cam.scale;
					let h =	c.h*cam.scale; 
					let x = pt.x-w/2, y = pt.y-h/2;
					
					//lay.add_button({te:claa.comp[k].name, x:x, y:y, dx:w, dy:h, ac:ac, type:"Compartment", col:c.col, p:p, cl:cl, i:k});
					lay.add_button({te:claa.comp[k].name, x:x, y:y, dx:w, dy:h, ac:ac, type:"CompGraph", value:Math.floor(value), col:col, p:p, cl:cl, i:k});
				}
				break;
			
			case "latlng":
				{
					let pt = trans_point(c.x,c.y,cam,lay);
					let w =	2*latlng_radius*cam.scale*cam.ruler, h = w; 
					let x = pt.x-w/2, y = pt.y-h/2;
					lay.add_button({te:claa.comp[k].name, x:x, y:y, dx:w, dy:h, ac:ac, type:"CompLatLng", col:c.col, p:p, cl:cl, i:k});
				}
				break;
				
			case "boundary":
				{
					let box = c.feature.box;
				
					let pmin = trans_point(box.xmin,box.ymin,cam,lay);
					let pmax = trans_point(box.xmax,box.ymax,cam,lay);
					
					lay.add_button({x:pmin.x, y:pmin.y, dx:pmax.x-pmin.x, dy:pmax.y-pmin.y, ac:ac, type:"CompMap", polygon:c.feature.polygon, mask:c.mask, col:c.col, p:p, cl:cl, i:k});
				}
				break;

			default: error("Option not recognised 82"); break;
			}
		}
	}


	/// Adds all the compartment buttons
	add_transition_buts(lay)
	{
		let p = this.op.p;
		let cl = this.op.cl;
		let claa = this.op.species[p].cla[cl];
		model.add_transition_buts2(p,cl,"No Click",claa,false,lay);
	}
	

	/// Adds buttons for 
	add_animcontrol_buts(lay)
	{
		let p = this.op.p;
		let cl = this.op.cl;
		
		lay.add_button({x:playbar_mar, y:0.4, dx:lay.dx-2*playbar_mar, dy:1.1, type:"Timebar", ac:"Timebar"});
			
		let play_r = 1;
		lay.add_button({x:lay.dx/2-play_r, y:lay.dy-2*play_r, dx:2*play_r, dy:2*play_r, type:"PlayButton", ac:"PlayButton"});
		
	
		let dx = 3;
		lay.add_button({x:lay.dx/2-dx-play_r, y:lay.dy-2*play_r, dx:2*play_r, dy:2*play_r, type:"PlayBackward", ac:"PlayBackward"});
			
		lay.add_button({x:lay.dx/2+dx-play_r, y:lay.dy-2*play_r, dx:2*play_r, dy:2*play_r, type:"PlayForward", ac:"PlayForward"});
		
		let xz = lay.dx-6, yz = lay.dy-2;
		lay.add_button({x:xz, y:yz, dx:1.8, dy:2, p:p, cl:cl, ac:"ZoomIn", type:"ZoomIn"});

		lay.add_button({x:xz+2, y:yz, dx:1.8, dy:1.8, p:p, cl:cl, ac:"ZoomOut", type:"ZoomOut"});
		
		lay.add_button({x:2.4, y:yz+0.2, dx:1.4, dy:1.4, ac:"Settings", type:"Settings"});
	}
	
	
	/// Setting for the distribution plots
	settings_dist_bubble(cont)
	{
		cont.dx = 9;
		bubble_addtitle(cont,"Settings",{});
		
		bubble_add_minititle(cont,"Distribution");
		
		let ds = inf_result.plot_filter.dist_settings;
		bubble_addradio(cont,0,"kde","KDE",ds.radio);
		
		if(ds.radio.value == "kde"){
			cont.y -= 1.5;
			bubble_adddropdown(cont,5,4,ds.sel_h,ds.pos_h);
			cont.y += 0.1;
		}
		
		cont.y += 0.2;
		
		bubble_addradio(cont,0,"bin","Binning",ds.radio);

		if(ds.radio.value == "bin"){
			cont.y -= 1.5;
			bubble_adddropdown(cont,5,4,ds.sel_bin,ds.pos_bin);
			cont.y += 0.1;
		}
	
		cont.y += 0.2;
		bubble_add_minititle(cont,"Annotation");
		
		bubble_addcheckbox(cont,-0.1,"Show mean",ds.show_mean);
		bubble_addcheckbox(cont,-0.1,"Show prior",ds.show_prior);
		
	}
		
	
	/// Setting for the trace plots
	settings_trace_bubble(cont)
	{
		cont.dx = 9;
		bubble_addtitle(cont,"Settings",{});
		bubble_input(cont,"Burn-In %",{type:"burnin"});
		add_end_button(cont,"OK","BurninOK",{});	
	}	
		
		
	/// Alls the user to change the speed of animations
	settings_speed_bubble(cont)
	{
		cont.dx = 4.5;
		bubble_addtitle(cont,"Speed",{});
		
		bubble_addradio(cont,0,"10","×10",this.animation.speed);
		bubble_addradio(cont,0,"5","×5",this.animation.speed);
		bubble_addradio(cont,0,"2","×2",this.animation.speed);
		bubble_addradio(cont,0,"1","×1",this.animation.speed);
		bubble_addradio(cont,0,"0.5","×0.5",this.animation.speed);
		bubble_addradio(cont,0,"0.2","×0.2",this.animation.speed);
		bubble_addradio(cont,0,"0.1","×0.1",this.animation.speed);
		add_end_button(cont,"OK","CloseBubble",{});	
	}
	
	// Creates the graph by adding layers for the content and axes
	create(x,y,w,h,lay) 
	{
		let vari = plot_variety(this.type);
		
		switch(vari){
		case "Line plot":
			{
				let mar = copy(graph_mar); mar.left += this.tick.wmax;
				
				add_layer("GraphContent",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
		
				add_layer("Axes",lay.x+x,lay.y+y,w,h,{mar:mar});
			}
			break;
			
		case "Scatter plot":
			{
				let mar = copy(graph_mar); mar.left += this.tick.wmax;
			
				add_layer("GraphContent",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
		
				add_layer("Axes",lay.x+x,lay.y+y,w,h,{mar:mar});
			}
			break;
			
		case "Comp plot":
			{
				let mar = { right:0, left:0, top:0, bottom:0};
				add_layer("GraphCompartments",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
				add_layer("GraphTransitions",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
		
				let hei = 4;
				add_layer("AnimControls",lay.x+x+mar.left,h+y-mar.top-hei-1,w-mar.right-mar.left,hei,{});
			}
			break;
			
		case "Individual plot":
			{
				let mar = { right:2, left:2, top:0, bottom:3.3};
	
				if(this.data.length == 0){
					let dx = 20;
					lay.add_paragraph("There are no individuals in the system.",dx,lay.dx/2-dx/2,lay.dy/2,BLACK,para_si,para_lh);	
				}
				else{
					add_layer("GraphContent",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left+scrollw,h-mar.top-mar.bottom,{});

					add_layer("TimeAxis",lay.x+x,lay.y+y,w,h,{mar:mar});
				}
			}
			break;
			
		case "Histogram plot":
			{
				let mar = copy(graph_mar); mar.left += this.tick.wmax;
			
				let fo = get_font(si_histo_label);
				let x_label = [];
				let wimax = 0;
				for(let i = 0; i < this.data.length; i++){
					let da = this.data[i];
					if(da.type == "Bar"){
						x_label.push({name:da.name, x:da.x});
						let wi = text_width(da.name,fo);
						if(wi > wimax) wimax = wi;
					}
				}
				
				let ra = this.range;
		
				let x_label_vert = false;
				let dx_label = (w-mar.right-mar.left)/(ra.xmax-ra.xmin);
				if(dx_label < wimax){
					x_label_vert = true;
					if(wimax > 7) wimax = 7;
					mar.bottom = 2.5+wimax;
				}
				
				add_layer("GraphContent",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
			
				add_layer("Axes",lay.x+x,lay.y+y,w,h,{mar:mar, x_label:x_label, x_label_vert:x_label_vert});
			}
			break;
			
		case "Matrix plot":
			{
				let mar = { right:2, left:3.3, top:2, bottom:3.3};
				let fo = get_font(si_histo_label);
				
				if(this.data.length != 1) error("Data not right");
				let da = this.data[0];
				if(da.type != "Matrix" && da.type != "Correlation") error("Data not right"); 
					
				let ra = this.range;
		
				let x_label = [];		
				let wimax = 0;
				for(let i = 0; i < da.xlab.length; i++){
					x_label.push({name:da.xlab[i], x:0.5+i});
					let wi = text_width(da.xlab[i],fo); if(wi > wimax) wimax = wi;
				}
			
				let x_label_vert = false;
				let dx_label = (w-mar.right-mar.left)/(ra.xmax-ra.xmin);
				if(dx_label < wimax){
					x_label_vert = true;
					if(wimax > 7) wimax = 7;
					mar.bottom = 2.5+wimax;
				}
				
				let y_label = [];		
				wimax = 0;
				for(let i = 0; i < da.ylab.length; i++){
					y_label.push({name:da.ylab[i], y:0.5+i});
					let wi = text_width(da.ylab[i],fo); if(wi > wimax) wimax = wi;
				}
			
				let y_label_hor = false;
				let dy_label = (h-mar.top-mar.bottom)/(ra.ymax-ra.ymin);
				if(dy_label < wimax){
					y_label_hor = true;
					if(wimax > 7) wimax = 7;
					mar.left = 2.5+wimax;
				}
				
				add_layer("GraphContent",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
			
				add_layer("Axes",lay.x+x,lay.y+y,w,h,{mar:mar, x_label:x_label, x_label_vert:x_label_vert, x_param:da.x_param, y_label:y_label, y_param:da.y_param, y_label_hor:y_label_hor});
			}
			break;
			
		case "Stat table plot":
			{
				let da = this.data[0];
				add_layer("TableContent",lay.x+x,lay.y+y,w,h,{table:da.table});
			}
			break;
			
		case "No graph plot":
			center_message(this.op.label,lay);
			break;
			
		default: error("Variety not recognised:"+vari); break;
		}
	}
	
	
	/// Exports an image of the graph
	export_image(filename)
	{
		/// Shifts scrollbars to the origin in tables (if necessary)	
		inter.export_image = true;
		
		inter.sca *= 2;
		generate_screen();
		
		let xmin = LARGE, xmax = -LARGE, ymin = LARGE, ymax = -LARGE;
	
		let right_bot_menu; 
	
		let list = [];
		for(let l = 0; l < inter.layer.length; l++){
			let lay = inter.layer[l];
			switch(lay.name){
			case "GraphContent": case "Axes": case "TimeAxis": case "RightBotMenu": case "TableContent":
			case "AnnotationMap": case "Annotation": case "Compartment": case "Transition": 			
				{
					if(lay.name == "RightBotMenu") right_bot_menu = lay;
					
					if(lay.name == "TableContent"){
						let box = get_but_box(lay.but);
					
						if(box.xmax > lay.dx){ lay.dx = box.xmax; lay.inner_dx = box.xmax; }
						if(box.ymax > lay.dy){ lay.dy = box.ymax; lay.inner_dy = box.ymax; }
					}
					
					if(lay.x < xmin) xmin = lay.x;
					if(lay.x + lay.dx > xmax) xmax = lay.x + lay.dx;
					
					if(lay.y < ymin) ymin = lay.y;
					if(lay.y + lay.dy > ymax) ymax = lay.y + lay.dy;
					
					list.push(l);
				}
				break;
			}
		}
		
		if(right_bot_menu){
			let lay = right_bot_menu;
			let box = get_but_box(lay.but);
			if(box.ymax > lay.dy){
				let y = lay.y - (box.ymax-lay.dy); if(y < ymin+1) y = ymin+1;
				lay.y = y;
				lay.dy = box.ymax;
				lay.inner_dy = lay.dy;
			}			
		}
		
		if(list.length == 0){ error("Could not find a source"); return;}
		
		let mar = 1;
		let w = ro(xmax-xmin+2*mar);
		let h = ro(ymax-ymin+2*mar);
		
		let outcan = document.createElement('canvas');
		outcan.width = w;
		outcan.height = h;
		let outcv = outcan.getContext('2d');
	
		outcv.beginPath();
		outcv.rect(0,0,w,h);
		outcv.fillStyle = WHITE;
		outcv.fill();

		for(let i = 0; i < list.length; i++){
			let lay = inter.layer[list[i]];
			lay.plot_buttons(["Settings","ZoomIn","ZoomOut"]);
			if(lay.can){
				outcv.drawImage(lay.can,ro(lay.x-xmin+mar),ro(lay.y-ymin+mar));
			}
		}

		inter.export_image = false;
		
		inter.sca /= 2;
		generate_screen();
		
		if(filename) save_image(outcan,filename);
		else print_image(outcan)
	}
}
