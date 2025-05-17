"use strict";
// Functions which create graphical plots

class Graph
{
	variety;                                         // The thing being plotted (e.g. "Population", "Transition")
	type;                                            // The type of graph (e.g. "Graph", "Compartment")
	data = [];                                       // Data for graph
	view;                                            // Shows how the data is viewed (graph or compartments)
	op;                                              // Options for the graph
	range;                                           // Range for axes 
	range_default;                                   // The default range
	animation = {speed:{te:"×1", value:1, noupdate:true}}; // Stores animation information
	colour_key;                                      // Information about colour key
	ind_sel;                                         // A selected individual
	mag;                                             // The magnitude.
	barh;                                            // The size of bar
	extra_flag;                                      // Set if extra space on right of line plot
	complink = [];                                   // Information about links between compartments
	density_info = {};                               // Information used to make density plots
	tick = {value:[], wmax:0, x:[], y:[]};           // Potential tick marks
	table_width = [];                                // Table width for statistics tables
	bf_store;                                        // Stores information about Bayes factor
	warn;                                            // Any warning to show why the graph could not be plotted
	
	constructor()
	{
		this.initialise_ticks();
	}
	
	
	/// Initialises tick marks used when drawing graphs
	initialise_ticks()           
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
		if(graph_dia){		
			error("GRAPH DIA: variety:"+variety+" type:"+type+" data:"); error(data); error("op:"); error(op);
		}
		
		this.extra_flag = false;
		
		if(data.length == 0 && variety != "No Graph"){
			let te = "No data";
			if(variety == "Transition") te = "No transitions data";	
			variety = "No Graph"; 
			type = "No Graph"; 
			op = {label:te, key:[]};
		}
	
		this.variety = variety;
		this.type = type;
		this.data = data;
		this.barh = 1;	
		this.op = op;
		this.bf_store = undefined;
		
		this.set_range();
	
		this.expand_boundaries();
		
		if(op.def_xrange){ this.range.xmin = op.def_xrange.min; this.range.xmax = op.def_xrange.max;}
		
		this.initialise_axes();
		this.range_default = copy(this.range);
	}
	
	
	/// Initialises the axes
	initialise_axes()
	{
		if(graph_dia) error("GRAPH DIA: initialise axes: "+ plot_variety(this.type));
	
		this.colour_key = undefined;
		
		switch(plot_variety(this.type)){
		case "Line plot": 
			this.set_axes(); 
			break;
		
		case "Scatter plot": 
			this.set_axes(); 
			break;
			
		case "Comp plot": 
			this.anim_init(); 
			this.set_colour_key(this.range.ymin,this.range.ymax); 
			this.comp_colour_init();
			break;
			
		case "Density plot": 
			this.anim_init(); 
			this.set_colour_key(this.range.ymin,this.range.ymax);		
			this.comp_colour_init();
			this.density_init(); 	
			break;
			
		case "Individual plot": case "TransTree plot": case "PhyloTree plot":
			this.set_ind_time(); 
			break;
		
		case "Histogram plot": 
			this.set_axes(); 
			break;
		
		case "HistoAnim plot": 
			this.anim_init(); 
			this.set_axes(); 
			break;
			
		case "Matrix plot":
			this.set_axes(); 
			this.set_matrix_colour(); 
			break;
			
		case "MatrixAnim plot":
			this.anim_init(); 
			this.set_axes(); 
			this.set_matrix_colour(); 
			break;
			
		case "Stat table plot":
			break;
		
		case "CompMatrix plot": 
			this.comp_matrix_init(); 
			break;
		
		case "CompMatrixAnim plot":
			this.anim_init(); 
			this.comp_matrix_init(); 
			break;
		
		case "CompVector plot":
			this.set_colour_key(this.range.ymin,this.range.ymax); 
			this.comp_colour_init();
			break;
			
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
				
			case "Individual": case "InfBar":
				if(da.tmin < xmin) xmin = da.tmin; if(da.tmax > xmax) xmax = da.tmax;
				break;
			
			case "Phylo":
				break;
				
			case "Transmission":
				break;
				
			case "Bar":
				if(da.x-0.5 < xmin) xmin = da.x-0.5; if(da.x+0.5 > xmax) xmax = da.x+0.5;
				if(da.y < ymin) ymin = da.y;
				if(da.y > ymax) ymax = da.y;
				break;
			
			case "Cross":
				if(da.x < xmin) xmin = da.x; if(da.x > xmax) xmax = da.x;
				if(da.y	< ymin) ymin = da.y; if(da.y > ymax) ymax = da.y;
				break;
			
			case "ErrorBar":
				if(da.x < xmin) xmin = da.x; if(da.x > xmax) xmax = da.x;
				if(da.ymin < ymin) ymin = da.ymin; if(da.ymax > ymax) ymax = da.ymax;
				break;
				
			case "Matrix": case "MatrixAnim": case "Correlation":
				xmin = 0; xmax = da.xlab.length;
				ymin = 0; ymax = da.ylab.length;
				break;
		
			case "CompMatrix": case "CompMatrixAnim":
				xmin = 0; xmax = 0;
				ymin = 0; ymax = 0;
				break;
				
			case "Table":
				break;
			
			case "VertLine": case "VertLine2": case "VertDashLine": 
				{
					let x = da.x;
					if(x < xmin) xmin = x; if(x > xmax) xmax = x;
				}
				break;
			
			case "PriorLine":
				break;
				
			case "BurninLine":
				for(let i = 0; i < da.point.length; i++){
					let po = da.point[i];		
					if(po.x < xmin) xmin = po.x; if(po.x > xmax) xmax = po.x;
				}
				break;
			
			case "SimX":
				{
					let x = da.x; if(x < xmin) xmin = x; if(x > xmax) xmax = x;
				}
				break;
				
			case "SimY":
				{
					let y = da.y; if(y < ymin) ymin = y; if(y > ymax) ymax = y;
				}
				break;
				
			case "Comp Video":
				{
					let y_vec = da.y_vec;
					for(let i = 0; i < y_vec.length; i++){
						let y = y_vec[i];
						if(y > ymax) ymax = y;
						if(y < ymin) ymin = y;
					}
				}
				break;
				
			default:
				for(let i = 0; i < da.point.length; i++){
					let po = da.point[i];		
					if(po.x < xmin) xmin = po.x; if(po.x > xmax) xmax = po.x;
					if(po.y < ymin) ymin = po.y; if(po.y > ymax) ymax = po.y;
					if(po.CImin != undefined){ if(po.CImin < ymin) ymin = po.CImin;}
					if(po.CImax != undefined){ if(po.CImax > ymax) ymax = po.CImax;}
				}	
				break;
			}
		}
	
		if(xmax < xmin+TINY) xmax = xmin;
		
		if(xmin == xmax){                              // Ensures thay xmin and xmax are not the same
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
		
			
		switch(this.variety){
		case "Histogram": case "HistoAnim":                   // On histogram place zero on one end
			if(!(ymin <= 0 && ymax >= 0)){
				if(ymin > 0) ymin = 0;
				else ymax = 0;
			}
			break;
		}
		
		if(ymax < ymin+TINY) ymax = ymin;
		
		if(ymin == ymax){                              // Ensures thay ymin and ymax are not the same
			if(ymax > 0){ ymin = 0; ymax *= 1.1;}
			else{
				if(ymin < 0){ ymax = 0; ymin *= 1.1;}
				else ymax = 1;
			}
		}
		else{
			if(this.variety != "Scatter"){
				// Puts lower bound at zero
				if(ymin > 0 && ymax > 0 && ymin < 0.8*ymax) ymin = 0; 
			}
		}
	
		if(this.op.yaxis == false) ymin = 0;
		
		if(this.type == "HistoAnim"){ xmin = 0; xmax = this.data.length;}
		
		switch(this.variety){
		case "TransBias": ymin = -1; ymax = 1; break;
		//case "TransP": ymin = 0; break;
		}
		
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
		if(this.variety == "Matrix" || this.variety == "MatrixAnim") return;
		
		let ra = this.range;
		if(ra.xmin != LARGE && this.variety != "Histogram" && this.type != "HistoAnim" ){
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
	
	
	/// For comp matrix initilises lines used to represent matrix
	comp_matrix_init()
	{
		let op = this.op;
		let comp = this.get_cla(op.p,op.cl).comp;
		
		let N = comp.length;
		let p=[];
		for(let i = 0; i < N; i++) p.push(model.comp_center(comp[i]));
		
		let val = this.data[0].value;
		if(val.length != N){ error("Value wrong size"); return;}
		if(val[0].length != N){ error("Value wrong size"); return;}
		
		// Determines if matrix is symetric
		let sym = true;
		let max = 0;
		for(let j = 0; j < N; j++){
			for(let i = 0; i < N; i++){
				if(j != i){
					if(this.variety == "CompMatrixAnim"){
						let vec = val[j][i];
					
						for(let k = 0; k < vec.length; k++){
							let va = Number(vec[k]);
							if(va != val[i][j][k]) sym = false; 
					
							if(va > 0){ if(va > max) max = va;}
							else{ if(-va > max) max = -va;}	
						}
					}
					else{
						let va = Number(val[j][i]);
						if(va != val[i][j]) sym = false; 
					
						if(va > 0){ if(va > max) max = va;}
						else{ if(-va > max) max = -va;}	
					}
				}					
			}
		}
		
		let box = get_model_box(op.p,op.cl); 
	
		let mar1 = 0.1, mar2 = 1-mar1;
	
		let loopmax = 1;
		if(this.variety == "CompMatrixAnim") loopmax = val[0][0].length;
	
		let list = [];
		
		let fac = Math.exp(2*inter.compmatrix_slider.value);
	
		for(let loop = 0; loop < loopmax; loop++){
			let clink = [];
		
			if(sym == true){
				for(let j = 0; j < N; j++){
					let x1 = p[j].x, y1 = p[j].y;
					for(let i = j+1; i < N; i++){			
						let va = Number(val[j][i]);
						if(this.variety == "CompMatrixAnim") va = va[loop];
						
						if(va != 0){
							let w = va/max;
							
							let thick = fac*w; if(thick < 0) thick = -thick;
							if(thick > 0.01){	
								let alpha = w; if(alpha < 0) alpha = -w;
							
								thick *= 0.5/box.scale;
								let x2 = p[i].x, y2 = p[i].y;
								let dx = x2-x1, dy = y2-y1;
								clink.push({alpha:0.1+0.9*alpha, val:va, thick:thick, x1:x1+mar1*dx, y1:y1+mar1*dy, x2:x1+mar2*dx, y2:y1+mar2*dy});
							}
						}
					}
				}
			}
			else{
				for(let j = 0; j < N; j++){
					let x1 = p[j].x, y1 = p[j].y;
					for(let i = 0; i < N; i++){
						let va = val[j][i];
						if(this.variety == "CompMatrixAnim") va = va[loop];
					
						if(j != i && va != 0){
							let w = va/max;
							let thick = fac*w; if(thick < 0) thick = -thick;
							
							let va_rev = val[i][j];
							if(this.variety == "CompMatrixAnim") va_rev = va_rev[loop];
					
							let w_rev = va_rev/max;
							let thick_rev = fac*w_rev; if(thick_rev < 0) thick_rev = -thick_rev;
							
							if(thick > 0.01){	
								let alpha = w; if(alpha < 0) alpha = -w;
							
								thick *= 0.25/box.scale;
								thick_rev *= 0.25/box.scale;
								
								let x2 = p[i].x, y2 = p[i].y;
								let dx = x2-x1, dy = y2-y1;
								let nx = -dy, ny = dx;
								let rr = Math.sqrt(nx*nx + ny*ny);
								
								let mar3 = mar2 - 3*thick/rr;
								let mar4 = mar3+0.01;
							
								if(mar3 < mar1) mar3 = mar1;
								 
								nx *= 0.8*thick/rr; ny *= 0.8*thick/rr;
								
								let shx = nx, shy = ny;
								if(thick_rev > thick){
									shx *= thick_rev/thick; shy *= thick_rev/thick; 
								}
							
								clink.push({alpha:0.1+0.9*alpha, val:va, thick:thick, arrow:true, x1:x1+mar1*dx+shx, y1:y1+mar1*dy+shy, x2:x1+mar4*dx+shx, y2:y1+mar4*dy+shy, x3:x1+mar2*dx+shx, y3:y1+mar2*dy+shy, x4:x1+mar3*dx+shx-nx, y4:y1+mar3*dy+shy-ny, x5:x1+mar3*dx+shx+nx, y5:y1+mar3*dy+shy+ny});
							}
						}
					}
				}
			}
			list.push(clink);
		}
	
		this.complink = list;
	}
		
		
	/// Initiailises a compartmental vector plot
	comp_colour_init()
	{	
		if(this.data.length == 0) return;
	
		let ckey = this.colour_key;
		let min = ckey.min;
		let max = ckey.max;
		let grad = ckey.grad;
		
		if(this.data[0].point){
			for(let k = 0; k < this.data.length; k++){
				let dat = this.data[k];
				for(let j = 0; j < dat.point.length; j++){
					dat.point[j].col = this.get_matrix_col(dat.point[j].y);
				}
			}
		}
		
		{
			let f = COLOUR_KEY_DIV/(max-min);
			for(let k = 0; k < this.data.length; k++){
				let dat = this.data[k];
				let col_k = [];
				if(dat.y_vec){
					let y_vec = dat.y_vec;	
					for(let j = 0; j < y_vec.length; j++){
						let k = Math.floor(f*(y_vec[j]-min));
						if(k < 0 || k > COLOUR_KEY_DIV) error("out of range");
						col_k.push(k);
					}
				}
				if(dat.point){
					let point = dat.point;	
					for(let j = 0; j < point.length; j++){
						let k = Math.floor(f*(point[j].y-min));
						if(k < 0 || k > COLOUR_KEY_DIV) error("out of range");
						col_k.push(k);
					}
				}
				dat.col_k = col_k;
			}
		}
	}
	
		
	/// Draws compartmental links
	draw_complink(lay)
	{
		let p = this.op.p;
		let cl = this.op.cl;
	
		let claa = this.get_cla(p,cl);
		
		let cam = claa.camera;
		
		let fr = this.animation.playframe;
		if(this.variety == "CompMatrix") fr = 0;
	
		for(let i = 0; i < this.complink[fr].length; i++){
			let clink = this.complink[fr][i];
			let p1 = trans_point(clink.x1,clink.y1,cam,lay);
			let p2 = trans_point(clink.x2,clink.y2,cam,lay);
		
			let thick = cam.scale*clink.thick*inter.sca;
			let col = BLUE; if(clink.val < 0) col = RED; 
		
			cv.globalAlpha = clink.alpha;
			draw_line(p1.x,p1.y,p2.x,p2.y,col,thick);
			if(clink.arrow == true){
				let polypoint=[];    
				polypoint.push(trans_point(clink.x3,clink.y3,cam,lay));
				polypoint.push(trans_point(clink.x4,clink.y4,cam,lay));
				polypoint.push(trans_point(clink.x5,clink.y5,cam,lay));
				draw_polygon(polypoint,col,col,1);
			}
		}
		cv.globalAlpha = 1;
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
				if(this.variety == "MatrixAnim"){
					for(let k = 0; k < mat[j][i].length; k++){
						let val = Number(mat[j][i][k]);
						if(val < min) min = val; if(val > max) max = val;
					}
				}
				else{
					let val = Number(mat[j][i]);
					if(val < min) min = val; if(val > max) max = val;
				}
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
		
		this.set_colour_key(min,max);
		
		for(let j = 0; j < mat.length; j++){
			for(let i = 0; i < mat[j].length; i++){
				if(this.variety == "MatrixAnim"){
					for(let k = 0; k < mat[j][i].length; k++){
						mat_col[j][i][k] = this.get_matrix_col(Number(mat[j][i][k]));
					}
				}
				else{
					mat_col[j][i] = this.get_matrix_col(Number(mat[j][i]));
				}
			}
		}
	}
	

	/// Sets the colour key based on a minimum and maximum
	set_colour_key(min,max)
	{
		this.mag = max; if(-min > this.mag) this.mag = -min;
	 
		let grad = [], grad_rgb = [];
		let imax = COLOUR_KEY_DIV;
		for(let i = 0; i <= imax; i++){
			let col = this.get_matrix_col_rgb(min+(max-min)*i/imax);
			grad_rgb.push(col);
			grad.push("#"+hex(col.r)+hex(col.g)+hex(col.b));
		}
		
		this.colour_key = {min:min, max:max, tick:this.get_tick_marks(min,max,3), grad:grad, grad_rgb:grad_rgb};
	}
	
	
	/// Gets a matrix colour from a value
	get_matrix_col(val)
	{
		let col = this.get_matrix_col_rgb(val);
		return "#"+hex(col.r)+hex(col.g)+hex(col.b);
	}
	
	
	/// Gets a matrix colour from a value uisng rgb
	get_matrix_col_rgb(val)
	{
		let ma = 240;
		
		if(val > 0){
			let f = val/this.mag;
			return {r:Math.floor(ma*(1-f)), g:Math.floor(ma*(1-f)), b:ma};
		}
		else{
			let f = -val/this.mag;
			return {r:ma, g:Math.floor(ma*(1-f)), b:ma*(1-f)};
		}
	}
	
				
	/// Initialises the view of the comparments
	anim_init()
	{
		if(this.variety == "CompVector") return;
		
		//let species = this.op.species;
		let anim = this.animation;
		
		anim.playframe = 0;
		anim.playing = false;
		
		let k = 0;
		while(k < this.data.length && this.data[k].type == "VertLine") k++;
		if(k == this.data.length){ error("Could not find data"); return;}
				
		let da = this.data[k];
		
		switch(this.variety){
		case "MatrixAnim": anim.playframe_max = da.mat[0][0].length-1; break;
		case "CompMatrixAnim": anim.playframe_max = da.value[0][0].length-1; break;
		case "Population":	
			if(da.y_vec) anim.playframe_max = da.y_vec.length-1;
			else anim.playframe_max = da.point.length-1;
			break;
		default: anim.playframe_max = da.point.length-1; break;
		}
	}
	
	
	/// Initialises the density view of the comparments
	density_init()
	{
		let p = this.op.p;
		let cl = this.op.cl;
	
		let claa = this.get_cla(p,cl);
		
		let po=[];
		let xmax = -LARGE, xmin = LARGE, ymax = -LARGE, ymin = LARGE; 
		for(let i = 0; i < claa.comp.length; i++){
			let co = claa.comp[i];
			let p = model.comp_center(co);
			let x = p.x, y = p.y;
			if(x > xmax) xmax = x; if(x < xmin) xmin = x;
			if(y > ymax) ymax = y; if(y < ymin) ymin = y;
			
			po.push({x:x, y:y});
		}

		if(xmin == xmax && ymin == ymax){
			xmin -= 10; xmax += 10;
			ymin -= 10; ymax += 10;
		}
		
		let dx = xmax-xmin, dy = ymax-ymin;
		let l = dx; if(dy > dx) l = dy;
	
		let info = inter.density_slider;
		
		let va = 1.0/Math.sqrt(claa.comp.length)
		if(va > 0.1) va = 0.1;
		let frac = va*Math.exp(info.value);
		let r = frac*l;
		xmin -= r; xmax += r;
		ymin -= r; ymax += r;
		dx += 2*r; dy += 2*r;
		l += 2*r;
		
		let rule = l/DEN_X;
		let DX = Math.round(dx/rule);
		let DY = Math.round(dy/rule);
		let sr = r/rule;
		if(sr < SR_MIN){
			rule = (l/DEN_X)*(sr/SR_MIN);
			DX = Math.round(dx/rule);
			DY = Math.round(dy/rule);
			sr = r/rule;
		}
		
		let mat=[], mat_dist=[];
		for(let j = 0; j < DY; j++){
			mat[j]=[]; mat_dist[j]=[];
		}
		
		for(let k = 0; k < po.length; k++){
			let sx = (po[k].x-xmin)/rule;
			let dxmin = Math.floor(sx-sr);
			let dxmax = Math.floor(1+sx+sr);
			
			//let sy = (po[k].y-ymin)/rule;
			let sy = (ymax-po[k].y)/rule;
			let dymin = Math.floor(sy-sr);
			let dymax = Math.floor(1+sy+sr);
			
			for(let j = dymin; j < dymax; j++){
				if(j >= 0 && j < DY){
					for(let i = dxmin; i < dxmax; i++){
						if(i >= 0 && i < DX){
							let ddx = (i+0.5)-sx, ddy = (j+0.5)-sy;
							let dist = ddx*ddx+ddy*ddy;
							if(dist < sr*sr){
								let val = mat_dist[j][i];
								if(val == undefined){ mat[j][i] = k; mat_dist[j][i] = dist;}
								else{
									if(dist < val){ mat[j][i] = k; mat_dist[j][i] = dist;}
								}
							}
						}
					}
				}
			}
		}

		this.density_info = { x:xmin, y:ymax, dx:dx, dy:dy, mat:mat, DX:DX, DY:DY};
	}
	
	
	/// Activates when the play button is pressed
	press_play()
	{
		let anim = this.animation;
		
		if(anim.playing == false){
			anim.playing = true;
			if(anim.playframe == anim.playframe_max) anim.playframe = 0;
			anim.playframe_start = anim.playframe;
			anim.clock = clock();
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
		let anim = inter.graph.animation;
			
		let x = inter.mx-bu.x-lay.x;

		let frac = (x-bu.dy/2)/(bu.dx-bu.dy);
		
		anim.playing = false;
		anim.playframe = Math.round(frac*anim.playframe_max);
		if(anim.playframe < 0) anim.playframe = 0;
		if(anim.playframe > anim.playframe_max) anim.playframe = anim.playframe_max;
	}
	
	
	/// Activates when the slider is pressed
	press_slider(bu,update)
	{
		let info = bu.info;
	
		let lay = get_lay(info.lay);
		let anim = inter.graph.animation;
			
		let x = inter.mx-bu.x-lay.x;

		let frac = (x-slider_dx/2)/(bu.dx-slider_dx);
		if(frac < 0) frac = 0; if(frac > 1) frac = 1;
		
		info.value = info.min+frac*(info.max-info.min);
		
		if(update == true){
			switch(bu.info.update){
			case "density": inter.graph.density_init(); break;
			case "compmatrix": this.comp_matrix_init(); break;
			case "scale": model.update_pline(this.get_cla(bu.p,bu.cl)); break;
			default: error("update not recognised"); break;
			}
		}
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
			let frame_new = Math.round(anim.playframe_start + (anim.speed.value*(clock()-anim.clock)/4000)*anim.playframe_max);
			if(frame_new > anim.playframe_max) frame_new = anim.playframe_max;
			
			if(frame_new != anim.playframe){
				anim.playframe = frame_new;
				if(anim.playframe == anim.playframe_max) anim.playing = false;
				
				this.replot_anim();
			}
		
			if(anim.playframe < anim.playframe_max){
				setTimeout(function(){ inter.graph.playanim();}, 10);
			}
		}
	}
	
	
	/// Replots the animation
	replot_anim()
	{
		if(layer_exist("GraphCompartments")) replot_layer("GraphCompartments");
		if(layer_exist("GraphContent")) replot_layer("GraphContent");
		
		replot_layer("AnimControls");
		plot_screen();
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
		if(graph_dia) error("GRAPH DIA content: "+this.variety);
		
		let vari = plot_variety(this.type);
		
		switch(this.variety){
		case "PhyloTree":
			{
				let left_shift = lay.op.left_shift;
	
				let dy = 2*this.barh;
				
				let pnode = this.data[0].pnode;
				let ymax = 0;
				for(let i = 0; i < pnode.length; i++){ 	
					if(pnode[i].y > ymax) ymax = pnode[i].y;
				}
				let yend = (ymax+1)*dy;
				
				lay.add_button({x:left_shift, y:0, dx:lay.dx-scrollw-left_shift, dy:yend, barh:dy, pnode:pnode, type:"PhyloButton"});
				
				lay.add_button({x:0, y:0, dx:lay.dx-scrollw, dy:yend, ac:"TimelineGrab", type:"Nothing"});
			}
			break;
			
		case "TransTree":
			{
				let left_shift = lay.op.left_shift;
				let y = 0, dy = this.barh, gap = 0; 
				for(let i = 0; i < this.data.length; i++){
					let da = this.data[i];
					if(da.type ==  "InfBar"){
						lay.add_button({x:0, y:y, dx:lay.dx-scrollw, dy:dy, name:da.name, tmin:da.tmin, tmax:da.tmax, inf_line:da.inf_line, left_shift:left_shift, type:"InfBarBut"});
						y += dy+gap;
					}
				}
				
				let da = this.data[this.data.length-1];
				if(da.type != "Transmission") error("Should be transmission");
			
				lay.add_button({x:left_shift, y:0, dx:lay.dx-scrollw-left_shift, dy:y, barh:dy, transmission:da.transmission, type:"TransArrow"});
				
				lay.add_button({x:0, y:0, dx:lay.dx-scrollw, dy:y, ac:"TimelineGrab", type:"Nothing"});
			}
			break;
			
		case "Individual":
			{
				let y = 0, dy = 2.5*this.barh, gap = 0.4;
				for(let i = 0; i < this.data.length; i++){
					let da = this.data[i];
				
					lay.add_button({x:0, y:y, dx:lay.dx-scrollw, dy:dy, name:da.name, info:da.info, col_timeline:da.col_timeline, obs:da.obs, type:"IndTimelineBut"});
					
					y += dy+gap;
				}
				lay.add_button({x:0, y:0, dx:lay.dx-scrollw, dy:y, ac:"TimelineGrab", type:"Nothing"});
			
				{ // Adds links to individuals
 					let si = dy/3;
					let fo = get_font(si);
			
					let y = 0;
					for(let i = 0; i < this.data.length; i++){
						let da = this.data[i];
					
						lay.add_button({te:da.name, x:0.3, y:y, dx:text_width(da.name,fo)+1, dy:si, type:"Link", ac:"SelectInd", si:si, font:fo});
					
						y += dy+gap;
					}
				}
			}
			break;
			
		case "Histogram": case "HistoAnim":
			lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, type:"GraphContentBut"});
			break;
			
		case "Matrix": case "MatrixAnim":
			{
				let da = this.data[0];
			
				let fr = this.animation.playframe;
				
				let labels = {ylab:da.ylab, xlab:da.xlab};
				
				for(let j = 0; j < da.ylab.length; j++){
					let fracj = j/da.ylab.length;
					for(let i = 0; i < da.xlab.length; i++){	
						let fraci = i/da.xlab.length; 
						
						let dx = lay.dx/da.xlab.length;
						let dy = lay.dy/da.ylab.length;
						let marx = 0.03*dx;
						let mary = 0.03*dy;
						
						let stat; if(da.mat_stat) stat = da.mat_stat[j][i];
					
						let col, value, CImin, CImax, te;
						if(this.variety == "MatrixAnim"){
							col = da.mat_col[j][i][fr];
							value = precision(da.mat[j][i][fr]);
							if(da.CImin) CImin = precision(da.CImin[j][i][fr]);
							if(da.CImax) CImax = precision(da.CImax[j][i][fr]);
							te = da.mat[j][i][fr];
						}
						else{
							col = da.mat_col[j][i];
							value = precision(da.mat[j][i]);
							if(da.CImin) CImin = precision(da.CImin[j][i]);
							if(da.CImax) CImax = precision(da.CImax[j][i]);
							te = da.mat[j][i];
						}							
						
						lay.add_button({te:precision(te), x:fraci*lay.dx+marx, y:(1-fracj)*lay.dy-dy+mary, dx:dx-2*marx, dy:dy-2*mary, value:value, CImin:CImin, CImax:CImax, i:i, j:j, labels:labels, stat:stat, col:col, type:"MatrixEleBut", ac:"MatrixEleBut"});
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
		if(graph_dia){
			error("GRAPH DIA draw content button: data:");
			error(this.data);
		}
		
		fill_rectangle(x,y,dx,dy,WHITE); 
		if(this.extra_flag) dx -= GRAPH_EXTRA;
		
		let ra = this.range;
	
		for(let li = 0; li < this.data.length; li++){
			let da = this.data[li];
		
			switch(da.type){
			case "Line": case "BurninLine": this.draw_graph_line(x,y,dx,dy,ra,da); break;
			case "PriorLine": this.draw_graph_line(x,y,dx,dy,ra,da); break;
			case "Line CI": this.draw_graph_CI(x,y,dx,dy,ra,da); this.draw_graph_line(x,y,dx,dy,ra,da); break;
			case "Bar": this.draw_bar(x,y,dx,dy,ra,da); break;
			case "SplineAnim": this.draw_splinebar(x,y,dx,dy,ra,da,li); break;
			case "ErrorBar": this.draw_errorbar(x,y,dx,dy,ra,da); break;
			case "Cross": this.draw_data_cross(x,y,dx,dy,ra,da); break;
			case "Distribution": this.draw_distribution(x,y,dx,dy,ra,da); break;
			case "VertLine": this.draw_vert_line(x,y,dx,dy,ra,da); break;
			case "VertLine2": this.draw_vert_line2(x,y,dx,dy,ra,da); break;
			case "VertDashLine": this.draw_vert_dash_line(x,y,dx,dy,ra,da); break;
			case "Points": this.draw_point(x,y,dx,dy,ra,da); break;
			case "SimX": this.draw_vert(x,y,dx,dy,ra,da); break;
			case "SimY": this.draw_hor(x,y,dx,dy,ra,da); break;
			default: error("Option not recongnised: "+da.type); break;
			}
		}

		if(this.bf_store) this.draw_bf(x,y,dx,dy,ra);
			
		let title = this.op.title;
	
		if(title != undefined){
			center_text(title,x+dx/2,y+1,get_font(1,"bold"),BLACK,dx-1); 
		}
	}
	
	
	/// Draws an individual timeline
	draw_ind_timeline_button(x,y,dx,dy,name,info,col_timeline,obs)
	{
		let bar = 1;
		let si = dy/3.5;
		let yte = y+dy/2-0.15*dy;
		
		right_text(info,x+dx-0.2,yte,get_font(si),BLACK,0.6*dx); 
		
		let ra = this.range;
				
		if(col_timeline.length == 0){
			fill_rectangle(x,y+dy/2,dx,dy/2,"#eeeeee");
			center_text("Unobserved",x+dx/2,y+0.87*dy,get_font(0.8),BLACK,dx);
		}
		else{	
			let y1 = ro_down(y+dy/2), y2 = ro_down(y+dy);
		
			let pos = [];
			for(let i = 0; i < col_timeline.length; i++){
				pos.push(ro_down(x+dx*((col_timeline[i].t-ra.xmin)/(ra.xmax-ra.xmin))));
			}
		
			for(let i = 0; i < col_timeline.length-1; i++){
				cv.globalAlpha = col_timeline[i].alpha;
				cv.beginPath();
				cv.rect(pos[i],y1,pos[i+1]-pos[i],y2-y1);
				cv.fillStyle = "rgb("+col_timeline[i].col+")";
				cv.fill();
			}
			cv.globalAlpha = 1;
		
			let fac = dy/2.5;
			let r = 0.32*fac, si = 0.34*fac, si2 = 0.4*fac;
			let yy = y+0.74*dy;
			
			for(let i = 0; i < obs.length; i++){
				let ob = obs[i];
				let xx = x+dx*((ob.t-ra.xmin)/(ra.xmax-ra.xmin));
	
				switch(ob.type){
				case "AddObs":
					{
						let siw = 0.5, sih = 0.7;
						multi_triangle(xx,yy-sih/2,siw,sih,ob.col_list,black_not_list(ob.col_list),NORMLINE);
					}
					break;
				
				case "RemObs":
					{
						let siw = 0.5, sih = 0.7;
						rem_triangle(xx-siw,yy-sih/2,siw,sih,BLACK);
					}
					break;
				
				case "MoveObs":
					if(ob.col_list.length == 1){
						fill_semicircle(xx,yy,r,ob.col_list[0].col,black_not_list(ob.col_list),NORMLINE); 
					}
					break;
					
				case "TransObs":
					if(ob.col.length == 2){
						draw_trans_obs(xx-si2/2,yy-si2/2,si2,si2,ob.col[0],ob.col[1]);
					}
					else{
						draw_line(xx,yy-si2/2,xx,yy+si2/2,BLACK,NORMLINE);
						draw_rectangle(xx-si2/2,yy-si2/2,si2,si2,BLACK,NORMLINE); 
					}
					break;
					
				case "CompObs":
					if(ob.col_list.length == 1){
						let co = ob.col_list[0].col;
						if(co == "notalive"){
							fill_circle_cross(xx,yy,r,WHITE,BLACK,NORMLINE); 
						}
						else{
							fill_circle(xx,yy,r,co,black_not_list(ob.col_list),NORMLINE);
						}
					}						
					else multi_circle(xx,yy,r,ob.col_list,black_not_list(ob.col_list),NORMLINE); 
					break;
				
				case "DiagObs":
					if(ob.res == true) draw_rect(xx-si/2,yy-si/2,si,si,DDGREY,BLACK,NORMLINE); 
					else draw_rect(xx-si/2,yy-si/2,si,si,WHITE,BLACK,NORMLINE);  
					break;
				
				case "GeneticObs":
					draw_X(xx,yy,1.2*si);
					break;
				}
			}
		}
	}
	
	
	/// Draws infection timeline (used for trans tree)
	draw_ind_infbar_button(x,y,dx,dy,name,inf_line,tmin,tmax,left_shift)
	{
		fill_rectangle(x,y,dx,dy,WHITE);
	
		let bar = 0.5;
		
		let ra = this.range;
			
		let dt = (tmax-tmin)/inf_line.length;
	
		let pos = [];
		for(let i = 0; i <= inf_line.length; i++){
			pos.push(ro_down(x+left_shift+(dx-left_shift)*((tmin+i*dt-ra.xmin)/(ra.xmax-ra.xmin))));
		}

		let y1 = ro_up(y+0.3*dy), y2 = ro_down(y+0.7*dy);

		cv.fillStyle = RED;
		for(let i = 0; i < inf_line.length; i++){
			cv.globalAlpha = inf_line[i];
			cv.beginPath();
			cv.rect(pos[i],y1,pos[i+1]-pos[i],y2-y1);
			cv.fill();
		}
		cv.globalAlpha = 1;
		
		fill_rectangle(x,y,left_shift,dy,WHITE);
		
		plot_text(name,x+0.2,y+dy/2+0.3,get_font(si_transtree),BLACK); 
		
		let xx = x+left_shift;
		draw_line(xx,y,xx,y+dy,BLACK,NORMLINE);
	}


	/// Draws arrows to indicate transmission of infection
	draw_transarrow_button(x,y,dx,dy,barh,transmission)
	{
		let ra = this.range;
		
		let fac = dx/(ra.xmax-ra.xmin);
		for(let j = 0; j < transmission.length; j++){
			let tr = transmission[j];
			
			let xx = x+fac*(tr.t-ra.xmin);
			
			if(xx > x && xx < x+dx){
				{
					let r = 0.25*barh*Math.sqrt(tr.frac);
					let yy = y+barh*(tr.k+0.5);
					fill_circle(xx,yy,r,RED,RED);  
					
					r = 0.20*barh*Math.sqrt(tr.frac_out);
					fill_circle(xx,yy,r,WHITE,WHITE);  
				}
				
				if(tr.min != undefined){
					let y1 = y+barh*(tr.min+0.5), y2 = y+barh*(tr.max+0.5);
					if((y1 > y && y1 < y+dy) || (y2 > y && y2 < y+dy) || 
						(y1 < y || y2 > y+dy) || (y2 < y || y1 > y+dy)){
							
						draw_line(xx,y1,xx,y2,BLACK,NORMLINE);
										
						for(let j = 0; j < tr.inf_list.length; j++){
							let ilist = tr.inf_list[j];
						
							let r = 0.25*barh*Math.sqrt(ilist.frac);
							let yy = y+barh*(ilist.k+0.5);
							fill_circle(xx,yy,r,BLACK,BLACK);  
						}
					}
				}
			}				
		}
	}		
	
	
	/// Converts number to reduce precision
	num_text(num)
	{
		if(num > 5) return Math.round(num);
		else{
			return Math.round(num*10)/10;
		}
	}
	
	
	/// Draws phylogentic tree
	draw_phylo_button(x,y,dx,dy,barh,pnode)
	{
		let ra = this.range;
		
		let fac = dx/(ra.xmax-ra.xmin);
		
		// Draw lines
		for(let i = 0; i < pnode.length; i++){
			let po = pnode[i];
			
			let xx = x+fac*(po.t-ra.xmin);
			let xx2;
			if(po.branch != undefined) xx2 = x+fac*(po.branch.t-ra.xmin);
			else{
				if(po.obs.length > 0) xx2 = x+fac*(po.obs[po.obs.length-1].t-ra.xmin);
			}
			
			let yy = y+barh*po.y;

			if(po.ty == "ROOT"){
				fill_circle(xx,yy,0.05*barh,BLACK,BLACK);  
			}

			draw_line(xx,yy,xx2,yy,BLACK,NORMLINE);
			
			if(po.branch != undefined){
				let br = po.branch.br;
				switch(br.length){
				case 1:
					{
						let d = 0.1*barh;
						draw_line(xx2,yy-d,xx2,yy+d,BLACK,NORMLINE);
					}
					break;
					
				case 2:
					{
						let yy1 = y+barh*pnode[br[0]].y;
						let yy2 = y+barh*pnode[br[1]].y;
						draw_line(xx2,yy1,xx2,yy2,BLACK,NORMLINE);
					}
					break;
					
				default: error("wrong num"); break;
				}
			}
		}
		
		// Draws text
		let si; if(barh > 1) si = Math.sqrt(barh); else si = barh;
		
		let fo = get_font(0.4*si);
		let fo2 = get_font(0.4*si);
		
		for(let i = 0; i < pnode.length; i++){
			let po = pnode[i];
			let xx = x+fac*(po.t-ra.xmin);
			let yy = y+barh*po.y;
		
			{
				let xx2;
				if(po.branch != undefined) xx2 = x+fac*(po.branch.t-ra.xmin);
				else{
					if(po.obs.length > 0) xx2 = x+fac*(po.obs[po.obs.length-1].t-ra.xmin);
				}
				
				let dx = xx2-xx;
				if(dx > 2){
					center_text(po.label,(xx+xx2)/2,yy-0.2*si,fo2,BLACK,dx-1); 
				}
			}
			
			if(po.ty == "ROOT"){
				right_text(this.num_text(po.mnum),xx-barh*0.1,yy+0.2*si,fo,BLUE); 
			}
			
			for(let j = 0; j < po.obs.length; j++){
				let ob = po.obs[j];
				let xx2 = x+fac*(ob.t-ra.xmin);
				let dx = xx2-xx;
				if(dx > 2){
					center_text(this.num_text(ob.mnum),(xx+xx2)/2,yy+0.4*si,fo2,BLUE,dx-1); 
				}
			
				xx = xx2;
			}
			
			if(po.branch != undefined){
				let xx2 = x+fac*(po.branch.t-ra.xmin);
		
				let dx = xx2-xx;
				if(dx > 2){	
					center_text(this.num_text(po.branch.mnum),(xx+xx2)/2,yy+0.4*si,fo2,BLUE,dx-1); 
				}
				
				let te = Math.floor(po.branch.frac*100)+"%";
				
				plot_text(te,xx2+barh*0.1,yy+0.12*si,fo,DGREEN); 
			}
		}
		
		// Draw observations
		for(let i = 0; i < pnode.length; i++){
			let po = pnode[i];
			let yy = y+barh*po.y;
			for(let j = 0; j < po.obs.length; j++){
				let xx = x+fac*(po.obs[j].t-ra.xmin);
				draw_X(xx,yy,0.2*si);	
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
		
		let w = 2;
		if(da.thick) w = da.thick;
		if(da.view == "Graph (lines)") w = 1;
		if(inter.printing) w *= print_line_factor;
		cv.lineWidth = w;
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
	
	
	/// Draws a vertical line
	draw_vert(x,y,dx,dy,ra,da)
	{
		let pox = []; pox.push({x:da.x,y:ra.ymin}); pox.push({x:da.x,y:ra.ymax});
		let datx = {point:pox, col:da.col, type:"Line", dash:da.dash};
		
		this.draw_graph_line(x,y,dx,dy,ra,datx);
	}
			
		
	/// Draws a horizontal line
	draw_hor(x,y,dx,dy,ra,da)
	{	
		let poy = []; poy.push({x:ra.xmin, y:da.y}); poy.push({x:ra.xmax,y:da.y});
		let daty = {point:poy, col:da.col, type:"Line", dash:da.dash};
		
		this.draw_graph_line(x,y,dx,dy,ra,daty);
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
			cv.fillStyle = da.col;
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
			cv.fillStyle = da.col;
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
		let w = MEDIUMLINE;
		if(inter.printing) w *= print_line_factor;
		cv.lineWidth = w;
		cv.strokeStyle = da.col;
		cv.stroke();	
	}
	
	
	
	/// Draws Bayes factor values
	draw_bf(x,y,dx,dy,ra)
	{
		let bf = this.bf_store;
		
		let xx = x+dx*((bf.val-ra.xmin)/(ra.xmax-ra.xmin));
		let y_post = y+dy-dy*((bf.post_int-ra.ymin)/(ra.ymax-ra.ymin));
		let y_prior = y+dy-dy*((bf.prior_int-ra.ymin)/(ra.ymax-ra.ymin));
		let y_base = y+dy-dy*((0-ra.ymin)/(ra.ymax-ra.ymin));
		let r = 0.3;
		let yy = y_post; if(y_prior < yy) yy = y_prior;
		
		draw_line(xx,y_base,xx,yy,BLACK,THICKLINE);    
		fill_circle(xx,y_post,r,RED,RED);  
		fill_circle(xx,y_prior,r,BLUE,BLUE); 
	}
	
	
	/// Draws a verticle line
	draw_vert_line(x,y,dx,dy,ra,da)
	{
		let xx = x+dx*((da.x-ra.xmin)/(ra.xmax-ra.xmin));
		let xp = ro(xx);
		cv.beginPath();
		
		cv.moveTo(xp,ro(y));
		cv.lineTo(xp,ro(y+dy));
		
		let w = da.thick; if(inter.printing) w *= print_line_factor;
		cv.lineWidth = w;
		cv.strokeStyle = da.col;
		cv.stroke();	
		
		let fo = get_font(0.8);
		let wid = text_width(da.te,fo);
		vert_text(da.te,xx-0.3,y+wid+0.2,fo,da.col,10); 
	}
	
	/// Draws a verticle line with text the other side
	draw_vert_line2(x,y,dx,dy,ra,da)
	{
		let xx = x+dx*((da.x-ra.xmin)/(ra.xmax-ra.xmin));
		let xp = ro(xx);
		cv.beginPath();
		
		cv.moveTo(xp,ro(y));
		cv.lineTo(xp,ro(y+dy));
		let w = da.thick; if(inter.printing) w *= print_line_factor;
		cv.lineWidth = w;
		cv.strokeStyle = da.col;
		cv.stroke();	
		
		let fo = get_font(0.8);
		let wid = text_width(da.te,fo);
		vert_text(da.te,xx+0.9,y+wid+0.2,fo,da.col,10); 
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

		if(y1 < y2){ let yt = y1; y1 = y2; y2 = yt;}

		fill_rectangle(x1,y1,x2-x1,y2-y1,da.col,dark_colour(da.col),NORMLINE);
	}
	
	
	/// Draws a bar on a histogram
	draw_splinebar(x,y,dx,dy,ra,da,i)
	{
		let fr = this.animation.playframe;
		let po = da.point[fr];
	
		let x1 = x+dx*((i+(1-bar_thick)/2-ra.xmin)/(ra.xmax-ra.xmin));
		let x2 = x+dx*((i+1-(1-bar_thick)/2-ra.xmin)/(ra.xmax-ra.xmin));
		
		let y1 = y+dy-dy*((po.y-ra.ymin)/(ra.ymax-ra.ymin));
		let y2 = y+dy-dy*((0-ra.ymin)/(ra.ymax-ra.ymin));
			
		let col = po.col;
		fill_rectangle(x1,y1,x2-x1,y2-y1,col);
		draw_rectangle(x1,y1,x2-x1,y2-y1,dark_colour(col),NORMLINE);
		
		if(po.CImin != undefined){
			let da2 = { x:i+0.5, y:po.y, ymin:po.CImin, ymax:po.CImax, col:BLACK};
			this.draw_errorbar(x,y,dx,dy,ra,da2);
		}
	}
	
	
	/// Draws an errorbar 
	draw_errorbar(x,y,dx,dy,ra,da)
	{
		let d = 0.3;
		let xp = x+dx*((da.x-ra.xmin)/(ra.xmax-ra.xmin));
		let yp = y+dy-dy*((da.y-ra.ymin)/(ra.ymax-ra.ymin));
		let ymin = y+dy-dy*((da.ymin-ra.ymin)/(ra.ymax-ra.ymin));
		let ymax = y+dy-dy*((da.ymax-ra.ymin)/(ra.ymax-ra.ymin));

		let thick = MEDIUMLINE;
		draw_line(xp,ymin,xp,ymax,da.col,thick);
		draw_line(xp-d,ymin,xp+d,ymin,da.col,thick);
		draw_line(xp-d/2,yp,xp+d/2,yp,da.col,thick);
		draw_line(xp-d,ymax,xp+d,ymax,da.col,thick);
	}
	
	
	/// Draws a data cross
	draw_data_cross(x,y,dx,dy,ra,da)
	{
		let d = 0.3;
		let xp = x+dx*((da.x-ra.xmin)/(ra.xmax-ra.xmin));
		let yp = y+dy-dy*((da.y-ra.ymin)/(ra.ymax-ra.ymin));
		
		let thick = NORMLINE;
		draw_line(xp-d,yp-d,xp+d,yp+d,da.col,thick);
		draw_line(xp-d,yp+d,xp+d,yp-d,da.col,thick);
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
	draw_xtick_label(x,y,dx,dy,mar,x_lab,x_vert,x_param,si,ov)
	{
		fill_rectangle(x,y,dx-1,dy,WHITE); 

		let ra = this.range;
		
		let col = BLACK; if(ov) col = DGREY;
		if(x_lab){
			let wid = (dx-mar.left-mar.right)/x_lab.length;
			
			for(let i = 0; i < x_lab.length; i++){
				let xx = x+mar.left + (dx-mar.left-mar.right)*(x_lab[i].x-ra.xmin)/(ra.xmax-ra.xmin);
				
				if(x_param == true){
					if(x_vert == true) top_vert_param_text(x_lab[i].name,xx+0.3*si,y+0.2,si,col,dy-0.2);
					else center_param_text(x_lab[i].name,xx,y+0.6,si,col,wid);
				}	
				else{
					if(x_vert == true){
						let tsa = text_sup_anno(x_lab[i].name,si,dy,"arial");
						vert_text_tsa(tsa,xx+0.3*si,y+0.2,si,col);
					}
					else{
						let tsa = text_sup_anno(x_lab[i].name,si,wid,"arial");
						center_text_tsa(tsa,xx,y+0.6,col);
					}
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
	draw_ytick_label(x,y,dx,dy,mar,y_lab,y_hor,y_param,si,ov)
	{
		fill_rectangle(x,y,dx,dy,WHITE); 
		
		let ra = this.range;
		
		let col = BLACK; if(ov) col = DGREY;
	
		if(y_lab){
			let fo = get_font(si);
			for(let i = 0; i < y_lab.length; i++){
				let yy = y + dy -yaxisgap - (dy-yaxisgap-mar.top)*(y_lab[i].y-ra.ymin)/(ra.ymax-ra.ymin);
			
				if(y_param == true){
					if(y_hor == true) right_param_text(y_lab[i].name,x+dx,yy+0.3*si,si,col,mar.left-2,20);
					else center_vert_param_text(y_lab[i].name,x+dx,yy,si,col,20); 
				}
				else{
					if(y_hor == true) right_text(y_lab[i].name,x+dx,yy+0.3*si,fo,col,mar.left-2);
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
		if(graph_dia){ error("GRAPH_DIA: axes: op:"); error(lay.op); error(this.op);}
		
		lay.add_button({x:mar.left, y:lay.dy-mar.bottom, dx:lay.dx-mar.left, dy:0, type:"x-axis"});
		
		let x_lab = lay.op.x_label;
		let x_vert = lay.op.x_label_vert;
		
		let x_param = lay.op.x_param;
			
		let yl = lay.dy-mar.bottom;
		if(x_lab){
			let dy = 0.8; 
			if(x_vert == true){ yl += 0.3; dy = mar.bottom-2.3;} 
			else yl += 0.8;
		
			lay.add_button({x:0, y:yl, dx:lay.dx, dy:dy, mar:mar, x_lab:x_lab, x_vert:x_vert, x_param:x_param, x_si:lay.op.x_si, type:"x-tick-label"});
			yl += 1;
		}
		else{
			lay.add_button({x:mar.left, y:yl, dx:lay.dx-mar.left-mar.right, dy:tick_si, type:"x-tick"});
			yl += tick_si+0.3;
			lay.add_button({x:0, y:yl, dx:lay.dx, dy:0.8, mar:mar, x_si:si_histo_label, type:"x-tick-label", ac:"X Tick"});
			yl += 1;
			
			let xz = lay.dx-4, yz = lay.dy-1.7, f= 0.8;
			lay.add_button({x:xz, y:yz+0.15, dx:1.8*f, dy:2*f, ac:"ZoomInGraph", type:"ZoomIn"});

			lay.add_button({x:xz+2*f, y:yz+0.15, dx:1.8*f, dy:2*f, ac:"ZoomOutGraph", type:"ZoomOut"});
			
			if(this.variety == "Distribution"){
				lay.add_button({x:mar.left, y:yz+0.2, dx:1.4, dy:1.4, ac:"GraphSettings", type:"Settings"});
				if(tab_name() == "Inference" && subsubtab_name() == "Parameters"){
					lay.add_button({x:mar.left+2, y:yz+0.2, dx:1.4, dy:1.4, ac:"BayesFactor", type:"BayesFactor"});
				}
			}
			
			if(this.variety == "Trace"){
				lay.add_button({x:mar.left, y:yz, dx:1.4, dy:1.4, ac:"TraceSettings", type:"Settings"});
			}
			
			if(this.variety == "Scatter"){
				lay.add_button({x:mar.left-2.4, y:yz, dx:1.4, dy:1.4, ac:"ScatterSettings", type:"Settings"});
			}
		}
		
		lay.add_button({te:this.op.x_label, x:mar.left, y:lay.dy-si_graph_label-0.21, dx:lay.dx-mar.left-mar.right, dy:si_graph_label+0.1, mar:mar, param:this.op.x_param, italic:this.op.italic, type:"x-label"});
	
		lay.add_button({x:mar.left, y:0, dx:0, dy:lay.dy-mar.bottom, param:this.op.param, italic:this.op.italic, type:"y-axis"});
		
		let y_lab = lay.op.y_label;
		let y_hor = lay.op.y_label_hor;
		let y_param = lay.op.y_param;
	
		if(y_lab){
			let xl = 1;
			lay.add_button({x:xl, y:0, dx:mar.left-1.5, dy:lay.dy-mar.bottom+yaxisgap, mar:mar, y_lab:y_lab, y_hor:y_hor, y_param:y_param, y_si:lay.op.y_si, type:"y-tick-label"});
			xl += 1;
		}
		else{
			let xl = mar.left-0.2-this.tick.wmax;
			xl -= tick_si;
			lay.add_button({x:mar.left-tick_si, y:mar.top, dx:tick_si, dy:lay.dy-mar.top-mar.bottom, type:"y-tick"});
			lay.add_button({x:xl, y:0, dx:this.tick.wmax, dy:lay.dy-mar.bottom+yaxisgap, mar:mar, y_si:si_histo_label, type:"y-tick-label", ac:"Y Tick"});
			xl -= si_graph_label+0.4;
		}
		
		lay.add_button({te:this.op.y_label, x:0, y:mar.top, dx:si_graph_label+0.2, dy:lay.dy-mar.top-mar.bottom, mar:mar, param:this.op.y_param, italic:this.op.italic, type:"y-label"});
	}
	
	
	/// Plots the time axis
	time_axis(lay)
	{
		let mar = lay.op.mar;
		
		lay.add_button({x:mar.left, y:lay.dy-mar.bottom, dx:lay.dx-mar.left, dy:0, type:"x-axis"});
		
		let yl = lay.dy-mar.bottom;
		lay.add_button({x:mar.left, y:yl, dx:lay.dx-mar.left-mar.right, dy:tick_si, type:"x-tick"});
		yl += tick_si+0.3;
		lay.add_button({x:0, y:yl, dx:lay.dx, dy:0.8, mar:mar, x_si:si_histo_label, type:"x-tick-label", ac:"X Tick"});
		yl += 1;
		lay.add_button({te:this.op.x_label, x:mar.left, y:yl, dx:lay.dx-mar.left-mar.right, dy:si_graph_label+0.1, mar:mar, param:this.op.x_paramr, italic:this.op.italic, type:"x-label"});
		
		let xz = lay.dx-6, yz = lay.dy-1.7, f= 0.8;
		lay.add_button({x:xz, y:yz, dx:1.8*f, dy:2*f, ac:"ZoomInTimeline", type:"ZoomIn"});

		lay.add_button({x:xz+2*f, y:yz, dx:1.8*f, dy:2*f, ac:"ZoomOutTimeline", type:"ZoomOut"});
	}
	
	
	/// Zooms into the timeline
	zoom_timeline_factor(fac,fx,fy)
	{
		if(fx == undefined) fx = 0.5;
		if(fy == undefined) fy = 0.5;
			
		let ra = this.range;
		let mean = ra.xmax*fx+ra.xmin*(1-fx);
		
		ra.xmin = mean - (mean-ra.xmin)/fac;
		ra.xmax = mean + (ra.xmax-mean)/fac;
		
		this.set_ind_time();
		
		this.barh *= fac;
	
		let fr_bef = fy;
		let lay_bef = get_lay("Yscroll");
		if(lay_bef){
			let spos = inter.scroll_position[lay_bef.scroll_ref];
			fr_bef = (spos.shift + fy*lay_bef.inner_dy)/spos.max;
		}
	
		generate_screen();
		
		let lay = get_lay("Yscroll");
		if(lay){
			let spos = inter.scroll_position[lay.scroll_ref];	
			let shift_new = fr_bef*spos.max - fy*lay.inner_dy
			change_scroll(shift_new,lay.but[0],"page_set");
		}
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
			
		case "Individual": case "TransTree": case "PhyloTree":
			this.set_ind_time();
			replot_layer("TimeAxis");
			break;
		}
		
		replot_layer("GraphContent");
		plot_screen();
	}
	
	
	/// Gets the classiciation
	get_cla(p,cl)
	{
		if(subtab_name() == "Results"){
			switch(tab_name()){
			case "Simulation": return model_sim.species[p].cla[cl];
			case "Inference": return model_inf.species[p].cla[cl];
			case "Post. Simulation": return model_inf.species[p].cla[cl];
			default: error("Not found"); break;
			}
		}
		return model.species[p].cla[cl];
	}
	
	
	/// Adds all the annotation buttons
	add_annotation_buts(lay)
	{
		let p = this.op.p;
		let cl = this.op.cl;
		
		let claa = this.get_cla(p,cl);  
		
		let cam = claa.camera;    
		
		for(let k = 0; k < claa.annotation.length; k++){
			let an = claa.annotation[k];
			switch(an.type){
				case "map":
					{
						let ms = find_map_store(an.map_ref);
						let box = ms.feature.box;
			
						for(let i = 0; i < ms.feature.length; i++){
							
							let fea = ms.feature[i];
							let box = fea.box;
						
							let p1 = trans_point(box.xmin,box.ymax,cam,lay);
							let p2 = trans_point(box.xmax,box.ymin,cam,lay);
							
							lay.add_button({x:p1.x, y:p1.y, dx:p2.x-p1.x, dy:p2.y-p1.y, type:"Feature", polygon:fea.polygon});
						}
					}
					break;

				default: break;
			}
		}
	}
	
	
	/// Adds all the compartment buttons
	add_compartment_buts(lay) 
	{
		if(graph_dia) error("GRAPH DIA: add_compartment_buts");
		
		let vari = this.variety;
		
		let ckey = this.colour_key;
	
		let grad;
		if(ckey != undefined) grad = ckey.grad;
		
		let p = this.op.p;
		let cl = this.op.cl;
	
		lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, p:p, cl:cl, ac:"ClassGraphBack", type:"Nothing"});
	
		let claa = this.get_cla(p,cl);  	
		let cam = claa.camera;
		
		if(this.type == "Density"){
			let di = this.density_info;
			
			let pt = trans_point(di.x,di.y,cam,lay);
			
			/// Calculates the colour matrix
			let fr = this.animation.playframe;
			if(vari == "CompVector") fr = 0;
			
			let col_val = [];

			for(let k = 0; k < claa.ncomp; k++){
				col_val[k] = LLGREY_CODE;
			}

			for(let k = 0; k < this.data.length; k++){	
				let dat = this.data[k];
				col_val[dat.c] = dat.col_k[fr];
			}
		
			let DX = di.DX, DY = di.DY;
			let mat = di.mat;
			
			let col_k = [];
			for(let j = 0; j < DY; j++){
				col_k[j] = [];
				for(let i = 0; i < DX; i++){
					let m = mat[j][i];
					if(m != undefined) col_k[j][i] = col_val[m];
				}
			}
			di.col_k = col_k;

			lay.add_button({x:pt.x, y:pt.y, dx:di.dx*cam.scale, dy:di.dy*cam.scale, type:"DensityBut"});
		}
		else{
			let ac;
			
			if(vari == "CompVector") ac = "CompVectorClick";
			if(vari == "Population") ac = "CompVectorClick";
				
			let co_graph = "CompGraph", co_latlnggraph = "CompLatLngGraph", co_map = "CompMapGraph";

			if(vari == "Population" && this.op.number_comp == true){
				co_graph = "CompGraph2"; co_latlnggraph = "CompLatLngGraph2";
			}
			
			let fr = this.animation.playframe;
			if(vari == "CompVector") fr = 0; 
			
			for(let k = 0; k < claa.ncomp; k++){
				let c = claa.comp[k];
			
				let value, CImin, CImax;
				let col = c.col;
				let name = claa.comp[k].name;
				
				switch(this.type){
				case "CompMatrix": case "CompMatrixAnim": 
					break;
				default:
					{
						let j = k;
						if(j >= this.data.length || this.data[j].c != k){
							j = find(this.data,"c",k);
						}
						if(j == undefined){
							value = "";
							CImin = "";
							CImax = "";
							col = WHITE;
						}
						else{
							let dat = this.data[j];
							if(dat.y_vec){
								value = dat.y_vec[fr];
								CImin = "";
								CImax = "";
								col = grad[dat.col_k[fr]];
							}
							else{
								let dp = dat.point[fr];
								value = dp.y;
								CImin = dp.CImin;
								CImax = dp.CImax;
								col = dp.col; 
							}
						}
					}
					break;		
				}
				
				if(!isNaN(value)){
					if(Math.floor(value) != value){ 	
						value = precision(value,4)
					}
				}
				
				switch(c.type){
				case "box":
					{
						let pt = trans_point(c.x,c.y,cam,lay);
						let w =	c.w*cam.scale*Math.exp(cam.slider.value);
						let h =	c.h*cam.scale*Math.exp(cam.slider.value); 
						let x = pt.x-w/2, y = pt.y-h/2;
					
						lay.add_button({te:name, x:x, y:y, dx:w, dy:h, ac:ac, type:co_graph, value:value, CImin:CImin, CImax:CImax, col:col, col_dark:dark_colour(col), p:p, cl:cl, i:k, show:true});
					}
					break;
				
				case "latlng":
					{
						let si = cam.scale*cam.ruler*Math.exp(cam.slider.value);
						if(si > si_limit_circle){
							let pt = trans_point(c.x,c.y,cam,lay);
							let w =	2*latlng_radius*si, h = w; 
							let x = pt.x-w/2, y = pt.y-h/2;
							lay.add_button({te:name, x:x, y:y, dx:w, dy:h, ac:ac, type:co_latlnggraph, value:value, CImin:CImin, CImax:CImax, col:col, col_dark:dark_colour(col), p:p, cl:cl, i:k});
						}
					}
					break;
					
				case "boundary":
					{
						let ms = find_map_store(c.map_ref);
						let box = ms.feature.box;
			
						let pmin = trans_point(box.xmin,box.ymax,cam,lay);
						let pmax = trans_point(box.xmax,box.ymin,cam,lay);
						
						lay.add_button({te:name, value:value, CImin:CImin, CImax:CImax, x:pmin.x, y:pmin.y, dx:pmax.x-pmin.x, dy:pmax.y-pmin.y, ac:ac, type:co_map, polygon:ms.feature.polygon, mask:ms.mask, col:col, p:p, cl:cl, i:k});
					}
					break;

				default: error("Option not recognised 82"); break;
				}
			}
		}
		
		switch(this.type){
		case "CompMatrix": case "CompMatrixAnim":
			lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, lay:lay, type:"CompLink"});
			break;
		}
		
		let timepoint = this.op.timepoint;
		let anim = inter.graph.animation;
		if(timepoint && anim.playing == false){
			let si = 1.2;
			lay.add_button({te:"t="+timepoint[anim.playframe], x:lay.dx/2-3, y:0.3, dx:6, dy:si, type:"CenterText", font:get_font(si)}); 
		}
	}


	/// Adds all the compartment buttons
	add_transition_buts(lay)
	{
		let p = this.op.p;
		let cl = this.op.cl;
		let claa = this.get_cla(p,cl);  
		model.add_transition_buts2(p,cl,"No Click",claa,false,lay);
	}

	
	/// Makes a density plot
	density_plot(x,y,dx,dy)
	{
		let di = this.density_info;
		let DX = di.DX, DY = di.DY;
		let col_k = di.col_k;
	
		let dencan = document.createElement('canvas');
		dencan.width = DX;
		dencan.height = DY;
		let dencv = dencan.getContext('2d');

		let ckey = this.colour_key;
		let rgb = ckey.grad_rgb;

		let imageData = new ImageData(DX,DY);
		let data = imageData.data;
	
		for(let j = 0; j < DY; j++){
			for(let i = 0; i < DX; i++){
				let k = col_k[j][i];
				if(k != undefined){
					let p = 4*(j*DX+i);
					if(k == LLGREY_CODE){
						data[p] = 240; data[p+1] = 240; data[p+2] = 240; data[p+3] = 255;
					}
					else{					
						let co = rgb[k];
						
						data[p] = co.r; data[p+1] = co.g; data[p+2] = co.b; data[p+3] = 255;
					}
				}
			}
		}
		dencv.putImageData(imageData,0,0);
	
		let x1 = ro(x), y1 = ro(y), x2 = ro(x+dx), y2 = ro(y+dy);

		cv.drawImage(dencan,x1,y1,x2-x1,y2-y1);
	}

	
	/// Adds buttons for 
	add_animcontrol_buts(lay)
	{
		let xz = lay.dx-6, yz = lay.dy-2;

		let xslid = xz-6;
	
		switch(this.variety){
		case "CompMatrix": case "CompVector":
			break;
		
		default:
			{	
				lay.add_button({x:playbar_mar, y:0.4, dx:lay.dx-2*playbar_mar, dy:1.1, type:"Timebar", ac:"Timebar"});
					
				let play_r = 1;
				lay.add_button({x:lay.dx/2-play_r, y:lay.dy-2*play_r, dx:2*play_r, dy:2*play_r, type:"PlayButton", ac:"PlayButton"});
				
				let dx = 3;
				lay.add_button({x:lay.dx/2-dx-play_r, y:lay.dy-2*play_r, dx:2*play_r, dy:2*play_r, type:"PlayBackward", ac:"PlayBackward"});
					
				lay.add_button({x:lay.dx/2+dx-play_r, y:lay.dy-2*play_r, dx:2*play_r, dy:2*play_r, type:"PlayForward", ac:"PlayForward"});
				
				lay.add_button({x:2.4, y:yz+0.2, dx:1.4, dy:1.4, ac:"Settings", type:"Settings"});
			}
			break;
		}
		
		switch(this.type){
		case "CompMatrix": case "CompMatrixAnim":
			lay.add_button({x:xslid, y:yz+0.5, dx:5, dy:1, info:inter.compmatrix_slider, ac:"Slider", type:"Slider"});
			break;
			
		case "Density": 
			lay.add_button({x:xslid, y:yz+0.5, dx:5, dy:1, info:inter.density_slider, ac:"Slider", type:"Slider"});
			break;
			
		case "Compartment":
			{
				let op = this.op;
				let claa = this.get_cla(op.p,op.cl);
				let cam = claa.camera;
			
				cam.slider.lay = lay.name;
				lay.add_button({x:xslid, y:yz+0.5, dx:5, dy:1, info:cam.slider, p:op.p, cl:op.cl, ac:"Slider", type:"Slider"});
			}
		}
		
		switch(this.type){
		case "HistoAnim": break;
		
		default:
			{
				let p = this.op.p;
				let cl = this.op.cl;
		
				lay.add_button({x:xz, y:yz, dx:1.8, dy:2, p:p, cl:cl, ac:"ZoomIn", type:"ZoomIn"});

				lay.add_button({x:xz+2, y:yz, dx:1.8, dy:1.8, p:p, cl:cl, ac:"ZoomOut", type:"ZoomOut"});
			}
			break;
		}
	}

	
	/// Setting for the distribution plots
	settings_dist_bubble(cont)
	{
		cont.dx = 10;
		bubble_addtitle(cont,"Settings",{});
		
		bubble_add_minititle(cont,"Distribution");
	
		let ds = get_inf_res().plot_filter.dist_settings;
		
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
	
		if(this.op.param){
			bubble_addcheckbox(cont,-0.1,"Show prior",ds.show_prior);
		
			let ds2 = get_inf_res().plot_filter.sim_val;	
			bubble_addcheckbox(cont,-0.1,"Show simulated values",ds2);
		}
		
		add_end_button(cont,"OK","DistributionPlotOK",{});	
	}
		
	
	/// Setting for the trace plots
	settings_trace_bubble(cont)
	{
		cont.dx = 10;
		bubble_addtitle(cont,"Settings",{});
		bubble_input(cont,"Burn-In %",{type:"burnin"});
		
		let ds = get_inf_res().plot_filter.sim_val;
		
		bubble_addcheckbox(cont,0,"Show simulated values",ds);
		
		add_end_button(cont,"OK","BurninOK",{});	
	}	
		
		
	/// Setting for the scatter plots
	settings_scatter_bubble(cont)
	{
		cont.dx = 10;
		bubble_addtitle(cont,"Settings",{});
		
		if(subsubtab_name() == "Individuals"){
			let seb = get_inf_res().plot_filter.scatter_settings.show_eb;
			bubble_addcheckbox(cont,0,"Show error bars",seb);
		}
		else{
			let ds = get_inf_res().plot_filter.sim_val;	
			bubble_addcheckbox(cont,0,"Show simulated values",ds);
		}
		
		add_end_button(cont,"OK","ScatterOK",{});	
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
	
	
	/// Adds a graph warning
	warning(lay)
	{
		lay.add_button({te:lay.op.warn, x:0, y:0, dx:lay.dx, dy:lay.dy, type:"BannerWarn"});
	}
	
	
	/// Creates the graph by adding layers for the content and axes
	create(x,y,w,h,lay) 
	{
		let vari = plot_variety(this.type);
		if(graph_dia) error("GRAPH DIA: createL  vari:"+vari);
		
		let hei = 4;
	
		let warn = "";
		let warn_upper = false;
		if(this.op.line_max) warn = "Too many lines to plot";
		if(this.op.ind_max) warn = "Too many individuals";
		if(this.op.point_max) warn = "Too many points";
		if(this.op.param_max){ warn = "Too many parameters to calculate"; warn_upper = true;}
		
		switch(vari){
		case "Line plot":
			{
				let mar = copy(graph_mar); mar.left += this.tick.wmax;
				
				this.extra_flag = true;
				add_layer("GraphContent",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left+GRAPH_EXTRA,h-mar.top-mar.bottom,{});
		
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
			
		case "CompMatrix plot": case "CompMatrixAnim plot": case "CompVector plot":
			{
				let mar = { right:0, left:0, top:0, bottom:0};
				add_layer("GraphAnnotations",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
				add_layer("GraphCompartments",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
		
				add_layer("AnimControls",lay.x+anim_mar,lay.y+lay.dy-anim_mar_bot,lay.dx-right_menu_width-anim_mar,hei,{});		
			}
			break;
			
		case "Comp plot":
			{
				let mar = { right:0, left:0, top:0, bottom:0};	
				let op = this.op;
				let cam = this.get_cla(op.p,op.cl).camera;
				let xx = lay.x+x+mar.left;
				let yy = lay.y+y+mar.top;
				let ww = w-mar.right-mar.left;
				let hh = h-mar.top-mar.bottom;
				
				model.ensure_lng(cam,ww);
				add_layer("GraphAnnotations",xx,yy,ww,hh,{});
				add_layer("GraphCompartments",xx,yy,ww,hh,{});
				
				if(this.variety != "CompVector"){
					add_layer("GraphTransitions",xx,yy,ww,hh,{});
				}
				
				add_layer("AnimControls",lay.x+anim_mar,lay.y+lay.dy-anim_mar_bot,lay.dx-right_menu_width-anim_mar,hei,{});
			}
			break;
			
		case "Density plot":
			{
				let mar = { right:0, left:0, top:0, bottom:0};
				add_layer("GraphAnnotations",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
				add_layer("GraphCompartments",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
			
				let hei = 4;

				add_layer("AnimControls",lay.x+anim_mar,lay.y+lay.dy-anim_mar_bot,lay.dx-right_menu_width-anim_mar,hei,{});
			}
			break;
			
		case "Individual plot": case "TransTree plot": case "PhyloTree plot":
			{
				let mar = { right:2, left:2, top:0, bottom:3.3};
				let left_shift = 0;
				if(vari == "TransTree plot" || vari == "PhyloTree plot"){
					let fo = get_font(si_transtree);
					for(let i = 0; i < this.data.length; i++){
						let da = this.data[i];
						if(da.type == "InfBar"){
							let w = text_width(da.name,fo);
							if(w > left_shift) left_shift = w;
						}
					}
					left_shift += 0.4;
				}
				mar.left += left_shift;
				
				if(this.data.length == 0){
					let dx = 20;
					lay.add_paragraph("There are no individuals in the system.",dx,lay.dx/2-dx/2,lay.dy/2,BLACK,para_si,para_lh);	
				}
				else{
					add_layer("GraphContent",lay.x+x+mar.left-left_shift,lay.y+y+mar.top,w-mar.right-mar.left+scrollw+left_shift,h-mar.top-mar.bottom,{left_shift:left_shift,j:"here"});

					add_layer("TimeAxis",lay.x+x,lay.y+y,w,h,{mar:mar});
				}
			}
			break;
			
		case "Histogram plot": case "HistoAnim plot":
			{
				let mar = copy(graph_mar); mar.left += this.tick.wmax;
			
				let x_param = this.op.x_param;
				
				if(vari == "HistoAnim plot") h -= 4;
				
				let lab=[];
				
				if(vari == "HistoAnim plot"){
					let op = this.op;
					let comp = this.get_cla(op.p,op.cl).comp;
					for(let c = 0; c < comp.length; c++){
						lab.push(comp[c].name);
					}
				}
				else{
					for(let i = 0; i < this.data.length; i++){
						let da = this.data[i];
					
						if(da.type == "Bar"){
							let name = da.name;
							lab.push(da.name);
						}
					}
				}
				x_param = false;
				
				let lab_info = this.get_axis_label_info(lab,"x",w-mar.right-mar.left,mar,x_param);
				let x_label = lab_info.label;
				let x_label_vert = lab_info.turn;
				let x_si = lab_info.si;
				if(x_label_vert == true) mar.bottom = 2.5+lab_info.wimax;
				else mar.bottom = 3.5;
			
				add_layer("GraphContent",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});
		
				add_layer("Axes",lay.x+x,lay.y+y,w,h,{mar:mar, x_label:x_label, x_param:x_param, x_label_vert:x_label_vert, x_si:x_si});
				
				if(vari == "HistoAnim plot"){
					let hei = 4;
					
					add_layer("AnimControls",lay.x+anim_mar,lay.y+lay.dy-anim_mar_bot,lay.dx-right_menu_width-anim_mar,hei,{});
				}
			}
			break;
		
		case "Matrix plot": case "MatrixAnim plot":
			{
				if(vari == "MatrixAnim plot") h -= 4;
			
				let mar = { right:2, left:3.3, top:2, bottom:3.3};
			
				if(this.data.length != 1) error("Data not right");
				let da = this.data[0];
				if(da.type != "Matrix" && da.type != "Correlation") error("Data not right"); 
				
				//let x_param = this.op.x_param;
				let x_param = da.x_param;
				
				let lab_info = this.get_axis_label_info(da.xlab,"x",w-mar.right-mar.left,mar,x_param);
	
				let x_label = lab_info.label;
				let x_label_vert = lab_info.turn;
				let x_si = lab_info.si;
				if(x_label_vert == true) mar.bottom = 2.5+lab_info.wimax;
				else mar.bottom = 3.5;
				
				//let y_param = this.op.y_param;
				let y_param = da.y_param;
						
				lab_info = this.get_axis_label_info(da.ylab,"y",h-mar.top-mar.bottom,mar,y_param);
				let y_label = lab_info.label;
				let y_label_hor = lab_info.turn;
				let y_si = lab_info.si;
			
				if(y_label_hor) mar.left = 2+lab_info.wimax;
				else mar.left = 3;
					
				add_layer("GraphContent",lay.x+x+mar.left,lay.y+y+mar.top,w-mar.right-mar.left,h-mar.top-mar.bottom,{});

				add_layer("Axes",lay.x+x,lay.y+y,w,h,{mar:mar, x_label:x_label, x_label_vert:x_label_vert, x_param:da.x_param, x_si:x_si, y_label:y_label, y_param:da.y_param, y_label_hor:y_label_hor, y_si:y_si});
				
				if(vari == "MatrixAnim plot"){
					let hei = 4;
					add_layer("AnimControls",lay.x+anim_mar,lay.y+lay.dy-anim_mar_bot,lay.dx-right_menu_width-anim_mar,hei,{});
				}
			}
			break;
			
		case "Stat table plot":
			{
				let da = this.data[0];
				add_layer("TableContent",lay.x+x,lay.y+y,w,h,{table:da.table});
			}
			break;
			
		case "No graph plot":
			center_message_box(x,y,w,h,this.op.label,lay);
			break;
			
		default: error("Variety not recognised:"+vari); break;
		}
		
		if(warn != ""){
			let yy = lay.y+y-1; if(yy < 3.7) yy = 3.7;
			if(warn_upper) yy = 1.8;
			add_layer("GraphWarn",lay.x+x,yy,w,si_graph_label+0.1,{warn:warn});
		}	
	}
	
	
	/// Works out how to fit labels around graph
	get_axis_label_info(lab,dir,room,mar,param)
	{
		let si = si_histo_label;
		let fo = get_font(si);
		let ra = this.range;
		
		let label = [];		
		let wimax = 0;
		for(let i = 0; i < lab.length; i++){
			if(dir == "x") label.push({name:lab[i], x:0.5+i});
			else label.push({name:lab[i], y:0.5+i});
			let wi;
			if(param == true) wi = plot_param_info(lab[i],si).w;
			else wi = text_width(lab[i],fo); 
		
			if(wi > wimax) wimax = wi;
		}
	
		let turn = false;
		let d_label = room/lab.length;
		
		if(d_label < wimax){
			turn = true;
			
			let frac = 0.9;                               // Ensures font size fits in gap
			if(d_label*frac < si){ wimax *= d_label*frac/si; si = d_label*frac;}
			
			if(wimax > 7){                              // Ensures axes are not too big
				si *= 7/wimax; wimax = 7;
			}
		}
		
		return { turn:turn, si:si, wimax:wimax, label:label};
	}


	/// Creates qa message when there is no graph
	no_graph(x,y,w,h,lay,te) 
	{
		if(te == undefined) te = "No data";
		center_message_box(x,y,w,h,te,lay);
	}
	
	
	/// Exports an image of the graph
	export_image(filename)
	{
		// Shifts scrollbars to the origin in tables (if necessary)	
		inter.export_image = true;
		if(subsubtab_name() == "Individuals" && this.variety == "Individual"){	
			inter.export_image = false;
		}
		
		inter.printing = true;
		generate_screen();
		
		let xmin = LARGE, xmax = -LARGE, ymin = LARGE, ymax = -LARGE;
	
		let right_bot_menu; 
	
		let list = [];
		for(let l = 0; l < inter.layer.length; l++){
			let lay = inter.layer[l];
			let name = lay.name;
			
			if(name == "RightBotMenu"){
				let fl = false;
				for(let i = 0; i < lay.but.length; i++){
					if(find_in(graph_but_not_print,lay.but[i].type) == undefined) fl = true;
				}
				if(fl == true) right_bot_menu = lay;
				else name = ""; 
			}
					
			switch(name){
			case "GraphContent": case "Axes": case "TimeAxis": case "RightBotMenu": case "TableContent":
			case "AnnotationMap": case "Annotation": case "PointLabel": 
			case "Compartment": case "Transition":	
			case "GraphAnnotations": case "GraphCompartments":
			case "GraphTransitions":
				{
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
			
			lay.plot_buttons(graph_but_not_print);
			if(lay.can){
				outcv.drawImage(lay.can,ro(lay.x-xmin+mar),ro(lay.y-ymin+mar));
			}
		}

		inter.export_image = false;
		
		inter.printing = false;
		generate_screen();
		
		if(filename) save_image(outcan,filename);
		else print_image(outcan)
	}
	
	
	/// Exports a figure
	export_figure(filename)
	{
		//inter.export_image = true;	
		inter.printing = true;
		generate_screen();
		
		save_image(inter.canvas,filename);
		
		//inter.export_image = false;
		inter.printing = false;
		generate_screen();
	}
	
	
	/// Gets the intersection between the distribuion and the value
	get_intersection(val,point)
	{
		if(point.length == 0) return 0;
		if(val < point[0].x) return 0;
		if(val > point[point.length-1].x) return 0;
	
		let i = 0; while(i < point.length-1 && point[i+1].x <= val) i++;
		let f = (val-point[i].x)/(point[i+1].x-point[i].x);
	
		return point[i].y*(1-f) + point[i+1].y*f;
	}
	
	
	/// Calculate the Bayes factor
	calculate_BF()
	{
		let data = this.data;
		let i = 0; while(i < data.length && data[i].type != "Distribution") i++;
		if(i == data.length){ alert_help("The distribution is not found"); return;}
		
		let j = 0; while(j < data.length && data[j].type != "PriorLine") j++;
		if(j == data.length){ alert_help("The prior is not found"); return;}
		
		let val = Number(inter.bubble.BF_val);
	
		let post_int = this.get_intersection(val,data[i].point);
		let prior_int = this.get_intersection(val,data[j].point);
		
		if(prior_int == 0){ alert_help("The value is outside the prior range"); return;}
		
		let te;
		if(post_int == 0) te = "The Bayes factor is large in favour of the model without the fixed parameter.";
		else{
			let bf = post_int/prior_int;
			if(bf > 1) te = "The Bayes factor is '"+bf.toFixed(1)+"' in favour of the model with the fixed parameter.";
			else te = "The Bayes factor is '"+(1/bf).toFixed(1)+"' in favour of the model without the fixed parameter.";
		}		
		
		this.bf_store = { te:te, val:val, post_int:post_int, prior_int:prior_int};
		
		change_bubble_mode("bf");
	}
}
