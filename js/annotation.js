"use strict";
// Functions which deal with annotations on the model

/// Adds all buttons associated with annotation
function add_annotation_buts(lay)
{
	let p = model.get_p();
	let cl = model.get_cl();
	if(model.species.length == 0 || cl == model.species[p].ncla) return;
	let claa = model.species[p].cla[cl];
	
	let cam = claa.camera;   
	
	if(cam.grid == "on"){
		lay.add_button({x:0, y:0, dx:lay.dx, dy:lay.dy, p:p, cl:cl, type:"Grid"});
	}
	
	for(let k = 0; k < claa.annotation.length; k++){
		let an = claa.annotation[k];
		switch(an.type){
			case "text":
				{
					let x = lay.dx/2 + (an.x-cam.x)*cam.scale; 
					let y = lay.dy/2 - (an.y-cam.y)*cam.scale;	

					let si = si_anno*(an.size/size_annotation_default)*cam.scale*cam.ruler;
					let fo = get_font(si);
					let w = text_width(an.te,fo);
					
					lay.add_button({te:an.te, x:x-w/2, y:y-0.5*si, dx:w+1, dy:si, p:p, cl:cl, i:k, si:si, font:fo, type:"LabelText", col:an.color, ac:"LabelText"});
				}
				break;
			
			case "box":
				{
					let xmin = LARGE, xmax = -LARGE;
					let ymin = LARGE, ymax = -LARGE;
					
					for(let k = 0; k < an.comps.length; k++){
						let na = an.comps[k];
						let c = hash_find(claa.hash_comp,na);
						let co = claa.comp[c];
					
						switch(co.type){
						case "box": 
							{
								let x = co.x-co.w/2; if(x < xmin) xmin = x;
								let y = co.y-co.h/2; if(y < ymin) ymin = y;
						
								x = co.x+co.w/2; if(x > xmax) xmax = x;
								y = co.y+co.h/2; if(y > ymax) ymax = y;
							}
							break;
							
						case "latlng": 
							{
								let x = co.x; if(x < xmin) xmin = x; if(x > xmax) xmax = x;
								let y = co.y; if(y < ymin) ymin = y; if(y > ymax) ymax = y;
							}
							break;
								
						case "boundary": 	
							{
								let x = co.xmid; if(x < xmin) xmin = x; if(x > xmax) xmax = x;
								let y = co.ymid; if(y < ymin) ymin = y; if(y > ymax) ymax = y;
							}
							break;

						default: error("Option not recognised 2"); break;
						}
					}
					
					let mar = cam.ruler;
					if(xmin == xmax || ymin == ymax) mar *= 2;
					
					let xc = xmin - mar;
					let yc = ymax + mar;
					let wc = xmax-xmin+2*mar; 
					let hc = ymax-ymin+2*mar;
					
					let si = si_anno*(an.size/size_annotation_default)*cam.scale*cam.ruler;
					
					let x = lay.dx/2 + (xc-cam.x)*cam.scale; 
					let y = lay.dy/2 - (yc-cam.y)*cam.scale - si;	
	
					let w = wc*cam.scale; 
					let h = hc*cam.scale+si; 
				
					lay.add_button({te:an.te, x:x, y:y, dx:w, dy:h, type:"Box", p:p, cl:cl, i:k, col:an.color, ac:"Box", si:si});
				}
				break;
			
			case "map": break;

			default: error("Option not recognised 1"+an.type); break;
		}
	}
}

/// Adds point labels in latlng view
function add_pointlabel_buts(lay)
{
	let p = model.get_p();
	let cl = model.get_cl();
	if(model.species.length == 0 || cl == model.species[p].ncla) return;
	let claa = model.species[p].cla[cl];
	let cam = claa.camera;   
	
	let si = cam.scale*cam.ruler*Math.exp(cam.slider.value);

	let fo = get_font(si);
	if(si > si_limit_label || claa.ncomp < 100){
		for(let k = 0; k < claa.ncomp; k++){
			let c = claa.comp[k];
			if(c.type == "latlng"){
				let p = trans_point(c.x,c.y,cam,lay);
				let te = claa.comp[k].name
		
				lay.add_button({te:te, x:p.x, y:p.y+1*si, dx:0, dy:0, si:si, type:"LatLngLabel", col:BLACK});
			}
		}
	}
}

	
/// Adds all the buttons associated with annotation maps
function add_annotation_map_buts(lay)
{
	let p = model.get_p();
	let cl = model.get_cl();
	
	if(model.species.length == 0 || cl == model.species[p].ncla) return;

	let claa = model.species[p].cla[cl];

	let cam = claa.camera; 

	for(let k = 0; k < claa.annotation.length; k++){
		let an = claa.annotation[k];
		switch(an.type){
			case "map":
				{
					let feature = find_feature(an.map_ref);
					
					for(let i = 0; i < feature.length; i++){
						let fea = feature[i];
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


/// Finds a feature from a reference
function find_feature(name)
{	
	let k = find(map_store,"name",name); 
	if(k == undefined) error("Cannot find feature");
	return map_store[k].feature;
}	


/// Finds a feature from a reference
function find_map_store(name)
{	
	let k = find(map_store,"name",name); 
	if(k == undefined) error("Cannot find feature");
	return map_store[k];
}	


/// Transforms a lng lat position into an x,y position 
function transform_latlng(lng,lat)
{
	let x = lng*Math.PI/180;
	let sign = 1;
	let phi_p = lat*Math.PI/180; if(phi_p < 0){ phi_p = -phi_p; sign = -1;}
	let y = sign*Math.log(Math.tan(0.5*phi_p + 0.25*Math.PI));
	
	return {x:x, y:y};
}


/// Transforms a x,y position into an lng lat position 
function transform_latlng_inv(x,y)
{
	let lng = x*180/Math.PI;
	while(lng > 180) lng -= 360;
	while(lng < -180) lng += 360;
	
	let sign = 1;
	if(y < 0){ sign = -1; y = -y;}
	
	let phi_p = (Math.atan(Math.exp(y)) - 0.25*Math.PI)/0.5;
	
	let lat = sign*phi_p*180/Math.PI;
	
	return {lng:lng, lat:lat};
}


/// Transforms a point based on camera
function trans_point(x,y,cam,lay)
{
	return {x:lay.dx/2 + (x-cam.x)*cam.scale, y:lay.dy/2 - (y-cam.y)*cam.scale};
}


/// Does the reverse function (i.e. gets position in model based on mouse location
function trans_point_rev(xx,yy,cam,lay)
{
	return {x:cam.x + (xx - lay.dx/2)/cam.scale, y:cam.y - (yy - lay.dy/2)/cam.scale};
}
	
	
/// Loads up a map into the annotations
function load_annotation_map(te)
{			
	let data = JSON.parse(te);

	let feature = get_feature(data);
	if(feature.length == 0){
		alert_help("Error loading file","This map has no features");
		return;
	}	
	let box = get_map_bound_box(feature);
	
	return {type:"map", feature:feature, box:box}
}


/// Loads up a map into compartments
function load_compartment_map(te)
{	
	let da = JSON.parse(te);

	let feature = get_feature(da);
	if(feature.length == 0){
		alert_help("Error loading file","This map has no features");
		return;
	}	
		
	let tab = get_feature_table(da,feature);
	
	post({type:"Load Comp Map", tab:tab});
}


/// Adds a compartment map
function add_compartment_map(col,so)
{
	let p = so.p, cl = so.cl;
	clear_classification(p,cl);
	
	map_store = [];
	
	let tab = so.table;
	let ele = so.table.ele;
	for(let r = 0; r < tab.nrow; r++){
		let fea = tab.ele[r][0];
		let res = model.add_map_compartment(tab.ele[r][1],p,cl,fea,get_mask(fea),col,false);
		output_help(res);
	}
}


/// Clears a classification such that it can be reloaded
function clear_classification(p,cl,op)
{
	let claa = model.species[p].cla[cl];
	claa.comp=[]; claa.ncomp=0; claa.hash_comp=[];
	claa.tra=[]; claa.ntra = 0; claa.hash_tra=[];
	claa.nsource = 0; claa.nsink = 0;
	
	if(op == undefined){
		let list = model.find_clones(p,cl);	
		for(let j = 0; j < list.length; j++){
			clear_classification(list.p,list.cl,"no");
		}
	}
}


/// Adds an individual compartmental bounardy based on file
function add_individual_compartment_boundary(name,p,cl,file,color,infected)
{	
	let tab = import_geojson(file);
	if(tab == undefined) return err("Error importing");
	
	for(let c = 1; c < tab.ncol; c++){
		for(let r = 0; r < tab.nrow; r++){
			if(tab.ele[r][c] == name){
				let fea = tab.ele[r][0];
				return model.add_map_compartment(name,p,cl,fea,get_mask(fea),color,infected);
			}
		}
	}
 
	return err(in_file_text(file)+", '"+name+"' could not be found");
}

	
/// Gets a mask of the feature so that mouse over works
function get_mask(fea)
{
	const can = new OffscreenCanvas(mask_size,mask_size);
	
  let cv = can.getContext('2d');
	
	can.width = mask_size;
	can.height = mask_size;
	
	let poly = fea.polygon;
	for(let p = 0; p < poly.length; p++){
		let cor = poly[p].cor;
		
		cv.beginPath();
		for(let i = 0; i < cor.length; i++){
			let xx = mask_size*cor[i][0];
			let yy = mask_size*cor[i][1];
			if(i == 0) cv.moveTo(xx,yy);
			else cv.lineTo(xx,yy);
		}
		cv.closePath();
		cv.fillStyle = WHITE;
		cv.fill();
	}

	let data = cv.getImageData(0,0,mask_size,mask_size).data;
	let mask = [];
	for(let i = 0; i < mask_size*mask_size; i++){
		if(data[4*i+3] < 128) mask.push(0); else mask.push(1);
	}
	
	return mask;
}


/// Works out if mouse is over the compartment map
function mouse_over_comp_map(frx,fry,mask)
{
	if(mask[Math.floor(fry*mask_size)*mask_size+Math.floor(frx*mask_size)] == 1){
		return true;
	}		

	return false;	
}


/// Sets the camera to fit in all the objects 
function set_camera(p,cl)
{
	let claa = model.species[p].cla[cl];
	
	let cam = claa.camera;

	let box = get_model_box(p,cl); if(box.unset) return;
	
	cam.x = 0.5*(box.xmin+box.xmax); cam.y = 0.5*(box.ymin+box.ymax);
	cam.scale = box.scale;
	
	set_ruler(p,cl);
}


/// Sets the camera ruler (this is used to detemine text size in map view)
function set_ruler(p,cl)
{
	let box = get_model_box(p,cl);
	
	let claa = model.species[p].cla[cl];
	let cam = claa.camera;

	switch(cam.coord){
	case "cartesian":
		cam.ruler = 1; 
		break;
		
	case "latlng":
		{
			// This sets default ruller size to 2 unit on final page
		
			let num = claa.ncomp; if(num == 0) num = 1;
			let area_per_comp = (page_char_wid-menu_width)*page_char_hei/num;
			
			let dist = 0.1*Math.sqrt(area_per_comp);
			
			cam.ruler = dist/box.scale; 
			
			/*
			// This works out the average minimum size to nearest compartment
			let xmin = LARGE, xmax = -LARGE;
			let ymin = LARGE, ymax = -LARGE;
			for(let c = 0; c < claa.ncomp; c++){
				let p = comp_pos(claa.comp[c]);
				if(p.x < xmin) xmin = p.x; if(p.x > xmax) xmax = p.x;
				if(p.y < ymin) ymin = p.y; if(p.y > ymax) ymax = p.y;
			}
			
			let N = 20;
			let grid=[];
			for(let k = 0; k < N*N; k++) grid[k]=[];
			
			for(let c = 0; c < claa.ncomp; c++){
				let co = claa.comp[c];
				let p = comp_pos(claa.comp[c]);
				
				let i = Math.floor(ALMOST_ONE*N*(p.x-xmin)/(xmax-xmin+TINY));
				let j = Math.floor(ALMOST_ONE*N*(p.y-ymin)/(ymax-ymin+TINY));
				grid[j*N+i].push(c);
			}
			
			let av = 0, nav = 0;
			for(let k = 0; k < N*N; k++){
				let gr = grid[k];
				if(gr.length > 1){
					for(let i = 0; i < gr.length; i++){
						let min = LARGE;
						let pi = comp_pos(claa.comp[gr[i]]);
				
						for(let j = 0; j < gr.length; j++){
							if(j != i){
								let pj = comp_pos(claa.comp[gr[j]]);
								let d = (pi.x-pj.x)*(pi.x-pj.x) + (pi.y-pj.y)*(pi.y-pj.y);
								if(d < min) min = d;
							}
						}
						//av += Math.sqrt(min); nav++;
						av += Math.sqrt(min); nav++;
					}
				}
			}
			
			if(nav > 0){
				av /= nav;
				let rule = 0.1*av;
				if(rule < cam.ruler) cam.ruler = rule;	
			}
			*/
			/*
			let dx = xmax-xmin, dy = ymax-ymin;
			let d = dx; if(dy > d) d = dy;
			
			if(d == 0){
				if(cam.ruler < 0.05) cam.ruler = 0.05;
			}
			else{
				if(cam.ruler > 0.1*d){
					cam.ruler = 0.1*d;
				}
			}
			*/
		}
		break;

	default: error("Option not recognised 3"); break;
	}
}


/// Gets the position of a compartment
function comp_pos(co)
{
	if(co.type == "boundary") return {x:co.xmid, y:co.ymid};
	return {x:co.x, y:co.y};
}


/// Gets a box surrounding the complete model
function get_model_box(p,cl)
{
	let claa = model.species[p].cla[cl];
	let cam = claa.camera;
	
	let xmin = LARGE, xmax = -LARGE;
	let ymin = LARGE, ymax = -LARGE;

	for(let c = 0; c < claa.ncomp; c++){
		let co = claa.comp[c];

		switch(co.type){
		case "boundary":
			{
				let feature = find_feature(co.map_ref);
				let box = feature.box;
					
				if(box.xmin < xmin) xmin = box.xmin; 
				if(box.xmax > xmax) xmax = box.xmax; 
				if(box.ymin < ymin) ymin = box.ymin; 
				if(box.ymax > ymax) ymax = box.ymax; 
			}
			break;
			
		case "box": case "latlng":
			if(co.x < xmin) xmin = co.x; 
			if(co.x > xmax) xmax = co.x;
			if(co.y < ymin) ymin = co.y; 
			if(co.y > ymax) ymax = co.y;
			break;

		default: error("Option not recognised 4"); break;
		}
	}
			
	for(let i = 0; i < claa.ntra; i++){
		let tr = claa.tra[i];
		for(let j = 0; j < tr.midp.length; j++){
			let po = tr.midp[j];
			if(po.x < xmin) xmin = po.x;
			if(po.x > xmax) xmax = po.x;
			if(po.y < ymin) ymin = po.y;
			if(po.y > ymax) ymax = po.y;			
		}
	}
	
	for(let i = 0; i < claa.annotation.length; i++){
		let anno = claa.annotation[i];

		let ms;
		if(anno.type == "map"){
			let i = find(map_store,"name",anno.map_ref);
			if(i == undefined) error("Cannot find map");
			ms = map_store[i];
		}

		if(!(ms && ms.feature.length == 248)){	
			switch(anno.type){
			case "text":
				if(anno.x < xmin) xmin = anno.x; 
				if(anno.x > xmax) xmax = anno.x;
				if(anno.y < ymin) ymin = anno.y; 
				if(anno.y > ymax) ymax = anno.y;
				break;
				
			case "box":
				break;
				
			case "map":
				let box = ms.box;
				if(box.xmin < xmin) xmin = box.xmin;
				if(box.xmax > xmax) xmax = box.xmax;
				if(box.ymin < ymin) ymin = box.ymin;
				if(box.ymax > ymax) ymax = box.ymax;
				break;

			default: error("Option not recognised 5"+anno.type); break;
			}		
		}
	}

	let unset = false;

	let scale;
	if(xmin == LARGE || xmax == xmin){
		scale = cam.scale;
		unset = true;
	}
	else{
		let scalex = 0.7*(page_char_wid-menu_width)/(xmax-xmin); 
		if(cam.coord == "cartesian" && scalex > 1) scalex = 1;
		
		let scaley = 0.7*page_char_hei/(ymax-ymin);
		if(cam.coord == "cartesian" && scaley > 1) scaley = 1;

		if(scalex < scaley) scale = scalex;
		else scale = scaley;
	}
	
		
	return {unset:unset, xmin:xmin, xmax:xmax, ymin:ymin, ymax:ymax, scale:scale};
}


/// Gets the bounding box for a series of polygons
function get_map_bound_box(feature)
{
	let xmin = LARGE, xmax = -LARGE;
	let ymin = LARGE, ymax = -LARGE;
	for(let i = 0; i < feature.length; i++){
		let box = feature[i].box;
		if(box.xmin < xmin) xmin = box.xmin;
		if(box.xmax > xmax) xmax = box.xmax;
		if(box.ymin < ymin) ymin = box.ymin;
		if(box.ymax > ymax) ymax = box.ymax;
	}
	return {xmin:xmin, xmax:xmax, ymin:ymin, ymax:ymax}; 
}

	
/// Gets polygons from the data
function get_feature(data) 
{
	let feature = [];
	
	for(let i = 0; i < data.features.length; i++){
		let geo = data.features[i].geometry;

		let polygon = [];
	
		let cor = geo.coordinates;
		
		switch(geo.type){
		case "Polygon":	
			for(let j = 0; j < cor.length; j++){
				polygon.push({cor:cor[j], foreground:true}); 
			}
			break;
			
		case "MultiPolygon":
			for(let pn = 0; pn < cor.length; pn++){
				let cor2 = cor[pn];
				for(let j = 0; j < cor2.length; j++){
					polygon.push({cor:cor2[j], foreground:true}); 
				}
			}
			break;
			
		default: error(geo[j].type+" Type not syported"); break;
		}
	
		// Transforms lng lat to x y 
		for(let i = 0; i < polygon.length; i++){
			let points = polygon[i].cor;
			for(let j = 0; j < points.length; j++){
				let p = points[j];
				let pt = transform_latlng(p[0],p[1]);
				p[0] = pt.x; p[1] = pt.y; 
			}
		}

		// Gets bounding box
		let xmin = LARGE, xmax = -LARGE;
		let ymin = LARGE, ymax = -LARGE;
		
		for(let i = 0; i < polygon.length; i++){
			let points = polygon[i].cor;
			for(let j = 0; j < points.length; j++){
				let p = points[j];
				if(p[0] < xmin) xmin = p[0];
				if(p[0] > xmax) xmax = p[0];
				if(p[1] < ymin) ymin = p[1];
				if(p[1] > ymax) ymax = p[1];
			} 
		}
			
		// Scales points within the bounding box
		if(xmin != xmax && ymin != ymax){	
			for(let i = 0; i < polygon.length; i++){
				let points = polygon[i].cor;
				for(let j = 0; j < points.length; j++){
					let p = points[j];
					p[0] = (p[0] - xmin)/(xmax - xmin);
					p[1] = (p[1] - ymin)/(ymax - ymin);
				} 
			}
		}
		
		let box = {xmin:xmin, xmax:xmax, ymin:ymin, ymax:ymax};
		feature.push({ name:"", polygon:polygon, box:box});
	}
	
	return feature;
}


/// Gets names for all the features
function get_feature_table(da,feature) 
{
	let nrow = da.features.length;
	let geo = da.features[0].properties;
	
	let keys = Object.keys(geo);

	let ele=[];
	let heading=[];
	heading.push("Boundary");
	for(let r = 0; r < nrow; r++){		
		ele[r]=[];
		ele[r].push(feature[r]);
	}
		
	for(let i = 0; i < keys.length; i++){
		let prop = keys[i];
		let r = 0; while(r < nrow && da.features[r].properties[prop] != undefined) r++;
		if(r == nrow){
			heading.push(prop);
			for(let r = 0; r < nrow; r++){
				ele[r].push(da.features[r].properties[prop]);
			}
		}
	}

	return {filename:"file", heading:heading, col_used:[], ele:ele, ncol:heading.length, nrow:nrow, edit:false};
}

	
/// Sets up removal of annotation bubble
function remove_annotation_init()
{
	let p = model.get_p();
	let cl = model.get_cl();
	let anno = model.species[p].cla[cl].annotation;
	let list = [];
	let wmax = 0;
	
	for(let i = 0; i < anno.length; i++){			
		let an = anno[i];
		let te;
		switch(an.type){
		case "map": te = "Map with "+find_feature(an.map_ref).length+" features";	break;
		case "text": te = "Label '"+an.te+"'"; break;
		case "box": te = "Box with text '"+an.te+"'"; break;
		default: error("Annotation type not recognised"); break;
		}
	
		let w = text_width(te,get_font(si_radio));
		if(w > wmax) wmax = w;
		list.push({te:te, box:{check:false}});
	}
	wmax += 3.5;
	if(wmax < 11.4) wmax = 11.4;
	
	inter.bubble.rem_annotation = { w:wmax, p:p, cl:cl, list:list};
	change_bubble_mode("remove");
}
	
	
/// Adds a scrollable list to allow for removal of annotations
function annotation_scrollable(lay)
{
	let list = lay.op.list;

	let cy = 0;
	for(let i = 0; i < list.length; i++){
		lay.add_checkbox(0,cy,list[i].te,list[i].te,list[i].box);
		cy += 1.2;
	}
	return cy;
}


/// Deletes selected annotations
function delete_annotation()
{
	let rem = inter.bubble.rem_annotation;
	let anno = model.species[rem.p].cla[rem.cl].annotation;

	for(let j = anno.length-1; j >= 0; j--){
		if(rem.list[j].box.check == true){
			let ann = anno[j];
			if(ann.type == "map") delete_map_store(anno[j].map_ref);
		
			anno.splice(j,1);
		}
	}		
}


/// Deletes an element in map_list
function delete_map_store(na)
{
	let k = find(map_store,"name",na);
	if(k == undefined) error("cannot delete map store");
	map_store.splice(k,1);
}


/// Zooms in and out
function zoom_factor(p,cl,fac)
{
	let show = model.get_show_model();
			
	let cam;

	if(show == true) cam = model.species[p].cla[cl].camera;
	else cam = inter.graph.get_cla(p,cl).camera;
	
	cam.scale *= fac;
	
	model.ensure_lng(cam);
	
	if(show == true) clone_camera(p,cl);
}


/// Zooms into a location selected by the mouse
function zoom_double_click(p,cl)
{
	let mx = inter.mx, my = inter.my;
	let lay;
	
	let show = model.get_show_model();
	
	let cam;
	if(show == true){
		cam = model.species[p].cla[cl].camera;
		lay = get_lay("Main");
	}
	else{
		cam = inter.graph.get_cla(p,cl).camera;
		lay = get_lay("GraphCompartments");
	}

	let x = mx-lay.x, y = my-lay.y;
	
	let f = 0.25;
	let xx = f*x + (1-f)*lay.dx/2;
	let yy = f*y + (1-f)*lay.dy/2;
	
	let rev = trans_point_rev(xx,yy,cam,lay);
	cam.x = rev.x; 
	cam.y = rev.y; 
	
	zoom_factor(p,cl,1.3);
}


/// Adds a layer which allows the user to select compartments
function add_select_buts(lay)
{
	lay.add_button({lay_x:lay.x, lay_y:lay.y, x:0, y:0, dx:lay.dx, dy:lay.dy, ac:"Select", type:"Select"});
}


/// Adds a box to the model
function add_box()
{
	let p = model.get_p();
	let cl = model.get_cl();
	
	let claa = model.species[p].cla[cl];
	let comps = [];
	let lay = get_lay("Compartment");

	let mo = inter.mode;
	let x1 = mo.x1, x2 = mo.x2, y1 = mo.y1, y2 = mo.y2;
	if(x1 > x2){ let t = x1; x1 = x2; x2 = t;}
	if(y1 > y2){ let t = y1; y1 = y2; y2 = t;}
	
	for(let i = 0; i < lay.but.length; i++){
		let bu = lay.but[i];
	
		switch(bu.type){
		case "Compartment": case "CompLatLng":
			let x = bu.x+bu.dx; 
			let y = bu.y+bu.dy;
			if(x > x1 && x < x2 && y > y1 && y < y2){
				comps.push(bu.te);
			}
			break;
		}
	}

	claa.annotation.push({type:"box", te:"", size:size_annotation_default, comps:comps, color:annotation_col_default});
}


/// Creates a GEOJSON feature from a feature (used to create files from loaded boundaries)
function generate_JSON_feature(name,feature)
{
	let cor = [];
	let poly = feature.polygon;
	let box = feature.box;
	for(let i = 0; i < poly.length; i++){
		cor[i] = [];
		cor[i].push(copy(poly[i].cor));
	}
	
	for(let i = 0; i < cor.length; i++){
		for(let j = 0; j < cor[i].length; j++){
			for(let k = 0; k < cor[i][j].length; k++){
				let x = box.xmin + cor[i][j][k][0]*(box.xmax-box.xmin);
				let y = box.ymin + cor[i][j][k][1]*(box.ymax-box.ymin);
				let p = transform_latlng_inv(x,y);
				cor[i][j][k][0] = Number(precision(p.lng,6));
				cor[i][j][k][1] = Number(precision(p.lat,6));
			}
		}
	}
	
	let fea = { type:"Feature", geometry:{type:"MultiPolygon", coordinates:cor}, properties:{name:name}};

	return fea;
}


/// Loads a map from worker reply
function load_map(ans)
{
	let p = ans.info.p, cl = ans.info.cl; 
	
	let name = "file"+Math.random();
	map_store.push({name:name, feature:ans.feature, box:ans.box});
	
	let claa = model.species[p].cla[cl];
	claa.annotation.push({type:"map", map_ref:name, default_map:ans.info.default_map});
	
	set_camera(p,cl);
	
	if(ans.info.default_map == true){
		claa.default_map = true;
		
		let cam = claa.camera;
		let po = transform_latlng(0,40);
		cam.x = po.x;
		cam.y = po.y;
		cam.scale = (page_char_wid-menu_width)/(2*Math.PI);
	}
}


/// Gets the number of steps
function get_step(v1,v2,min)
{
	let ideal = 20;
	let dif = v2-v1;

	let step = min; 
	while(dif/step > ideal){
		if(dif/(2*step) < ideal){ step *= 2; break;} 
		if(dif/(5*step) < ideal){ step *= 5; break;} 
		step *= 10;
	}
	return step;
}


/// Gets equally spaced division between two values
function get_div(v1,v2,step)
{
	let list=[];
	for(let i = Math.floor(v1/step)+1; i <= Math.floor(v2/step); i++){
		list.push(round_small(step*i));
	}
	return list;
}


/// Determines if the circle is used for the origin
function circle_origin(divx,divy)
{
	if(divx[0] <= 0 && divx[divx.length-1] > 0){
		if(divy[0] <= 0 && divy[divy.length-1] > 0) return true;
	}
	return false;
}


/// Gets the position of the origin
function get_origin(div)
{
	let i = 0;
	if(div[0] <= 0 && div[div.length-1] > 0){
		while(i < div.length && div[i] < 0) i++;
	}
	else{
		if(div[0] <= 0) i = div.length-2;
	}
	return i;
}

			
/// Plots the grid behind the model
function plot_grid(p,cl,dx,dy)
{
	let claa = model.species[p].cla[cl];
	
	let cam = claa.camera;
	
	let lay = get_lay("AnnotationMap");
	
	let col = MLLBLUE;
	let col_text = LBLUE;
	
	let xst=[], yst=[];
	let divx, divy;
	
	switch(cam.coord){
	case "cartesian":
		{	
			let p1 = trans_point_rev(0,dy,cam,lay);
			let p2 = trans_point_rev(dx,0,cam,lay);
			
			let step = get_step(p1.x,p2.x,1);
			
			divx = get_div(p1.x,p2.x,step);
			divy = get_div(p1.y,p2.y,step);
			
			for(let i = 0; i < divx.length; i++){
				let tp = trans_point(divx[i],0,cam,lay);
				xst.push(tp.x);
			}
			
			for(let i = 0; i < divy.length; i++){
				let tp = trans_point(0,divy[i],cam,lay);
				yst.push(tp.y);
			}
		}
		break;
	
	case "latlng":
		{	
			let p1 = trans_point_rev(0,dy,cam,lay);
			let p2 = trans_point_rev(dx,0,cam,lay);
			
			let pg1 = transform_latlng_inv(p1.x,p1.y,cam,lay);
			let pg2 = transform_latlng_inv(p2.x,p2.y,cam,lay);
			
			let step = get_step(pg1.lng,pg2.lng,0.1);
			
			divx = get_div(pg1.lng,pg2.lng,step);
			divy = get_div(pg1.lat,pg2.lat,step);
		
			for(let i = 0; i < divx.length; i++){
				let tp = transform_latlng(divx[i],0);
				let p = trans_point(tp.x,tp.y,cam,lay);
				xst.push(p.x);
			}
			
			for(let i = 0; i < divy.length; i++){
				let tp = transform_latlng(0,divy[i]);
				let p = trans_point(tp.x,tp.y,cam,lay);
				yst.push(p.y);
			}
		}
		break;
	}
	
	let xi = get_origin(divx);
	let xmid = 0.95*xst[xi]+0.05*xst[xi+1];
	
	let yi = get_origin(divy);
	let ymid = 0.95*yst[yi]+0.05*yst[yi+1];
	
	let ori = circle_origin(divx,divy);
	
	let fo = get_font(0.9);
	for(let i = 0; i < divx.length; i++){
		let xx = xst[i];
		draw_line(xx,0,xx,dy,col,NORMLINE);
		if(!(divx[i]==0 && ori) && !(i == xi && !ori)){
			plot_vert_text_rev(divx[i],xx+0.2,ymid,fo,col_text);  
		}			
	}
	
	for(let i = 0; i < divy.length; i++){
		let yy = yst[i];
		draw_line(0,yy,dx,yy,col,NORMLINE);
		if(!(divy[i]==0 && ori)) plot_text(divy[i],xmid,yy-0.2,fo,col_text); 
	}
	
	if(ori) fill_circle(xmid+0.3,ymid-0.3,0.3,col,col_text,NORMLINE);
}
