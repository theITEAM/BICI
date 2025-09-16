"use strict";
// Functions related to pop-up bubbles

/// Adds buttons associated with a speech bubble
function add_bubble_buts(lay)
{
	let bub = inter.bubble;
	
	if(bub.i == undefined) return;

	let bu_lay = get_lay(bub.lay_name);
	if(bu_lay == undefined) return;
	
	let bu = bu_lay.but[bub.i];
	
	let cont = {y:0, input_list:[], drop_list:[], lay:lay, bu:bu, bu_lay:bu_lay, end_button:[]};
	
	bub.warning = false;

	switch(bu.type){	
	case "NameLink":
		cont.dx = 14;
		bubble_addtitle(cont,"Data source name",{te:data_source_text});
		bubble_input(cont,"Name:",{type:"dataname"});
		add_end_button(cont,"Done","Done",{});	
		break;
		
	case "LinkPara":
		switch(bu.ac){
		case "LoadAMatrix":
			cont.dx = 9;
			bubble_addtitle(cont,"Load options",{te:load_options_text});
			bubble_addradio(cont,0,"A","A matrix",bub.type_radio);
			bubble_addradio(cont,0,"Ainv","Inverse of A matrix",bub.type_radio);
			bubble_addradio(cont,0,"pedigree","Pedigree",bub.type_radio); 
			add_end_button(cont,"Done","LoadAMatrix2",{p:bu.p,i:bu.i});	
			break;
		}
		break;
		
	case "BayesFactor":
		cont.dx = 9;
		bubble_addtitle(cont,"Bayes Factor",{te:bayesfactor_text});
		switch(bub.mode){
		case "bf":
			bubble_addparagraph(cont,inter.graph.bf_store.te,0,cont.dx);
			add_end_button(cont,"Done","Done",{});	
			break;
			
		default:
			
			bubble_input(cont,"Value:",{type:"BF_val"});
			add_end_button(cont,"Calculate","CalculateBF",{});	
		}
		break;
		
	case "ShowModel":
		{
			let ex = bu.mod_but;
			
			cont.dx = 14;
			bubble_addtitle(cont,"Model used",{});
			let pic = find_pic(ex.pic);
		
			bubble_addpicture(cont,pic);
		
			let si = 0.9, lh = 1.1, ma = 0.5;  
			let text_anno = text_convert_annotation("<b>"+ex.mod+": "+ex.te+"</b>",si,lh,cont.dx-2*ma,"center",DDBLUE);
			let word = text_anno.word;
			
			let dy = text_anno.height;
		
			cont.lay.add_button({word:word, x:0, y:cont.y, dx:cont.dx, dy:dy, type:"PlotWord"}); 
			cont.y += dy;
			
			add_end_button(cont,"Done","Done",{});	
		}
		break;
	
	case "ModType":
		{
			let title, te;
			switch(bu.variety){
			case "IBM":
				title = "Individual-based model"; 
				te = "This example uses an individual-based model (IBM). This means that individuals are explicitly modelled with each having it's own timeline.";
				break;
				
			case "POP":
				title = "Population-based model"; 
				te = "This example uses a population-based model (PBM). This means that only the compartmental populations are modelled. Such a model cannot make used of individual-level data or individual-level variation (e.g. through fixed or individual effects).";
				break;
				
			case "IBM/POP": 
				title = "Combined population and individual-based model"; 
				te = "This example uses both an individual (IBM) and a population-based model (PBM). This is often useful for modelling different aspects of the model.";
				break;
			}
			cont.dx = 17;
			bubble_addtitle(cont,title,{});
			bubble_addbigparagraph(cont,te,0,cont.dx);
			add_end_button(cont,"Done","Done",{});	
		}
		break;
		
	case "CopyPic":
		cont.dx = 8;
		bubble_addtitle(cont,"Copy",{});
		bubble_input(cont,"Suffix:",{type:"suffix"});
		add_end_button(cont,"Done","DoneCopyPic",{});	
		break;
		
	case "CompAlpha": case "MultiCompAlpha":
		cont.dx = 10;
		bubble_addtitle(cont,"Edit α",{});
		bubble_input(cont,"Value:",{type:"alpha_val"});
		add_end_button(cont,"Done","DoneAlpha",{});	
		break;
		
	case "ProbEqn":
		cont.dx = 10;
		bubble_addtitle(cont,"Detection prob.",{});
		add_end_button(cont,"Done","Done",{});	
		break;
		
	case "CompGraph": case "CompGraph2": case "CompLatLngGraph": case "CompMapGraph": 
		{
			cont.dx = 9.3;
			bubble_addtitle(cont,"Compartment",{});
			bubble_addparagraph(cont,"<b>Name:</b> "+bu.te,0,cont.dx);
			if(bu.CImin){	
				bubble_addparagraph(cont,"<b>Mean:</b> "+precision(bu.value,5),0,cont.dx);
				bubble_addparagraph(cont,"<b>CI min:</b> "+precision(bu.CImin,5),0,cont.dx);
				bubble_addparagraph(cont,"<b>CI max:</b> "+precision(bu.CImax,5),0,cont.dx);
			}
			else{
				bubble_addparagraph(cont,"<b>Value:</b> "+precision(bu.value,5),0,cont.dx);
			}		
			add_end_button(cont,"Done","Done",{});	
		}
		break;
	
	case "SliceTime":
		{
			cont.dx = 9.3;
			bubble_addtitle(cont,"Slice Time",{});
			bubble_input(cont,"Time:",{type:"slice_time"});
			add_end_button(cont,"Done","DoneSliceTime",{});	
		}
		break;
		
	case "PopFilt":
		{
			let name = subsubtab_name();
			let popf = bub.popfilt;
			let filt = popf.filter;
			let cl_sel = popf.rpf.species[filt.p].sel_class.cl;
			
			let te;
			switch(name){
			case "Populations":
				if(filt.cl == cl_sel) te = filter_pop_incl_text;
				else te = filter_pop_notincl_text;
				break;
				
			case "Transitions":
				if(filt.cl == cl_sel) te = filter_trans_incl_text;
				else te = filter_trans_notincl_text;
				break;
				
			case "Individuals":
				{
					let ind_sel = inter.graph.ind_sel;
					if(ind_sel != undefined){
						if(filt.cl == cl_sel) te = filter_ind_sel_incl_text;
						else te = filter_ind_sel_notincl_text;
					}
					else{
						if(filt.cl == cl_sel) te = filter_ind_incl_text;
						else te = filter_ind_notincl_text;
					}
				}
				break;
			}
				
			cont.dx = 14.3;
			bubble_addtitle(cont,"Filter",{te:te});
			
			if(filt.cl != cl_sel){
				bubble_addradio(cont,0,"select","Select",bub.popfilt.filter.radio); cont.y -= 1.4;
				bubble_addradio(cont,5,"single","Single",bub.popfilt.filter.radio);
			}
			
			cont.y += 0.2;
			bubble_addscrollable(cont,{type:"filterpos", ymax:10});
		
			if(name == "Populations"){		
				cont.y += 0.3;
				bubble_addcheckbox(cont,0,"Calculate fraction",bub.popfilt.filter.fraction);
				cont.y += 0.1;
			}
			
			add_end_button(cont,"Delete","RemoveFilter",{num:bu.num, p: bu.filter.p, rpf:bub.popfilt.rpf});	
			add_end_button(cont,"Done","DoneFilter",{});	
		}
		break;
		
	case "MatrixEleBut":	
		{
			cont.dx = 8.5;
			bubble_addtitle(cont,"Matrix element",{});
			if(bu.stat){
				bubble_addparagraph(cont,"<b>Mean:</b> "+bu.stat.mean,0,cont.dx);
				bubble_addparagraph(cont,"<b>CI min:</b> "+bu.stat.CImin,0,cont.dx);
				bubble_addparagraph(cont,"<b>CI max:</b> "+bu.stat.CImax,0,cont.dx);
				bubble_addparagraph(cont,"<b>ESS:</b> "+bu.stat.ESS,0,cont.dx);
				bubble_addparagraph(cont,"<b>GR:</b> "+bu.stat.GR,0,cont.dx);
			}
			else{
				if(bu.CImin){
					bubble_addparagraph(cont,"<b>Mean:</b> "+bu.value,0,cont.dx);
					bubble_addparagraph(cont,"<b>CI min:</b> "+bu.CImin,0,cont.dx);
					bubble_addparagraph(cont,"<b>CI max:</b> "+bu.CImax,0,cont.dx);
				}
				else{
					bubble_addparagraph(cont,"<b>Value:</b> "+bu.value,0,cont.dx);
				}
			}
			cont.y += 0.2;
			bubble_addparagraph(cont,"<b>x:</b> <e>"+bu.labels.xlab[bu.i]+"</e>",0,cont.dx);
			cont.y += 0.2;
			bubble_addparagraph(cont,"<b>y:</b> <e>"+bu.labels.ylab[bu.j]+"</e>",0,cont.dx);
		}
		break;
		
	case "Link":
		switch(bu.ac){
		case "AddConstantTime":
			cont.dx = 11;
			bubble_addtitle(cont,"Set column time",{});
			bubble_addradio(cont,0,"start","Start time",bub.radio);
			bubble_addradio(cont,0,"fixed","Fixed time",bub.radio);
		
			if(bub.radio.value == "fixed") bubble_input(cont,"Fixed Time:",{type:"fixed_time"});
			else bubble_input(cont,"",{hidden:true, type:"fixed_time"});
			add_end_button(cont,"OK","AddFixedtimeColumn",{});	
			break;
		
		case "AddConstantValue":
			cont.dx = 11;
			bubble_addtitle(cont,"Add constant value",{});
			bubble_input(cont,"Value:",{type:"const_column"});
			add_end_button(cont,"OK","AddConstantValue2",{});	
		}
		break;
		
		
	case "Settings":
		{
			switch(bu.ac){
			case "GraphSettings": inter.graph.settings_dist_bubble(cont); break;
			case "TraceSettings": inter.graph.settings_trace_bubble(cont); break;
			case "ScatterSettings": inter.graph.settings_scatter_bubble(cont); break;
			case "HistogramSettings": inter.graph.settings_histogram_bubble(cont); break;
			default: inter.graph.settings_speed_bubble(cont); break;
			}
		}
		break;
		
	case "AddFilter":
		filter_bubble(bu,cont);
		break;
		
	case "x-tick-label":
		cont.dx = 11;
		bubble_addtitle(cont,"X-axis range",{});
		bubble_input(cont,"Minimum:",{type:"min"});
		bubble_input(cont,"Maximum:",{type:"max"});
		add_end_button(cont,"Default","Defaultaxis",{});	
		add_end_button(cont,"OK","xaxisOK",{});	
		break;
		
	case "y-tick-label":
		cont.dx = 11;
		bubble_addtitle(cont,"Y-axis range",{});
		bubble_input(cont,"Minimum:",{type:"min"});
		bubble_input(cont,"Maximum:",{type:"max"});
		add_end_button(cont,"Default","Defaultaxis",{});	
		add_end_button(cont,"OK","yaxisOK",{});	
		break;
	
	case "RedInvalid":
		cont.dx = 11;
			
		bubble_addtitle(cont,"Invalid",{});
		bubble_addparagraph(cont,"Because of changes to the model this data source is no longer valid. Please delete.",0,11);
		break;
		
	case "Compartment": case "CompMap": case "CompLatLng":
		{
			cont.dx = 11.3;
			
			bubble_addtitle(cont,"Compartment",{te:compartment_text2});
			bubble_input(cont,"Name:",{type:"compartment", p:bu.p, cl:bu.cl, i:bu.i});
			
			bubble_colour(cont);
			cont.y += 0.3;
			
			let claa = model.species[bu.p].cla[bu.cl];
		
			let co = claa.comp[bu.i];
	
			if(claa.camera.grid == "on" && (co.type == "box" || co.type == "latlng")){
				cont.lay.add_button({te:"Coordinates:", x:0, y:cont.y, dx:cont.dx, dy:0.8, type:"InputBoxName"});
			
				let te;	
				if(co.type == "box"){
					te = "x: "+co.x+",  y: "+co.y;
				}
				else{
					let p = transform_latlng_inv(co.x,co.y,claa.camera,lay);
					te = "latitude: "+p.lat.toPrecision(5)+",  longitude: "+p.lng.toPrecision(5);
				}
				let tab = 5;
				bubble_addparagraph(cont,te,tab,cont.dx-tab);
				cont.y += 0.3;
			}
			
			if(co.choose_branch == true){
				bubble_addcheckbox(cont,-0.1,"Add branching probability",bub.checkbox);
			}
			
			if(co.branch == true){
				bubble_addcheckbox(cont,-0.1,"Use branching factors",bub.checkbox2);
			}
			
			if(co.type != "boundary"){
				bubble_addcheckbox(cont,-0.1,"Position fixed",co.fixed);
			}
			
			let sp = model.species[bu.p];
			if(sp.trans_tree.check == true){
				if(sp.infection_cl.te == sp.cla[bu.cl].name){
					bubble_addcheckbox(cont,-0.1,"Infected state",co.infected);
				}
			}
			
			if(co.type != "boundary"){
				add_end_button(cont,"Copy","CopyComp",{p:bu.p, cl:bu.cl, i:bu.i});		
			}
			add_end_button(cont,"Delete","DeleteComp",{p:bu.p, cl:bu.cl, i:bu.i});		
			add_end_button(cont,"OK","ChangeComp",{p:bu.p, cl:bu.cl, i:bu.i});		
		}
		break;
	
	case "Transition": case "TransitionPoint":
		{
			let tr = bu.tr;

			cont.dx = 12;
			
			if(tr.i == SOURCE){
				bubble_addtitle(cont,"Source",{te:source_text2});
			}
			else{
				if(tr.f == SINK){
					bubble_addtitle(cont,"Sink",{te:sink_text2});
				}
				else{
					bubble_addtitle(cont,"Transition",{te:trans_text2});
				}
			}
			
			cont.y += 0.1;
		
			if(tr.branch == true){
				if(tr.branch_select == true || tr.all_branches == true){
					let si = 1.1;
					if(tr.all_branches == true){
						bubble_input(cont,"Branching factor:",{type:"trans_bp", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
					}
					else{
						cont.lay.add_button({te:"Hel", x:11, y:cont.y+1.2, dx:si, dy:si, type:"Delete", back_col:BUBBLE_COL,  ac:"DeleteBranchProb", p:bu.p, cl:bu.cl, i:bu.i}); 
					
						bubble_input(cont,"Branching probability:",{type:"trans_bp", w:10.8, eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
					}
				}
				else{
					cont.lay.add_button({te:"Branching probability:", x:0, y:cont.y, dx:6, dy:0.8, type:"InputBoxName"});
					cont.y += 1.3;
					
					bubble_addparagraph(cont,"Set by other transitions.",0.5,cont.dx-0.5);
					
					cont.y += 0.4;
				}
				cont.y += 0.3;
			}
			cont.y += 0.3;
			
			let dist = get_dist_pos();
			if(tr.i == SOURCE) dist = source_dist_pos;
			
			bubble_addtext(cont,"Distribution:");
			cont.y -= 1.7;
			
			let pos=[]; for(let i = 0; i < dist.length; i++) pos.push({te:dist[i], trans_dist:true});
			
			bubble_adddropdown(cont,4.6,7.2,inter.bubble.trans_type,pos);

			cont.y += 0.5;
			
			let mean_fl = false, rate_fl = false, shape_fl = false;
			let cv_fl = false, scale_fl = false, shape_erlang_fl = false; 
			
			switch(inter.bubble.trans_type.te){
			case "exp(rate)":
				rate_fl = true;
				bubble_input(cont,"Rate:",{type:"trans_rate", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
				if(tr.i == SOURCE){ cv_fl = true; shape_fl = true; scale_fl = true; shape_erlang_fl = true;}
				break;
				
			case "exp(mean)":
				mean_fl = true;
				bubble_input(cont,"Mean:",{type:"trans_mean", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
				if(tr.i == SOURCE){ cv_fl = true; shape_fl = true; scale_fl = true; shape_erlang_fl = true;}
				break;
				
			case "gamma":
				mean_fl = true;
				bubble_input(cont,"Mean:",{type:"trans_mean", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
				cv_fl = true;
				bubble_input(cont,"CV:",{type:"trans_cv", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
				break;
				
			case "erlang":
				mean_fl = true;
				bubble_input(cont,"Mean:",{type:"trans_mean", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
				shape_erlang_fl = true;
				bubble_input(cont,"Shape:",{type:"trans_shape_erlang", p:bu.p, cl:bu.cl, i:bu.i});
				break;
				
			case "log-normal":
				mean_fl = true;
				bubble_input(cont,"Mean:",{type:"trans_mean", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
				cv_fl = true;
				bubble_input(cont,"CV:",{type:"trans_cv", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
				break;
				
			case "weibull":
				scale_fl = true;
				bubble_input(cont,"Scale:",{type:"trans_scale", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
				shape_fl = true;
				bubble_input(cont,"Shape:",{type:"trans_shape", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
				break;

			case "period":
				mean_fl = true;
				bubble_input(cont,"Time:",{type:"trans_period", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
				break;
			}
		
			if(mean_fl == false){
				bubble_input(cont,"",{hidden:true, type:"trans_mean", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
			}

			if(rate_fl == false){
				bubble_input(cont,"",{hidden:true, type:"trans_rate", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
			}
			
			if(cv_fl == false){
				bubble_input(cont,"",{hidden:true, type:"trans_cv", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});	
			}
			
			if(shape_fl == false){
				bubble_input(cont,"",{hidden:true, type:"trans_shape", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
			}
			
			if(scale_fl == false){
				bubble_input(cont,"",{hidden:true, type:"trans_scale", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
			}
			
			if(shape_erlang_fl == false){
				bubble_input(cont,"",{hidden:true, type:"trans_shape_erlang", p:bu.p, cl:bu.cl, i:bu.i});
			}
			
			add_end_button(cont,"Delete","DeleteTrans",{p:bu.p, cl:bu.cl, i:bu.i});			
			add_end_button(cont,"OK","ChangeTransValue",{p:bu.p, cl:bu.cl, i:bu.i});		
		}
		break;
		
	case "ClassTab":
		{
			cont.dx = 13.5;
			bubble_addtitle(cont,"Classification",{te:class_text});
		
			bubble_input(cont,"Name:",{type:"classification_name",p:bu.p, cl:bu.cl, w:10});
				
			cont.y -= 1.8+1.2;//1.8;
			bubble_addtext(cont,"Index:",10.7);
					
			cont.y -= 0.4;
			bubble_adddropdown(cont,10.7,3, bub.drop,bub.drop.pos);
				
			cont.y += 0.7;

			add_end_button(cont,"Delete","DeleteClass");
			add_end_button(cont,"Done","ChangeClassName");
		}
		break;
		
	case "AddButton":
		{	
			switch(bu.ac){
			case "AddClassification": 
				switch(bub.mode){
				case "clone":
					{
						if(bub.drop.pos.length == 0){
							cont.dx = 13;
							bubble_addtitle(cont,"Clone Classification",{});
							bubble_addparagraph(cont,"There are no potential classifications which can be cloned.",0,cont.dx);
						}
						else{
							cont.dx = bub.drop.w+0.2;
						
							bubble_addtitle(cont,"Clone Classification",{});
						
							bubble_addparagraph(cont,"Select classification to be cloned:",0,cont.dx);
							bubble_adddropdown(cont,0,bub.drop.w,bub.drop,bub.drop.pos);
						
							cont.y += 0.3;
						
							let ac; if(bub.drop.te != select_str) ac = "CloneClass";
							add_end_button(cont,"Clone",ac);
						}							
					}
					break;
					
				default:
					if(bub.op.type == "Popup"){
						cont.dx = 10;
						bubble_addparagraph(cont,"Click here to add a classification.",0,8);
					}
					else{
						cont.dx = 13.5;
						bubble_addtitle(cont,"Add classification",{te:class_text});
			
						let y_st = cont.y;
						
						bubble_addtext(cont,"Index:",10.7);
						
						cont.y -= 0.4;
						bubble_adddropdown(cont,10.7,3,bub.drop,bub.drop.pos);
						
						cont.y = y_st;
				
						bubble_input(cont,"Name:",{type:"add_classification_name",p:bu.op.p, w:10});
					
						bubble_addparagraph(cont,"Select coordinate system:",0,11);
						cont.y += 0.2;
						bubble_addradio(cont,0,"cartesian","Cartesian (x/y)",bub.radio);
						bubble_addradio(cont,0,"latlng","Geographic (Lat/Long)",bub.radio);
						if(bub.radio.value == "latlng"){
							bubble_addcheckbox(cont,1,"Load default map",bub.checkbox);
						}
						
						cont.y += 0.3;
						bubble_addlink(cont,"Clone classification","CloneClassInit",{p:bu.op.p});
						cont.y -= 1.3;
						
						add_end_button(cont,"Add","AddClassificationOK");				
					}
					break;
				}
				break;
				
			case "AddCompartment":
				if(bub.op.type == "Popup"){
					cont.dx = 10;
					bubble_addparagraph(cont,"Click here to add a compartment.",0,8);
				}
				break;
			
			case "AddTransition":
				if(bub.op.type == "Popup"){
					cont.dx = 10;
					bubble_addparagraph(cont,"Click here to add a transition.",0,8);
				}
				break;
	
			case "ImportModel":
				{
					switch(bub.mode){
					case "LoadComps":
						{
							cont.dx = 11.5;
							bubble_addtitle(cont,"Load compartments",{te:loadcomp_text});
							
							bubble_addradio(cont,0,"SelectColour","Select colour",bub.radio);
							bubble_colour(cont);
							cont.y += 0.3;

							bubble_addradio(cont,0,"FileColour","Load colour from file",bub.radio);
							
							cont.y += 0.5;
							
							let op = bub.bu.op;
							if(model.species[op.p].cla[op.cl].camera.coord == "cartesian"){
								bubble_addcheckbox(cont,0,"Load position data",bub.pos_check);
							}

							add_end_button(cont,"Back","MenuBack");		
							add_end_button(cont,"Next","LoadCompartmentNext");		
						}
						break;
						
					case "LoadCompMaps":
						{
							cont.dx = 10.8;
							bubble_addtitle(cont,"Load comp. map",{te:loadcomp_text});
							
							bubble_colour(cont,0);
							cont.y += 0.5;
							
							add_end_button(cont,"Back","MenuBack");		
							add_end_button(cont,"Next","ImportCompMap2");
						}			
						break;
					
					
					case "LoadTrans":
						{
							cont.dx = 12;
							bubble_addtitle(cont,"Load transitions",{te:loadtrans_text});
								
							bubble_addparagraph(cont,"Select if position data is included in data file:",0,cont.dx);
							cont.y += 0.3;
							
							bubble_addcheckbox(cont,0,"Load position data",bub.pos_check);
							
							add_end_button(cont,"Back","MenuBack");		
							add_end_button(cont,"Load","ImportTrans2");		
						}
						break;

					default:
						cont.dx = 10;
						bubble_addtitle(cont,"Import",{te:label_text});
						bubble_addradio(cont,0,"Compartments","Compartments",inter.bubble.radio); 
						let claa = model.get_cla();
						if(claa.camera.coord == "latlng"){
							bubble_addradio(cont,0,"Comp. Map","Comp. Map",inter.bubble.radio); 
						}
						bubble_addradio(cont,0,"Transitions","Transitions",inter.bubble.radio); 
					
						add_end_button(cont,"Next","ImportModelNext");	
					}
				}
					
			/*
				let sub = [];
				let active = cont.bu.show_model;
				sub.push({te:"Compartments", ac:"ImportComp", active:active});
				let active2 = active; if(active2 == true && cont.bu.coord != "latlng") active2 = false;
				sub.push({te:"Comp. Map", ac:"ImportCompMap", active:active2});			
				*/			
			
				break;
			
			case "AddAnnotation":
				cont.dx = 12;
				switch(bub.mode){
				case "add_label":
					bubble_addtitle(cont,"Add label",{te:label_text});
					bubble_input(cont,"Text:",{type:"label"});
					bubble_input(cont,"Text size:",{type:"label_size"});
					bubble_colour(cont);
					cont.y += 0.3;
			
					add_end_button(cont,"Add","AddLabelOK");		
					break;
				
				case "remove":
					let rem = bub.rem_annotation;
					cont.dx = rem.w;
					bubble_addtitle(cont,"Edit annotations",{te:edit_anno_text});
					if(rem.list.length == 0){
						bubble_addparagraph(cont,"There are currently no annotations.",0,cont.dx);
					}
					else{
						bubble_addparagraph(cont,"Delete selected annotations:",0,cont.dx);
						cont.y += 0.2;
						bubble_addscrollable(cont,{type:"annotation", list:rem.list, ymax:10});
						cont.y += 0.4;
						add_end_button(cont,"Next","DeleteAnnotation");
					}
					break;
				
				default:
					bubble_addtitle(cont,"Annotation",{te:anno_text});
					bubble_addradio(cont,0,"Label","Add label",bub.radio);
					bubble_addradio(cont,0,"Box","Add bounding box",bub.radio);
					if(model.get_coord() == "latlng"){
						bubble_addradio(cont,0,"Map","Add map (.geoJSON)",bub.radio);
					}
					//bubble_addradio(cont,0,"Remove","Edit annotations",bub.radio);
					add_end_button(cont,"Next","AddAnnotationOK");				
					break;
				}
				break;
				
			case "AddInitialPopulation": initpop_data_bubble(cont,"add"); break;
			case "AddInitPopPrior": initpopprior_data_bubble(cont,"add"); break;
			case "MoveIndividuals": move_data_bubble(cont,"add"); break;
			case "CompData": comp_data_bubble(cont,"add"); break;
			case "TransData": trans_data_bubble(cont,"add","TransVar"); break;
			case "SourceData": trans_data_bubble(cont,"add","SourceVar"); break;
			case "SinkData": trans_data_bubble(cont,"add","SinkVar"); break;
			case "DiagTestData": diagtest_data_bubble(cont,"add"); break;
			case "PopulationData": population_data_bubble(cont,"add"); break;
			case "PopTransData": poptrans_data_bubble(cont,"add"); break;
			case "SeqData": sequence_data_bubble(cont,"add"); break;
			case "IndEffData": ind_eff_data_bubble(cont,"add"); break;
			case "IndGroupData": ind_group_data_bubble(cont,"add"); break;
			case "SetConstant": set_constant_bubble(cont); break;
			case "SetFactor": set_factor_bubble(cont); break;
			case "SetReparam": set_reparam_bubble(cont); break;
			case "SetDistribution": set_distribution_bubble(cont); break;
			case "SetDerived": set_derived_bubble(cont); break;
			case "AddParamMult": add_param_mult_bubble(cont); break;
			default: error("Cannot find source2"); break;
			}
		}
		break;
	
	case "GreyView":
		if(bu.ac == "CombineIEBubble") combine_ie_bubble(cont);
		else{	
			switch(edit_source.type){
			case "Init. Pop.": initpop_data_bubble(cont,"view"); break;
			case "Move Ind.": move_data_bubble(cont,"view"); break;
			case "Compartment": comp_data_bubble(cont,"view"); break;
			case "Transition": trans_data_bubble(cont,"view","TransVar"); break;
			case "Diag. Test": diagtest_data_bubble(cont,"view"); break;
			case "Ind. Eff.": ind_eff_data_bubble(cont,"view"); break;
			case "Ind. Group": ind_group_data_bubble(cont,"view"); break;
			case "Genetic": sequence_data_bubble(cont,"view"); break;
			case "Population": population_data_bubble(cont,"view"); break;
			case "Pop. Trans.": poptrans_data_bubble(cont,"view"); break;
			default: error("Cannot find source: "+edit_source.type); break;
			}
		}
		break;
	
	case 	"AddSpecies":
		if(bub.op.type == "Popup"){
			cont.dx = 9;
			bubble_addparagraph(cont,"Click here to add a species to the model.",0,8);
		}
		else{		
			cont.dx = 12;
			bubble_addtitle(cont,"Add species",{te:species_text});
			bubble_input(cont,"Name:",{type:"add_species_name"});
			add_species_bubble(cont);
			add_end_button(cont,"Add","AddSpeciesOK");	
		}			
		break;
	
	case "Menu":
		drop_menu(cont,lay); 
		break;
			
	case "PageSubSub":
		cont.dx = 12;
		bubble_addtitle(cont,"Species",{te:species_text});
		bubble_input(cont,"Name:",{type:"species_name",p:bu.val});
		add_species_bubble(cont);
		add_end_button(cont,"Delete","DeleteSpecies");
		add_end_button(cont,"Done","ChangeSpeciesName");		
		break;
	
	case "CompPop":
		cont.dx = 10;
		switch(bu.ac){
		case "EditCompPop":
			bubble_addtitle(cont,"Set population",{te:setpop_text});
			bubble_input(cont,"Population:",{type:"init_population", cl:bu.cl, c:bu.c});
			add_end_button(cont,"Done","Done");
			break;
			
		case "EditCompPopPercent":
			bubble_addtitle(cont,"Set percentage",{te:setpopper_text});
			bubble_input(cont,"Percent:",{type:"init_per", cl:bu.cl, c:bu.c});
			add_end_button(cont,"Done","Done");
			break;

		default: error("Option not recognised 8"); break;
		}
		break;
	
	case "MultiComp":
		cont.dx = 10;
		bubble_addtitle(cont,"Set population",{te:setmultipop_text});
		bubble_input(cont,"Population:",{type:"init_globpopulation", c:bu.c});
		add_end_button(cont,"Done","Done");
		break;
	
	case "ParamElement":
		cont.dx = 10;
		bubble_addtitle(cont,"Edit value",{te:editconstparam_text});
		bubble_input(cont,"Value:",{type:"element_param", pindex:bu.pindex});
		add_end_button(cont,"Done","Done");	
		break;
		
	case "ParamElementConst":
		cont.dx = 10;
		bubble_addtitle(cont,"Edit constant",{te:editconstparam_text});
		bubble_input(cont,"Value:",{type:"element_param_const", pindex:bu.pindex});
		add_end_button(cont,"Done","Done");	
		break;
		
	case "ParamFactorConst":
		cont.dx = 10;
		bubble_addtitle(cont,"Edit factor",{te:editfactor_text});
		bubble_input(cont,"Value:",{type:"element_factor_const", pindex:bu.pindex});
		add_end_button(cont,"Done","Done");	
		break;
		
	case "ParamWeightConst":
		cont.dx = 10;
		bubble_addtitle(cont,"Edit weight",{te:editweight_text});
		bubble_input(cont,"Value:",{type:"element_weight_const", pindex:bu.pindex});
		add_end_button(cont,"Done","Done");	
		break;
		
	case "AmatrixElement":
		cont.dx = 10;
		bubble_addtitle(cont,"Edit value",{te:editamatrix_text});
		bubble_input(cont,"Value:",{type:"element_Amatrix", i:bu.i, j:bu.j});
		add_end_button(cont,"Done","Done");	
		break;
		
	case "XvectorElement":
		cont.dx = 10;
		bubble_addtitle(cont,"Edit value",{te:editxvector_text});
		bubble_input(cont,"Value:",{type:"element_Xvector", i:bu.i, j:bu.j});
		add_end_button(cont,"Done","Done");	
		break;
	
	case "ReparamElement":
		cont.dx = 15;
		bubble_addtitle(cont,"Edit equation",{te:editeqnparam_text});	
		bubble_input(cont,"Equation:",{type:"reparam_eqn", eqn:true});
		add_end_button(cont,"Done","DoneReparamEquation");	
		break;
		
	case "ReparamEqn":
		cont.dx = 15;
		bubble_addtitle(cont,"Edit equation",{te:editeqnparam_text});	
		bubble_input(cont,"Equation:",{type:"reparam_equation", eqn:true});
		add_end_button(cont,"Done","DoneReparamEquation");	
		break;
		
	case "ReparamTableElement":
		cont.dx = 15;
		bubble_addtitle(cont,"Edit equation",{te:editeqnparam_text});	
		bubble_input(cont,"Equation:",{type:"element_eqn", eqn:true, pindex:bu.pindex});
		add_end_button(cont,"Done","Done");	
		break;
	
	case "ParamSimElement": case "DistSimElement":
		cont.dx = 10;
		bubble_addtitle(cont,"Edit value",{te:editsimparam_text});
		bubble_input(cont,"Value:",{type:"param_val", val:bu.i});
		add_end_button(cont,"Done","Done");	
		break;
		
	case "PriorElement": case "DistElement": 
	case "PriorSplitElement": case "DistSplitElement":
	case "CompPrior":
		{
			let ac;
			let ae = false;                              // Determines if allow eqn
		
			cont.dx = 8;
			let pp;
			if(bu.type == "CompPrior"){
				bubble_addtitle(cont,"Pop. prior",{title:"Population prior", te:editpopprior_text});
						
				pp = prior_pos_positive;
				ac = "DoneCompPrior";
			}
			else{
				switch(bu.type){
				case "PriorElement": case "PriorSplitElement":
					{
						ac = "DonePrior";
						let te = editprior_text;
					
						bubble_addtitle(cont,"Edit prior",{te:te});
						
						if(bu.type == "PriorSplitElement") ac = "DonePriorSplit";
					}
					break;
					
				case "DistElement": case "DistSplitElement":
					ac = "DoneDist";
					ae = true;
					bubble_addtitle(cont,"Edit dist.",{title:"Edit Distribution", te:distribution_text});
					if(bu.type == "DistSplitElement") ac = "DoneDistSplit";
					break;
				
				default: error("Option not recognised 1A"); break;
				}

				pp = model.param[bu.i].pri_pos;
			}
			
			let pos=[]; 
			for(let i = 0; i < pp.length; i++) pos.push({te:pp[i]});
		
			cont.y += 0.2;
			bubble_addtext(cont,"Distribution:");
			cont.y -= 0.5;
		
			bubble_adddropdown(cont,0,8,inter.bubble.prior.type,pos);
		
			let str = "prior_"; if(ae == true) str += "dist_";
			
			let min_fl = false, max_fl = false, mean_fl = false;
			let sd_fl = false, cv_fl = false, alpha_fl = false, beta_fl = false; 
		
			cont.y += 0.2;
			switch(inter.bubble.prior.type.te){
			case "fix":
				bubble_input(cont,"Value:",{type:str+"mean", val:bu.i, eqn:ae});
				mean_fl = true; 
				break;
				
			case "uniform":		
				bubble_input(cont,"Min:",{type:str+"min", val:bu.i, eqn:ae});
				bubble_input(cont,"Max:",{type:str+"max", val:bu.i, eqn:ae});
				min_fl = true; max_fl = true; 
				break;
				
			case "exp":
				bubble_input(cont,"Mean:",{type:str+"mean", val:bu.i, eqn:ae});
				mean_fl = true; 
				break;
			
			case "normal":
				bubble_input(cont,"Mean:",{type:str+"mean", val:bu.i, eqn:ae});
				bubble_input(cont,"sd:",{type:str+"sd", val:bu.i, eqn:ae});
				mean_fl = true; sd_fl = true;
				break;
				
			case "gamma":
				bubble_input(cont,"Mean:",{type:str+"mean", val:bu.i, eqn:ae});
				bubble_input(cont,"CV:",{type:str+"cv", val:bu.i, eqn:ae});
				mean_fl = true; cv_fl = true;
				break;
				
			case "log-normal":
				bubble_input(cont,"Mean:",{type:str+"mean", val:bu.i, eqn:ae});
				bubble_input(cont,"CV:",{type:str+"cv", val:bu.i, eqn:ae});
				mean_fl = true; cv_fl = true;
				break;
			
			case "beta":
				bubble_input(cont,"alpha:",{type:str+"alpha", val:bu.i, eqn:ae});
				bubble_input(cont,"beta:",{type:str+"beta", val:bu.i, eqn:ae});
				alpha_fl = true; beta_fl = true;
				break;
				
			case "bernoulli":
				bubble_input(cont,"Mean:",{type:str+"mean", val:bu.i, eqn:ae});
				mean_fl = true; 
				break;
			
			case "dirichlet":
				bubble_input(cont,"Alpha:",{type:str+"alpha", val:bu.i, eqn:ae});
				alpha_fl = true; 
				break;
			
			case "mdir":
				bubble_input(cont,"Sigma:",{type:str+"sigma", val:bu.i, eqn:ae});
				alpha_fl = true; 
				break;
				
			case "flat":
				break;
				
			default:
				cont.y += 0.4;
				break;
			}
		
			if(min_fl == false){
				bubble_input(cont,"",{hidden:true, type:str+"min", val:bu.i, eqn:ae});
			}
			
			if(max_fl == false){
				bubble_input(cont,"",{hidden:true, type:str+"max", val:bu.i, eqn:ae});
			}
			
			if(mean_fl == false){
				bubble_input(cont,"",{hidden:true, type:str+"mean", val:bu.i, eqn:ae});
			}
			
			if(cv_fl == false){
				bubble_input(cont,"",{hidden:true, type:str+"cv", val:bu.i, eqn:ae});
			}

			if(sd_fl == false){
				bubble_input(cont,"",{hidden:true, type:str+"sd", val:bu.i, eqn:ae});
			}
			
			if(alpha_fl == false){
				bubble_input(cont,"",{hidden:true, type:str+"alpha", val:bu.i, eqn:ae});
			}
			
			if(beta_fl == false){
				bubble_input(cont,"",{hidden:true, type:str+"beta", val:bu.i, eqn:ae});
			}
			
			add_end_button(cont,"Done",ac);	
		}
		break;
		
	case "SplineKnots":
		cont.dx = 25;
		bubble_addtitle(cont,"Spline knots",{te:knot_text});
		
		bubble_addradio(cont,0,"knot_time","Times",bub.knot_radio); cont.y -= 1.4;
		bubble_addradio(cont,7,"time_step","Time-step",bub.knot_radio); cont.y -= 1.4;
		bubble_addradio(cont,14,"load","Load",bub.knot_radio);
		
		switch(bub.knot_radio.value){
		case "knot_time":
			bubble_input(cont,"Knot times:",{type:"knot_times", val:bu.i, eqn:false});
			add_end_button(cont,"Done","DoneKnots");	
			break;
			
		case "time_step":
			{
				let ww = 7.7, gap = 0.7;
				bubble_input(cont,"",{hidden:true, type:"knot_times", val:bu.i});
				bubble_input(cont,"Start time:",{type:"knot_t_start", no_down:true, w:ww, eqn:false});
				bubble_input(cont,"End time:",{type:"knot_t_end", no_down:true, x:ww+gap, w:ww, eqn:false});
				bubble_input(cont,"Time-step:",{type:"knot_dt", x:2*(ww+gap), w:ww, eqn:false});
				add_end_button(cont,"Done","DoneKnots");	
			}
			break;
			
		case "load":
			cont.y += 0.2;
			bubble_addparagraph(cont,"Load from the column of a table",0,cont.dx);
			cont.y += 0.2;
			add_end_button(cont,"Upload","UploadKnots",{i:bu.i});	
			break;
			
		default: error("Option not recog"); break;
		}
		break;
		
	case "Element":
		switch(bu.ac){
		case "EditTableElement":
			if(bub.mode == "SearchResult"){
				cont.dx = 8;
				bubble_addtitle(cont,"Search",{te:search_text});
				bubble_addparagraph(cont,"Result "+(inter.bubble.search_select+1)+" / "+inter.bubble.row_find.length,0,8);
				cont.y += 0.3;
				add_end_button(cont,"Back","SearchBack");
				add_end_button(cont,"Next","SearchNext");
			}
			else{
				let te = edit_source.table.ele[bu.r][bu.c];
				if(typeof(te) == "object"){
					cont.dx = 8;
					bubble_addtitle(cont,"Boundary",{te:boundary_text});
					bubble_addparagraph(cont,	"This stores polygon boundary data which cannot be editted.",0,8);
					add_end_button(cont,"Close","CloseBubble");
				}
				else{
					cont.dx = 10;
					bubble_addtitle(cont,"Edit element",{te:editele_text});
					bubble_input(cont,"Value:",{type:"element", c:bu.c, r:bu.r});
					add_end_button(cont,"Done","Done");
				}
			}
			break;
			
		case "EditTableHead":
			switch(bub.mode){
			case "Search":
				cont.dx = 10;
				bubble_addtitle(cont,"Search",{te:search_text});
				bubble_input(cont,"Find:",{type:"find"});
				add_end_button(cont,"Search","DoneSearch");
				break;
			
			case "Replace":
				cont.dx = 10;
				bubble_addtitle(cont,"Replace",{te:replace_text});
				bubble_input(cont,"Find:",{type:"find"});
				bubble_input(cont,"Replace:",{type:"replace"});
				add_end_button(cont,"Replace","DoneReplace");
				break;
			
			case "ReplaceResult":
				{
					cont.dx = 10;
					bubble_addtitle(cont,"Replace",{te:replace_text});
					let te = "Replaced "+inter.bubble.num+" occurrences";
					if(inter.bubble.num > 1) te += ".";
					bubble_addparagraph(cont,te,0,10);
				}
				break;
			
			case "Order":
				cont.dx = 10;
				bubble_addtitle(cont,"Order column",{te:order_text});
				bubble_addradio(cont,0,"A-Z","Alphabetically (A-Z)",inter.bubble.order_radio);			
				bubble_addradio(cont,0,"Z-A","Alphabetically (Z-A)",inter.bubble.order_radio);
				bubble_addradio(cont,0,"Low-High","Numerically (Low-High)",inter.bubble.order_radio);
				bubble_addradio(cont,0,"High-Low","Numerically (High-Low)",inter.bubble.order_radio);
				add_end_button(cont,"Order","DoneOrder");	
				break;
			
			case "Delete":
				cont.dx = 10;
				bubble_addtitle(cont,"Delete rows",{te:deleterows_text});
				bubble_input(cont,"Find:",{type:"find"});
				add_end_button(cont,"Delete","DoneDeleteRows");
				break;
		
			case "DeleteRowsResult":
				{
					cont.dx = 10;
					bubble_addtitle(cont,"Delete rows",{te:deleterows_text});
					let te = inter.bubble.num+" row"; if(inter.bubble.num > 1) te += "s";
					te += " deleted"; 
					bubble_addparagraph(cont,te,0,10);
				}
				break;

			default:
				let dx = 4, dy = 1.2, gap = 0.4;
				cont.dx = 2*dx+gap;

				bubble_addtitle(cont,"Edit column",{te:editcol_text});

				lay.add_button({te:"Search", x:0, y:cont.y, dx:dx, dy:dy, ac:"Search", type:"BubbleEndBut"});
			
				lay.add_button({te:"Replace", x:dx+gap, y:cont.y, dx:dx, dy:dy, ac:"Replace", type:"BubbleEndBut"});
				
				cont.y += 1.6;
				
				lay.add_button({te:"Order", x:0, y:cont.y, dx:dx, dy:dy, ac:"Order", type:"BubbleEndBut"});
				
				lay.add_button({te:"Delete", x:dx+gap, y:cont.y, dx:dx, dy:dy, ac:"DeleteRows", type:"BubbleEndBut"});
				break;
			}
		}
		break;
	
	case "RowNumber":
		cont.dx = 8;
		bubble_addtitle(cont,"Table row");
		bubble_addparagraph(cont,"Delete this row from the data table?",0,8);
		cont.y += 0.3;
		add_end_button(cont,"Delete","DeleteRow");
		break;

	case "LabelText":
		cont.dx = 12;
		bubble_addtitle(cont,"Label",{te:label_text});
		bubble_input(cont,"Text:",{type:"label_anno",p:bu.p, cl:bu.cl, i:bu.i});
		bubble_input(cont,"Text size:",{type:"label_anno_size",p:bu.p, cl:bu.cl, i:bu.i});
		bubble_colour(cont);	
		cont.y += 0.2;
		add_end_button(cont,"Delete","DeleteAnno",{p:bu.p, cl:bu.cl, i:bu.i});		
		add_end_button(cont,"Done","LabelOK");		
		break;
		
	case "Box":
		cont.dx = 10;
		bubble_addtitle(cont,"Box",{te:box_text});
		bubble_input(cont,"Text:",{type:"label_anno", p:bu.p, cl:bu.cl, i:bu.i});
		bubble_input(cont,"Text size:",{type:"label_anno_size",p:bu.p, cl:bu.cl, i:bu.i});
		bubble_colour(cont);	
		cont.y += 0.2;
		add_end_button(cont,"Delete","DeleteAnno",{p:bu.p, cl:bu.cl, i:bu.i});		
		add_end_button(cont,"OK","LabelOK",{p:bu.p, cl:bu.cl, i:bu.i});	
		break;
	
	case "PlotEquation":
		cont.dx = 15;
		bubble_addtitle(cont,"Set derived equation",{te:derived_text});
		
		bubble_addparagraph(cont,"Set how parameter depends on others:",0,cont.dx);
		cont.y += 0.5;
	
		bubble_input(cont,"",{type:"derive_eqn", val:bu.val, eqn:true});
		add_end_button(cont,"OK","DeriveOK",{});	
		break;
	
	case "Derive": alter_derived_bubble(cont,bu); break;
	
	case "SmoothValue":
		cont.dx = 12;
		bubble_addtitle(cont,"Smoothing",{te:smoothing_text});
		bubble_input(cont,"Value:",{type:"smooth_value", i:bu.i});
		
		add_end_button(cont,"Done","SmoothValueDone",{});	
		break;
		
	case "ParamLabel2":
		set_paramselect_bubble(cont,bu.eqn_appear);
		break;
	
	case "PopSlice":
		pop_slice_bubble();
		break;
		
	default: error(bu.type+": Bubble Not recognised"); break;
	}

	setup_bubble_back(cont);
	
	if(bub.find_focus == true){
		set_focus_first();
	
		bub.find_focus = false;
	}
}


/// Displays the species bubble 
function add_species_bubble(cont)
{
	let bub = inter.bubble;

	bubble_addtext(cont,"Type:");
	cont.y -= 1.4;
	
	bubble_addradio(cont,2.3,"Population","Population-based",bub.radio);
	bubble_addradio(cont,2.3,"Individual","Individual-based",bub.radio);
	cont.y += 0.2;
	
	if(bub.radio.value == "Individual"){
		bubble_addcheckbox(cont,-0.1,"Calculate transmission tree",bub.checkbox);
	
		if(bub.checkbox.check == true && bub.infection_cl != undefined && bub.infection_cl_pos.length > 0){
			bubble_addtext(cont,"Infection classification:");	
			cont.y -= 0.4;
			bubble_adddropdown(cont,1,cont.dx-2,bub.infection_cl,bub.infection_cl_pos);
			
			cont.y += 0.4;
		}
	}
}

			
/// Determines if a bubble is on
function bubble_on()
{
	if(inter.bubble.lay_name != undefined) return true;
	return false;
}


/// Determines if a button is selected
function selected(bu)
{
	let bb = inter.bubble.bu;
	if(bb != undefined){
		if(bb.x == bu.x && bb.y == bu.y && bb.dx == bu.dx && bb.dy == bu.dy) return true;
	}
	return false;
}


/// Determines if a button is in selected column
function selected_column(bu)
{
	let bb = inter.bubble.bu;
	if(bb != undefined){
		if(bb.x == bu.x && bb.dx == bu.dx) return true;
	}
	return false;
}

			
/// Selects the button the mouse is over and adds a speech bubble 
function select_bubble_over()
{
	let l = inter.over.layer;
	let i = inter.over.i;
	
	if(l == undefined || i == undefined){ error("Select bubble"); return;}

	select_bubble(inter.layer[l].name,i);
}


/// Selects a bubble based on a layer name and a button property
function select_bubble_but_prop(lay_name,type,info,value,op)
{
	close_bubble();
	generate_screen();
		
	let lay = get_lay(lay_name);

	for(let i = 0; i < lay.but.length; i++){
		let bu = lay.but[i];
		if(bu.type == type){
			let j = 0;
			while(j < info.length && bu[info[j]] != undefined){ bu = bu[info[j]]; j++;}
			if(j == info.length){
				if(bu == value){
					select_bubble(lay_name,i,op);
					return;
				}
			}
		}
	}
	
	error("Could not find");
}


/// If there is a bubble OK button then press it
function press_bubble_OK()
{
	if(inter.bubble.lay_name == undefined) return;
	let lay = get_lay("Bubble");
	
	for(let i = 0; i < lay.but.length; i++){
		let bu = lay.but[i];
		if(bu.type == "BubbleEndBut" && (bu.te == "OK" || bu.te == "Done")){
			activate_button(lay,i);
			return;
		}
	}
	
	error("Could not find");
}


/// Selects a bubble based on a layer name and a button property
function press_button_prop(lay_name,type,info,value,op)
{
	close_bubble();
	generate_screen();
	
	if(loading_sym()) return;
	
	let lay = get_lay(lay_name);
	if(lay == undefined){
		error("Could not find layer '"+lay_name+"'");
		print_layer();
		return;
	}
	
	for(let i = 0; i < lay.but.length; i++){
		let bu = lay.but[i];
		if(bu.type == type){
			let j = 0;
			while(j < info.length && bu[info[j]] != undefined){ bu = bu[info[j]]; j++;}
			if(j == info.length){
				if(bu == value){ activate_button(lay,i); return;}
			}
		}
	}
	
	error("Could not find '"+type+"'");
}
	
	
/// Selects a button and adds a speech bubble 
function select_bubble(lay_name,i,op)
{ 
	if(op == undefined) op = {};

	let lay = get_lay(lay_name);
	let bu = lay.but[i];
	
	scroll_to_view(lay,bu);
	
	inter.bubble = { lay_name:lay_name, bu:bu, i:i, show_warning:false, warning:false, op:op, find_focus:true, check_radio_press:false};
	
	switch(bu.type){
	case "Transition": case "TransitionPoint":
		inter.bubble.value = bu.tr.value;
		inter.bubble.trans_type = {te:bu.tr.type};
		break;
	}
	
	if(subtab_name() != "Description") reset_text_box();
}


/// Changes the mode for a bubble
function change_bubble_mode(mode)
{
	inter.bubble.mode = mode;
	reset_text_box();
	inter.bubble.find_focus = true;
}


/// Closes a speech bubble
function close_bubble()
{
	if(inter.help.title != undefined && inter.equation.te != undefined){
		error("Cannot close bubble");
	}

	turn_off_cursor();

	inter.bubble = {};
}


/// Adds buttons to the back layer of a bubble
function add_bubble_back_buts(lay)
{
	lay.add_button({type:"BubbleBack"});
	lay.add_button({ac:"Nothing", type:"Nothing"});
	lay.add_button({dx:1, dy:1, ac:"CloseBubble", type:"CloseBubble"});
}


/// Adds a button to the end of the bubble 
function add_end_button(cont,te,ac,op)
{
	cont.end_button.push({te:te, ac:ac, op:op});	
}


/// Adds a termination point to a bubble
function add_bubble_end(cont)
{
	cont.lay.add_button({x:0,y:cont.y, type:"Nothing"});
}


/// Adds a mini title 
function bubble_add_minititle(cont,te,xx)
{
	let w = cont.dx;
	let x = 0; if(xx != undefined){ x = xx; w -= x;}
	cont.lay.add_button({te:te, x:x, y:cont.y, dx:w, dy:0.8, type:"InputBoxName"});
	cont.y += 1.4;	
}


/// Adds an input text box to the bubble		
function bubble_input(cont,te,op)
{	
	let yst = cont.y;
	
	let w = cont.dx;
	let x = 0; if(op.x != undefined){ x = op.x; w -= x;}
	if(op.w != undefined) w = op.w;
	 
	add_ref(op);
	
	if(op.hidden == true){
		cont.input_list.push({x:x, y:-1000, dx:cont.dx, dy:0.8, op:op});
	}
	else{
		if(te != ""){
			cont.lay.add_button({te:te, x:x, y:cont.y, dx:w, dy:0.8, type:"InputBoxName"});
			cont.y += 1.0;
		}
		
		cont.input_list.push({x:x+0.3, y:cont.y, dx:w-0.6, op:op});
		
		cont.y += 1.8;
		
		let sto = inter.textbox_store;

		let warn;
		let k = find(sto,"ref",op.ref);
		if(k != undefined) warn = sto[k].warning;
	
		if(warn == undefined && inter.bubble.error_warning != undefined){
			warn = inter.bubble.error_warning;
		}
		
		if(warn != "" && warn != undefined){
			cont.y = cont.lay.add_paragraph(warn,w,x,cont.y-0.2,RED,warn_si,warn_lh,undefined,"center");
			cont.y += 0.2;
		}				
	}
	
	if(op.no_down == true) cont.y = yst;
}


/// Adds an box which shows a previously input value	
function bubble_set_input(cont,te,te2)
{	
	let w = cont.dx;
	cont.lay.add_button({te:te, x:0, y:cont.y, dx:w, dy:0.8, type:"InputBoxName"});
	
	let x = text_width(te,get_font(0.7,"Bold"))+0.3;
	
	let si = 0.7;
	cont.lay.add_button({te:te2, x:x, y:cont.y, dx:w-x, dy:0.8, si:si, font:get_font(si), type:"Text", col:BLACK});
	
	cont.y += 1.2;
}


/// Adds two input text box to the bubble		
function bubble_double_input(cont,te,op,te2,op2)
{	
	if(op.hidden == true){
		cont.input_list.push({x:x, y:-1000, dx:cont.dx, dy:0.8, op:op});
		cont.input_list.push({x:x, y:-1000, dx:cont.dx, dy:0.8, op:op2});
	}
	else{
		let warn = false, ywarn = -LARGE;
		let x = 0, gap = 1;
		let w = (cont.dx-gap)/2;
		for(let i = 0; i < 2; i++){
			let tex = te, opp = op; if(i == 1){ tex = te2; opp = op2;}
			cont.lay.add_button({te:tex, x:x, y:cont.y, dx:w, dy:0.8, type:"InputBoxName"});
			cont.input_list.push({x:x+0.3, y:cont.y+1.2, dx:w-0.6, op:opp});
		
			let k = cont.input_list.length-1;
			let sto = inter.textbox_store;
	
			if(k < sto.length && sto[k].warning != undefined && sto[k].warning != ""){
				let y = cont.lay.add_paragraph(sto[k].warning,w,x,cont.y+1.4+1.8-0.3,RED,warn_si,warn_lh)+0.1;
				if(y > ywarn) ywarn = y;
				warn = true;
			}
			x += w+gap;
		}	
		
		cont.y += 1.4+1.8;		
		if(warn == true){ cont.y = ywarn;}
	}
}


/// Adds a picture to the bubble
function bubble_addpicture(cont,pic)
{
	let dx = cont.dx;
	let dy = pic.height*dx/pic.width;
	cont.lay.add_button({pic:pic, x:0, y:cont.y, dx:dx, dy:dy, type:"BubblePic"});
	cont.y += dy+0.2;
}


/// Adds a title to the bubble
function bubble_addtitle(cont,te,op)
{
	cont.lay.add_button({te:te, x:0, y:cont.y, dx:cont.dx, dy:1, type:"BubbleTitle"});
	
	if(op && op.te){
		let title = te; if(op.title) title = op.title;
		let w = text_width(te,get_font(si_bubble_title,"bold"));
		cont.lay.add_help_button(w+0.1,cont.y+1.3,{title:title, te:op.te, back_col:BUBBLE_COL});
	}
		
	cont.y += 1.4;
}


/// Adds a dropdown menu to the bubble
function bubble_adddropdown(cont,tab,w,source,pos)
{
	cont.drop_list.push({x:tab, y:cont.y, dx:w, source:source, pos:pos});
	cont.y += 1.4;
}


/// Adds a scrollable box for putting content in
function bubble_addscrollable(cont,op)
{
	cont.scrollable_y = cont.y; 

	add_layer("BubbleScrollable",0,0,cont.dx,0,op);
	let lay = get_lay("BubbleScrollable");

	cont.y += lay.dy;
}


/// Adds potentially scrollable region to layer
function add_bubble_scrollable_buts(lay)
{
	let cy;
	
	lay.background = BUBBLE_COL;

	switch(lay.op.type){
	case "filterpos": cy = filterpos_scrollable(lay); break;
	case "comp list": cy = diagtest_scrollable(lay); break;
	case "ind list": cy = individual_scrollable(lay); break;
	case "pop list": cy = population_scrollable(lay); break;
	case "poptrans list": cy = poptrans_scrollable(lay); break;
	case "annotation": cy = annotation_scrollable(lay); break;
	case "const param sel": cy = param_sel_scrollable(lay,"const"); break;
	case "reparam param sel": cy = param_sel_scrollable(lay,"reparam"); break;
	case "dist param sel": cy = param_sel_scrollable(lay,"dist"); break;
	case "fac param sel": cy = param_sel_scrollable(lay,"fac"); break;
	case "trans param sel": cy = param_sel_scrollable(lay,"trans"); break;
	case "param details": cy = param_details_scrollable(lay); break;
	case "combine IE": cy = combineIE_scrollable(lay); break;
	default: error("Option not recognised 9"); break;
	}

	let box = get_but_box(lay.but);
	
	cy = box.ymax+0.8;

	if(cy > lay.op.ymax){
		cy = lay.op.ymax;
		lay.background = BUBBLE_SCROLL_COL;
		for(let i = 0; i < lay.but.length; i++){
			if(lay.but[i].back_col == BUBBLE_COL) lay.but[i].back_col = BUBBLE_SCROLL_COL;
		}
	}

	lay.dy = cy;
}


/// Adds a dropdown menu to the bubble
function bubble_adddropdown_text(cont,te,tab,w,source,pos)
{
	cont.lay.add_button({te:te, x:0, y:cont.y+0.3, dx:tab, dy:0.8, type:"InputBoxName"});
	bubble_adddropdown(cont,tab,w,source,pos);
}


/// Adds a link to the bubble
function bubble_addlink(cont,te,ac,op)
{
	let si = 0.8;
	let fo = get_font(si);
	let w = text_width(te,fo)+0.1;
	
	cont.lay.add_button({te:te, x:0.5, y:cont.y, dx:w, si:si, dy:1, font:fo, ac:ac, op:op, type:"BubbleLink"});
	cont.y += 1.4;
}


/// Adds a radio button to the bubble 
function bubble_addradio(cont,tab,value,te,source,disable)
{
	let ac; if(disable != true) ac = "RadioButton";
	
	let xx = nearest_pixel(tab+0.5), yy = nearest_pixel(cont.y);
	
	cont.lay.add_button({value:value, x:xx, y:yy, dx:1.2, dy:1.2, ac:ac, type:"RadioButton", source:source});
	
	let w = text_width(te,get_font(si_radio));
	cont.lay.add_button({te:te, x:xx+1.2, y:yy+0.1, dx:w+0.4, dy:1, type:"RadioButtonText", source:source, disable:disable, col:BLACK});

	cont.y += 1.4;
}


/// Adds a checkbox to the bubble 
function bubble_addcheckbox(cont,tab,te,source)
{
	let ac = "CheckboxButton";
	cont.lay.add_button({x:tab+0.5, y:cont.y, dx:1, dy:1, ac:ac, type:"CheckboxButton", source:source});
	cont.lay.add_button({te:te, x:tab+1.6, y:cont.y, dx:cont.dx-1.5-tab, dy:1, type:"CheckboxButtonText", source:source, col:BLACK});
	
	cont.y += 1.4;
}


/// Adds text to the bubble
function bubble_addtext(cont,te,tab)
{
	if(tab == undefined) tab = 0;
	cont.lay.add_button({te:te, x:tab, y:cont.y, dx:cont.dx-tab, dy:0.8, type:"InputBoxName"});
	cont.y += 1.4;
}


/// Adds simple text to the bubble
function bubble_simple_text(cont,te,si,col,tab)
{
	if(tab == undefined) tab = 0;
	cont.lay.add_button({te:te, x:tab, y:cont.y, dy:si, type:"Text", si:si, font:get_font(si),col:col});
	
	cont.y += si+0.3;
}


/// Adds a paragraph to the bubble
function bubble_addparagraph(cont,te,x,dx)
{
	cont.y = cont.lay.add_paragraph(te,dx,x,cont.y,BLACK,bub_si,bub_lh);
}


/// Adds a big paragraph to the bubble
function bubble_addbigparagraph(cont,te,x,dx)
{
	cont.y = cont.lay.add_paragraph(te,dx,x,cont.y,BLACK,bubbig_si,bubbig_lh);
}


/// Adds button to the bubble
function bubble_addbutton(cont,te,x,y,dx,dy,type,ac,op)
{	
	cont.lay.add_button({te:te, x:x, y:y, dx:dx, dy:dy, ac:ac, type:type, op:op});
}


/// Checks for error in the bubble (activated when done is pressed)
function bubble_check_error()
{
	let bub = inter.bubble;
	
	switch(bub.check){
	case "checkbox":
		{
			let cb = edit_source.spec.check_box.value; 
			let i = 0; while(i < cb.length && cb[i].check == false) i++;
			if(i == cb.length){
				bub.check_warning = "At least one must be selected";
				return true;
			}				
			delete bub.check_warning;
		}
		break;
			
	case "trans_checkbox":
		{
			let cb = edit_source.spec.filter.tra; 
			let obs_mod = edit_source.spec.filter.trans_obs_model.value;
			switch(obs_mod){
			case "on": 
				{
					let i = 0; while(i < cb.length && cb[i].prob_eqn.te == "0") i++;
					if(i == cb.length){
						bub.check_warning = "At least one probability must be non-zero";
						return true;
					}
				}
				break;
			
			case "off":
				{
					let i = 0; while(i < cb.length && cb[i].check == false) i++;
					if(i == cb.length){
						bub.check_warning = "At least one transition must be selected";
						return true;
					}
				}
				break;
			}
			
			let source_fl = false, sink_fl = false, tra_fl = false;
			for(let i = 0; i < cb.length; i++){
				let cbb = cb[i];
				
				let fl = false;
				switch(obs_mod){
				case "on":{ if(cbb.prob_eqn.te != "0") fl = true;}
				case "off":{ if(cbb.check == true) fl = true;}
				}
				
				if(fl == true){
					if(cbb.i == SOURCE) source_fl = true;
					else{
						if(cbb.f == SINK) sink_fl = true;
						else tra_fl = true;
					}
				}
			}

			if(source_fl && tra_fl){
				bub.check_warning = "Cannot include source and non-source transitions";
				return true;
			}
			
			if(sink_fl && tra_fl){
				bub.check_warning = "Cannot include sink and non-sink transitions";
				return true;
			}
			
			if(source_fl && sink_fl){
				bub.check_warning = "Cannot include source and sink transitions";
				return true;
			}
			
			delete bub.check_warning;
		}
		break;
		
	case "pop_checkbox":
		{
			let sp = model.get_sp();
			let clz = edit_source.spec.filter.cla;
		
			for(let cl = 0; cl < clz.length; cl++){
				if(clz[cl].radio.value == "Comp"){
					let cb = clz[cl].comp; 
					let i = 0; while(i < cb.length && cb[i].check == false) i++;
					if(i == cb.length){
						inter.bubble.check_warning = "At least one compartment must be selected";
						return true;
					}	
				}
			}
		}
		break;
	
	default: break;
	}

	if(check_error_textbox() == false && inter.bubble.warning == false) return false;
	
	inter.bubble.show_warning = true;
	return true;
}


/// Adds a clour table to bubble
function bubble_colour(cont,tab)
{
	let dby = 1.2, gap = 0, mar = 0.;
	if(tab == undefined) tab = 0;
	
	let dbx = (cont.dx-tab+gap*6-2*mar)/8;

	cont.lay.add_button({te:"Colour:", x:tab, y:cont.y, dx:cont.dx-tab, dy:0.8, type:"InputBoxName"});
	cont.y += 1.0;
	
	for(let j = 0; j < 4; j++){
		for(let i = 0; i < 8; i++){
			cont.lay.add_button({x:tab+mar+(dbx+gap)*i, y:cont.y+(dby+gap)*j, dx:dbx, dy:dby, ac:"ColourSelect", type:"ColourSelect", col:collist[j*8+i], sel_bu:cont.bu});
		}
	}
	
	cont.y += (dby+gap)*4;
	
	let col2 = cont.bu.col;
	if(cont.bu.type == "AddButton"){
		col2 = inter.bubble.radio.col;
		if(inter.bubble.radio.value != "SelectColour") col2 = BLACK;
	}
			
	if(col2 == WHITE){
		inter.bubble.warning = true;
		if(inter.bubble.show_warning == true){
			cont.y = cont.lay.add_paragraph("A colour must be selected!",cont.dx,tab,cont.y+0.1,RED,warn_si,warn_lh,undefined,"center");
			cont.y -= 0.1;
		}
	}
}


/// Shifts the bubble to ensure it is within view
function shift_bubble(dx,dy,lay,pa,pb,pc)
{
	lay.x += dx; pa.x += dx; pb.x += dx; pc.x += dx;
	lay.y += dy; pa.y += dy; pb.y += dy; pc.y += dy;
}


/// Sets up the positioning of layers for the bubble 
function setup_bubble_back(cont)
{
	let lay = cont.lay;
	let bu = cont.bu;

	let bu_lay = cont.bu_lay;

	let lay_but = get_lay(inter.bubble.lay_name);

	let bux = nearest_pixel(bu.x-lay_but.x_shift);
	let buy = nearest_pixel(bu.y-lay_but.y_shift);
	let budx = nearest_pixel(bu.dx), budy = nearest_pixel(bu.dy);
	
	if(bu.type == "Transition"){ 
		bux = nearest_pixel(bu.center.x); buy = nearest_pixel(bu.center.y);
		budx = 0; budy = 0;
	}
	
	if(bux+budx < 0 || buy+budy < 0 || bux > bu_lay.dx || buy > bu_lay.dy){
		inter.bubble.out_of_range = true;
	}
	
	let n = cont.end_button.length;
	
	if(n > 0){
		let box = get_but_box(lay.but);
		
		let dy = nearest_pixel(1.2), gap = nearest_pixel(0.4);
		
		let xx = box.xmax;
		let yy = cont.y;
		

		for(let i = n-1; i >= 0; i--){
			let endb = cont.end_button[i];
			let dx = nearest_pixel(3.5);
			if(endb.te == "Calculate") dx = nearest_pixel(4.5);
			lay.add_button({te:endb.te, x:xx-dx, y:yy, dx:dx, dy:dy, ac:endb.ac, type:"BubbleEndBut", op:endb.op});
			if(i == n-1){
				inter.bubble.final_button = lay.but[lay.but.length-1];
			}
			xx -= dx+gap;
		}
	}
	
	let back_lay = inter.layer[lay.index-1];
	
	let box = get_but_box(lay.but);

	lay.dx = box.xmax+0.1;
	lay.dy = box.ymax+0.1;
	lay.inner_dx = lay.dx;
	lay.inner_dy = lay.dy;
	
	let bx = bu_lay.x + bux, by = bu_lay.y + buy;
	let bw = budx, bh = budy;
	
	let gapx = nearest_pixel(1.5), gapy = nearest_pixel(1.6), gap = nearest_pixel(1);
	let marx = nearest_pixel(0.6), mary = nearest_pixel(0.6);
	
	let w_right = distance(bu_lay.x+bux+budx+gap+0.5*lay.dx+marx,bu_lay.y+buy+0.5*budy);

	let w_left = distance(bu_lay.x+bux-gap-0.5*lay.dx-marx,bu_lay.y+buy+0.5*budy);
	
	let w_top = distance(bu_lay.x+bux+0.5*budx,bu_lay.y+buy-gap-0.5*lay.dy-mary);
	let w_bottom = distance(bu_lay.x+bux+0.5*budx,bu_lay.y+buy+budy+gap+0.5*lay.dy+mary);
		
	let overlap_red = 30;
	
	let bubwid = lay.dx+2*marx;
	let bubhei = lay.dy+2*mary;
	
	// Penalises if bubble is outside page
	if(bu_lay.x+bux + bw + bubwid+gap > page_char_wid) w_right += overlap_red;
	if(bu_lay.x+bux - bubwid-gap < 0) w_left += overlap_red;
	if(bu_lay.y+buy + bh + bubhei+gap > page_char_hei) w_bottom += overlap_red;
	if(bu_lay.y+buy - bubhei-gap < 0) w_top += overlap_red;

	// Penalises non-right if in menu
	if(bu_lay.x+bux < menu_width){ w_top += overlap_red; w_bottom += overlap_red;}

	let w = w_right, orient = "right";
	if(w_left < w){ w = w_left; orient = "left";}
	if(w_top < w){ w = w_top; orient = "top";}
	if(w_bottom < w){ w = w_bottom; orient = "bottom";}
	
	let pa={}, pb={}, pc={};

	switch(bu.type){
	case "BayesFactor":  orient = "top"; break;
	case "AddSpecies": orient = "right"; break;
	case "ClassTab": orient = "bottom"; break;
	case "Settings": orient = "top"; break;
		
	case "AddButton":
		if(bu.ac == "AddClassification") orient = "bottom";
		else orient = "top"; 
		break;

	case "GreyView": case "AddFilter": 
		orient = "left"; 
		break;

	default: break;
	}
	
	switch(orient){
	case "right":
		{
			pa.x = bx+bw;
			pa.y = by+0.5*bh;

			pb.x = pa.x + gapx;
			pb.y = pa.y + gapy;
			
			pc.x = pb.x;
			pc.y = pb.y-gap;
		
			lay.x = pb.x+marx;
			lay.y = pb.y - 1*lay.dy+mary;
			if(lay.y < 1.5) lay.y = 1.5;
			
			if(pb.y-lay.y < 1.3) lay.y = pb.y-1.3;
		}
		break;
		
	case "left":
		{
			pa.x = bx;
			pa.y = by+0.5*bh;

			pb.x = pa.x - gapx;
			pb.y = pa.y + gapy;
			
			pc.x = pb.x;
			pc.y = pb.y-gap;
		
			lay.x = pb.x-lay.dx-marx;
			lay.y = pb.y - 1*lay.dy+mary;
			if(lay.y < 1.5) lay.y = 1.5;
			
			if(pb.y-lay.y < 1.3) lay.y = pb.y-1.3;
		}
		break;
		
	case "bottom":
		{
			pa.x = bx+0.5*bw;
			pa.y = by+bh;

			pb.x = pa.x + gapy;
			pb.y = pa.y + gapx;
			
			pc.x = pb.x-gap;
			pc.y = pb.y;
		
			lay.x = pb.x - 1*lay.dx+marx;
			lay.y = pb.y+mary;
			
			if(lay.x < back_lay.x+1.5) lay.x = back_lay.x+1.5;
			
			if(pb.x-lay.x < 1.) lay.x = pb.x-1.;
		}
		break;
		
	case "top":
		{
			pa.x = bx+0.5*bw;
			pa.y = by;

			pb.x = pa.x + gapy;
			pb.y = pa.y - gapx;
			
			pc.x = pb.x-gap;
			pc.y = pb.y;
		
			lay.x = pb.x - 1*lay.dx+marx;
			lay.y = pb.y-lay.dy-mary;
			
			if(lay.x < back_lay.x+1.5) lay.x = back_lay.x+1.5;
			
			if(pb.x-lay.x < 1.) lay.x = pb.x-1.;
		}
		break;

	default: error("Option not recognised 12"); break;
	}
	
	let mar_l = 1.5, mar_r = 0.5;
	let mar_u = 1.3, mar_d = 1.3;
	if(lay.x-mar_l < menu_width) shift_bubble(menu_width-(lay.x-mar_l),0,lay,pa,pb,pc);
	if(lay.x+lay.dx+mar_r > page_char_wid) shift_bubble(page_char_wid-(lay.x+lay.dx+mar_r),0,lay,pa,pb,pc);
	if(lay.y-mar_u < 0) shift_bubble(0,-(lay.y-mar_u),lay,pa,pb,pc);
	if(lay.y+lay.dy+mar_d > page_char_hei) shift_bubble(0,page_char_hei-(lay.y+lay.dy+mar_d),lay,pa,pb,pc);

	let xref = back_lay.x, yref = back_lay.y;
	
	let back_bu = back_lay.but[0];

	back_bu.tag = {x0:pa.x-xref, y0:pa.y-yref, x1:pb.x-xref, y1:pb.y-yref, x2:pc.x-xref, y2:pc.y-yref};
	
	back_bu.x = lay.x - marx - xref;
	back_bu.y = lay.y - mary - yref;
	
	back_bu.dx = lay.dx+2*marx;
	back_bu.dy = lay.dy+2*mary;

	let nothing_bu = back_lay.but[1];	
	nothing_bu.x = back_bu.x; nothing_bu.y = back_bu.y; 
	nothing_bu.dx = back_bu.dx; nothing_bu.dy = back_bu.dy; 

	let close_bu = back_lay.but[2];
	
	close_bu.x = back_bu.x+back_bu.dx-1.5;
	close_bu.y = back_bu.y+0.5;
	
	let sc_y = cont.scrollable_y;
	if(sc_y != undefined){
		shift_layer("BubbleScrollable",lay.x,lay.y+sc_y);
	}

	back_lay.plot_buttons();
	for(let i = 0; i < cont.input_list.length; i++){
		let ip = cont.input_list[i];
		if(ip.op == undefined) error("Problem with input list");

		lay.add_input(ip.x,ip.y,ip.dx,ip.op);
	}
	
	for(let i = 0; i < cont.drop_list.length; i++){
		let ip = cont.drop_list[i];
		lay.add_dropdown(ip.x,ip.y,ip.dx,10,ip.source,ip.pos);
	}
}


/// Returns the distance to the center of the screen (used for orient)
function distance(x,y)
{
	return Math.sqrt((x-page_char_wid/2)*(x-page_char_wid/2) + (y-page_char_hei/2)*(y-page_char_hei/2));
}


/// Selects a given transition in the model
function select_bubble_transition(p,cl,i)
{
	change_page({pa:"Model",su:"Compartments",susu:p,sususu:cl});
	
	let lay_name = "Transition";
	let lay = get_lay(lay_name);
	let bu = lay.but;
	
	let j = 0;
	while(j < bu.length && !(bu[j].type == "Transition" && bu[j].i == i)) j++;
	if(j == bu.length){ error("Problem selecting1"); return;}
	
	let mar = 2;
	let x = bu[j].center.x, y = bu[j].center.y;
	
	// Recenters camera if outside range
	if(x < mar || x > lay.dx-mar || y < mar || y > lay.dy-mar){  
		let cam = model.species[p].cla[cl].camera;
		let tra = model.species[p].cla[cl].tra[i];
		cam.x = tra.center.x; cam.y = tra.center.y;
	}
	
	select_bubble(lay_name,j);
	generate_screen();
}


/// Selects the specification bubble for a data source
function select_bubble_data_spec(p,i)
{
	change_page({pa:"Inference",su:"Data",susu:p});
	
	let lay_name = "TableContent";
	let lay = get_lay(lay_name);
	
	let bu = lay.but;

	let j = 0;
	while(j < bu.length && !(bu[j].te == "View" && bu[j].info.i == i)) j++;
	if(j == bu.length){ error("Problem selecting2"); return;}
	
	activate_button(lay,j);   
}	


/// Selects the data table element
function select_bubble_data_element(p,i,r,c)
{
	change_page({pa:"Inference",su:"Data",susu:p});
	
	let lay_name = "TableContent";
	let lay = get_lay(lay_name);
	
	let bu = lay.but;

	let j = 0;
	while(j < bu.length && !(bu[j].te == "Edit" && bu[j].info.i == i)) j++;
	if(j == bu.length){ error("Problem selecting3"); return;}
	
	activate_button(lay,j); 
	
	select_table_elelent(r,c);
}	


/// Selects a reparameterisation element with parameter name and index
function select_reparam_element(par_name,index)
{
	goto_param_page();
	
	let lay_name = "ModelParamContent";
	let lay = get_lay(lay_name);
	
	let bu = lay.but;

	let i = find(model.param,"name",par_name);
	if(i == undefined){ error("Problem selecting6"); return;}

	let j = 0;
	while(j < bu.length && !(bu[j].ac	== "EditReparamValue" && bu[j].i == i)) j++;
	if(j == bu.length){ error("Problem selecting5"); return;}
	
	activate_button(lay,j); 
	
	if(index != undefined) select_param_element(index);
}


/// Selects a distribution element with parameter name and index
function select_dist_element(par_name,index)
{
	goto_param_page();
	
	let lay_name = "ModelParamContent";
	let lay = get_lay(lay_name);
	
	let bu = lay.but;

	let i = find(model.param,"name",par_name);
	if(i == undefined){ error("Problem selecting6"); return;}

	let j = 0;
	while(j < bu.length && !((bu[j].ac	== "EditPriorElement" || bu[j].ac	== "EditDistSplitValue") && bu[j].i == i)) j++;
	if(j == bu.length){ error("Problem selecting5"); return;}
	
	activate_button(lay,j);
	if(index != undefined) select_param_element(index);
}


/// Goes to the parameter page but ensures that warning are turned off
function goto_param_page()
{
	change_page({pa:"Model",su:"Parameters"});
}


/// Selects the bubble for a derived quantitiy
function select_bubble_derived(i)
{
	goto_param_page();

	let lay_name = "ModelParamContent";
	let lay = get_lay(lay_name);
	
	let bu = lay.but;

	let j = 0;
	while(j < bu.length && !(bu[j].ac	== "EditDerive" && bu[j].val == i)) j++;
	if(j == bu.length){ error("Problem selecting4"); return;}
	
	activate_button(lay,j); 
}	


/// Selects the bubble for a derived quantitiy
function select_bubble_param(par)
{
	goto_param_page();

	let lay_name = "ModelParamContent";
	let lay = get_lay(lay_name);
	
	let bu = lay.but;

	let i = find(model.param,"name",par.name);
	
	let ac, ac2;
	switch(par.variety){
	case "const": ac = "EditSimValue"; break;
	case "reparam": ac = "EditReparamValue"; break;
	case "dist": ac = "EditPriorElement"; ac2 = "EditDistSplitValue"; break;
	}		
	
	let j = 0;
	while(j < bu.length && !((bu[j].ac	== ac || bu[j].ac	== ac2) && bu[j].i == i)) j++;
	if(j == bu.length){ error("Problem selecting7"); return;}
	
	activate_button(lay,j); 
}	


/// Selects a given transition in the model
function select_bubble_compartment(p,cl,i)
{
	change_page({pa:"Model",su:"Compartments",susu:p,sususu:cl});
	
	let lay_name = "Compartment";
	let lay = get_lay(lay_name);
	
	let bu = lay.but;
	
	let j = 0;
	while(j < bu.length && !((bu[j].type == "Compartment" || bu[j].type == "CompLatLng") && bu[j].i == i)) j++;
	if(j == bu.length){ error("Problem selecting"); return;}
	
	let x = bu[j].x, y = bu[j].y;
	
	let mar = 2;
	
	// Recenters camera if outside range
	if(x < mar || x > lay.dx-mar || y < mar || y > lay.dy-mar){  
		let cam = model.species[p].cla[cl].camera;
		let co = model.species[p].cla[cl].comp[i];
		
		cam.x = co.x; cam.y = co.y;
	}
	
	select_bubble(lay_name,j);
	generate_screen();
}


/// Selects a given transition in the model
function select_bubble_box(p,cl,i)
{
	change_page({pa:"Model",su:"Compartments",susu:p,sususu:cl});
	
	let lay_name = "Annotation";
	let lay = get_lay(lay_name);
	
	let bu = lay.but;

	let j = 0;
	while(j < bu.length && !(bu[j].type == "Box" && bu[j].i == i)) j++;
	if(j == bu.length){ error("Problem selecting"); return;}
	
	let x = bu[j].x, y = bu[j].y;
	
	let mar = 2;
	
	// Recenters camera if outside range
	if(x < mar || x > lay.dx-mar || y < mar || y > lay.dy-mar){  
		let cam = model.species[p].cla[cl].camera;
		let co = model.species[p].cla[cl].comp[i];
		
		cam.x = co.x; cam.y = co.y;
	}
	
	select_bubble(lay_name,j);
	generate_screen();
}


/// Gets the posible types of distributions
function get_dist_pos()
{
	let dist = dist_pos;
	let p = model.get_p();
	if(model.species[p].type == "Population") dist = exp_dist_pos;
	return dist;
}
