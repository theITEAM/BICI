"use strict";

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
	
	//pr(bu.type+" Bubble type");
	switch(bu.type){		
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
			default: inter.graph.settings_speed_bubble(cont); break;
			}
		}
		break;
		
	case "AddFilter":
		filter_bubble(bu,cont);
		break;
		
	case "CombineIE":
		combine_ie_bubble(cont);
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
			cont.dx = 11;
			
			bubble_addtitle(cont,"Compartment",{te:compartment_text2});
			bubble_input(cont,"Name:",{type:"compartment", p:bu.p, cl:bu.cl, i:bu.i});
			
			bubble_colour(cont);
			cont.y += 0.3;
			
			let co = model.species[bu.p].cla[bu.cl].comp[bu.i];
		
			if(co.choose_branch == true){
				bubble_addcheckbox(cont,-0.1,"Add branching probability",bub.checkbox);
			}
			
			if(co.type != "boundary"){
				bubble_addcheckbox(cont,-0.1,"Position fixed",co.fixed);
			}
			
			add_end_button(cont,"Delete","DeleteComp",{p:bu.p, cl:bu.cl, i:bu.i});		
			add_end_button(cont,"OK","ChangeComp",{p:bu.p, cl:bu.cl, i:bu.i});		
		}
		break;
	
	case "Transition": case "TransitionPoint":
		let tr = bu.tr;

		cont.dx = 12;
		
		if(tr.i == "Source"){
			bubble_addtitle(cont,"Source",{te:source_text2});
		}
		else{
			if(tr.f == "Sink"){
				bubble_addtitle(cont,"Sink",{te:sink_text2});
			}
			else{
				bubble_addtitle(cont,"Transition",{te:trans_text2});
			}
		}
		
		cont.y += 0.4;
		
		let dist = get_dist_pos();
		if(tr.i == "Source") dist = exp_dist_pos;
		
		bubble_addtext(cont,"Distribution:");
		cont.y -= 1.7;
		
		let pos=[]; for(let i = 0; i < dist.length; i++) pos.push({te:dist[i]});
		
		bubble_adddropdown(cont,4.6,6.5,inter.bubble.trans_type,pos);

		cont.y += 0.5;
		
		if(tr.branch == true){
			bubble_input(cont,"Branching probability:",{type:"trans_bp", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
		}
		
		let mean_fl = false, rate_fl = false, shape_fl = false;
		let cv_fl = false, scale_fl = false, shape_erlang_fl = false; 
		
		switch(inter.bubble.trans_type.te){
		case "exp(mean)":
			mean_fl = true;
			bubble_input(cont,"Mean:",{type:"trans_mean", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
			if(tr.i == "Source"){ cv_fl = true; shape_fl = true; scale_fl = true; shape_erlang_fl = true;}
			break;
			
		case "exp(rate)":
			rate_fl = true;
			bubble_input(cont,"Rate:",{type:"trans_rate", eqn:true, p:bu.p, cl:bu.cl, i:bu.i});
			if(tr.i == "Source"){ cv_fl = true; shape_fl = true; scale_fl = true; shape_erlang_fl = true;}
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
					
						//cont.y += 0.3;
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
	
			case "AddAnnotation":
				cont.dx = 10;
				switch(bub.mode){
				case "add_label":
					bubble_addtitle(cont,"Add label",{te:label_text});
					bubble_input(cont,"Text:",{type:"label"});
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
					bubble_addradio(cont,0,"Remove","Edit annotations",bub.radio);
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
			case "SetConstant": set_constant_bubble(cont); break;
			case "SetReparam": set_reparam_bubble(cont); break;
			case "SetDistribution": set_distribution_bubble(cont); break;
			case "SetDerived": set_derived_bubble(cont); break;
			case "SetTrapsData": settraps_data_bubble(cont,"add"); break;
			default: error("Cannot find source"); break;
			}
		}
		break;
	
	case "GreyView":
		switch(inter.edit_source.type){
		case "Init. Pop.": initpop_data_bubble(cont,"view"); break;
		case "Move Ind.": move_data_bubble(cont,"view"); break;
		case "Compartment": comp_data_bubble(cont,"view"); break;
		case "Transition": trans_data_bubble(cont,"view","TransVar"); break;
		case "Source": trans_data_bubble(cont,"view","SourceVar"); break;
		case "Sink": trans_data_bubble(cont,"view","SinkVar"); break;
		case "Diag. Test": diagtest_data_bubble(cont,"view"); break;
		case "Population": population_data_bubble(cont,"view"); break;
		case "Pop. Trans.": poptrans_data_bubble(cont,"view"); break;
		case "Set Traps": settraps_data_bubble(cont,"view"); break;
		default: error("Cannot find source: "+inter.edit_source.type); break;
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
			
			bubble_addtext(cont,"Type:");
			
			cont.y -= 1.2;
			
			bubble_addradio(cont,2.3,"Individual","Individual-based",bub.radio);
			bubble_addradio(cont,2.3,"Population","Population-based",bub.radio);
			cont.y += 0.2;
			
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
		
		bubble_addtext(cont,"Type:");
			
		cont.y -= 1.2;
		
		bubble_addradio(cont,2.3,"Individual","Individual-based",bub.radio);
		bubble_addradio(cont,2.3,"Population","Population-based",bub.radio);
		cont.y += 0.2;
		
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
		cont.dx = 10;
		bubble_addtitle(cont,"Edit equation",{te:editeqnparam_text});	
		bubble_input(cont,"Equation:",{type:"reparam_eqn", eqn:true});
		add_end_button(cont,"Done","DoneEquation");	
		break;
		
	case "ReparamTableElement":
		cont.dx = 10;
		bubble_addtitle(cont,"Edit equation",{te:editeqnparam_text});	
		bubble_input(cont,"Equation:",{type:"element_eqn", eqn:true, pindex:bu.pindex});
		add_end_button(cont,"Done","Done");	
		break;
	
	case "ParamSimElement":
		cont.dx = 10;
		bubble_addtitle(cont,"Edit value",{te:editsimparam_text});
		bubble_input(cont,"Value:",{type:"param_val", val:bu.i});
		add_end_button(cont,"Done","Done");	
		break;
		
	case "PriorElement": case "DistElement": 
	case "PriorSplitElement": case "DistSplitElement":
		{
			let ac;
			let ae = false; // determines if allow eqn
		
			cont.dx = 8;
			switch(bu.type){
			case "PriorElement": case "PriorSplitElement":
				{
					ac = "DonePrior";
					let te = editprior_text;
					if(model.param[bu.i].pri_pos == bp_prior_pos) te = editpriorbp_text;
					if(model.param[bu.i].pri_pos == zeroone_prior_pos) te = editpriorzeroone_text;
			
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

			let pp = model.param[bu.i].pri_pos;
			
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
		bubble_input(cont,"Knot times:",{type:"knot_times", val:bu.i, eqn:false});
			
		add_end_button(cont,"Done","DoneKnots");	
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
				let te = inter.edit_source.table.ele[bu.r][bu.c];
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
					let te = "Replaced "+inter.bubble.num+" occurence";
					if(inter.bubble.num > 1) te += ".";
					bubble_addparagraph(cont,te,0,10);
				}
				break;
			
			case "Order":
				cont.dx = 10;
				bubble_addtitle(cont,"Order column",{te:order_text});
				bubble_addradio(cont,0,"A-Z","Alphabetically (A-Z)",inter.bubble.order_radio);			bubble_addradio(cont,0,"Z-A","Alphabetically (Z-A)",inter.bubble.order_radio);
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
		cont.dx = 10;
		bubble_addtitle(cont,"Label",{te:label_text});
		bubble_input(cont,"Text:",{type:"label_anno",p:bu.p, cl:bu.cl, i:bu.i});
		add_end_button(cont,"Delete","DeleteAnno",{p:bu.p, cl:bu.cl, i:bu.i});		
		add_end_button(cont,"Done","LabelOK");		
		break;
		
	case "Box":
		cont.dx = 10;
		bubble_addtitle(cont,"Box",{te:box_text});
		bubble_input(cont,"Text:",{type:"label_anno", p:bu.p, cl:bu.cl, i:bu.i});
			
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
	
	default: error(bu.type+": Bubble Not recognised"); break;
	}

	setup_bubble_back(cont);
	
	if(bub.find_focus == true){
		set_focus_first();
	
		bub.find_focus = false;
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
	
	inter.bubble = { lay_name:lay_name, bu:bu, i:i, show_warning:false, warning:false, op:op, find_focus:true};
	
	switch(bu.type){
	case "Transition": case "TransitionPoint":
		inter.bubble.value = bu.tr.value;
		inter.bubble.trans_type = {te:bu.tr.type};
		break;
	}
	reset_text_box();
}


/// Changes the mode for a bubble
function change_bubble_mode(mode)
{
	inter.bubble.mode = mode;
	reset_text_box();
	inter.bubble.find_focus = true;
	//	let over_new = get_button_over(x,y);
}


/// Closes a speech bubble
function close_bubble()
{
	if(inter.help.title != undefined && inter.equation.te != undefined){
		error("Cannot close bubble");
	}

	turn_off_cursor();
	inter.bubble = {};
	//generate_screen();
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


///Adds a termination point to a bubble
function add_bubble_end(cont)
{
	cont.lay.add_button({x:0,y:cont.y, type:"Nothing"});
}


/// Adds a mini title 
function bubble_add_minititle(cont,te)
{
	cont.lay.add_button({te:te, x:0, y:cont.y, dx:cont.dx, dy:0.8, type:"InputBoxName"});
	cont.y += 1.4;	
}


/// Adds an input text box to the bubble		
function bubble_input(cont,te,op)
{	
	let x = 0; if(op.x != undefined) x = op.x;
	let w = cont.dx; if(op.w != undefined) w = op.w;
	
	add_ref(op);
	
	if(op.hidden == true){
		cont.input_list.push({x:x, y:-1000, dx:cont.dx, dy:0.8, op:op});
	}
	else{
		if(te != ""){
			cont.lay.add_button({te:te, x:x, y:cont.y, dx:w, dy:0.8, type:"InputBoxName"});
			cont.y += 1.2;
		}
		
		cont.input_list.push({x:x+0.3, y:cont.y, dx:w-0.6, op:op});
		
		if(op.no_down == true) cont.y -= 1.4;
		else cont.y += 1.8;
		
		//let k = cont.input_list.length-1;
		let sto = inter.textbox_store;

		let warn;
		let k = find(sto,"ref",op.ref);
		if(k != undefined) warn = sto[k].warning;
		//let k = 0; while(k < sto.length && sto[k].ref != op.ref) k++;
		//if(k < sto.length) warn = sto[k].warning;
	
		if(warn == undefined && inter.bubble.error_warning != undefined){
			warn = inter.bubble.error_warning;
		}
		
		if(warn != "" && warn != undefined){
			cont.y = cont.lay.add_paragraph(warn,w,x,cont.y-0.2,RED,warn_si,warn_lh,undefined,"center");
			cont.y += 0.2;
		}				
	}
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

	add_layer("BubbleScrollable",0,0,cont.dx+0.4,0,op);
	let lay = get_lay("BubbleScrollable");

	cont.y += lay.dy-0.2;
}


function add_bubble_scrollable_buts(lay)
{
	let cy;
	
	lay.background = BUBBLE_COL;

	switch(lay.op.type){
	case "comp list": cy = diagtest_scrollable(lay); break;
	case "pop list": cy = population_scrollable(lay); break;
	case "poptrans list": cy = poptrans_scrollable(lay); break;
	case "annotation": cy = annotation_scrollable(lay); break;
	case "param sel": cy = param_sel_scrollable(lay); break;
	case "param details": cy = param_details_scrollable(lay); break;
	case "combine IE": cy = combineIE_scrollable(lay); break;
	default: error("Option not recognised 9"); break;
	}

	let box = get_but_box(lay.but);
	cy = box.ymax;
	
	if(cy > lay.op.ymax){
		cy = lay.op.ymax;
		lay.background = WHITE;
		for(let i = 0; i < lay.but.length; i++){
			if(lay.but[i].back_col == BUBBLE_COL) lay.but[i].back_col = WHITE;
		}
	}

	lay.dy = cy;

	//lay.but[0].dy = box.ymax;
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
	let w = text_width(te,fo);
	cont.lay.add_button({te:te, x:0.5, y:cont.y, dx:w, si:si, dy:1, font:fo, ac:ac, op:op, type:"BubbleLink"});
	cont.y += 1.4;
}


/// Adds a radio button to the bubble 
function bubble_addradio(cont,tab,value,te,source,disable)
{
	let ac; if(disable != true) ac = "RadioButton";
	
	cont.lay.add_button({value:value, x:tab+0.5, y:cont.y, dx:1, dy:1, ac:ac, type:"RadioButton", source:source});
	
	cont.lay.add_button({te:te, x:tab+1.6, y:cont.y+0.05, dx:cont.dx-(tab+1.6), dy:1, type:"RadioButtonText", source:source, disable:disable, col:BLACK});

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
			let cb = inter.edit_source.spec.check_box.value; 
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
			let cb = inter.edit_source.spec.filter.tra; 
			let i = 0; while(i < cb.length && cb[i].check == false) i++;
			if(i == cb.length){
				bub.check_warning = "At least one transition must be selected";
				return true;
			}				
			delete bub.check_warning;
		}
		break;
		
	case "pop_checkbox":
		{
			let sp = model.get_sp();
			let clz = inter.edit_source.spec.filter.cla;
		
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
	let dbx = (cont.dx-gap*6-2*mar)/7;

	if(tab == undefined) tab = 0;
	
	cont.lay.add_button({te:"Colour:", x:tab, y:cont.y, dx:cont.dx, dy:0.8, type:"InputBoxName"});
	cont.y += 1.0;
	
	for(let j = 0; j < 4; j++){
		for(let i = 0; i < 7; i++){
			cont.lay.add_button({x:tab+mar+(dbx+gap)*i, y:cont.y+(dby+gap)*j, dx:dbx, dy:dby, ac:"ColourSelect", type:"ColourSelect", col:collist[j*7+i], sel_bu:cont.bu});
		}
	}
	
	cont.y += (dby+gap)*4;
	
	let col2 = cont.bu.col;
	if(cont.bu.type == "Menu"){
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


/// Sets up the positioning of layers for the bubble 
function setup_bubble_back(cont)
{
	let lay = cont.lay;
	let bu = cont.bu;

	let bu_lay = cont.bu_lay;

	let lay_but = get_lay(inter.bubble.lay_name);

	let bux = bu.x-lay_but.x_shift, buy = bu.y-lay_but.y_shift, budx = bu.dx, budy = bu.dy;
	
	if(bu.type == "Transition"){ bux = bu.center.x; buy = bu.center.y; budx = 0; budy = 0;}
	
	//pr(bux+budx
	if(bux+budx < 0 || buy+budy < 0 || bux > bu_lay.dx || buy > bu_lay.dy){
		inter.bubble.out_of_range = true;
	}
	
	let n = cont.end_button.length;
	
	if(n > 0){
		let box = get_but_box(lay.but);
		
		let dx = 3.5, dy = 1.2, gap = 0.4;
		
		let xx = box.xmax - n*dx - gap*(n-1);
		let yy = cont.y;
		
		for(let i = 0; i < n; i++){
			let endb = cont.end_button[i];
			lay.add_button({te:endb.te, x:xx, y:yy, dx:dx, dy:dy, ac:endb.ac, type:"BubbleEndBut", op:endb.op});
			xx += dx+gap;
		}
		
		inter.bubble.final_button = lay.but[lay.but.length-1];
	}
	
	let back_lay = inter.layer[lay.index-1];
	
	let box = get_but_box(lay.but);

	lay.dx = box.xmax+0.1;
	lay.dy = box.ymax+0.1;
	lay.inner_dx = lay.dx;
	lay.inner_dy = lay.dy;
	
	let bx = bu_lay.x + bux, by = bu_lay.y + buy;
	let bw = budx, bh = budy;
	
	let gapx = 1.5, gapy = 1.6, gap = 1;
	let marx = 0.6, mary = 0.6;
	
	let w_right = distance(bux+budx+gap+0.5*lay.dx+marx,buy+0.5*budy,back_lay);

	let w_left = distance(bux-gap-0.5*lay.dx-marx,buy+0.5*budy,back_lay);
	
	let w_top = distance(bux+0.5*budx,buy-gap-0.5*lay.dy-mary,back_lay);
	let w_bottom = distance(bux+0.5*budx,buy+budy+gap+0.5*lay.dy+mary,back_lay);
		
	let overlap_red = 30;
	if(bux + bw > bu_lay.dx){	
		w_right += overlap_red;
		w_top += overlap_red;
		w_bottom += overlap_red;
	}
	
	/*
	if(bux < 0){	
		w_left -= overlap_red;
		w_top -= overlap_red;
		w_bottom -= overlap_red;
	}
	
	if(bux+bw+lay.dx > back_lay.dx) w_right = 0;
	if(bux-lay.dx < 0) w_left = 0;
	
	if(buy+bh+lay.dy > back_lay.dy) w_bottom = 0;
	if(buy-lay.dy < 0) w_top = 0;
	let ma = 4;
	
	let w_hor = 0; if(w_top < ma || w_bottom < ma) w_hor = -20;
	let w_vert = 0; //if(w_right < ma || w_left < ma) w_vert = -2000;

	w_right += w_hor; w_left += w_hor;
	w_top += w_vert; w_bottom += w_vert;
	*/
	
	let w = w_right, orient = "right";
	if(w_left < w){ w = w_left; orient = "left";}
	if(w_top < w){ w = w_top; orient = "top";}
	if(w_bottom < w){ w = w_bottom; orient = "bottom";}
	
	let pa={}, pb={}, pc={};
	
	//pr(bu.type+"ty");
	switch(bu.type){
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
	
	let xref = back_lay.x, yref = back_lay.y;
	
	let back_bu = back_lay.but[0];

	back_bu.tag = {x0:nearest_pixel(pa.x-xref), y0:nearest_pixel(pa.y-yref), x1:nearest_pixel(pb.x-xref), y1:nearest_pixel(pb.y-yref), x2:nearest_pixel(pc.x-xref), y2:nearest_pixel(pc.y-yref)};
	
	back_bu.x = nearest_pixel(lay.x - marx - xref);
	back_bu.y = nearest_pixel(lay.y - mary - yref);
	
	back_bu.dx = nearest_pixel(lay.dx+2*marx);
	back_bu.dy = nearest_pixel(lay.dy+2*mary);

	let nothing_bu = back_lay.but[1];	
	nothing_bu.x = back_bu.x; nothing_bu.y = back_bu.y; nothing_bu.dx = back_bu.dx; nothing_bu.dy = back_bu.dy; 

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

/*
/// Selects a given transition equation
function select_transition_equation(p,cl,i,type)
{
	select_bubble_transition(p,cl,i);
	
	for(let l = 0; l < inter.layer.length; l++){
		let lay = inter.layer[l];
		if(lay.name == "Input"){
			if(lay.op.source.type == type){
				let bs = lay.get_text_box_store();
				start_equation(bs.te,bs.eqn,lay.op.source,0);
				//equation_done();
				return;
			}
		}
	}
	
	error("Could not select transition equation");
}
*/


/// Returns the distance to the center of the screen (used for orient)
function distance(x,y,back)
{
	return Math.sqrt((x-back.dx/2)*(x-back.dx/2) + (y-back.dy/2)*(y-back.dy/2));
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
	
	if(x < mar || x > lay.dx-mar || y < mar || y > lay.dy-mar){  // Recenters camera if outside range
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
	
	select_param_element(index);
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
	if(model.warn.length > 0){
		model.warn.length = 0; 
		generate_screen();
	}
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
	if(x < mar || x > lay.dx-mar || y < mar || y > lay.dy-mar){  // Recenters camera if outside range
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
