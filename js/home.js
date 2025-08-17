"use strict";
// Functions for generating the home page

/// Adss buttons for home page
function add_home_page_buts(lay) 
{
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("About",cx,cy);
	
	let dx = 42;
	
	lay.add_button({x:cx+dx, y:cy+0.5, dx:4.5, dy:4.5, ac:"NewModel", type:"NewModel"});
	
	cy = lay.add_paragraph("BICI stands for 'Bayesian individual-based compartmental inference'. Either start a new analysis (right), load up a previous analysis from the menu (top right) or try one of the examples below.\nClick on the [?] icons for help.",dx,cx,cy,BLACK,para_si,para_lh);

	cy += 1;
	
	cy = lay.add_title("Examples",cx,cy,{te:example_text});

	let x1 = lay.x+cx;
	let x2 = lay.x+lay.dx-cx;
	let y1 = lay.y+cy+0.5;
	let y2 = lay.y+lay.dy-1;

	add_layer("Examples",x1,y1,x2-x1,y2-y1,{});
}


/// Adds buttons for example section
function add_examples_buts(lay) 
{
	let mod_mess = "\nEach example can be simulated from to see the dynamics it exhibits (using different visualisations). These simulated results can be used to generate data onto which inference can then be performed.";
	
	let ex_mod=[];
	
	//ex_mod.push({mod:"M", te:"Temp",file:"temp",mod_type:"IBM", pic:"SI"});
	
	if(true){
		if(false){
			ex_mod.push({te:"M) Science Review",help:"Illustrates some simple examples for the Science Review."});
			ex_mod.push({mod:"SR1", te:"SIR model (IBM)",file:"EX_SR1",mod_type:"IBM", pic:"SIR"});
			
			ex_mod.push({mod:"SR1", te:"SIR model (IBM)",file:"EX_B1",mod_type:"IBM", pic:"SIR"});
			
			ex_mod.push({mod:"SR2", te:"COVID-19 model",file:"EX_SR2",mod_type:"POP", pic:"covid"});
		
			ex_mod.push({mod:"SR3", te:"Metapopulation model using geographical regions",file:"EX_SR3",mod_type:"POP", pic:"scotland"});
		}
	
	
		//ex_mod.push({mod:"M5.2", te:"Multiple contact groups, investigating infectivity",file:"EX_E1",mod_type:"IBM", pic:"group"});
		/*
		ex_mod.push({mod:"M2.6", te:"SIR model (PBM) with demographic stratification",file:"temp2",mod_type:"POP", pic:"SIR_sex2"});
		
		ex_mod.push({mod:"M", te:"Temp",file:"temp",mod_type:"IBM", pic:"SI"});
	
	
		ex_mod.push({mod:"M", te:"Geno model",file:"geno",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"M", te:"Geno model with cattle/badgers",file:"geno_bad_cat",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"M", te:"Geno farm",file:"geno_farm",mod_type:"IBM", pic:"SI"});
		*/

		ex_mod.push({te:"M1) Simple epidemiological models",help:"This section gives some examples of simple epidemiological models, illustrating both individual and population-based approaches."+mod_mess});
		
		ex_mod.push({mod:"M1.1", te:"SI population-based model (PBM)",file:"EX_M1-1",mod_type:"POP", pic:"SI"});
		ex_mod.push({mod:"M1.2", te:"SI individual-based model (IBM)",file:"EX_M1-2",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"M1.3", te:"SIR model (IBM)",file:"EX_M1-3",mod_type:"IBM", pic:"SIR"});
		ex_mod.push({mod:"M1.4", te:"SIR model (PBM) with Erlang distribution",file:"EX_M1-4",mod_type:"POP", pic:"SIRErlang"});
		ex_mod.push({mod:"M1.5", te:"SIR model (IBM) with gamma distribution",file:"EX_M1-5",mod_type:"IBM", pic:"SIRgamma"});
		ex_mod.push({mod:"M1.6", te:"SEIR model (IBM) with exposed period",file:"EX_M1-6",mod_type:"IBM", pic:"SEIR"});
	
		ex_mod.push({te:"M2) Additional epidemiological models",help:"These examples look at some additional features, such as branching, reverse transitions and stratification."+mod_mess});
		
		ex_mod.push({mod:"M2.1", te:"SIRD model (PBM) with branching using transition rates",file:"EX_M2-1",mod_type:"POP", pic:"SIRD"});
		ex_mod.push({mod:"M2.2", te:"SIRD model (IBM) with branching probability",file:"EX_M2-2",mod_type:"IBM", pic:"SIRD_IBM"});
		ex_mod.push({mod:"M2.3", te:"SIRD model (IBM) with branching factors",file:"EX_M2-3",mod_type:"IBM", pic:"SIRD_IBM2"});

		ex_mod.push({mod:"M2.4", te:"SIS model (PBM)",file:"EX_M2-4",mod_type:"POP", pic:"SIS"});
		ex_mod.push({mod:"M2.5", te:"SIRS model (IBM) with waning immunity",file:"EX_M2-5",mod_type:"IBM", pic:"SIRwaning"});
		ex_mod.push({mod:"M2.6", te:"SIR model (PBM) with demographic stratification",file:"EX_M2-6",mod_type:"POP", pic:"SIR_sex2"});
		ex_mod.push({mod:"M2.7", te:"SIR model (PBM) with differential infectivity and demographic stratification",file:"EX_M2-7",mod_type:"POP", pic:"SIR_sex"});
		ex_mod.push({mod:"M2.8", te:"SIR model (IBM) with demographic stratification",file:"EX_M2-8",mod_type:"POP", pic:"SIR_sex2"});
	
		ex_mod.push({te:"M3) Spatial epidemiological models",help:"These models incorporate spatial stratification in different ways."+mod_mess});
		
		ex_mod.push({mod:"M3.1", te:"Metapopulation model using geographical regions",file:"EX_M3-1",mod_type:"POP", pic:"scotland"});
		ex_mod.push({mod:"M3.2", te:"Metapopulation model using geographical points",file:"EX_M3-2",mod_type:"POP", pic:"world"});
		ex_mod.push({mod:"M3.3", te:"Metapopulation model using a distance kernel",file:"EX_M3-3",mod_type:"POP", pic:"uk"});
		ex_mod.push({mod:"M3.4", te:"Farm-based model using a distance kernel",file:"EX_M3-4",mod_type:"IBM", pic:"farm"});
	
		ex_mod.push({te:"M4) Ecological models",help:"Compartmental models are often applied to epidemiological applications, but they can equally be used in other settings. Here we look at some examples from ecology."+mod_mess});
		ex_mod.push({mod:"M4.1", te:"Logistic growth population model",file:"EX_M4-1",mod_type:"IBM", pic:"carry"});
		ex_mod.push({mod:"M4.2", te:"A predator-prey model",file:"EX_M4-2",mod_type:"POP", pic:"predprey"});
		ex_mod.push({mod:"M4.3", te:"A spatial diffusion model",file:"EX_M4-3",mod_type:"IBM", pic:"grid"});
		ex_mod.push({mod:"M4.4", te:"A species presence/absence distribution model",file:"EX_M4-4",mod_type:"IBM", pic:"UK_grid"});
		ex_mod.push({mod:"M4.5", te:"Logistic growth model with age dependent mortality",file:"EX_M4-5",mod_type:"IBM", pic:"carryage"});
		
		ex_mod.push({te:"M5) Disease transmission experiments",help:"Disease transmission experiments can be used to identify factors affecting disease transmission. These can either be fixed categorical effects, such as vaccination status, covariates, such as weight, or quantitative genetic traits."+mod_mess});
		
		ex_mod.push({mod:"M5.1", te:"Single contact group, investigating susceptibility",file:"EX_M5-1",mod_type:"IBM", pic:"dte"});

		ex_mod.push({mod:"M5.2", te:"Multiple contact groups, investigating infectivity",file:"EX_M5-2",mod_type:"IBM", pic:"group"});
		
		ex_mod.push({mod:"M5.3", te:"Quantitative genetics model for susc./inf.",file:"EX_M5-3",mod_type:"IBM", pic:"groupQG"});

		ex_mod.push({mod:"M5.4", te:"Environmental pathogen accumulation 	model",file:"EX_M5-4",mod_type:"IBM/POP", pic:"envpathogen"});
		
		ex_mod.push({te:"M6) COVID-19 models",help:"These examples look at a variety of different models aimed at understanding a real-world epidemic outbreak."+mod_mess});
		
		ex_mod.push({mod:"M6.1", te:"Simple",file:"EX_M6-1",mod_type:"POP", pic:"covid"});
		ex_mod.push({mod:"M6.2", te:"Age-structured model",file:"EX_M6-2",mod_type:"POP", pic:"covidage"});
	}
	else{
		/*
		ex_mod.push({mod:"M1.1", te:"SI population-based model (PBM)",file:"EX_M1-1",mod_type:"POP", pic:"SI"});
		ex_mod.push({mod:"M1.2", te:"SI individual-based model (IBM)",file:"EX_M1-2",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"M2.6", te:"SIRD model (PBM) with demographic stratification",file:"EX_M2-6",mod_type:"POP", pic:"SIR_sex"});
		ex_mod.push({mod:"M2.8", te:"SIR model (IBM) with demographic stratification",file:"EX_M2-8",mod_type:"POP", pic:"SIR_sex2"});
		ex_mod.push({mod:"M1.3", te:"SIR model (IBM)",file:"EX_M1-3",mod_type:"IBM", pic:"SIR"});
		*/
	}
	
	let gap = 0.3;
	
	ex_mod.push({te:"A) Simulation features and initial conditions",help:"These examples look at various simulation features and ways to implement initial conditions."+mod_mess, gap:gap});
	ex_mod.push({te:"A1: Multiple simulations", mod:"M1.1", link:true, mod_type:"POP"});
	ex_mod.push({te:"A2: Uncertain initial conditions for PBM — single classification", mod:"M1.1", link:true, mod_type:"POP"});
	ex_mod.push({te:"A3: Uncertain initial conditions for PBM — multiple classifications — focal selected", mod:"M2.6", link:true, mod_type:"POP"});
	ex_mod.push({te:"A4: Uncertain initial conditions for PBM — multiple classifications — total population selected", mod:"M2.6",link:true, mod_type:"POP"});
	ex_mod.push({te:"A5: Uncertain initial conditions for IBM using individual state", mod:"M1.2",link:true, mod_type:"IBM"});
	ex_mod.push({te:"A6: Uncertain initial conditions for IBM using population distribution", mod:"M1.2",link:true, mod_type:"IBM"});
	ex_mod.push({te:"A7: Add / remove individuals from PBM", mod:"M2.6", link:true, mod_type:"POP"});
	ex_mod.push({te:"A8: Add / move / remove individuals from IBM",link:true, mod_type:"IBM"});
	
	ex_mod.push({te:"B) Population-level data types",help:"These examples look at different types of population data."+mod_mess, gap:gap});
	ex_mod.push({te:"B1: Time series population observations", mod:"M1.1", link:true, mod_type:"POP"});
	ex_mod.push({te:"B2: Time series population-level transition observations for PBM", mod:"M1.1", link:true, mod_type:"POP"});
	ex_mod.push({te:"B3: Stratified time series population observations for PBM", mod:"M2.6",link:true, mod_type:"POP"});
	ex_mod.push({te:"B4: Population observations from multiple compartments", mod:"M2.6", link:true, mod_type:"POP"});
	ex_mod.push({te:"B5: Combined population-based data sources in a COVID-19 model", mod:"M6.1", link:true, mod_type:"POP"});
	ex_mod.push({te:"B6: Time series population observations with IBM", mod:"M1.2", link:true, mod_type:"IBM"}); 

	ex_mod.push({te:"C) Individual-level data types",help:"These examples look at different types of individual-level data."+mod_mess, gap:gap});
	ex_mod.push({te:"C1: Known transition events — infection and recovery", mod:"M1.3", link:true, mod_type:"IBM"});

	ex_mod.push({te:"C2: Incomplete transition events - recovery only", mod:"M1.3", link:true,  mod_type:"IBM"}); 
	
	ex_mod.push({te:"C3: Compartmental observations", mod:"M1.3", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"C4: Disease diagnostic test results", mod:"M1.3", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"C5: A partially observed transition", mod:"M1.3", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"C6: A transition observed over a time window", mod:"M1.3", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"C7: A transition observed in a demographic category", mod:"M2.8", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"C8: Uncertain compartmental observations", mod:"M1.3", link:true,  mod_type:"IBM"}); 
	
	ex_mod.push({te:"D) Time variation",help:"These examples look at different ways in which time variation can be added into the model"+mod_mess, gap:gap});

	ex_mod.push({te:"D1: Time variation in transmission rate", mod:"M1.1", link:true, mod_type:"POP"}); 
	
	ex_mod.push({te:"D2: Time variation in transmission rate using a trigonometric function", mod:"M1.1", link:true, mod_type:"POP"}); 
	
	ex_mod.push({te:"D3: Time variation in transmission rate through a covariate", mod:"M1.1", link:true, mod_type:"POP"}); 
	
	ex_mod.push({te:"D4: Time variation in population-level transition observation probability", mod:"M1.1", link:true, mod_type:"POP"}); 

	ex_mod.push({te:"D5: Time variation in individual transition observation probability", mod:"M1.2", link:true, mod_type:"IBM"}); 
		
	ex_mod.push({te:"D6: Time variation in population observation probability", link:true,  mod_type:"POP"}); 
	
	ex_mod.push({te:"D7: Time variation in branching probability", mod:"M2.1", link:true,  mod_type:"POP"}); 
	
	ex_mod.push({te:"D8: Time-varying covariate affecting branching probability", mod:"M2.2", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"E) Individual-based variation",help:"These examples look at how individual-based variation can be incorporated into the model"+mod_mess, gap:gap});
	ex_mod.push({te:"E1: Individual fixed effect applied to a transition", mod:"M1.2", link:true, mod_type:"IBM"}); 
	ex_mod.push({te:"E2: Individual fixed effect applied to a population", link:true, mod_type:"IBM"}); 
	ex_mod.push({te:"E3: Individual effect applied to a transition", mod:"M1.2", link:true, mod_type:"IBM"}); 
	ex_mod.push({te:"E4: Correlated individual effect applied to a transition", mod:"M1.2", link:true, mod_type:"IBM"});
	ex_mod.push({te:"E5: Correlated individual effect applied to a population", mod:"M5.3", link:true, mod_type:"IBM"});
	ex_mod.push({te:"E6: Individual fixed effect applied to a branching probability", mod:"M2.3", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"E7: Individual effect applied to a branching probability", mod:"M2.3", link:true, mod_type:"IBM"}); 
	ex_mod.push({te:"E8: Correlated individual effect applied to a transition with pedigree", mod:"M1.2", link:true, mod_type:"IBM"});

	ex_mod.push({te:"F) Parameter definitions",help:"These examples look at how parameters can be reparameterised, set to distributions or derived from other quantities"+mod_mess, gap:gap});
	ex_mod.push({te:"F1: Reparameterisation", mod:"M3.4", link:true, mod_type:"IBM"}); 

	ex_mod.push({te:"F2: Parameter distribution", mod:"M5.2", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"F3: Derived quantities", mod:"M1.3", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"G) Incorporating pathogen genetics",help:"", gap:gap});
	ex_mod.push({te:"G1: Matrix of genetic differences", mod:"M1.2", link:true, mod_type:"IBM"}); 
	
	// Adds file names
	for(let i = 0; i < ex_mod.length; i++){
		let em = ex_mod[i];
		if(em.link){
			let spl = em.te.split(":");
			if(spl.length != 2) error("Does not split in two");
			em.file = "EX_"+spl[0];
		}
	}
	
	let ncol = 3;
	let cx = 0;
	let cy = 0.5;
	let gapx = 2;
	let gapy = 1;
	
	let W = lay.dx-2;
	
	let dx = (W-2*cx-(ncol-1)*gapx)/ncol;
	let dy = (0.8*dx);

	let ibm = find_pic("IBM");
	let pop = find_pic("POP");

	for(let j = 0; j < ex_mod.length; j++){ 
		let ex = ex_mod[j];
		
		let mod_but;
	
		if(ex.mod){
			for(let k = 0; k < j; k++){
				if(ex_mod[k].mod == ex.mod){ mod_but = ex_mod[k]; break;}
			}
		}
		
		let mod_ty, mod_ty2;
		if(ex.mod_type){
			switch(ex.mod_type){
			case "IBM":
				mod_ty = ibm; 
				break;
				
			case "POP":
				mod_ty = pop; 
				break;
				
			case "IBM/POP": 
				mod_ty = ibm; mod_ty2 = pop;
				break;
			}
		}
			
		let sel = false; if(ex.file == model.example) sel = true;
		
		if(!ex.pic){
			if(ex.link){
				cy = lay.add_example_link(ex.te,ex.file,cx+1,cy,mod_ty,mod_ty2,ex.mod_type,mod_but,sel);
			}
			else{
				cy += 1;		
				
				if(cx != 0){ cx = 0; cy += dy+gapy;}
				
				let te = ex.te;
				let si = 1;
				let font = get_font(si,"bold");
				let w = text_width(te,font)+0.3;
				
				let x = 0;
				
				lay.add_button({te:te, x:x, y:cy-0.2, dx:w, dy:si, type:"LeftText", font:font}); 
			
				lay.add_help_button(x+w,cy+si-0.2+0.2,{title:te, te:ex.help});
				
				cy += 0.8;
				if(ex.gap) cy += ex.gap;
			}
		}
		else{		
			let pic = find_pic(ex.pic);
		
			let si = 0.9, lh = 1.1, ma = 0.5;  
			let text_anno = text_convert_annotation("<b>"+ex.mod+": "+ex.te+"</b>",si,lh,dx-2*ma,"center",DDBLUE);
			let word = text_anno.word;
			
			for(let i = 0; i < word.length; i++){
				word[i].x += ma, word[i].y += dy-0.5*text_anno.height-2;
			}
			
			lay.add_button({word:word, x:cx, y:cy, dx:dx, dy:dy, ac:"ExampleModel", type:"ExampleModel", mod_ty:mod_ty, mod_ty2:mod_ty2, pic:pic, sel:sel, file:ex.file}); 
			cx += dx+gapx; if(cx > W){ cx = 0; cy += dy+gapy;}
		}
	}
	
	if(cx != 0){ cx = 0; cy += dy+gapy;}
}


/// Closes the edit description view
function close_description()
{
	if(model.description) model.description.edit = false;
}
