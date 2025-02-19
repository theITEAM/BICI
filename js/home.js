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
function add_examples_buts2(lay) 
{
	let mod_mess = "\nEach example can be simulated from to see the dynamics it exhibits (using different visualisations). These simulated results can be used to generate data onto which inference can then be performed.";
	
	let ex_mod=[];
	
	if(1 == 1){
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
		
		ex_mod.push({te:"M5) Disease transmission experiments",help:"Disease transmission experiments can be used to identify factors affecting disease transmission. These can either be fixed categorical effects, such as vaccination status, covariates, such as weight, or quantitative genetic traits."+mod_mess});
		
		ex_mod.push({mod:"M5.1", te:"Single contact group, investigating susceptibility",file:"EX_M5-1",mod_type:"IBM", pic:"dte"});

		ex_mod.push({mod:"M5.2", te:"Multiple contact groups, investigating infectivity",file:"EX_M5-2",mod_type:"IBM", pic:"group"});
		
		ex_mod.push({mod:"M5.3", te:"Quantitative genetics model for susceptibility/infectivity",file:"EX_M5-3",mod_type:"IBM", pic:"groupQG"});

		ex_mod.push({mod:"M5.4", te:"Environmental pathogen accumulation model",file:"EX_M5-4",mod_type:"IBM/POP", pic:"envpathogen"});

		ex_mod.push({te:"M6) Covid-19 models",help:"These examples look at a variety of different models aimed at understanding a real-world epidemic outbreak."+mod_mess});
		
		ex_mod.push({mod:"M6.1", te:"Simple",file:"EX_M6-1",mod_type:"POP", pic:"covid"});
		ex_mod.push({mod:"M6.2", te:"Age-structured model",file:"EX_M6-2",mod_type:"POP", pic:"covidage"});
	}
	else{
		ex_mod.push({mod:"M1.1", te:"SI population-based model (PBM)",file:"EX_M1-1",mod_type:"POP", pic:"SI"});
		ex_mod.push({mod:"M1.2", te:"SI individual-based model (IBM)",file:"EX_M1-2",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"M2.6", te:"SIRD model (PBM) with demographic stratification",file:"EX_M2-6",mod_type:"POP", pic:"SIR_sex"});
		ex_mod.push({mod:"M2.8", te:"SIR model (IBM) with demographic stratification",file:"EX_M2-8",mod_type:"POP", pic:"SIR_sex2"});
		ex_mod.push({mod:"M1.3", te:"SIR model (IBM)",file:"EX_M1-3",mod_type:"IBM", pic:"SIR"});
	}
	
	let gap = 0.3;
	
	
	ex_mod.push({te:"A) Simulation features",help:"These examples look at various simulation features."+mod_mess, gap:gap});
	ex_mod.push({te:"A1: Multiple simulations", mod:"M1.1", link:true, mod_type:"POP"});
	ex_mod.push({te:"A2: Uncertain initial conditions — single classification", mod:"M1.1", link:true, mod_type:"POP"});
	ex_mod.push({te:"A3: Uncertain initial conditions — multiple classifications — focal selected", mod:"M2.6", link:true, mod_type:"POP"});
	ex_mod.push({te:"A4: Uncertain initial conditions — multiple classifications — total population selected", mod:"M2.6",link:true, mod_type:"POP"});
	
	ex_mod.push({te:"A5: Add / remove individuals from PBM", mod:"M2.6", link:true, mod_type:"POP"});
		
	ex_mod.push({te:"A6: Add / move / remove individuals from IBM",link:true, mod_type:"IBM"});
	
	ex_mod.push({te:"B) Population-level data types",help:"These examples look at different types of population data."+mod_mess, gap:gap});
	ex_mod.push({te:"B1: Time series population observations", mod:"M1.1", link:true, mod_type:"POP"});
	ex_mod.push({te:"B2: Time series population-level transition observations for PBM", mod:"M1.1", link:true, mod_type:"POP"});
	ex_mod.push({te:"B3: Stratified time series population observations for PBM", mod:"M2.6",link:true, mod_type:"POP"});
	ex_mod.push({te:"B4: Population observations from multiple compartments", mod:"M2.6", link:true, mod_type:"POP"});
	ex_mod.push({te:"B5: Combined population-based data sources in a Covid-19 model", mod:"M6.1", link:true, mod_type:"POP"});
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
	ex_mod.push({te:"E2: Fixed-effect applied to a population", link:true, mod_type:"IBM"}); 
	ex_mod.push({te:"E3: Individual effect applied to a transition", mod:"M1.2", link:true, mod_type:"IBM"}); 
	ex_mod.push({te:"E4: Correlated individual effect applied to a transition", mod:"M1.2", link:true, mod_type:"IBM"});
	ex_mod.push({te:"E5: Correlated individual effect applied to a population", mod:"M5.3", link:true, mod_type:"IBM"});
	ex_mod.push({te:"E6: Individual fixed effect applied to a branching probability", mod:"M2.3", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"E7: Individual effect applied to a branching probability", mod:"M2.3", link:true, mod_type:"IBM"}); 
	

	ex_mod.push({te:"F) Parameter definitions",help:"These examples look at how parameters can be reparameterised, set to distributions or derived from other quantities"+mod_mess, gap:gap});
	ex_mod.push({te:"F1: Reparameterisation", mod:"M3.4", link:true, mod_type:"IBM"}); 

	ex_mod.push({te:"F2: Parameter distribution", mod:"M5.2", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"F3: Derived quantities", mod:"M1.3", link:true, mod_type:"IBM"}); 
	
	ex_mod.push({te:"G) Incorporating pathogen genetics",help:"", gap:gap});
	ex_mod.push({te:"G1: Matrix of genetic differences", mod:"M1.2", link:true, mod_type:"IBM"}); 
		
	
	/*	ex_mod.push({te:"M1) Simple epidemiological models",help:"This gives some examples of simple epidemiolgoical model, illustrating different ways in which BICI can be used."+mod_mess});
		
		ex_mod.push({mod:"M1.1", te:"SI population-based model (PBM)",file:"EX_M1-1",mod_type:"POP", pic:"SI"});
		ex_mod.push({mod:"M1.2", te:"SI individual-based model (IBM)",file:"EX_M1-2",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"M1.3", te:"SIR model (IBM)",file:"EX_M1-3",mod_type:"IBM", pic:"SIR"});
		ex_mod.push({mod:"M1.4", te:"SIR model (PBM) with Erlang distribution",file:"EX_M1-4",mod_type:"POP", pic:"SIRErlang"});
		ex_mod.push({mod:"M1.5", te:"SIR model (IBM) with gamma distribution",file:"EX_M1-5",mod_type:"IBM", pic:"SIRgamma"});
		ex_mod.push({mod:"M1.6", te:"SEIR model (IBM) with exposed period",file:"EX_M1-6",mod_type:"IBM", pic:"SEIR"});
	
		ex_mod.push({te:"M2) Additional epidemiological models",help:"This looks at some additional features, such as branching, reverse transitions and stratification."+mod_mess});
		
		ex_mod.push({mod:"M2.1", te:"SIRD model (PBM) with branching rate",file:"EX_M2-1",mod_type:"POP", pic:"SIRD"});
		ex_mod.push({mod:"M2.2", te:"SIRD model (IBM) with branching probability",file:"EX_M2-2",mod_type:"IBM", pic:"SIRD_IBM"});
		ex_mod.push({mod:"M2.3", te:"SIRD model (IBM) with branching factors",file:"EX_M2-3",mod_type:"IBM", pic:"SIRD_IBM2"});

		ex_mod.push({mod:"M2.4", te:"SIS model (PBM)",file:"EX_M2-4",mod_type:"POP", pic:"SIS"});
		ex_mod.push({mod:"M2.5", te:"SIRS model (IBM) with waning immunity",file:"EX_M2-5",mod_type:"IBM", pic:"SIRwaning"});
		ex_mod.push({mod:"M2.6", te:"SIR model (PBM) with demographic stratification",file:"EX_M2-6",mod_type:"POP", pic:"SIR_sex2"});
		ex_mod.push({mod:"M2.7", te:"SIR mode (PBM) with differential infectivity and demographic stratification",file:"EX_M2-7",mod_type:"POP", pic:"SIR_sex"});
		ex_mod.push({mod:"M2.8", te:"SIR model (IBM) with demographic stratification",file:"EX_M2-8",mod_type:"POP", pic:"SIR_sex2"});
		
		ex_mod.push({te:"M3) Spatial epidemiological models",help:"The models incorporate spatial stratification in different ways."+mod_mess});
		
		ex_mod.push({mod:"M3.1", te:"Metapopulation model using geographical regions",file:"EX_M3-1",mod_type:"POP", pic:"scotland"});
		ex_mod.push({mod:"M3.2", te:"Metapopulation model using geographical points",file:"EX_M3-2",mod_type:"POP", pic:"world"});
		ex_mod.push({mod:"M3.3", te:"Metapopulation model using a distance kernel",file:"EX_M3-3",mod_type:"POP", pic:"uk"});
		ex_mod.push({mod:"M3.4", te:"Farm-based model using a distance kernel",file:"EX_M3-4",mod_type:"IBM", pic:"farm"});
		
		ex_mod.push({te:"M4) Ecological models",help:"Comaprtmental model are often applied to epidemilogical applications, but they can equally be used in other settings. Here we look at some examples from ecology."+mod_mess});
		ex_mod.push({mod:"M4.1", te:"Logistic growth population model",file:"EX_M4-1",mod_type:"IBM", pic:"carry"});
		ex_mod.push({mod:"M4.2", te:"Predator-prey model",file:"EX_M4-2",mod_type:"POP", pic:"predprey"});
		ex_mod.push({mod:"M4.3", te:"Spatial diffusion model",file:"EX_M4-3",mod_type:"IBM", pic:"grid"});
		
		ex_mod.push({te:"M5) Disease transmission experiments",help:"Disease transmission experiments can be used to identify factors which affect disease transmission. These can either be fixed categorical effects, such as vaccination status, covariated, such as weight, or estimation of quantitative genetic traits."+mod_mess});
		
		ex_mod.push({mod:"M5.1", te:"Single contact group investigating susceptibility",file:"EX_M5-1",mod_type:"IBM", pic:"dte"});

		ex_mod.push({mod:"M5.2", te:"Multiple contact groups investigating infectivity",file:"EX_M5-2",mod_type:"IBM", pic:"group"});
		
		ex_mod.push({mod:"M5.3", te:"Quantitative genetics model for susc./infect.",file:"EX_M5-3",mod_type:"IBM", pic:"groupQG"});

		ex_mod.push({mod:"M5.4", te:"Environmental pathogen accumulation model",file:"EX_M5-4",mod_type:"IBM/POP", pic:"envpathogen"});

		ex_mod.push({te:"M6) Covid-19 models",help:"This looks at a variety of different models."+mod_mess});
		
		ex_mod.push({mod:"M6.1", te:"Simple",file:"EX_M6-1",mod_type:"POP", pic:"covid"});
		ex_mod.push({mod:"M6.2", te:"Age structured model",file:"EX_M6-2",mod_type:"POP", pic:"covidage"});
		*/
	
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

/*
	
	


	
	{
		let ex=[];
		ex.push({te:"SIR model with specified initial population",file:"Example 1"});
		ex.push({te:"stratified model with classification initial population",file:"Example 2"});
		ex.push({te:"stratified model with classification initial population",file:"Example 3"});
		ex.push({te:"multiple simulations",file:"Example 4"});
		ex.push({te:"movement and removal data",file:"Example 5"});
		
		example.push({model:"Simulation features", ex:ex, help_title:" "});
	}
	
	
	
	{
		let ex=[];
		ex.push({te:"population observations",file:"Example 1"});
		ex.push({te:"population-level transition observations",file:"Example 2"});
		
		example.push({model:"population data types", ex:ex, help_title:" "});
	}
	
	{
		let ex=[];
		ex.push({te:"temporal covariate",file:"Example 1"});
		ex.push({te:"spatial covariate",file:"Example 2"});
		
		example.push({model:"estimating covariate", ex:ex, help_title:" "});
	}
	

	{
		let ex=[];
		ex.push({te:"Single classification",file:"Example 1"});
		ex.push({te:"Stratified model wirh focal classifcation",file:"Example 2"});
		ex.push({te:"Stratified model wirhout focal classifcation",file:"Example 3"});
	
		example.push({model:"Inferring the initial state", ex:ex, help_title:" "});
	}


	{
		let ex=[];
		ex.push({te:"mark-recapture data",file:"Example 1"});
		ex.push({te:"geographic-diffusion model",file:"Example 2"});
	
		example.push({model:"Further applications", ex:ex, help_title:" "});
	}
	
	{
		let ex=[];
		ex.push({te:"predicting the future",file:"Example 1"});
		ex.push({te:"conterfactual",file:"Example 2"});
		ex.push({te:"scenario analysis",file:"Example 3"});
		
		example.push({model:"posterior simulation", ex:ex, help_title:" "});
	}
	
	{
		let ex=[];
		ex.push({te:"adding annotation to the model",file:"Example 1"});
		ex.push({te:"",file:"Example 2"});
		ex.push({te:"",file:"Example 3"});
		ex.push({te:"",file:"Example 4"});
		ex.push({te:"",file:"Example 5"});
		ex.push({te:"",file:"Example 6"});
		ex.push({te:"",file:"Example 7"});
		
		example.push({model:"Miscelaneous", ex:ex, help_title:" "});
	}


	let cx = 0, cy = 1;
	for(let i = 0; i < example.length; i++){
		let exam = example[i];
	
		for(let j = 0; j < exam.ex.length; j++){
			let e = exam.ex[j];
			cy = lay.add_example_link(e.te,e.file,cx,cy);
		}
	}
	*/
}

/// Adds buttons for example section
function add_examples_buts(lay) 
{
	let example=[];
		
	{
		let ex=[];
			ex.push({te:"Test",file:"test"});
	
		ex.push({te:"Individual-based model (IBM) with compartmental obsevations.",file:"SIR_ind"});
		ex.push({te:"Population-based model (PBM) with infected population obs.",file:"SIR_pop"});
		ex.push({te:"IBM with infection and recovery time obsevations.",file:"SIR_times"});
		ex.push({te:"IBM with disease diagnostic test obsevations.",file:"SIR_diagtest"});
		ex.push({te:"IBM with susceptibility individual variation from a fixed effect.",file:"SIR_fix"});
		ex.push({te:"IBM with susceptibility individual variation from a random effect.",file:"SIR_rand"});
		ex.push({te:"IBM with susceptibility quantitative genetics model.",file:"SIR_qg"});
		
		
		//ex.push({te:"EX 1: Complete knowledge of events",file:"EX1"});
		//ex.push({te:"EX 2: Known initial state and recoveries",file:"EX2"});
		//ex.push({te:"EX 3: Known recoveries only",file:"EX3"});
		//ex.push({te:"EX 4: Periodic disease status data",file:"EX4"});
		//ex.push({te:"EX 5: Diagnostic test results",file:"EX5"});
		//ex.push({te:"EX 6: Future and past predition",file:"EX6"});
		//ex.push({te:"EX 7: Time classification",file:"EX7"});
		
		example.push({model:"SIR model", dx:14, ex:ex, pic:"SIR", help_title:"The SIR model", help_text:"This is a simple model used to describe epidemic behaviour. Individuals are classified as being either susceptible to infection (S), infected and infectious (I), or recovered / removed / dead (R).\nA susceptible individual has a probability per unit time of becoming infected as a result of other infected individuals. For those individuals that do become infected, a recovery rate determines the duration over which they are infectious."}); 
	}

	{
		let ex=[];
		ex.push({te:"IBM with sex stratification of transition rate.",file:"SIR_sex"});
		ex.push({te:"Time variation.",file:"time_vary"});
		
		example.push({model:"SIR stratified model", dx:12, ex:ex, pic:"SIRsex", help_title:"The stratified SIR model", help_text:"SIR models can be stratified by individual attributes. Here we consider splitting individuals into males and female which have a different transmission rate."}); 
	}
	
	{
		let ex=[];
		
		ex.push({te:"Non-Markovian latent period",file:"SEIR"});
		
		//ex.push({te:"EX 8: Non-Markovian latent period",file:"EX8"});
		//ex.push({te:"EX 9: Disease diagnostic tests",file:"EX9"});
		//ex.push({te:"EX 10: Environmental and diagnostic tests",file:"EX10"});
		//ex.push({te:"EX 11: Estimating test Se and Sp",file:"EX11"});

		example.push({model:"SEIR model", dx:14, ex:ex, pic:"SEIR", help_title:"The SEIR model", help_text:"This is a simple model used to describe epidemic behaviour. Individuals are classified as being either susceptible to infection (S), exposed (which means infected but not infectious) (E), infectious (I), or recovered / removed / dead (R).\nA susceptible individual has a probability per unit time of becoming infected as a result of other infected individuals. For those individuals that do become infected, a latent period determines when they become infectious and a recovery rate determines the duration over which they are infectious."}); 
	}
	
	{
		let ex=[];
		
		ex.push({te:"A disease transition experiment",file:"TransExp"});
		ex.push({te:"Incorporating individual-based variation",file:"TransExpIndVar"});
		ex.push({te:"Transmission through the environment",file:"SI_Env"});
		
		//ex.push({te:"EX 12: Known infection times and initial state",file:"EX12"});
		//ex.push({te:"EX 13: Diagnostic test results",file:"EX13"});
		//ex.push({te:"EX 14: Initial / final disease status",file:"EX14"});
		//ex.push({te:"EX 15: Incorporating vaccination status",file:"EX15"});
		//ex.push({te:"EX 16: Incorporating group effect",file:"EX16"});
		//ex.push({te:"EX 17: SNP effects on susceptiblity / infectivity",file:"EX17"});
		//ex.push({te:"EX 18: SIR with non-Markovian recovery",file:"EX18"});
	
		example.push({model:"Disease transmission experiment", dx:14, ex:ex, pic:"transexp", help_title:"Disease transmission experiment", help_text:"These are used to discover the rate at which individuals become infected and factors affecting this rate (e.g. the role of fixed effects such as vaccination status). \nTypically populations of initially infected and uninfected individuals are placed in closed 'contact groups'. Data is them collected on the ensuing epidemics (in a variety of ways) to estimate model parameters."}); 
	}
	
	{
		let ex=[];
		ex.push({te:"Population model with carrying capacity",file:"SimplePop"});
		
		//ex.push({te:"EX 19: Periodic population estimates",file:"EX19"});
		//ex.push({te:"EX 20: Death times",file:"EX20"});
		//ex.push({te:"EX 21: Captures",file:"EX21"});
		//ex.push({te:"EX 22: Age dependent mortality",file:"EX22"});
		//ex.push({te:"EX 23: Captures and a fraction of death times",file:"EX23"});
	
		example.push({model:"Population with births and deaths", dx:14, ex:ex, pic:"Simple population", help_title:"Population with births and deaths", help_text:"This simple model can be used to capture stochastic variation in the population number in a wildlife setting."}); 
	}
	
	{
		let ex=[];
		ex.push({te:"Population model with carrying capacity",file:"SI_source_sink"});
		
		//ex.push({te:"EX 24: Periodic disease status measurements",file:"EX24"});
		//ex.push({te:"EX 25: Disease status from captures",file:"EX25"});
		//ex.push({te:"EX 26: Disease diagnostic tests from captures",file:"EX26"});
		//ex.push({te:"EX 27: Disease induced mortality",file:"EX27"});
		
		example.push({model:"SIS model with births and deaths", dx:14, ex:ex, pic:"SIbd", help_title:"SIS model with births and deaths", help_text:"This model not only captures stochastic variation in population numbers in a wildlife setting, but also accounts for changes in individual disease status. As well as the usual infection transition from 'S' to 'I', this model includes the reverse transition to account for the fact that individuals do not gain lifelong immunity (which is true for some diseases, such as influenza)."}); 
	}
	
	
	{
		let ex=[];
		
		ex.push({te:"Dynamic from the stochastic Lotka-Volterra model",file:"PredPrey"});
	
		example.push({model:"Predator-prey model", dx:14, ex:ex, pic:"LV", help_title:"Predator-prey model", help_text:"These stochastic models predict the dynamic variation in predator and prey populations. Typically the food supply for prey (e.g. rabbits) is plentiful. Rather than being limited by these local resources, the number of prey is regulated by being consumed by predators (e.g. foxes). On the other hand, predators rely on prey for food, which leads to unstable dynamic behaviour. Here a simple Lotka-Volterra model is assumed."}); 
	}
	
	
	{
		let ex=[];
		ex.push({te:"A local authority model of Scotland",file:"Scotland_Region"});
		example.push({model:"Metapopulations", dx:10, ex:ex, pic:"scotland_region", help_title:"Spatial model", help_text:""}); 
	}
	
	/*
	{
		let ex=[];
		ex.push({te:"A UK meta population model",file:"metapopUK"});
		example.push({model:"Metapopulations", dx:10, ex:ex, pic:"metapopUK", help_title:"Spatial model", help_text:""}); 
	}
	*/
	
	
	/*
	{
		let ex=[];
		ex.push({te:"EX 30: State measurements at all locations",file:"EX30"});
		ex.push({te:"EX 31: Captures at some locations",file:"EX31"});
		example.push({model:"Spatial diffusion model", ex:ex, pic:"grid2", help_title:"Spatial diffusion model", help_text:"This simple model captures the movement of individuals across a landscape."}); 
	}
	
	{
		let ex=[];
		ex.push({te:"EX 32: Location and disease status",file:"EX32"});
		ex.push({te:"EX 33: Captures with disease diagnostic tests",file:"EX33"});
		example.push({model:"Spatial disease model", ex:ex, pic:"gridSI", help_title:"Spatial disease model", help_text:"This model captures not only the movement of individuals across a landscape, but also their disease status (which is represented by a simple SI model)."}); 
	}
	
	{
		let ex=[];
		ex.push({te:"EX 34: Power-law spatial kernel",file:"EX34"});
		ex.push({te:"EX 35: Exponential spatial kernel",file:"EX35"});
		example.push({model:"Location-based spatial disease model", ex:ex, pic:"grid2", help_title:"Location-based spatial disease model", help_text:"This model imagines that the 'individuals' themselves are fixed, e.g. they could represent farms. Each farm is classified as being either susceptible to infection (S), infected (I), or recovered (R). Disease is spread spatially by means of a transition kernel."}); 
	}
	
	{
		let ex=[];
		ex.push({te:"EX 36: A model of badger social groups",file:"EX34"});
		example.push({model:"Spatial disease model", ex:ex, pic:"badger", help_title:"Spatial disease model with births and deaths", help_text:"This model assumes spatially distributed groups of individuals (e.g. badger social groups). These individuals are born, die, can move from location to location and are classified by disease status (through an SI model)."}); 
	}
	*/
	
	let cx_link = lay.dx*0.42;
		
	let cy = 1;
	let cx_bracket = cx_link - lay.dx*0.03;
	
	let cx_mid = lay.dx*0.18;
	
	for(let i = 0; i < example.length; i++){
		let exam = example[i];
	
		let pic = find_pic(exam.pic);

		let picx = pic.width;
		let picy = pic.height;
		
		let ratio = picy/picx;
		
		let dxx = exam.dx;
		
		let dyy = dxx*ratio;
		
		let picheight = dyy+2;
		let padding = 0;
		if(picheight > exam.ex.length*example_space) padding = (picheight-(exam.ex.length*example_space))/2;
		
		cy += padding;
			
		let cyst = cy;
		
		for(let j = 0; j < exam.ex.length; j++){
			let e = exam.ex[j];
			cy = lay.add_example_link(e.te,e.file,cx_link,cy);
		}
		
		lay.add_button({x:cx_bracket, y:cyst-0.2, dx:0.6, dy:(cy-cyst+0.4), type:"Bracket"}); 
	
		lay.example_picture(exam.model,cx_mid-dxx/2,(cy+cyst)/2-dyy/2,dxx,dyy,pic,{title:exam.help_title, te:exam.help_text});
		cy += padding;
	
		cy += 1.3;
	}
}	
	
	
/// Closes the edit description view
function close_description()
{
	if(model.description) model.description.edit = false;
}
