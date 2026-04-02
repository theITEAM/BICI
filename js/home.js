"use strict";
// Functions for generating the home page

/// Adss buttons for home page
function add_home_page_buts(lay) 
{
	let cx = corner.x;
	let cy = corner.y;
	
	let sel = inter.section_sel;
	if(sel == undefined){
		cy = lay.add_title("Welcome to BICI!",cx,cy,{te:welcome_text});
	}
	else{
		let sec = inter.section[sel];
		cy = lay.add_title(sec.desc+" ("+sec.name+")",cx,cy,{te:sec.help});
	}
	
	lay.add_button({x:39, y:cy-2.2, dx:7, dy:1.5, ac:"NewModel", type:"NewModel"});
		
	let x1 = lay.x+cx;
	let x2 = lay.x+lay.dx-cx;
	let y1 = lay.y+cy+0.5;
	let y2 = lay.y+lay.dy-1;
		
	if(sel == undefined){
		x2 += 1;
		add_layer("Sections",x1,y1,x2-x1,y2-y1,{});
	}
	else{
		y2 -= 3;
		add_layer("Examples",x1,y1,x2-x1,y2-y1,{});
		
		lay.add_corner_button([["Back","Grey","BackSection"]],{x:lay.dx-button_margin.dx, y:lay.dy-button_margin.dy});
	}
}


/// Creates a tree of sections
function initialise_section()
{
	let mod_mess = "\nEach example can be simulated from to see the dynamics it exhibits (using different visualisations). These simulated results can be used to generate data onto which inference can then be performed.";
	
	let gap = 0.3;
		
	{
		let ex_mod=[];
		ex_mod.push({mod:"SE1.1", te:"SI population-based model (PBM)",file:"EX_SE1-1",mod_type:"POP", pic:"SI"});
		ex_mod.push({mod:"SE1.2", te:"SI individual-based model (IBM)",file:"EX_SE1-2",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"SE1.3", te:"SIR model (IBM)",file:"EX_SE1-3",mod_type:"IBM", pic:"SIR"});
		ex_mod.push({mod:"1.4", te:"SIR model (PBM) with Erlang distribution",file:"EX_SE1-4",mod_type:"POP", pic:"SIRErlang"});
		ex_mod.push({mod:"SE1.5", te:"SIR model (IBM) with gamma distribution",file:"EX_SE1-5",mod_type:"IBM", pic:"SIRgamma"});
		ex_mod.push({mod:"SE1.6", te:"SEIR model (IBM) with exposed period",file:"EX_SE1-6",mod_type:"IBM", pic:"SEIR"});
		
		let sec = { name:"simp-epi", desc:"Simple epidemiological", pic:"SI", ex_mod:ex_mod, help:simp_epi_text+mod_mess};
		
		inter.section.push(sec);
	}
	
	{
		let ex_mod=[];
	
		ex_mod.push({te:"1) Branching models",help:"Branching points differentiate the long-term behaviour of individuals."+mod_mess});
	
		ex_mod.push({mod:"CE1.1", te:"SIRD model (PBM) with branching using transition rates",file:"EX_CE1-1",mod_type:"POP", pic:"SIRD"});
		ex_mod.push({mod:"CE1.2", te:"SIRD model (IBM) with branching probability",file:"EX_CE1-2",mod_type:"IBM", pic:"SIRD_IBM"});
		ex_mod.push({mod:"CE1.3", te:"SIRD model (IBM) with branching factors",file:"EX_CE1-3",mod_type:"IBM", pic:"SIRD_IBM2"});

		ex_mod.push({te:"2) Waning immunity",help:"Reverse transitions allow for individuals to become reinfected due to waning immunity."+mod_mess});
	
		ex_mod.push({mod:"CE2.1", te:"SIS model (PBM)",file:"EX_CE2-1",mod_type:"POP", pic:"SIS"});
		ex_mod.push({mod:"CE2.2", te:"SIRS model (IBM) with waning immunity",file:"EX_CE2-2",mod_type:"IBM", pic:"SIRwaning"});
		
		
		ex_mod.push({te:"3) Demographic stratification",help:"This splits the population into subgroups based on shared characteristics like age and/or gender."+mod_mess});
		
		ex_mod.push({mod:"CE3.1", te:"SIR model (PBM) with demographic stratification",file:"EX_CE3-1",mod_type:"POP", pic:"SIR_sex2"});
		ex_mod.push({mod:"CE3.2", te:"SIR model (PBM) with differential infectivity and demographic stratification",file:"EX_CE3-2",mod_type:"POP", pic:"SIR_sex"});
		ex_mod.push({mod:"CE3.3", te:"SIR model (IBM) with demographic stratification",file:"EX_CE3-3",mod_type:"IBM", pic:"SIR_sex2"});
		
		ex_mod.push({te:"4) COVID-19 models",help:"These examples look at a variety of different models aimed at understanding a real-world epidemic outbreak."+mod_mess});
		
		ex_mod.push({mod:"CE4.1", te:"Simple COVID-19 model",file:"EX_CE4-1",mod_type:"POP", pic:"covid"});
		
		ex_mod.push({mod:"CE4.2", te:"Age-structured COVID-19 model",file:"EX_CE4-2",mod_type:"POP", pic:"covidage"});
		
	
		let sec = { name:"comp-epi", desc:"Complex epidemiological",  pic:"SIRD_IBM2", ex_mod:ex_mod, help:comp_epi_text+mod_mess};
		
		inter.section.push(sec);
	}

	{
		let ex_mod=[];
		ex_mod.push({te:"1) Metapopulation models",help:"Metapopulation models partition a population into discrete sub-populations, or 'patches' (such as cities or households), connected by movement or contact. Unlike simple compartmental models that assume everyone interacts equally, metapopulation models capture spatial heterogeneity and the role of mobility in spreading disease."+mod_mess});
	
		ex_mod.push({mod:"SPE1.1", te:"Metapopulation model using geographical regions",file:"EX_SPE1-1",mod_type:"POP", pic:"scotland"});
		ex_mod.push({mod:"SPE1.2", te:"Metapopulation model using geographical points",file:"EX_SPE1-2",mod_type:"POP", pic:"world"});
		ex_mod.push({mod:"SPE1.3", te:"Metapopulation model using a distance kernel",file:"EX_SPE1-3",mod_type:"POP", pic:"uk"});
		
		ex_mod.push({te:"2) Farm-based model",help:"These models work at the level of the status of the farm being the fundamental unit."+mod_mess});
	
		ex_mod.push({mod:"SPE2.1", te:"Farm-based model using a distance kernel",file:"EX_SPE2-1",mod_type:"IBM", pic:"farm"});
		ex_mod.push({mod:"SPE2.2", te:"Farm-based model with density dependency",file:"EX_SPE2-2",mod_type:"IBM", pic:"farm2"});
		
		let sec = { name:"spatial-epi", desc:"Spatial epidemiological", pic:"world", ex_mod:ex_mod, help:spatial_epi_text+mod_mess};
		
		inter.section.push(sec);
	}
	
	{
		let ex_mod=[];
		
		ex_mod.push({mod:"EE1.1", te:"Single contact group, investigating susceptibility",file:"EX_EE1-1",mod_type:"IBM", pic:"dte"});

		ex_mod.push({mod:"EE1.2", te:"Multiple groups, investigating infectivity",file:"EX_EE1-2",mod_type:"IBM", pic:"group"});
		
		ex_mod.push({mod:"EE1.3", te:"Environmental pathogen accumulation model",file:"EX_EE1-3",mod_type:"IBM/POP", pic:"envpathogen"});
		
		let sec = { name:"exp-epi", desc:"Disease transmission experiments", pic:"dte", ex_mod:ex_mod, help:exp_epi_text+mod_mess};
		
		inter.section.push(sec);
	}
	
	
	{
		let ex_mod=[];

		ex_mod.push({te:"Applications for individual effect",help:"These examples look at different ways in which individual effects can be incorporated into models."+mod_mess});

		ex_mod.push({mod:"EQG1.1", te:"Individual effect applied to a transition",file:"EX_EQG1-1",mod_type:"IBM", pic:"SIind"});
	
		ex_mod.push({mod:"EQG1.2", te:"Correlated individual effect applied to a transition",file:"EX_EQG1-2",mod_type:"IBM", pic:"SIind"});
		
		ex_mod.push({mod:"EQG1.3", te:"Correlated individual effect applied to a transition with pedigree",file:"EX_EQG1-3",mod_type:"IBM", pic:"SIind"});
			
		ex_mod.push({mod:"EQG1.4", te:"Correlated individual effect applied to a population",file:"EX_EQG1-4",mod_type:"IBM", pic:"SIindgroup"});
		
		ex_mod.push({mod:"EQG1.5", te:"Individual effect applied to a branching probability ",file:"EX_EQG1-5",mod_type:"IBM", pic:"branchie"});
		
		ex_mod.push({te:"Complete models",help:"These examples incorporate multiple individual/fixed effects, as would typically be done in a real analysis."+mod_mess});

		ex_mod.push({mod:"EQG2.1", te:"Quantitative genetics model for susc./inf.",file:"EX_EQG2-1",mod_type:"IBM", pic:"groupQG"});

		ex_mod.push({mod:"EQG2.2", te:"Quantitative genetics model for susc./inf./rec.",file:"EX_EQG2-2",mod_type:"IBM", pic:"qgSEIR"});
		
		let sec = { name:"epiqG", desc:"Quantitative genetics in epidemiology", pic:"groupQG", ex_mod:ex_mod, help:epiqg_text+mod_mess};
		
		inter.section.push(sec);
	}
	
	{
		
		let ex_mod=[];
		
		ex_mod.push({mod:"GE1.1", te:"Matrix of genetic differences",file:"EX_GE1-1", mod_type:"IBM", pic:"SI"});
	
		let sec = { name:"gen-epi", desc:"Pathogen genetics in epidemiology", pic:"SI", ex_mod:ex_mod, help:gen_epi_text+mod_mess};
		
		inter.section.push(sec);
	}
	
	{
		let ex_mod=[];
		ex_mod.push({mod:"SEC1.1", te:"Logistic growth population model",file:"EX_SEC1-1",mod_type:"IBM", pic:"carry"});
		ex_mod.push({mod:"SEC1.2", te:"A predator-prey model",file:"EX_SEC1-2",mod_type:"POP", pic:"predprey"});
		ex_mod.push({mod:"SEC1.3", te:"Logistic growth model with age dependent mortality",file:"EX_SEC1-3",mod_type:"IBM", pic:"carryage"});
		
		let sec = { name:"simp-eco", desc:"Simple ecological", pic:"predprey", ex_mod:ex_mod, help:simp_eco_text+mod_mess};
		
		inter.section.push(sec);
	}
	
	{
		let ex_mod=[];
		ex_mod.push({mod:"SPEC1.1", te:"A spatial diffusion model",file:"EX_SPEC1-1",mod_type:"IBM", pic:"grid"});
		ex_mod.push({mod:"SPEC1.2", te:"A species presence/absence distribution model",file:"EX_SPEC1-2",mod_type:"IBM", pic:"UK_grid"});
	
		let sec = { name:"spatial-eco", desc:"Spatial ecological", pic:"UK_grid", ex_mod:ex_mod, help:spatial_eco_text+mod_mess};
		
		inter.section.push(sec);
	}
	
	{
		let ex_mod=[];
		ex_mod.push({te:"1) Stochastic behaviour",help:"The underlying equations are inherently stochastic."+mod_mess});
	
		ex_mod.push({te:"S1.1: Multiple simulations", mod:"SE1.1", link:true, mod_type:"POP"});
		
		ex_mod.push({te:"2) Uncertainty in initial conditions",help:"Usually the initial populations are entirely fixed, but BICI allow for uncertainty in these quantities."+mod_mess});
	
		ex_mod.push({te:"S2.1: Uncertain initial conditions for PBM — single classification", mod:"SE1.1", link:true, mod_type:"POP"});
		ex_mod.push({te:"S2.2: Uncertain initial conditions for PBM — multiple classifications — focal selected", mod:"CE3.1", link:true, mod_type:"POP"});
		ex_mod.push({te:"S2.3: Uncertain initial conditions for PBM — multiple classifications — total population selected", mod:"CE3.1",link:true, mod_type:"POP"});
		ex_mod.push({te:"S2.4: Uncertain initial conditions for IBM using individual state", mod:"SE1.2",link:true, mod_type:"IBM"});
		ex_mod.push({te:"S2.5: Uncertain initial conditions for IBM using population distribution", mod:"SE1.2",link:true, mod_type:"IBM"});
		
		ex_mod.push({te:"3) Enfored movements",help:"Here enforce movement of individuals/populations are considered."+mod_mess});
	
		ex_mod.push({te:"S3.1: Add / remove individuals from PBM", mod:"CE3.1", link:true, mod_type:"POP"});
		ex_mod.push({te:"S3.2: Add / move / remove individuals from IBM",link:true, mod_type:"IBM"});
	
		let sec = { name:"simu-scen", desc:"Simulation scenarios", pic:"sim", shade:true, ex_mod:ex_mod, help:sim_scen_text+mod_mess};
		
		inter.section.push(sec);
	}
	
	{
		let ex_mod=[];
		
		
		ex_mod.push({te:"1) No data",help:"This example looks at the simplest possible scenario, that in which when there is no data (only the initial conditions are specified). Here the inferred parameters are expected to follow the prior and the state maps out the same distribution as if directly simulating from the model (albeit in a much less computationally efficient way due to correlations inherent in MCMC). "+mod_mess, gap:gap});
		
		ex_mod.push({te:"D1.1: No data (only initial conditions)", mod:"SE1.1", link:true, mod_type:"POP"});
		
		ex_mod.push({te:"2) Population-level data types",help:"These examples look at different types of population data."+mod_mess, gap:gap});
		ex_mod.push({te:"D2.1: Time series population observations", mod:"SE1.1", link:true, mod_type:"POP"});
		ex_mod.push({te:"D2.2: Time series population-level transition observations for PBM", mod:"SE1.1", link:true, mod_type:"POP"});
		ex_mod.push({te:"D2.3: Stratified time series population observations for PBM", mod:"CE3.1",link:true, mod_type:"POP"});
		ex_mod.push({te:"D2.4: Population observations from multiple compartments", mod:"CE3.1", link:true, mod_type:"POP"});
		ex_mod.push({te:"D2.5: Combined population-based data sources in a COVID-19 model", mod:"CE4.1", link:true, mod_type:"POP"});
		ex_mod.push({te:"D2.6: Time series population observations with IBM", mod:"SE1.2", link:true, mod_type:"IBM"}); 
		
		ex_mod.push({te:"3) Individual-level data types",help:"These examples look at different types of individual-level data."+mod_mess, gap:gap});
		ex_mod.push({te:"D3.1: Known transition events — infection and recovery", mod:"SE1.3", link:true, mod_type:"IBM"});

		ex_mod.push({te:"D3.2: Incomplete transition events - recovery only", mod:"SE1.3", link:true,  mod_type:"IBM"}); 
		
		ex_mod.push({te:"D3.3: Compartmental observations", mod:"SE1.3", link:true, mod_type:"IBM"}); 
		
		ex_mod.push({te:"D3.4: Disease diagnostic test results", mod:"SE1.3", link:true, mod_type:"IBM"}); 
		
		ex_mod.push({te:"D3.5: A partially observed transition", mod:"SE1.3", link:true, mod_type:"IBM"}); 
		
		ex_mod.push({te:"D3.6: A transition observed over a time window", mod:"SE1.3", link:true, mod_type:"IBM"}); 
		
		ex_mod.push({te:"D3.7: A transition observed in a demographic category", mod:"CE3.3", link:true, mod_type:"IBM"}); 
		
		ex_mod.push({te:"D3.8: Uncertain compartmental observations", mod:"SE1.3", link:true,  mod_type:"IBM"}); 	

		let sec = { name:"data-scen", desc:"Data scenarios", pic:"data", shade:true, ex_mod:ex_mod, help:data_scen_text+mod_mess};
		
		inter.section.push(sec);
	}
	
	{
		let ex_mod=[];
		
		ex_mod.push({te:"1) Time variation",help:"These examples look at different ways in which time variation can be added into the model"+mod_mess, gap:gap});

		ex_mod.push({te:"F1.1: Time variation in transmission rate", mod:"SE1.1", link:true, mod_type:"POP"}); 
	
		ex_mod.push({te:"F1.2: Time variation in transmission rate using a trigonometric function", mod:"SE1.1", link:true, mod_type:"POP"}); 
	
		ex_mod.push({te:"F1.3: Time variation in transmission rate through a covariate", mod:"SE1.1", link:true, mod_type:"POP"}); 
		
		ex_mod.push({te:"F1.4: Time variation in population-level transition observation probability", mod:"SE1.1", link:true, mod_type:"POP"}); 

		ex_mod.push({te:"F1.5: Time variation in individual transition observation probability", mod:"SE1.2", link:true, mod_type:"IBM"}); 
			
		ex_mod.push({te:"F1.6: Time variation in population observation probability", link:true,  mod_type:"POP"}); 
		
		ex_mod.push({te:"F1.7: Time variation in branching probability", mod:"CE1.1", link:true,  mod_type:"POP"}); 
		
		ex_mod.push({te:"F1.8: Time-varying covariate affecting branching probability", mod:"CE1.2", link:true, mod_type:"IBM"}); 
	
		ex_mod.push({te:"F1.9: Time-varying disease intervention measure", mod:"SE1.3", link:true, mod_type:"POP"}); 
	
		ex_mod.push({te:"2) Individual-based variation",help:"These examples look at how individual-based variation can be incorporated into the model"+mod_mess, gap:gap});
		ex_mod.push({te:"F2.1: Individual fixed effect applied to a transition", mod:"SE1.2", link:true, mod_type:"IBM"}); 
		ex_mod.push({te:"F2.2: Individual fixed effect applied to a population", link:true, mod_type:"IBM"}); 
		ex_mod.push({te:"F2.3: Individual effect applied to a transition", mod:"EQG1.1", link:true, mod_type:"IBM"}); 
		ex_mod.push({te:"F2.4: Correlated individual effect applied to a transition", mod:"EQG1.2", link:true, mod_type:"IBM"});
		ex_mod.push({te:"F2.5: Correlated individual effect applied to a population", mod:"EQG1.4", link:true, mod_type:"IBM"});
		ex_mod.push({te:"F2.6: Individual fixed effect applied to a branching probability", mod:"CE1.3", link:true, mod_type:"IBM"}); 
		
		ex_mod.push({te:"F2.7: Individual effect applied to a branching probability", mod:"EQG1.5", link:true, mod_type:"IBM"}); 
		ex_mod.push({te:"F2.8: Correlated individual effect applied to a transition with pedigree", mod:"EQG1.2", link:true, mod_type:"IBM"});
	
		ex_mod.push({te:"3) Parameter definitions",help:"These examples look at how parameters can be reparameterised, set to distributions or derived from other quantities"+mod_mess, gap:gap});
		ex_mod.push({te:"F3.1: Reparameterisation", mod:"SPE2.1", link:true, mod_type:"IBM"}); 
		ex_mod.push({te:"F3.2: Parameter distribution", mod:"EE1.2", link:true, mod_type:"IBM"}); 
		ex_mod.push({te:"F3.3: Derived quantities", mod:"SE1.3", link:true, mod_type:"IBM"}); 
		ex_mod.push({te:"F3.4: Factor", mod:"CE3.1", link:true, mod_type:"PBP"}); 
		ex_mod.push({te:"F3.5: Spline reparameterisation", mod:"SPE2.1", link:true, mod_type:"IBM"}); 
		
		let sec = { name:"features", desc:"BICI features", pic:"features", shade:true, ex_mod:ex_mod, help:features_text+mod_mess};
		
		inter.section.push(sec);
	}
	
	
	if(false){
		let ex_mod=[];
		
		if(false){
			ex_mod.push({te:"M) Science Review",help:"Illustrates some simple examples for the Science Review."});
			ex_mod.push({mod:"SR1", te:"SIR model (IBM)",file:"EX_SR1",mod_type:"IBM", pic:"SIR"});
			
			ex_mod.push({mod:"SR1", te:"SIR model (IBM)",file:"EX_B1",mod_type:"IBM", pic:"SIR"});
			
			ex_mod.push({mod:"SR2", te:"COVID-19 model",file:"EX_SR2",mod_type:"POP", pic:"covid"});
		
			ex_mod.push({mod:"SR3", te:"Metapopulation model using geographical regions",file:"EX_SR3",mod_type:"POP", pic:"scotland"});
		}
	
	
		//ex_mod.push({mod:"M5.2", te:"Multiple contact groups, investigating infectivity",file:"EX_E1",mod_type:"IBM", pic:"group"});
		/*
		ex_mod.push({mod:"CE3.1", te:"SIR model (PBM) with demographic stratification",file:"temp2",mod_type:"POP", pic:"SIR_sex2"});	
		*/
		
		ex_mod.push({te:"B1: Time series population observations", mod:"SE1.1", link:true, mod_type:"POP"});
	
		ex_mod.push({mod:"M", te:"Temp",file:"temp",mod_type:"IBM", pic:"SI"});

		ex_mod.push({mod:"M", te:"Test",file:"test",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"M", te:"Geno model",file:"geno",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"M", te:"Geno model with cattle/badgers",file:"geno_bad_cat",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"M", te:"Geno farm",file:"geno_farm",mod_type:"IBM", pic:"SI"});
		
		ex_mod.push({mod:"SE1.1", te:"SI population-based model (PBM)",file:"EX_M1-1",mod_type:"POP", pic:"SI"});
		ex_mod.push({mod:"SE1.2", te:"SI individual-based model (IBM)",file:"EX_M1-2",mod_type:"IBM", pic:"SI"});
		ex_mod.push({mod:"CE3.1", te:"SIRD model (PBM) with demographic stratification",file:"EX_M2-6",mod_type:"POP", pic:"SIR_sex"});
		ex_mod.push({mod:"CE3.3", te:"SIR model (IBM) with demographic stratification",file:"EX_M2-8",mod_type:"POP", pic:"SIR_sex2"});
		ex_mod.push({mod:"SE1.3", te:"SIR model (IBM)",file:"EX_M1-3",mod_type:"IBM", pic:"SIR"});
	
		let sec = { name:"test", desc:"Testing", pic:"data", shade:true, ex_mod:ex_mod};
		
		inter.section.push(sec);
	}
	
	
	// Puts in folder
	for(let i = 0; i < inter.section.length; i++){
		let sec = inter.section[i];
		for(let j = 0; j < sec.ex_mod.length; j++){
			let ex = sec.ex_mod[j];
	
			ex.file = sec.name+"/"+ex.file;
		}
	}
}


/// Adds clickable sections to the home screen
function add_sections_buts(lay)
{
	let cx = 0, cy = 0;

	//let dx = 11, dy = 13.5, gap = 1.2;
	let dx = 11, dy = 10.5, gap = 1.2;

	for(let i = 0; i < inter.section.length; i++){
		let sec = inter.section[i];
		
		let si = 0.9, lh = 1.1, ma = 0.5;  
	
		let text_anno = text_convert_annotation(sec.desc,si,lh,dx-2*ma,"center",BLACK);
			
		let word = text_anno.word;
			
		for(let i = 0; i < word.length; i++){
			word[i].x += ma, word[i].y += dy-0.5*text_anno.height-1.5;
		}
		
		let pic; if(sec.pic != undefined) pic = find_pic(sec.pic);
			
		lay.add_button({ word:word, te:sec.name, desc:sec.desc, pic:pic, x:cx, y:cy, dx:dx, dy:dy, ac:"Section", shade:sec.shade, i:i, type:"Section"});
		
		cx += dx+gap;
		if(cx+dx > lay.dx){ cx = 0; cy += dy+gap;}
	}
}


/// Adds buttons for example section
function add_examples_buts(lay) 
{
	let sec = inter.section[inter.section_sel];
	let ex_mod = sec.ex_mod;
	
	let gap = 0.3;
		
	// Adds file names
	for(let i = 0; i < ex_mod.length; i++){
		let em = ex_mod[i];
		if(em.link){
			let spl = em.te.split(":");
			if(spl.length != 2) error("Does not split in two");
			let fi = spl[0].replace(".","-");
			em.file = sec.name+"/EX_"+fi;
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
			for(let se = 0; se <= inter.section_sel; se++){
				let ex_mod2 = inter.section[se].ex_mod;
				for(let k = 0; k < ex_mod2.length; k++){
					if(ex_mod2[k].mod == ex.mod){ mod_but = ex_mod2[k]; break;}
				}
				if(mod_but != undefined) break;
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
				if(cy != 0.5) cy += 1;		
				
				if(cx != 0){ cx = 0; cy += dy+gapy;}
				
				let te = ex.te;
				let si = 1;
				let font = get_font(si,"bold");
				let w = text_width(te,font)+0.3;
				
				let x = 0;
				
				lay.add_button({te:te, x:x, y:cy-0.2, dx:w, dy:si, type:"LeftText", font:font}); 
			
				lay.add_help_button(x+w,cy+si-0.2+0.2,{title:te, te:ex.help});
				
				cy += 1;
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
