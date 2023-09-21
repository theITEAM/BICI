"use strict";

/// Adss buttons for home page
function add_home_page_buts(lay) 
{
	let cx = corner.x;
	let cy = corner.y;
	
	cy = lay.add_title("About",cx,cy);
	
	let dx = 45;
	
	lay.add_button({x:cx+dx-1.3, y:cy+0.5, dx:4.5, dy:4.5, ac:"NewModel", type:"NewModel"});
	
	//lay.add_paragraph("A description of how to use this software is provided in the attached manual. Further details are given in a paper hosted on bioRxiv.",18,cx+dx+5,cy,BLACK,para_si,para_lh);
	
	//lay.add_button({x:lay.dx-6, y:cy+0.5, dx:3.5, dy:3.5, ac:"Pdf", type:"Pdf"});
	//lay.add_button({x:lay.dx-6, y:cy+4.5, dx:3.5, dy:3.5, ac:"Pdf2", type:"Pdf2"});

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
	let example=[];
	
	/*
	{
		let ex=[];
		ex.push({te:"Simple SIR model",file:"SIR"});
		ex.push({te:"COVID-19 model",file:"covid"});
		ex.push({te:"Marek's disease model",file:"marek"});
		ex.push({te:"Farm",file:"farm"});
		ex.push({te:"Badger and Cattle",file:"badger"});
	
		example.push({model:"Miscellaneous", ex:ex, pic:"SIR", help_title:"The SIR model", help_text:"This is a simple model used to describe epidemic behaviour. Individuals are classified as being either susceptible to infection (S), infected and infectious (I), or recovered / removed / dead (R).\nA susceptible individual has a probability per unit time of becoming infected as a result of other infected individuals. For those individuals that do become infected, a recovery rate determines the duration over which they are infectious."}); 
	}
	*/
	
	{
		let ex=[];
		ex.push({te:"Individual-based model (IBM) with compartmental obsevations.",file:"SIR_ind"});
		ex.push({te:"Population-based model (PBM) with infected population obs.",file:"SIR_pop"});
		ex.push({te:"IBM with infection and recovery time obsevations.",file:"SIR_times"});
		ex.push({te:"IBM with disease diagnostic test obsevations.",file:"SIR_diagtest"});
		ex.push({te:"IBM with susceptibility individual variation from a fixed effect.",file:"SIR_fix"});
		ex.push({te:"IBM with susceptibility individual variation from a random effect.",file:"SIR_rand"});
		ex.push({te:"IBM with susceptibility quantitative genetics model.",file:"SIR_qg"});
		
		/*
		ex.push({te:"EX 1: Complete knowledge of events",file:"EX1"});
		ex.push({te:"EX 2: Known initial state and recoveries",file:"EX2"});
		ex.push({te:"EX 3: Known recoveries only",file:"EX3"});
		ex.push({te:"EX 4: Periodic disease status data",file:"EX4"});
		ex.push({te:"EX 5: Diagnostic test results",file:"EX5"});
		ex.push({te:"EX 6: Future and past predition",file:"EX6"});
		ex.push({te:"EX 7: Time classification",file:"EX7"});
		*/
		
		example.push({model:"SIR model", dx:14, ex:ex, pic:"SIR", help_title:"The SIR model", help_text:"This is a simple model used to describe epidemic behaviour. Individuals are classified as being either susceptible to infection (S), infected and infectious (I), or recovered / removed / dead (R).\nA susceptible individual has a probability per unit time of becoming infected as a result of other infected individuals. For those individuals that do become infected, a recovery rate determines the duration over which they are infectious."}); 
	}

	{
		let ex=[];
		ex.push({te:"IBM with sex stratification of transition rate.",file:"SIR_sex"});
		
		example.push({model:"SIR stratified model", dx:12, ex:ex, pic:"SIRsex", help_title:"The stratified SIR model", help_text:"SIR models can be stratified by individual attributes. Here we consider splitting individuals into males and female which have a different transmission rate."}); 
	}
	
	{
		let ex=[];
		
		ex.push({te:"Non-Markovian latent period",file:"SEIR"});
	
		/*
		ex.push({te:"EX 8: Non-Markovian latent period",file:"EX8"});
		ex.push({te:"EX 9: Disease diagnostic tests",file:"EX9"});
		ex.push({te:"EX 10: Environmental and diagnostic tests",file:"EX10"});
		ex.push({te:"EX 11: Estimating test Se and Sp",file:"EX11"});
		*/

		example.push({model:"SEIR model", dx:14, ex:ex, pic:"SEIR", help_title:"The SEIR model", help_text:"This is a simple model used to describe epidemic behaviour. Individuals are classified as being either susceptible to infection (S), exposed (which means infected but not infectious) (E), infectious (I), or recovered / removed / dead (R).\nA susceptible individual has a probability per unit time of becoming infected as a result of other infected individuals. For those individuals that do become infected, a latent period determines when they become infectious and a recovery rate determines the duration over which they are infectious."}); 
	}
	
	{
		let ex=[];
		
		ex.push({te:"A disease transition experiment",file:"TransExp"});
		ex.push({te:"Incorporating individual-based variation",file:"TransExpIndVar"});
		ex.push({te:"Transmission through the environment",file:"SI_Env"});
		
		
		/*
		ex.push({te:"EX 12: Known infection times and initial state",file:"EX12"});
		ex.push({te:"EX 13: Diagnostic test results",file:"EX13"});
		ex.push({te:"EX 14: Initial / final disease status",file:"EX14"});
		ex.push({te:"EX 15: Incorporating vaccination status",file:"EX15"});
		ex.push({te:"EX 16: Incorporating group effect",file:"EX16"});
		ex.push({te:"EX 17: SNP effects on susceptiblity / infectivity",file:"EX17"});
		ex.push({te:"EX 18: SIR with non-Markovian recovery",file:"EX18"});
		*/
		example.push({model:"Disease transmission experiment", dx:14, ex:ex, pic:"transexp", help_title:"Disease transmission experiment", help_text:"These are used to discover the rate at which individuals become infected and factors affecting this rate (e.g. the role of fixed effects such as vaccination status). \nTypically populations of initially infected and uninfected individuals are placed in closed 'contact groups'. Data is them collected on the ensuing epidemics (in a variety of ways) to estimate model parameters."}); 
	}
	
	{
		let ex=[];
		ex.push({te:"Population model with carrying capacity",file:"SimplePop"});
		/*
		ex.push({te:"EX 19: Periodic population estimates",file:"EX19"});
		ex.push({te:"EX 20: Death times",file:"EX20"});
		ex.push({te:"EX 21: Captures",file:"EX21"});
		ex.push({te:"EX 22: Age dependent mortality",file:"EX22"});
		ex.push({te:"EX 23: Captures and a fraction of death times",file:"EX23"});
		*/
		example.push({model:"Population with births and deaths", dx:14, ex:ex, pic:"Simple population", help_title:"Population with births and deaths", help_text:"This simple model can be used to capture stochastic variation in the population number in a wildlife setting."}); 
	}
	
	{
		let ex=[];
		ex.push({te:"Population model with carrying capacity",file:"SI_source_sink"});
		/*	
		ex.push({te:"EX 24: Periodic disease status measurements",file:"EX24"});
		ex.push({te:"EX 25: Disease status from captures",file:"EX25"});
		ex.push({te:"EX 26: Disease diagnostic tests from captures",file:"EX26"});
		ex.push({te:"EX 27: Disease induced mortality",file:"EX27"});
		*/
		example.push({model:"SIS model with births and deaths", dx:14, ex:ex, pic:"SIbd", help_title:"SIS model with births and deaths", help_text:"This model not only captures stochastic variation in population numbers in a wildlife setting, but also accounts for changes in individual disease status. As well as the usual infection transition from 'S' to 'I', this model includes the reverse transition to account for the fact that individuals do not gain lifelong immunity (which is true for some diseases, such as influenza)."}); 
	}
	
	/*
	{
		let ex=[];
		
		ex.push({te:"Dynamic from the stochastic Lotka-Volterra model",file:"PredPrey"});
	
		example.push({model:"Predator-prey model", dx:14, ex:ex, pic:"LV", help_title:"Predator-prey model", help_text:"These stochastic models predict the dynamic variation in predator and prey populations. Typically the food supply for prey (e.g. rabbits) is plentiful. Rather than being limited by these local resources, the number of prey is regulated by being consumed by predators (e.g. foxes). On the other hand, predators rely on prey for food, which leads to unstable dynamic behaviour. Here a simple Lotka-Volterra model is assumed."}); 
	}
	*/
	
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
	