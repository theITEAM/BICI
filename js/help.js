function helpbuts()                                        // Buttons for help page
{
	var wid, hei, title, text;

	wid = 550;
	switch(helptype){
	case 0:   // Importing
		title = "Importing";
		if(fileToLoadlist.length == 0){
			text = "Currently you have no files imported. Click on 'Load' to import the relevant file (in .csv or .txt format).";
		}
		else{
			wid = 570;
			text = "Either use a previously imported table or click on 'Load' to import the relevant file (in .csv or .txt format).";
		}
		break;
		
	case 1:
		title = "Examples";
		text = "This section introduces some simple examples to demonstrate the versatility of BICI applied to a variety of different model and data scenarios.\nThese examples can be altered or experimented on in any way. Reloading from this page will return them to their default settings.";
		break;
	
	case 2:
		title = "The SIR model";
		text = "This is a simple model used to describe epidemic behaviour. Individuals are classified as being either susceptible to infection (S), infected and infectious (I), or recovered / removed / dead (R).\nA susceptible individual has a probability per unit time of becoming infected as a result of other infected individuals. For those individuals that do become infected, a recovery rate determines the duration over which they are infectious.";
		break;
		
	case 3:
		title = "The SEIR model";
		text = "This is a simple model used to describe epidemic behaviour. Individuals are classified as being either susceptible to infection (S), exposed (which means infected but not infectious) (E), infectious (I), or recovered / removed / dead (R).\nA susceptible individual has a probability per unit time of becoming infected as a result of other infected individuals. For those individuals that do become infected, a latent period determines when they become infectious and a recovery rate determines the duration over which they are infectious.";
		break;
		
	case 4:
		title = "Disease transmission experiment";
		text = "These are used to discover the rate at which individuals become infected and factors affecting this rate (e.g. the role of fixed effects such as vaccination status). Typically populations of initially infected and uninfected individuals are placed in closed 'contact groups'. Data is them collected on the ensuing epidemics (in a variety of ways) to estimate model parameters.";
		break;
		
	case 5:
		title = "Population with births and deaths";
		text = "This simple model can be used to capture stochastic variation in the population number in a wildlife setting.";
		break;
		
	case 6:
		title = "SIS model with births and deaths";
		text = "This model not only captures stochastic variation in population numbers in a wildlife setting, but also accounts for changes in individual disease status. As well as the usual infection transition from 'S' to 'I', this model includes the reverse transition to account for the fact that individuals do not gain lifelong immunity (which is true for some diseases, such as influenza).";
		break;
		
	case 7:
		title = "Spatial diffusion model";
		text = "This simple model captures the movement of individuals across a landscape.";
		break;
		
	case 8:
		title = "Predator-prey model";
		text = "These stochastic models predict the dynamic variation in predator and prey populations. Typically the food supply for prey (e.g. rabbits) is plentiful. Rather than being limited by these local resources, the number of prey is regulated by being consumed by predators (e.g. foxes). On the other hand, predators rely on prey for food, which leads to unstable dynamic behaviour. Here a simple Lotka-Volterra model is assumed.";		
		break;
		
	case 9:
		title = "Spatial disease model";
		text = "This model captures not only the movement of individuals across a landscape, but also their disease status (which is represented by a simple SI model).";
		break;
		
	case 10:
		title = "Loading and saving analyses";
		text = "BICI permits users to load and save analyses in a special “.bici” format. This is useful because it conveniently allows description, data and analysis to all be contained in a single file for future reference.\nWhen saving, two options are available: “With results” includes the posterior samples along with the model and data (so that inference does not need to be run again when the file is loaded), and “W/o results” which does not store the results from inference (leading to a much smaller file size, which can, for example, be emailed)."; 
		break;
	
	case 11:
		title = "Saving and exporting";
		text ="BICI permits users to load and save analyses in a special “.bici” format. This is useful because it conveniently allows description, data and analysis to all be contained in a single file for future reference. When saving, two options are available: “With results” includes the posterior samples along with the model and data (so that inference does not need to be run again when the file is loaded), and “W/o results” which does not store the results from inference (leading to a much smaller file size, which can, for example, be emailed).\nVarious outputs can be exported from BICI including graphs (such as trace, scatter and population plots), posterior samples for parameters and events (in text format) and MCMC diagnostic information."; 
		break;
		
	case 12:
		title = "Importing";
		text = "Constructing very complicated models using the point and click interface can be time consuming and cumbersome.\nIn these instances BICI allows the model (and/or data) to be imported as a simple script in a “.txt” text file.\nThis script has a special syntax (see manual for details). For example “transition from='S' to='I' type='exponential' rate='[β]{I}'” would add an exponential transition between the compartments S and I.";
		break;
	
	case 13:
		title = "Location-based spatial disease model";
		text = "This model imagines that the 'individuals' themselves are fixed, e.g. they could represent farms. Each farm is classified as being either susceptible to infection (S), infected (I), or recovered (R). Disease is spread spatially by means of a transition kernel.";
		break;
	
	case 14:
		title = "Spatial disease model with births and deaths";
		text = "This model assumes spatially distributed groups of individuals (e.g. badger social groups). These individuals are born, die, can move from location to location and are classified by disease status (through an SI model).";
		break;
		
	case 20:
		title = "Description";
		text = "BICI allows users to provide a brief description of the model and any analysis performed. This is not only useful to keep track for personal use, but also makes it easier and more transparent for others to understand what has been done.\nThe description can be edited by clicking on the 'Edit' button (note, bullet points are automatically generated for each carriage return in the editable text box).";
		break;
		
	case 21:
		title = "Priors";
		text = "Priors are specified for each of the model parameters.\nBICI supports the following prior specifications: flat, which relates to a uniform probability distribution across a range, and the gamma, normal, log-normal, exponential and Weibull distributions, as well as the possibility to fix parameters to specific known values.";
		break;
	
	case 22:
		title = "Data sources";
		text = "Rather than loading all data at once, BICI works by loading different sources of data one at a time (in any order). This is done by clicking on the data options at the bottom of the page. Data can be deleted by clicking on the corresponding red cross."; 
		break;
	
	case 23:
		title ="Setup";
		text = "This page provides several options which must be selected before inference can start. By clicking “advanced options” other possibilities can be selected.";
		break;
	
	case 24:
		title = "Advanced option";
		text = "The maximum number of stored parameter and event sequence samples is defined (when exceeded, BICI thins the existing samples by a factor of two and subsequently stores them at half the rate). The maximum number of individuals is also set (by default this is 100000).\nSuch limits are placed to avoid BICI crashing due to excessive memory usage. Note, increasing the sample maximums can lead to smoother posterior plots.\nThree option are available regarding inference termination: either it continues indefinately until manually stopped, a fixed number of MCMC updates are performed, or termination happens when convergence is acheived (as specified when every model parameter has an the effective sample size that exceeds a specified threshold and a Gelman-Rubin statistic that is below a  specified threshold).";
		break;
		
	case 25:
		title = "Set initial parameters";
		text = "Before the simulation can be run it is necessary to specify all of the model parameters.";
		break;
		
	case 26:
		title = "Set initial populations";
		text = "Before a simulation can be run it is necessary to specify the initial populations in each of the compartments. This can either be achieved by manually setting the population sizes, or by loading up a set of initial individuals.\nIn the former case, the number of individuals in a defined classification is given and the other classifications are randomly allocated based on specified percentages.\nIn the latter case, the file must consist of lines which contain the individual ID followed by a tab separator, and then the initial compartment (multiple classifications are tab separated). An example would be:\nInd. 1	  S   Male\nInd. 2	  I   Female\n...\n";
		break;
		
	case 28:
		title = "Generating data";
		text = "Hypothetical data can be generated based on simulated results obtained from the model. These datasets can then be used as inputs when performing inference. This not only provides a good test to show that inference is able to recover model parameters, but also is an invaluable tool in estimating the power of a given experimental setup.";
		break;
	
	case 29:
		title = "Distributions";
		text = "Rather than allowing parameters to take any value they can be assumed to be drawn from a distribution. In the case when this distribution is normal this is analogous to random effects in a mixed model.";
		break;	
	
	case 30:
		title = "Age classifications";
		text = "Different age classifications allow for model parameters to change as a function of age (e.g. older individuals might be more prone to infectious than younger individuals).\nTo incorporate this first add the ages at which transitions in classification occur and then modify the relevant model parameters (e.g. [ β_Age ]).";
		break;
		
	case 31:
		title = "Temporal variation";
		text = "Temporal variation in parameter values is incorporated by splitting up the timeline into discrete periods (in a similar way to any other classification). For example two time periods might relate to before and after a disease intervention strategy, or yearly intervals can account for variations in weather conditions.\nTo incorporate this first add the time points at which changes in time period occur and then modify the relevant model parameters (e.g. [ β_Time ]).";
		break;
	
	case 32:
		title = "Derived quantities";
		text = "These are quantities that are functionally dependant on the model parameters and/or populations (potentially through additional parameters). This is sometimes useful when the things actually being measured relate to these derived quantities rather than directly to anything in the model itself.";
		break;
	
	case 33:
		title = "Add a classification";
		text = "A classification is a discrete means of differentiating individuals. Typical classifications include disease status, location or sex. To add a new classification click on this button and give it a short but meaningful name (such as DS for disease status).";
		break;
		
	case 34:
		title = "Add derived quantity";
		text = "Click here to add a derived quantity. These are quantities that are functionally dependant on the model parameters and populations (potentially through additional parameters). This is sometimes useful when the things actually being measured relate to these derived quantities rather than directly to anything in the model itself.";
		break;
	
	case 35:
		title = "Add compartment";
		text = "Click here to add a new compartment. Compartments are possible values a classification can take. For example the classification 'Sex' would contain two compartments: 'Male' and 'Female'.";
		break;
	
	case 36:
		title = "Add transition";
		text = "Click here to add a new transition. Transitions are events in which an individual moves from one compartment to another. For example when an individual becomes infected it would move from a susceptible 'S' compartment to an infected 'I' compartment. \nDifferent distributions (specifically exponential, gamma and Weibull) can be used to model when the transitions actually happen. Furthermore, 'source' transitions allow for individuals to enter the model and sink transitions allow them to leave (e.g. these can be used to incorporate births and deaths).";
		break;
		
	case 37:
		title = "Add time transition";
		text = "Click here to add a new time transition. Time transitions divide up the timeline into different periods or classifications. This is useful because it allows for model parameters to change in time.";
		break;
		
	case 38:
		title = "Add age transition";
		text = "Click here to add a new age transition. Age transitions divide up the lifespan of an individual into different periods or classifications. This is useful because it allows for model parameters to change depending on how old an individual is.";
		break;
		
	case 39:
		title = "Markovian transition";
		text = "This transition has a certain probability per unit time (or rate) of occurring. The rate is expressed by a user defined equation which can contain model parameters, populations and be dependant on other classifications the individual may have.";
		break;
		
	case 40:
		title = "Gamma distributed transition";
		text = "If 'i' represents the initial compartment and 'f' represents the final compartment then the duration in time the individual spends in 'i' before moving to 'f' is taken to be gamma distributed. This distribution has a mean and shape parameter (which characterises the size of the standard deviation about the mean). Both of these quantities are user defined equations which can contain model parameters, populations and be dependant on other classifications the individual may have.";
		break;	
		
	case 41:
		title = "Weibull distributed transition";
		text = "If 'i' represents the initial compartment and 'f' represents the final compartment then the duration in time the individual spends in 'i' before moving to 'f' is taken to have a Weibull distribution. This distribution is characterised by two parameters, λ and k. Both of these quantities are user defined equations which can contain model parameters, populations and be dependant on other classifications the individual may have.";
		break;	
	
	case 42:
		title = "Source transition";
		text = "A source transition is one which causes individuals to enter the system (e.g. births). This transition is taken to be Markovian and has a certain probability per unit time (or rate) of occurring. The rate is expressed by a user defined equation which can contain model parameters, populations and be dependant on other classifications the individual may have.";
		break;	
		
	case 43:
		title = "Sink transition";
		text = "A sink transition is one which causes individuals to leave the system (e.g. deaths). This transition is taken to be Markovian and has a certain probability per unit time (or rate) of occurring. The rate is expressed by a user defined equation which can contain model parameters, populations and be dependant on other classifications the individual may have.";
		break;	
		
	case 44:
		title = "Trace plot";
		text = "MCMC works by successively drawing parameter samples (represented by the x-axis) from the posterior. Ideally these samples should be randomly distributed, but in reality they are correlated (which manifests itself by structure within these plots). Note, when the number of parameter samples exceeds the maximum value (as specified on the start page), samples are thinned by a factor of two and subsequently gathered at half the rate (this is implemented to ensure that computational memory is not exhausted). The vertical dashed red line represents the so-called burn-in period (before which samples are discarded). This period is dynamically shifted as more posterior sample are generated.\nMultiple variables can be selected by holding the control key.";
		break;
	
	case 45:
		title = "Probability distribution";
		text = "This graph shows the posterior probability distribution for a given model parameter (generated from the posterior samples using kernel density estimation with smoothing parameter adjusted by means of the slider in the bottom right hand corner).\nMultiple parameters can be selected by holding down the control key.";
		break;
		
	case 46:
		title = "Scatter plot";
		text = "Scatter plots enable the user to display the posterior samples of one variable against another. This is achieved by means of clicking the x-axis selecting the relevant variable and then doing the same for the y-axis. Scatter plots are a useful tool to investigate confounding between different model parameters.";
		break;
		
	case 47: 
		title = "Statistics";
		text = "BICI summarises the posterior probability distributions (specifically the mean and 95% credible intervals) for all the model parameters. The credible intervals for SNP or fixed effects are of particular importance because whether they bound zero can be used to establish if they are statistically significant or not.";
		break;
		
	case 48: 
		title = "Simulation statistics";
		text = "BICI summarises the parameter values used to run the simulations as well as other statistics.";
		break;
		
	case 49:
		title = "Correlation";
		text = "This page allows posterior correlations between different model parameters to be calculated. First the parameters of interest are selected by clicking the check-boxes on the right hand menu. The Pearson correlation matrix is then automatically calculated. Values near to one indicate a high degree of correlation (i.e. when one parameter is high then most likely the other one is also high). Values near to minus one indicate anti-corrleation (i.e. when one parameter is high then most likely the other one is low). Values near zero represent little correlation. Clicking one of the elements in the matrix displays a KDE or scatter plot for the corresponding variables.";
		break;
		
	case 50:
		title = "Simulated individual data";
		text = "This shows simulated transitions for each of the individuals in each of the classifications (colours, which can be changed, refer to the those specified in the compartmental model).";
		break;
		
	case 51:
		title = "Individual data";
		text = "This shows timelines for each of the classifications summarising data for each of the individuals in the population. Different type of symbol are used to represent different types of data:\nSquares with different colours in each half represent transitions.\nDashed lines with right and left arrows indicate the time range over which transitions are observed.\nCircles with a single colour represent precise state data\nCircles with multiple colours represent imprecise state data (such as a disease diagnostic test results).";
		break;
		
	case 52:
		title = "Individual posterior distribution";
		text = "This shows timelines for each of the classifications summarising the posterior distribution for each of the individuals in the population.\nThe state of individuals is indicated by the colour bar, which refers to the colour scheme used in the model.\nOn top of the posterior distribution is superimposed the data. Different type of symbol are used to represent different types of data:\nSquares with different colours in each half represent transitions.\nDashed lines with left and right and left arrows indicate the time range over which transitions are observed.\nCircles with a single colour represent precise state data\nCircles with multiple colours represent imprecise state data (such as a disease diagnostic test results).";
		break;
		
	case 60:
		title = "Population plots";
		text = "This shows how populations dynamically change as a function of time. This can either be viewed as a graph (where lines represent posterior means and the shaded areas give 95% credible intervals) or through an animation across the model itself (which shows population variation in greyscale).\nThe results can be filtered by classification, MCMC run, or sample number.";
		break;
		
	case 61:
		title = "Posterior predictive population plots";
		text = "This shows how populations dynamically change as a function of time when the model is simulated using posterior samples for the parameters and initial conditions (generating a posterior predictive check). This can either be viewed as a graph (where lines represent average behaviour and the shaded areas denote regions which contain 95% of the simulations) or through an animation across the model itself (which shows population variation in greyscale).\nThe results can be filtered by classification, MCMC run, or sample number.";
		break;
		
	case 62:
		title = "Simulated population plots";
		text = "This shows how populations dynamically change as a function of time when the model is simulated. This can either be viewed as a graph (when multiple simulations are performed the shaded areas denote regions which contain 95% of the simulations) or through an animation  across the model itself (which shows population variation in greyscale).\nThe results can be filtered by classification, MCMC run, or sample number.";
		break;
		
	case 63:
		title = "View";
		text = "Dynamic variation in populations can either be viewed as a graph or through an animation across the model itself (population variation shown in greyscale).";
		break;
		
	case 64:
		title = "Selecting variable type";
		text = "This drop-down menu allows the user to select different types of variable:\nFirst there is a list of all the classification names which contain parameters governing the transitions within these classifications.\n“Init. Prob.” gives parameters associated with the probability of an individual's initial state.\n“Trans.” gives the number of transitions that occur.\n“Hyper.” gives information about hyper-parameters.\n“Misc.” gives other quantities, e.g. likelihoods and prior.";
		break;
		
	case 65:
		title = "Select classification";
		text = "Select which classification in the model you would like to view.";
		break;
		
	case 67:
		title = "Select classification";
		text = "This selects which classification the populations are actually defined in. For the other classifications percentages are set and these define the probability an individual is in the corresponding states (which are randomly selected at the time of simulation).";
		break;
		
	case 68:
		title = "Style";
		text = "Selects whether the colour scheme from the model or greyscale is used to represent the population sizes.";
		break;
		
	case 69:
		title = "Label";
		text = "Determines if the compartmental name is used as a label or not.";
		break;
		
	case 70:
		title = "Individual view";
		text = "When individuals are being looked at this drop-down menu determines whether they are viewed on a timeline or transitions are superimposed on top of the model itself.";
		break;
		
	case 71:
		title = "Individual";
		text = "Select the individual which is to be viewed.";
		break;
		
	case 72:
		title = "Select classification filter";
		text = "Select which classification to filer by.";
		break;
		
	case 73:
		title = "Select simulation";
		text = "This drop-down menu selects results to be shown for a particular simulation or for all of them simultaneously.";
		break;
		
	case 74:
		title = "Select sample";
		text = "This drop-down menu selects results to be shown for a particular posterior sample or for all of them simultaneously.";
		break;
		
	case 75:
		title = "Label";
		text = "Determines if the compartmental name is used as a label, the population number is used as a label or no label is used.";
		break;
		
	case 76:
		title = "Bayes Factor";
		text = "A Bayes factor (BF) is the ratio of the likelihood of one particular hypothesis to the likelihood of another. The BF comparing the full model to one in which a particular parameter is fixed (usually to zero) can be calculated using this button. This is one way to determine which parameters in the model are redundant, so allowing the model to be simplified. A BF between 3 and 10 represent moderate evidence for one hypothesis over another and exceeding 10 is considered strong evidence.";
		break;
	
	case 77:
		title = "Bayes Factor";
		text = BFtext;
		break;
		
	case 78:
		title = "Select transition";
		text = "View the posterior distribution in the duration of a selected transition.";
		break;
		
	case 79:
		title = "Select individual";
		text = "Select which individuals to look at when getting transition durations.";
		break;
		
	case 80:
		title = "Style";
		text = "Selects whether the colour scheme from the model or greyscale is used to represent the population sizes. The 'rescale' option rescales the palette for each frame to help visualise variation.";
		break;
		
	case 81:
		title = "Posterior distribution";
		text = "This graph shows the posterior probability distribution for two selected model parameters (generated from the posterior samples using kernel density estimation with smoothing parameter adjusted by means of the slider in the bottom right hand corner).";
		break;
		
	case 82:
		title = "Maximum number of points";
		text = "This value sets the maximum number of samples to be plotted. If the number of posterior samples is less than this value, they are all plotted. If above, posterior samples are randomly selected (without replacement).";
		break;
		
	case 90:
		title = "State data";
		text = "Provide information regarding which state an individual is in at particular points in time. This data can either directly give the state (e.g. male or female) or associate certain probabilities to various possibilities (e.g. as is done in a disease diagnostic test with a sensitivity and specificity). The latter is incorporated using a so-called observation model.";
		break;
		
	case 91:
		title = "Capture data";
		text = "When the population (or a subpopulation) is sampled from this is referred to as a 'capture campaign' (in the ecologicial setting this corresponds to the case when traps are laid and animals are actually captured and measured as a means of studying wildlife groups).\nCapture data provides information about the timings of capture campaigns, the subpopulations being sampled from and the probability of detecting individuals.\nAny data recorded during a capture campaign is incorporated either as 'state' data or simply 'presence' data.\nIn reality a capture campaign may take a period of time to carry out. To be analysed using BICI, however, it is necessary to specify a single time point (i.e. at the middle of the capture campaign), which must also be used in any state data obtained from the capture.";
		break;
		
	case 92:
		title = "Transition data";
		text = "Provide information about the timings of individual transitions. These are events in which the compartment in one of the classifications changes (e.g. a susceptible individual becomes infected).";
		break;
		
	case 93:
		title = "Population data";
		text = "Provide data which estimates population or subpopulation sizes at given points in time. These population estimates have an associated uncertainty with them.";
		break;
		
	case 94:
		title = "Derived data";
		text = "Provide data which estimates any of the derived quantities in the model at given points in time. These estimates have an associated uncertainty with them.";
		break;
	
	case 95:
		title = "Presence data";
		text = "This is simply used to inform BICI that certain individuals are present in the population at given time points. Note, presence data is not required if state data is already available at those time points.";
		break;
		
	case 96:
		title = "Dependency";
		text = "In some instances model parameters are related through distributions with hyperparameters.";
		break;
		
	case 97:
		title = "Move data";
		text = "This is like transition data but here the change in state is not as a results of the model, but externally imposed (i.e. it has no likelihood associated with it). For example, it could be used to incorporate the movement of individuals from one location to another or the result of culling.";
		break;
		
	case 99:
		title = "Define model";
		text = "Before simulation or inference can be performed first the model must be defined. Click 'Done' when complete.";
		break;
		
	case 100:
		title = "An error has occurred!";
		text = errormsg;
		break;
		
	case 101:
		title = "Missing information";
		text = errormsg;
		break;

	case 102:
		title = "Inference has finished!";
		text = "The diagnostics suggest that the degree of convergence has been satisfied.";
		break;
		
	case 103:
		title = "Load example";
		text = "Are you sure you want to discard the currently loaded model?";
		break;
		
	case 104:
		title = "Please note...";
		text = errormsg;
		break;
	
	case 105:
		title = "Selecting run";
		text = "BICI can run multiple MCMC chains in parallel. This drop-down menu selects how these runs should be displayed. They can either be individually selected or the 'All Runs' options shows the results separately, or the 'Combine' option uses information from all runs simultaneously.";
		break;
		
	case 106:
		title = "Selecting run";
		text = "BICI can run multiple MCMC chains in parallel. This drop-down menu selects how these runs should be displayed. They can either be individually selected or the 'All Runs' options shows the results simultaneously.";
		break;
		
	case 107:
		title = "Selecting run";
		text = "BICI can run multiple MCMC chains in parallel. This drop-down menu selects how these runs should be displayed. They can either be individually selected or the 'Combine' option uses information from all runs simultaneously.";
		break;
		
	case 108:
		title = "Error importing file";
		text = errormsg;
		break;
		
	case 109:
		title = "Opps..";
		text = "For some unknown reason BICI crashed. We are very sorry!";
		break;
		
	case 110:
		title = "Stop inference?";
		text = "Inference must be stopped before simulation can begin.";
		break;
		
	case 111:
		title = "There was a problem starting the simulation!";
		text = errormsg;
		break;
		
	case 112:
		title = "Error with table";
		text = "The data table must have at least "+ncoldefmax+" columns."; 
		break;
		
	case 200:
		title = "Cannot add state data";
		text = "The model must have more than one state for state data to be added. If the data indicates that individuals simply exist, this can be added as 'presence' data.";
		break;
		
	case 201:
		title = "Click 'Next' to continue.";
		text = "This allows for the model to be setup.";
		break;
		
	case 202:
		title = "Cannot add transition data";
		text = "The model must have at least one transition for transition data to be added.";
		break;
		
	case 203:
		title = "Cannot add derived data";
		text = "The model does have any derived quantities.";
		break;
		
	case 204:
		title = "Cannot add dependency";
		text = "The model has no parameters.";
		break;
		
	case 205:
		title = "Cannot add to this classification";
		text = "Sink transitions cannot be added to more than one classification.";
		break;
		
	case 206:
		title = "Cannot add to this classification";
		text = "Source transitions cannot be added to more than one classification.";
		break;
		
	case 207:
		title = "Changes to the data";
		text = "The following changes have been made to the data:\n"; 
		for(i = 0; i < points.length; i++) text += "<bul>"+points[i]+"\n";
		break;
		
	case 208:
		title = "Import";
		text = "Model changes and/or data successfully incorporated"; 
		break;
		
	case 209:
		title = "Model specified!";
		text = "The model has successfully been specified, and now either simulation or inference can be performed."; 
		break;
		
	case 210:
		title = "Graph options";
		text = "These options specify the appearance of the graph."; 
		break;
	
	case 211:
		title = "Timelines";
		text = "Select which of the classification timelines are shown."; 
		break;
		
	case 212:
		title = "Edit observation model";
		text = "To make changes to the observation model click on the 'Edit' button."; 
		break;
		
	case 213:
		title = "Edit data table";
		text = "To make changes to the data table click on the 'Edit' button."; 
		break;
		
	case 214:
		title = "Capture ID data";
		text = "When captures are performed it is usually assumed that if an individual is observed at the same time then it is included in the capture. However in some circumstances capture campaigns are performed concurrently (e.g. at different locations). Here 'capture ID' data is used to specify which individuals are observed in which captures.";
		break;
		
	case 215:
		title = "Capture PD data";
		text = "Usually capture campaigns have associated with them a given detection probability (usually given by a parameter that is estimated). 'Capture PD' data allows for this detection probability to be specified seperately for each capture campaign.";
		break;
	
	case 216:
		title = "Edit detection probability";
		text = "To make changes to the detection probability click on the 'Edit' button."; 
		break;
		
	case 217:
		title = "Priors defined by model";
		text = "These are prior relationships defined by the model. These can be editted by going to the model section.";
		break;
	
	case 218:
		title = "Priors for model parameters";
		text = "These are priors on parameters used in the comapartmental model.";
		break;
	
	case 219:
		title = "Priors for observation model parameters";
		text = "These are priors on parameters used in the observation model. This is the part of the model which relates the data measurements made to the underlying system dynamics.";
		break;
	
	case 220:
		title = "Priors on the initial state";
		text = "The compartmental state for those individuals present at the inference start time is assumed to take a Dirichlet distribution. Here the α values set the relative probability of being in a given state. For example, α_S = 2, α_I = 1, α_R = 1 corresponds to an individual being twice as likely to initially be in the S state than in the other two states. By default, all α values are set to one representing a flat prior.";
		break;
		
	case 221:
		title = "Age smoothing";
		text = "For those parameters which depend on age it may be assumed that rather than each value being independent, they are 'smoothed' in some way. Two options can be selected: The first smooths the actual parameter values themselves (with a prior that penalises large differences between consecutive age classifications) and the second smooths the log of the values (appropriate only for positive quantities).";
		break;
	
	case 222:
		title = "Time smoothing";
		text = "For those parameters which depend on time it may be assumed that rather than each value being independent, they are 'smoothed' in some way. Two options can be selected: The first smooths the actual parameter values themselves (with a prior that penalises large differences at consecutive time classifications) and the second smooths the log of the values (appropriate only for positive quantities).";
		break;
		
	case 223:
		title = "Advance options";
		text = "The maximum number of individuals can be set (by default this is 100000). Such a limit is placed to avoid BICI crashing due to excessive memory usage, which typically only becomes a problem when the model specification leads to the number of individuals diverging.\nThe number of simulations can be specified. Usually this would just be one, but multiple simulations can be also be performed (this may be of use if the resulting 'samples' are exported from BICI for use elsewhere).";
		break;
	
	case 224:
		title = "Setup";
		text = "The time range over which the simulation is performed must be specified. By clicking “advanced options” other possibilities can be selected.";
		break;
		
	default:
		title = "Help page";
		text = "Help page could not be found!";
		break;
	}

	addbutton("",menux,0,width-menux,height,DONOTHINGBUT,DONOTHINGBUT,-1,-1);
	
	hei = 90;
	
	var d, div;
	div = text.split("\n");
	
	for(d = 0; d < div.length; d++){
		dd = 0; te = div[d]; if(te.substr(0,5) == "<bul>"){ dd = 20; te = te.substring(5);} 
		alignparagraph(te,wid-35-dd,HELPBUTFONT2);
		hei += nlines*25+10;
	}
	
	if(helptype == 0) hei += fileToLoadlist.length*30;
	
	x = menux+(width-menux)/2-wid/2; y = height/2-hei/2;

	addbutton(div,x,y,wid,hei,HELPBACKBUT,HELPBACKBUT,title,-1);
	
	switch(helptype){
	case 0:
		for(j = 0; j < fileToLoadlist.length; j++){
			te = fileToLoadlist[j].name;
			addbutton(te,x+40,y + nlines*25 + 70 + j*30,textwidth(te,HELPBUTFONT)+20,20,LOADFILEBUT,LOADFILEBUT,j,-1);
			addbutton("",x+40+textwidth(te,HELPBUTFONT)+20,y + nlines*25 + 70 + j*30+4,12,12,DELFILEBUT,DELFILEBUT,j,-1);
		}
		addbutton("Load",x+wid-90,y+hei-40,75,30,IMPBUT,IMPBUT,0,-1);
		break;
	
	case 103:
		addbutton("Cancel",x+wid-90-85,y+hei-40,75,30,IMPBUT,IMPBUT,-1,-1);
		addbutton("OK",x+wid-90,y+hei-40,75,30,LOADEXAC,IMPBUT,1,-1);
		break;	
		
	case 110:
		addbutton("Cancel",x+wid-90-85,y+hei-40,75,30,IMPBUT,IMPBUT,-1,-1);
		addbutton("OK",x+wid-90,y+hei-40,75,30,STOPSTARTAC,IMPBUT,1,-1);
		break;	
		
	case 209:
		addbutton("Simulation",x+wid-120-115,y+hei-40,105,30,SIMINFBUT,SIMINFBUT,0,-1);
		addbutton("Inference",x+wid-120,y+hei-40,105,30,SIMINFBUT,SIMINFBUT,1,-1);
		break;
		
	case 212: case 213: case 216:
		addbutton("Cancel",x+wid-90-85,y+hei-40,75,30,IMPBUT,IMPBUT,-1,-1);
		addbutton("Edit",x+wid-90,y+hei-40,75,30,DONEAC,SIMINFBUT,0,-1);
		break;
		
	default:
		addbutton("OK",x+wid-90,y+hei-40,75,30,IMPBUT,IMPBUT,-1,-1);
		break;
	}
	
	addbutton("",x+wid-40,y+10,25,25,HELPCLOSEBUT,HELPCLOSEBUT,-1,-1);
}

function addfullstop(st)                                  // Adds a fullstop to a string, if necessary
{
	if(st.length > 0){
		last = st.substr(st.length-1,1);
		if(last != "." && last != "!" && last != "?") st += ".";
	}
	return st;
}

function alertimp(st)                                     // Error message for imported file
{
	selectbub = -1;
	errormsg = addfullstop("On line "+(jsto+1)+": "+st);
	helptype = 108;
	load();
	buttoninit();
}

function alertimp2(st)                                    // Error message for imported file
{
	selectbub = -1;
	errormsg = addfullstop(st);
	helptype = 108;
	load();
	buttoninit();
}

function alertp(st)                                       // Error message
{
	errormsg = addfullstop(st);
	helptype = 100;
	buttoninit();
}

function alertp2(st)                                      // Error message
{
	errormsg = addfullstop(st);
	helptype = 101;
	buttoninit();
}

function alertp3(st)                                      // Error message
{
	errormsg = addfullstop(st);
	helptype = 104;
}

function alertsim(st)                                     // Error message for simulation
{
	errormsg = addfullstop(st);
	helptype = 111;
	simres.result = -1; loading = 0;
	buttoninit();
}

