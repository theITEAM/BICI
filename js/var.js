"use strict";
// Lists all the variables used in the graphical interface

let inter=                                         // Dynamic variables related to the interface
{
	layer:[],                                        // Stores output layers
	
	but:[],                                          // Stores information about buttons
	canbut:[],                                       // Stores information about canvas buttons

	ob:[],                                           // Stores objects (which get displayed on canvas)
	
	selectbub:{},                                    // Stores information about the object on which there is a selected sdpeech bubble popup
	
	over:{},                                         // Stores the button the mouse is over
											                  	
	density_slider:{update:"density",lay:"AnimControls", min:-3, max:3, value:0},      // Information about density slider
	compmatrix_slider:{update:"compmatrix",lay:"AnimControls", min:-1, max:1, value:0},// Slider used for matrix
	
	mp4quality_radio: { value:"Low"},                // The default quality for mp4    
	fps: 25,                                         // The frames per second for mp4
	
	mx:undefined,                                    // Position of the mouse
	my:undefined,
	
	mouse_pos: undefined,                            // The current position of the mouse
	
	mouse_time_down:undefined,                       // Time when the mouse button is pressed down
	
	pa:0,                                            // Page number
	page:[],                                         // Stores the page being viewed
	page_name:"",                                    // A string which represents the page name
		
	help:{},                                         // Stores information about help

	equation:{},                                     // Stores information when using equation editor
	
	cursor:{},                                       // Stores information about the cursor

	interval_store:[],                               // Stores any setInterval
	
	sca:undefined,                                   // Scaling factor from window units to pixels

	scroll_position:[],                              // Stores the position of scroll bars
	
	mode:{},                                         // Stores the mode (e.g. dragging)

	textbox_store:[],                                // Stores textbox values
	
	pics:[],                                         // Stores the pictures used in the interface

	arrow_icon:"",                                   // Determines if the grag icon is shown
	
	bubble:{},                                       // Stores information about speech bubble
	
	dropdown:{},                                     // Stores information about a dropdown
	
	tick:[],                                         // Stores potential tick marks
	
	printing:false,                                  // Determines if printing
	
	xi_sq:{},                                        // Stores xi-sq distribution
	
	loading_symbol:{},                               // Stores information about loading symbol
	
	file_type:"",                                    // The type of file to be loaded or saved
	
	imp:{},                                          // Store information about import file
	
	edit_source:{},                                  // Used when editing a data soure
	edit_param:{},                                   // Used when editing values for parameters
	edit_Amatrix:{},                                 // Used when editing values for an A matrix
	edit_Xvector:{},                                 // Used when editing design matrix for fixed effect
	view_graph:{},                                   // Used when viewing a specified graph

	inf_extend:"200%",                               // Default percentage increase when extending MCMC
	
	file_store:{filename:""},                        // Temporarily stores files before saving
	
	export_image:false,                              // Set to true when image is being exported
	
	options:false,                                   // Determines if options on the start inference page are shown
	
	corner_but:false,                                // Determines if there is a corner button
	
	graph:undefined,                                 // Used to plot graphs
	
	comp_select:{list:[]},                           // Stores information about selected compartments
	
	canvas:undefined,                                // Stores the main canvas for drawing
	canvas_cv:undefined,                             // Cv for main canvas
	
	running_status:false,                            // Determines if code is running
	
	child:[],                                        // Stores the c++ process 
	core:[],                                         // Stores information as it is proccessed
	
	figure:[],                                       // Stores letters used on figures           
	
	worker_mess:{active:false}                       // Used to send messages to the worker
}

let model = new Model();

let model_sim;                                     // The model used for simulation
let model_inf;                                     // The model used for inference

let map_store = [];

let data = new Data();

let edit_source={};                                   // Allow for a data source to be editted

let cv;                                            // Used for plotting to working  anvas       

let worker = new Worker('js/worker.js');           // Stores the worker
