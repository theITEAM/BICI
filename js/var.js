// Lists all the variables used in the graphical interface

"use strict";

let inter=                        // Dynamic variables related to the interface
{
	layer:[],
	
	but:[],                         // Stores information about buttons
	canbut:[],                      // Stores information about canvas buttons

	ob:[],                          // Stores objects (which get displayed on canvas)
	
	selectbub:{},                   // Stores information about the object on which there is a selected sdpeech bubble popup
	
	over:{},                        // Stores the button the mouse is over
	
	mx:undefined,                   // Position of the mouse
	my:undefined,
	
	mouse_pos: undefined,           // The current position of the mouse
	
	mouse_time_down:undefined,      // The time when the mouse button is pressed down
	
	pa:0,
	page:[],                        // Stores the page being viewed
	page_name:undefined,            // Gives an string which represents the page name
		
	help:{},                        // Stores information about help

	equation:{},                    // Stores information when using equation editor
	
	cursor:{},                      // Stores information about the cursor

	sca:undefined,                  // The scaling factor to go from window units to pixels

	scroll_position:[],             // Stores the position of scroll bars
	
	mode:{},                        // Stores the mode

	textbox_store:[],               // Stores textbox values
	
	pics:[],                        // Stores all the pictures used in the interface

	arrow_icon:"",                  // Determines if the grag icon is shown
	
	bubble:{},                      // Stores information about speech bubble
	
	dropdown:{},                    // Stores information about a dropdown
	
	tick:[],                        // Stores potential tick marks
	
	loading_symbol:{},              // Stores information about loading symbol
	
	file_type:"",                   // The type of file to be loaded or saved
	
	imp:{},                         // Store information about import file
	
	edit_source:{},                 // Used when editing a data soure
	edit_param:{},                  // Used when editing values for parameters
	edit_Amatrix:{},                // Used when editing values for an A matrix
	edit_Xvector:{},                // Used when editing design matrix for fixed effect
	view_graph:{},                  // Used when viewing a specified graph
	
	export_image:false,             // Set to true when image is being exported
	
	options:false,                  // Determines if options on the start inference page are shown
	
	graph:undefined,                // Used to plot graphs
	
	canvas:undefined,               // Stores the main canvas for drawing
	canvas_cv:undefined,            // Cv for main canvas
	
	running_status:false,           // Determines if code is running
	
	child:[],                       // Stores the c++ process 
	chain:[]                        // Stores information as it is proccessed
}

let model = new Model();

let data = new Data();

let sim_result = {};              // Stores results from simulation
let inf_result = {};              // Stores results from inference

let cv;

let xmlhttp = new XMLHttpRequest(); 

