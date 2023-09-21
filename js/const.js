// Lists all the constant values used within the interface

const make_file = 1;
const debug = false;           // Determines if debugger used

const check_clone_on = true;   // Turn on to check if cloning is working

"use strict";

const command_list = ["species","classification","class","set","camera","compartment","comp","transition","trans","source","sink","data-dir","description","desc","label","box","parameter","param","derived","der","init-pop","add-ind","remove-ind","move-ind","init-pop-sim","add-ind-sim","remove-ind-sim","move-ind-sim","init-pop-prior","comp-data","trans-data","source-data","sink-data","test-data","pop-data","pop-trans-data","set-traps-data","ind-trapped-data","genetic-data","simulation","sim","inference","inf","ind-effect","fixed-effect","sim-param","sim-state","inf-param","inf-state"];
//,"sim-const","inf-const"

const sim_alg_list = ["gillespie","tau"];
const inf_alg_list = ["DA-MCMC","MFA","ABC","ABC-SMC","ABC-MBP","ABC-PAS","PMCMC","HMC"];

// Possible types of equation
const eqn_types = [
{name:"Se", mode:"param only", obs_model:true},
{name:"Sp", mode:"param only", obs_model:true},
{name:"trap_prob", mode:"param only", obs_model:true},
{name:"comp_prob",mode:"param only", obs_model:true},
{name:"trans_bp", mode:"param with dep"},
{name:"trans_mean", mode:"all"},
{name:"trans_rate", mode:"all"},
{name:"trans_shape", mode:"all"},
{name:"trans_scale", mode:"all"},
{name:"trans_cv", mode:"all"},
{name:"reparam", mode:"param only"},
{name:"dist", mode:"param only"},
{name:"derive_param", mode:"derive_param"},
{name:"derive_eqn", mode:"derive"},
{name:"test", mode:"param only"}
];

const param_type = ["normal","const","dist","reparam"]; // Different varieties of parameter

const mask_size = 100;           // Size to make comp. map mask

const pre = 4;                   // The precision of statistics

const menu_width = 7.5;          // The width of the menu bar on the left
const right_menu_width = 8.5;    // The width of the menu bar on the right
const right_menu_top = 4;        // Distance from top to right menu
const graph_width = 50;          // The width of graphs
const wright = graph_width-6.4;    // The right edge reference used for plotting parameter definitions
const scrollw = 1.2;             // The width of the scroll bar (used in ploting individuals)
const but_width = 4.5;           // Size for a standard button
const but_height = 1.7;   
const but_width_small = 3.5;     // Size for a standard small button
const but_height_small = 1.4; 
const si_radio = 0.8;            // Font size for radio buttons
const si_title = 1.1;            // Font size for title
const si_bubble_title = 1;       // Font size for bubble title
const si_comp_text_frac = 0.85;  // Fraction of comp height for text size
const si_cogr_text_frac = 0.65;  // Fraction of comp height for text size (in graph mode)
const si_histo_label = 1;        // Font size for labels on histograms

const graph_mar = { right:2, left:2.6, top:2, bottom:3.5}; // The margins for a graph

const si_table = 1;              // The size of the font used in tables
const si_drop = 0.8;             // The font size used for drop down menus
const si_annotation = 1;         // Size of annotation text
const si_toolbut = 0.8;          // Size of text on toolbar buttons
const si_clatab = 0.8;           // Size of font in classification tab
const si_claindextab = 1;        // Size of font in classification tab for index

const si_graph_label = 1.2;      // Size of font for axis labels
const si_graph_param_label = 1.2;// Size of font for axis labels- (when parameter)
const si_graph_tick = 0.9;       // Size of font for ticks
const tick_si = 0.4;             // Size of tick marks
const yaxisgap = 0.5;            // the extra size on the y-axis tick marks

const si_big = 1.4;              // Used for displaying variable names
const si_sub = 0.7;
	
const dy_table = 1.3;            // The gap between lines of a table
const dy_table_param = 1.8;      // The gap between lines for parameter values

const expand_dx = 1;             // Determines the size of the expand textbox icon

const latlng_radius = 0.5;       // The size for geographic points

const import_scale_factor = 1;   // Works out how to scale coorinated given in import file

const outline_marx = 0.2;
const outline_mary = 0.1;
const inputbox_fontsi = 1.1;
const inputbox_linesi = 1.28;
const textbox_fontsi = 1.4;
const textbox_linesi = 1.8;
const textbox_space = 0.05;
const textbox_tab = 1.4;
const input_margin = 0.2;
const example_space = 1;

const warn_si = 0.8, warn_lh = 0.9; // Size and line height for warnings
const bub_si = 0.8, bub_lh = 0.9; // Size and line height for bubble paragraphs
const para_si = 1, para_lh = 1.4; // Size and line height for paragraphs
const para_eq_si = 1.2;

const dropdown_height = 1.4;     // The height for a dropdown menu
const dropdown_opheight = 1.2;   // The height for different options in a dropdown menu
const dropdown_dymax = 20;       // The maximum size for the dropdown menu 

const playbar_mar = 2;           // The margins around the playbar

const index_max = 6;             // The maximum number of indices for a tensor 

const page_char_wid = 60;//55;
const page_char_hei = 35;

const scroll_width = 1; 

const corner = {x:2, y:1.5};     // The position of the title button

const button_margin = {dx:2, dy:0.9}; // Position relative corner for buttons

const grid = 1;                  // The gridsize over which points are discretised

const compartment_height = 3;    // The height of compartments on this grid

// Lists parameter types not needed for simulation
const sim_param_not_needed = ["Se","Sp","trap_prob","comp_prob","derive_param"];
const inf_param_not_needed = ["derive_param"];

const dist_pos = ["exp(rate)","exp(mean)","gamma","erlang","log-normal","weibull"];
const exp_dist_pos = ["exp(rate)","exp(mean)","erlang"];

const prior_pos = ["uniform","exp","normal","gamma","log-normal","beta","bernoulli","fix"];
const prior_pos_all = ["uniform","exp","normal","gamma","log-normal","beta","bernoulli","fix","flat","dirichlet"];
const bp_prior_pos = ["flat","dirichlet","fix"];
const zeroone_prior_pos = ["uniform","beta","fix"];

const data_types = ["Init. Pop.", "Individual", "Population","Capture","Genetic"];

let chnotallowed = "\"£!$%&=;@~#\\?";
let notparam_list = "+-*×/0123456789.{<>}()|\n\r ";  // Sets not a parameter character
let paramend_list = "+-*×/0123456789.{<>}[]$\n\r ";	 // Sets characters which end parameter
	
const convert = [
		{command:"add-ind", type:"Add Ind."},
		{command:"remove-ind", type:"Remove Ind."},
		{command:"move-ind", type:"Move Ind."},
		{command:"fixed-eff", type:"Fixed Eff."},
		{command:"init-pop", type:"Init. Pop."},
		{command:"init-pop-prior", type:"Init. Pop. Prior"},
		{command:"comp-data", type:"Compartment"},
		{command:"trans-data", type:"Transition"},
		{command:"source-data", type:"Source"},
		{command:"sink-data", type:"Sink"},	
		{command:"test-data", type:"Diag. Test"},	
		{command:"pop-data", type:"Population"},	
		{command:"pop-trans-data", type:"Pop. Trans."},	
		{command:"set-traps-data", type:"Set Traps"},	
		{command:"ind-trapped-data", type:"Ind. Trapped"},	
		{command:"genetic-data", type:"Genetic"},	
	];
	
const data_template = [
// Used in initial population / add / move / remove
{type:"Init. Pop. Prior", cols:["cl_all","alpha"]},
{type:"Init. Pop.", cols:["cl_all","pop"]},
{type:"Add Ind.", cols:["ID","t,the time the individuals are added","cl_all"]},
{type:"Remove Ind.", cols:["ID","t,the time the individuals are removed"]},
{type:"Move Ind.", cols:["ID","t,the time the individuals are moved","to"]},

// Used in individual data 
{type:"Compartment", cols:["ID","t,the time the individual's state is measured","cl_prob"]},
{type:"Transition", cols:["ID","ttrans,the time the individuals undergo the transition"]},
{type:"Trans. in time range", cols:["ID","t,the time the individuals undergo the transition","start","end"]},
{type:"Source", cols:["ID","t,the time the individuals enter the system"]},
{type:"Source in time range", cols:["ID","t,the time the individuals enter the system","start","end"]},
{type:"Sink", cols:["ID","t,the time the individuals leave the system"]},
{type:"Sink in time range", cols:["ID","t,the time the individuals leave the system","start","end"]},
{type:"Diag. Test", cols:["ID","t,the time the individual's diagnostic test is taken","result"]},
{type:"Population", cols:["t,the time the population measurement is taken","filt_obspop"]},
{type:"Pop. Trans.", cols:["t,the time the individuals undergo the transition","filt_obspoptrans"]},

//Capture
{type:"Set Traps", cols:["setrap_name","t,the time the population measurement is taken","filt_settraps"]},
{type:"Ind. Trapped", cols:["ID","setrap_name"]},

// Genetic
{type:"Genetic", cols:["ID","t,the time the genetic measurement is taken","snp,the first SNP column"]},

// Loading directly into model
{type:"Fixed Effect", cols:["ID","value"]},
{type:"Comp File Pos", cols:["comp_name","comp_x","comp_y"]},
{type:"Comp File Pos Colour", cols:["comp_name","comp_x","comp_y","colour"]},
{type:"Comp File", cols:["comp_name"]},
{type:"Comp File Colour", cols:["comp_name","colour"]},
{type:"Trans exp(rate)", cols:["start_comp","end_comp","rate"]},
{type:"Trans exp(mean)", cols:["start_comp","end_comp","mean"]},
{type:"Trans gamma", cols:["start_comp","end_comp","mean","cv"]},
{type:"Trans log-normal", cols:["start_comp","end_comp","mean","cv"]},
{type:"Trans weibull", cols:["start_comp","end_comp","scale","shape"]},
{type:"Source exp(rate)", cols:["end_comp","rate"]},
{type:"Source exp(mean)", cols:["end_comp","mean"]},
{type:"Sink exp(rate)", cols:["start_comp","rate"]},
{type:"Sink exp(mean)", cols:["start_comp","mean"]},
{type:"Sink gamma", cols:["start_comp","mean","cv"]},
{type:"Sink log-normal", cols:["start_comp","mean","cv"]},
{type:"Sink weibull", cols:["start_comp","scale","shape"]},
{type:"Source Pos exp(rate)", cols:["end_comp","source_x","source_y","rate"]},
{type:"Source Pos exp(mean)", cols:["end_comp","source_x","source_y","mean"]},
{type:"Sink Pos exp(rate)", cols:["start_comp","sink_x","sink_y","rate"]},
{type:"Sink Pos exp(mean)", cols:["start_comp","sink_x","sink_y","mean"]},
{type:"Sink Pos gamma", cols:["start_comp","sink_x","sink_y","mean","cv"]},
{type:"Sink Pos log-normal", cols:["start_comp","sink_x","sink_y","mean","cv"]},
{type:"Sink Pos weibull", cols:["start_comp","sink_x","sink_y","scale","shape"]},

/// Load compartment map
{type:"CompMap", cols:["boundary","comp_name"]},

/// Loading tensor of values for a parameter
{type:"LoadTensor", cols:["dep","value"]},

/// Loading tensor of values for reparameterisation
{type:"LoadReparam", cols:["dep","eqn"]},

/// Loading tensor of values for priors
{type:"LoadPriorSplit", cols:["dep","prior"]},

/// Loading tensor of values for distributions
{type:"LoadDistSplit", cols:["dep","dist"]},

/// Loads the A matrix into and individual effect group
{type:"LoadAmatrix", cols:["A"]},

];
	
// Colours

const BLACK = "#000000", GREEN = "#22ff22", RED = "#ff2222", BLUE = "#4444ff", PURPLE = "#ff44ff";
const WHITE = "#ffffff", GREY = "#bbbbbb", CYAN = "#00bbbb"; 
const LGREEN = "#aaffaa", LRED = "#ffaaaa", LBLUE = "#aaaaff", LPURPLE = "#ffaaff", LGREY = "#dddddd", LCYAN = "#00ffff";
const DGREEN = "#009900", DRED = "#990000", DBLUE = "#000099", DPURPLE = "#aa00aa", DGREY = "#888888", DCYAN = "#009999";
const DDGREEN = "#005500", DDRED = "#550000", DDBLUE = "#000055", DDPURPLE = "#550055", DDGREY = "#444444", DDCYAN = "#005555";
const LLGREEN = "#ddffdd", LLRED = "#ffdddd", LLBLUE = "#ddddff", LLPURPLE = "#ffddff", LLGREY = "#eeeeee", LLCYAN = "#99ffff";
const LLLBLUE = "#eeeeff";
const LLLGREY = "#eeeeee";
const ORANGE = "#ff9900", DORANGE = "#cc6600", LORANGE = "#ffcc55",  DDORANGE = "#441100", LLORANGE = "#ffee99";
const BROWN = "#cb6a00", LBROWN = "#edbb99", DBROWN = "#974500", DDBROWN = "#431200";

const DDDGREEN = "#003300";
const BACKGROUND = "#000000";
const HELP_BLUE = "#444444";
const GREY_BACK = "#444444";
const BLUE_BACK = "#000033";
const MAP_DEFAULT = LGREY;

const BUBBLE_COL = LLLBLUE, BUBBLE_SURROUND = LBLUE, BUBBLE_TEXT_COL =DBLUE;
const INPUT_OUTLINE_COL = LBLUE;

const COMPARTMENT_CURVE = 0.25;                       // Determines the curvature put around compartments

const TRANS_OVER_RANGE = 0.5;                         // The width around which a transition is sensitive to mouse over

const SPLIT_SIZE = 0.5;                               // Size of mouse movement needed to split transition

const TRANS_POINT_R = 0.35;                           // The size of transitions points

const collist = [LGREEN,LRED,LBLUE,LPURPLE,LORANGE,LBROWN,LGREY,GREEN,RED,BLUE,PURPLE,ORANGE,BROWN,GREY,DGREEN,DRED,DBLUE,DPURPLE,DORANGE,DBROWN,DGREY,DDGREEN,DDRED,DDBLUE,DDPURPLE,DDORANGE,DDBROWN,BLACK];

//const auto_color = [GREEN,RED,BLUE,PURPLE,ORANGE,BROWN,GREY,DGREEN,DRED,DBLUE,DPURPLE,DORANGE,DBROWN,DGREY,LGREEN,LRED,LBLUE,LPURPLE,LORANGE,LBROWN,LGREY,DDGREEN,DDRED,DDBLUE,DDPURPLE,DDORANGE,DDBROWN,BLACK];
const auto_color = [BLUE,RED,GREEN,DRED,PURPLE,ORANGE,BROWN,GREY,GREEN,DRED,DBLUE,DPURPLE,DORANGE,DBROWN,DGREY,LGREEN,LRED,LBLUE,LPURPLE,LORANGE,LBROWN,LGREY,DDGREEN,DDRED,DDBLUE,DDPURPLE,DDORANGE,DDBROWN,BLACK];

const OK_butcol = "#444499", OK_butcol2 = "#333399";
const OK_butcol3 = "#6666cc", OK_butcol4 = "#6666cc";

const endl = '\n';

const set_str = "Please set";
const select_str = "Please Select";
const unset_type = "unset_type";

const functi = ["exp","cos","sin","log"];

const opbut = ["+","-","\u00d7","\u2215","^","(",")","Σ","'"];

const numbut = ["0","1","2","3","4","5","6","7","8","9","."];

const char_not_allowed = [":",",",".","!","[","]","{","}","〈","〉","|","*","_","×","-","+","<",">"]; 

const greek = ["\u03B1","\u03B2","\u03B3","\u03B4","\u03B5","\u03B6","\u03B7","\u03B8","\u03B9",
"\u03Ba","\u03Bb","\u03Bc","\u03Bd","\u03Be","\u03Bf","\u03C0","\u03C1","\u03C3","\u03C4","\u03C5","\u03C6","\u03C7","\u03C8","\u03C9"];
["Α","Β","Γ","Δ","Ε","Ζ","Η","Θ","Ι","Κ","Λ","Μ","Ν","Ξ","Ο","Π","Ρ","Σ","Τ","Υ","Φ","Χ","Ψ","Ω"];
const greek_latex = [["alpha","α"],["beta","β"],["gamma","γ"],["Gamma","Γ"],["delta","δ"],["Delta","Δ"],["epsilon","ε"],["zeta","ζ"],["eta","η"],["Eta","Η"],["theta","θ"],["Theta","Θ"],["Iota","ι"],["kappa","κ"],["lambda","λ"],["Lambda","Λ"],["mu","μ"],["nu","ν"],["xi","ξ"],["Xi","Ξ"],["omicron","ο"],["pi","π"],["Pi","Π"],["rho","ρ"],["sigma","σ"],["Sigma","Σ"],["tau","τ"],["upsilon","τ"],["phi","φ"],["Phi","Φ"],["chi","χ"],["psi","ψ"],["Psi","Ψ"],["omega","ω"],["Omega","Ω"],["sum","Σ"]];

const greek_capital = ["A","B","Γ","Δ","E","Z","H","Θ","I","K","Λ","M","N","Ξ","O","Π","P","Σ","T","Y","Φ","X","Ψ","Ω"];

const alphabet = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"];

const like_name = ["L^markov","L^non-markov","L^dist","L^ie","L^obs","Prior"];

// Numeric constants

const UNSET = 99999999;
const LARGE = 10000000;
const TINY = 0.000000001;
const ALMOST_ONE = 0.9999999;

const spawn = require('child_process').spawn;

// Different line thicknesses
const NOLINE = 0, THINLINE = 0.5, NORMLINE = 1, MEDIUMLINE = 1.5, THICKLINE = 2, MTHICKLINE = 3, VTHICKLINE = 4;
