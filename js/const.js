"use strict";
// Lists all the constant values used within the interface

// 53754 lines of code (16/11/24)
// 50018 lines of code (18/11/24)
// 59127 lines of code (16/15/25)
// 62667 lines of code (17/09/25)

//let ver="windows";         // Determines platform
let ver="linux";
//let ver="mac";

let win_linux = true;                            // When working on win running linux 

const try_on = false;//true;                              // Deterimines if try/catch is on (true)   
const turn_off_random_seed = false;               // Used in testing (false)
let debug = false;                                // Determines if debugger used (false)
let testing = false;                              // Used logo to load up results (false)
let load_map_fast = false;                        // If loads up world map from local computer
let make_one_chain = false;                       // Make into one chain (for diagnostic purposes)

if(win_linux || false){ testing = true; debug = true; load_map_fast = true;}
//if(win_linux){ testing = true; debug = true; load_map_fast = true;}
	
let make_file = false;                            // Determines if makes file or runs
let thick_line = false;                           // Used for making figures
let big_eqn = false;
const graph_dia = false;                          // Used to diagnose problems with graphs
const test_comment = false;                       // Used when checking help comments are correct 
const make_fig = false;                           // Sets editor to make figures for manual
// Here the leters A-H are pressed to add letters to the screen.
// "Home" is used to print the page

const set_default = false;                        // Used to set default when loading (turn off)

const check_clone_on = true;                      // Turn on to check if cloning is working

const command_list = [
{na:"species",ti:0},
{na:"classification",ti:0},
{na:"class",ti:0},
{na:"set",ti:0},
{na:"view",ti:0},
{na:"compartment",ti:0.01},
{na:"comp",ti:0.01},
{na:"comp-all",ti:0.1},
{na:"transition",ti:0.01},
{na:"trans",ti:0.01},
{na:"trans-all",ti:0.1},
{na:"data-dir",ti:0},
{na:"description",ti:0},
{na:"desc",ti:0},
{na:"label",ti:0},
{na:"box",ti:0},
{na:"parameter",ti:1},
{na:"param",ti:1},
{na:"derived",ti:0},
{na:"der",ti:0},
{na:"init-pop-inf",ti:0},
{na:"add-pop-inf",ti:0},
{na:"remove-pop-inf",ti:0},
{na:"add-ind-inf",ti:0},
{na:"remove-ind-inf",ti:0},
{na:"move-ind-inf",ti:0},
{na:"init-pop-sim",ti:0},
{na:"add-pop-sim",ti:0},
{na:"remove-pop-sim",ti:0},
{na:"add-ind-sim",ti:0},
{na:"remove-ind-sim",ti:0},
{na:"move-ind-sim",ti:0},
{na:"add-pop-post-sim",ti:0},
{na:"remove-pop-post-sim",ti:0},
{na:"add-ind-post-sim",ti:0},
{na:"remove-ind-post-sim",ti:0},
{na:"move-ind-post-sim",ti:0},
{na:"comp-data",ti:0},
{na:"trans-data",ti:0},
{na:"test-data",ti:0},
{na:"pop-data",ti:0},
{na:"pop-trans-data",ti:0},
{na:"ind-effect-data",ti:0},
{na:"ind-group-data",ti:0},
{na:"genetic-data",ti:0},
{na:"simulation",ti:0},
{na:"sim",ti:0},
{na:"inference",ti:0},
{na:"inf",ti:0},
{na:"post-sim",ti:0},
{na:"posterior-simulation",ti:0},
{na:"warning-post-sim",ti:0},
{na:"ind-effect",ti:0},
{na:"fixed-effect",ti:0},
{na:"param-sim",ti:0.5},
{na:"state-sim",ti:1},
{na:"warning-sim",ti:0},
{na:"param-inf",ti:5},
{na:"param-stats-inf",ti:1},
{na:"generation-inf",ti:1},
{na:"state-inf",ti:10},
{na:"warning-inf",ti:0},
{na:"trans-diag-inf",ti:1},
{na:"param-post-sim",ti:5},
{na:"state-post-sim",ti:10},
{na:"diagnostics-inf",ti:0.1},
{na:"map",ti:5},
{na:"post-sim",ti:0},
{na:"post-simulation",ti:0},
{na:"param-mult",ti:0}
];

// Lists all commands which need to load files
const data_command_list = ["init-pop-inf", "add-pop-inf", "remove-pop-inf", "add-ind-inf", "remove-ind-inf", "move-ind-inf", "init-pop-sim", "add-pop-sim", "remove-pop-sim","add-ind-sim", "remove-ind-sim", "move-ind-sim","add-pop-post-sim","remove-pop-post-sim","add-ind-post-sim", "remove-ind-post-sim","move-ind-post-sim", "comp-data", "trans-data", "test-data", "pop-data", "pop-trans-data", "ind-effect-data", "ind-group-data", "genetic-data","param-mult"];

const sim_alg_list = ["gillespie","tau"];
const inf_alg_list = ["DA-MCMC","PAS-MCMC","MFA","ABC","ABC-SMC","ABC-MBP","PMCMC","HMC"];

// Possible types of equation
const eqn_types = [
{name:"Se", mode:"param only", obs_model:true},
{name:"Sp", mode:"param only", obs_model:true},
{name:"mut_rate", mode:"param only", obs_model:true},
{name:"seq_var", mode:"param only", obs_model:true},
{name:"comp_prob",mode:"param only", obs_model:true},
{name:"trans_prob",mode:"param only", obs_model:true},
{name:"sim_comp_prob",mode:"param only"},
{name:"trans_bp", mode:"param with dep"},
{name:"trans_mean", mode:"all"},
{name:"trans_rate", mode:"all"},
{name:"trans_shape", mode:"all"},
{name:"trans_scale", mode:"all"},
{name:"trans_cv", mode:"all"},
{name:"reparam", mode:"param only"},
{name:"reparam_eqn", mode:"all"},
{name:"dist", mode:"param only"},
{name:"derive_param", mode:"derive_param"},
{name:"derive_eqn", mode:"derive"},
{name:"test", mode:"param only"}
];

const param_type = ["normal","const","dist","reparam"]; // Different varieties of parameter

const graph_but_not_print = ["Settings","ZoomIn",  // Buttons which are not printed
"ZoomOut","RadioButton","RadioButtonText"];
	
const mask_size = 100;                             // Size to make comp. map mask

const mem_state_sample_max = 1;                    // The memory (GB) allowed for state samples
const mem_param_sample_max = 0.25;                 // The memory (GB) allowed for parameter samples

const bytes_in_GB = 1073741824;                    // Bytes in a GB

const pre = 4;                                     // The precision of statistics

const menu_width = 7.5;                            // The width of the menu bar on the left
let right_menu_width = 8.5;                        // The width of the menu bar on the right
const right_menu_slider = 0.5;                     // The position of right slider
const right_menu_top = 4;                          // Distance from top to right menu
const graph_width = 50;                            // The width of graphs
const wright = graph_width-6.4;                    // The right edge reference used for plotting parameter definitions
const wright2 = graph_width-3.4;                   // The right edge reference used for plotting parameter definitions
const scrollw = 1.2;                               // The width of the scroll bar (used in ploting individuals)
const but_width = 4.5;                             // Width for a standard button
const but_height = 1.7;                            // Height for a standard button
const but_width_small = 3.5;                       // Width for a standard small button
const but_height_small = 1.4;                      // Height for a standard small button
const anim_mar = 2.5;                              // The margins around the time line
const anim_mar_bot = 5.5;                          // Margin at bottom 
const si_radio = 0.8;                              // Font size for radio buttons
const si_title = 1.0;                              // Font size for title
const si_transtree = 0.8;                          // Font size for transition tree labels
const si_bubble_title = 1;                         // Font size for bubble title
const si_comp_text_frac = 0.85;                    // Fraction of comp height for text size
const si_cogr_text_frac = 0.65;                    // Fraction of comp height for text size (in graph mode)
const si_histo_label = 1;                          // Default font size for labels on histograms
const si_sup_fac = 0.6;                            // Factor in font size for superscript
const bar_thick = 0.7;                             // The thickness of the histogram bars
const slider_dx = 0.5;                             // The width of slider button
const default_ruler = 0.362;                       // Sets the default ruler size

const graph_mar = { right:2, left:2.6, top:2, bottom:3.5}; // The margins for a graph

const si_table = 1;                                // The size of the font used in tables
const si_drop = 0.8;                               // The font size used for drop down menus
const si_toolbut = 0.9;                            // Size of text on toolbar buttons
const si_clatab = 0.8;                             // Size of font in classification tab
const si_claindextab = 1;                          // Size of font in classification tab for index
const si_graph_label = 1.2;                        // Size of font for axis labels
const si_graph_param_label = 1.2;                  // Size of font for axis labels- (when parameter)
const si_graph_tick = 0.9;                         // Size of font for ticks
const tick_si = 0.4;                               // Size of tick marks
const yaxisgap = 0.5;                              // the extra size on the y-axis tick marks
const width_table_max = 10;                        // The maximum width for table elements
const si_big = 1.4;                                // Used for displaying variable names
const si_sub = 0.7;                                // Size of subscript
const dy_spline_fac = 0.9;                         // Factor which multiplies line height for splines 
const dy_spline_fac2 = 0.7;                         // Factor which multiplies line height for splines 
const dy_table = 1.3;                              // The gap between lines of a table
const dy_table_param = 1.8;                        // The gap between lines for parameter values
const expand_dx = 1;                               // Determines the size of the expand textbox icon
const latlng_radius = 0.5;                         // The size for geographic points
const import_scale_factor = 1;                     // Works out how to scale coorinated given in import file
const outline_marx = 0.2;                          // Outline x margin for input box 
const outline_mary = 0.1;                          // Outline y margin for input box 
const inputbox_fontsi = 1.1;                       // Font size for input box
const inputbox_linesi = 1.28;                      // The line size for input box
const equation_fontsi = 1.9;                       // Font size for equation editor
const equation_brasi = 2.1;                        // Font size for bracket
const equation_linesi = 2.4;                       // Line size for equation editor
const textbox_fontsi = 1.1;                        // Font size for text box
const textbox_linesi = 1.4;                        // The line size for text box
const textbox_space = 0.05;                        // Size given to a space
const textbox_tab = 1.4;                           // Size of tab in text box
const input_margin = 0.2;                          // Margin in input box
const example_space = 1;                           // Used for examples spacing on home page
const bubblescroll_dymax = 16;                     // The maximum size of slider in bubbles
const warn_si = 0.8, warn_lh = 0.9;                // Size and line height for warnings
const bub_si = 0.8, bub_lh = 0.9;                  // Size and line height for bubble paragraphs
const bubbig_si = 0.9, bubbig_lh = 1.2;            // Size and line height for big bubble paragraphs
const para_si = 1, para_lh = 1.4;                  // Size and line height for paragraphs
const para_eq_si = 1.2;                            // Size of equation in paragraph
const table_lh = 1.4;                              // The line height used for tables 
const dropdown_height = 1.4;                       // The height for a dropdown menu
const dropdown_opheight = 1.2;                     // The height for different options in a dropdown menu
const dropdown_dymax = 17;                         // The maximum size for the dropdown menu 
const si_limit_label = 0.2;                        // Minimum font size for displaying on map
const si_limit_circle = 0.003;                     // Minimum size for displaying lng lat point on map
const playbar_mar = 2;                             // The margins around the playbar
const loading_si = 2.7;                            // The size of the loading button
const index_max = 6;                               // The maximum number of indices for a tensor 
const page_char_wid = 60;                          // The total page width in characters
const page_char_hei = 35;                          // The total page height in characters
const key_gap = 3;                                 // Gap from bottom of page for key
const scroll_width = 1;                            // The width of the scroll bar
const corner = {x:2, y:1.5};                       // The position of the title button
const button_margin = {dx:2, dy:0.9};              // Position relative corner for buttons
const grid = 1;                                    // The grid-size over which points are discretised
const compartment_height = 3;                      // The height of compartments on this grid
const compartment_width_min = 4;                   // The min width of compartments on this grid
const marnew = 0.1;                                // Margin used for new compartments added
const seed_max = 10000;                            // The maximum seed number
const letter_size = 2;                             // Letter size when making figures 

const image_scale_factor = 3;                      // Increase in screen size for print image
const mp4_factor_low = 1;                          // For low quality scale factor
const mp4_factor_medium = 2;                       // For medium quality scale factor 
const mp4_factor_high = 3;                         // For medium quality scale factor 
let print_scale_factor;                            // Increase is screen size for print
let print_line_factor;                             // Increase in line thickness print

// Lists parameter types not needed for simulation
const sim_param_not_needed = ["Se","Sp","mut_rate","seq_var","trans_prob","comp_prob","derive_param"];
const inf_param_not_needed = ["derive_param"];

const dist_pos = ["exp(rate)","exp(mean)","gamma","erlang","log-normal","weibull","period"];
const exp_dist_pos = ["exp(rate)","exp(mean)","erlang"];
const source_dist_pos = ["exp(rate)","exp(mean)"];
const prior_pos = ["inverse","uniform","power","exp","normal","gamma","log-normal","beta","bernoulli","fix"];
const prior_factor_pos = ["mdir"];
const prior_cv_pos = ["mvn-jeffreys","mvn-uniform"];
const prior_pos_positive = ["inverse","uniform","power","exp","gamma","log-normal","fix"];

const data_types = ["Init. Cond.", "Individual", "Population", "Additional"];

let chnotallowed = "\"$&=~#\\→";
let notparam_list = "+-*×/0123456789.{}<>〈〉()|\n\r Σ∫′→;";// Sets not a parameter character
let paramend_list = "+-*×/|.{}<>〈〉[]$\n\r Σ∫";	            // Sets characters which end parameter
let name_notallow = "|\"*×_ {}<>〈〉[]()=Σ∫′→;$&#\\";       // Character not allowed in species, class or comp names
let name_ch_max = 40;                                   // The maximum number of characters allowed for strings
let invalid_name = ["Compartment","Population","Alpha","Distribution","file"];

const convert = [
		{command:"add-pop-inf", type:"Add Pop."},
		{command:"remove-pop-inf", type:"Remove Pop."},
		{command:"add-ind-inf", type:"Add Ind."},
		{command:"remove-ind-inf", type:"Remove Ind."},
		{command:"move-ind-inf", type:"Move Ind."},
		{command:"fixed-eff", type:"Fixed Eff."},
		{command:"init-pop-inf", type:"Init. Pop."},
		{command:"comp-data", type:"Compartment"},
		{command:"trans-data", type:"Transition"},
		{command:"test-data", type:"Diag. Test"},	
		{command:"ind-effect-data", type:"Ind. Eff."},
		{command:"ind-group-data", type:"Ind. Group"},
		{command:"genetic-data", type:"Genetic"},	
		{command:"pop-data", type:"Population"},	
		{command:"pop-trans-data", type:"Pop. Trans."},	
		{command:"param-mult", type:"Parameter Mult."},
	];
	
const data_template = [
// Used in initial population / add / move / remove
{type:"Init. Pop.", is_data:true, title:"Initial population", help:load_initpop_text, cols:["cl_all","init_pop"]},

{type:"Add Pop.", is_data:true, title:"Add population", help:load_add_pop_text, cols:["t,the time the individuals are added","cl_all","add_pop"]},

{type:"Remove Pop.", is_data:true, title:"Remove population", help:load_remove_pop_text, cols:["t,the time the individuals are removed","cl_all","rem_pop"]},

{type:"Add Ind.", is_data:true, title:"Add individuals", help:add_ind_text, cols:["ID","t,the time the individuals are added","cl_all_prob"]},

{type:"Remove Ind.", is_data:true, title:"Remove individuals", help:rem_ind_text, cols:["ID","t,the time the individuals are removed"]},

{type:"Move Ind.", is_data:true, title:"Move individuals", help:move_ind_text, cols:["ID","t,the time the individuals are moved","to"]},

{type:"Ind. Eff.", is_data:true, title:"Individual effect values", help:ind_eff_data_text, cols:["ID","value"]},

{type:"Ind. Group", is_data:true, title:"Individual group", help:ind_group_data_text , cols:["ID"]},

// Used in individual data 
{type:"Compartment", is_data:true, title:"Compartmental data", help:load_compartment_text, cols:["ID","t,the time the individual's state is measured","cl_prob"]},

{type:"Transition", is_data:true, title:"Transition data", help:load_transition_text, cols:["ID","ttrans,the time the individuals undergo the transition"]},

{type:"Trans. in time range", is_data:true, title:"Transition in time range data",  help:load_transition_text, cols:["ID","t,the time the individuals undergo the transition","start","end"]},

{type:"Diag. Test", is_data:true, title:"Diagnostic test data", help:diag_test_data_text, cols:["ID","t,the time the individual's diagnostic test is taken","result"]},

{type:"Population", is_data:true, title:"Population data", help:load_population_text, cols:["t,the time the population measurement is taken","filt_obspop"]},

{type:"Pop. Trans.", is_data:true, title:"Population-level transition data", help:load_poptrans_text, cols:["tstart","tend","filt_obspoptrans"]},

// Genetic
{type:"Genetic", is_data:true, title:"Genetic data", help:seq_data_text, cols:["ID","t,the time the genetic measurement is taken","snp,the first SNP column"]},

// Loading directly into model
{type:"Fixed Effect", title:"Design matrix for fixed effect", help:fixed_eff_text, cols:["ID","value"]},

{type:"Comp File Pos", title:"Load compartments", help:load_comppos_text, cols:["comp_name","comp_x","comp_y"]},

{type:"Comp File Pos Colour", title:"Load compartments", help:load_compposcol_text, cols:["comp_name","comp_x","comp_y","colour"]},

{type:"Comp File", title:"Load compartments", help:load_comp_text, cols:["comp_name"]},

{type:"Comp File Colour", title:"Load compartments", help:load_compcol_text, cols:["comp_name","colour"]},

{type:"Trans File", title:"Load transitions", cols:["start_comp","end_comp","trans_value"]},
{type:"Trans File Pos", title:"Load transitions", cols:["start_comp","end_comp","source_x","source_y","trans_value"]},

// Load compartment map
{type:"CompMap", title:"Load compartment map", help:load_compmap_text, cols:["boundary","comp_name"]},

// Loading tensor of values for a parameter
{type:"LoadTensor", title:"Load tensor", cols:["dep","value"]},

// Loading tensor of values for reparameterisation
{type:"LoadReparam", title:"Load reparameterisation", help:load_reparam_text2, cols:["dep","eqn"]},

// Loading tensor of values for priors
{type:"LoadPriorSplit", title:"Load priors", help:load_priorsplit_text2, cols:["dep","prior"]},

// Loading tensor of values for distributions
{type:"LoadDistSplit", title:"Load distributions", help:load_distsplit_text2, cols:["dep","dist"]},

// Loads the A matrix into and individual effect group
{type:"LoadAmatrix", title:"Load relationship matrix", cols:["A"]},

// Loads the Ainv matrix into and individual effect group
{type:"LoadAinvmatrix", title:"Load relationship matrix", cols:["A"]},

// Loads the pedigree into and individual effect group
{type:"APed", title:"Pedigree", cols:["ID","sire","dam"]},


// Load compartment map
{type:"KnotTimes", title:"Load spline knot times", cols:["t,the knot times"]},
];
	
// Colours

const BLACK = "#000000", GREEN = "#22ff22", RED = "#ff2222", BLUE = "#4444ff", PURPLE = "#ff44ff";
const WHITE = "#ffffff", GREY = "#bbbbbb", CYAN = "#00bbbb", YELLOW = "#dddd00"; 
const LGREEN = "#aaffaa", LRED = "#ffaaaa", LBLUE = "#aaaaff", LPURPLE = "#ffaaff", LGREY = "#dddddd", LCYAN = "#00ffff", LYELLOW = "#ffff00";
const DGREEN = "#00bb00", DRED = "#bb0000", DBLUE = "#0000bb", DPURPLE = "#bb00bb", DGREY = "#888888", DCYAN = "#00bbbb", DYELLOW = "#bbbb00";
const DDGREEN = "#009900", DDRED = "#990000", DDBLUE = "#000099", DDPURPLE = "#990099", DDGREY = "#444444", DDCYAN = "#009999", DDYELLOW = "#999900";
const LLGREEN = "#ddffdd", LLRED = "#ffdddd", LLBLUE = "#ddddff", LLPURPLE = "#ffddff", LLGREY = "#eeeeee", LLCYAN = "#99ffff";
const LLLBLUE = "#eeeeff";
const MLLBLUE = "#bbbbff";
const LLLGREY = "#eeeeee";
const ORANGE = "#ff9900", DORANGE = "#dd6600", LORANGE = "#ffcc55",  DDORANGE = "#882200", LLORANGE = "#ffee99";
const BROWN = "#cb6a00", LBROWN = "#edbb99", DBROWN = "#974500", DDBROWN = "#431200";
const LLGREY_CODE = 99999;

const DDDGREEN = "#003300";
const BACKGROUND = "#000000";
const HELP_BLUE = "#444444";
const GREY_BACK = "#444444";
const BLUE_BACK = "#000033";
const MAP_DEFAULT = LGREY;

const BUBBLE_COL = LLLBLUE, BUBBLE_SURROUND = LBLUE, BUBBLE_TEXT_COL = BLACK;
const EDIT_MARGIN_COL = BLACK;
const BUBBLE_SCROLL_COL = BUBBLE_COL;
const INPUT_OUTLINE_COL = LBLUE;
const POP_COL = BLUE;
const IE_COL = DGREEN;
const FE_COL = DORANGE;
const PARAM_COL = DRED;
const SUM_COL = DPURPLE;

const TI_DIV_MAX = 10000;                          // Maximum number of time divisions
const DISPLAY_FILE_LEN_MAX = 100;
const KEY_LINE_MAX = 100;                          // Maximum number of key lines on a graph    
const GRAPH_LINE_MAX = 1000;                       // Maximum number of lines on graph    
const IND_PLOT_MAX = 20000;                        // Manimum number of individuals plotted
const MATRIX_PLOT_MAX = 40000;                     // The maximim number of nodes plotted
const MATRIX_COMP_MAX = 200;                       // The maximum comparments for MatrixAnim
const HISTO_PLOT_MAX = 100;                         // Maximum number of bars in histogram
const POINT_PLOT_MAX = 1000;                       // Maximum number of points in scatter
const PARAM_STAT_MAX = 1000;                       // Maximum number of parameter statistics
const NODE_PLOT_MAX = 1000;                        // Maximum number nodes is phylogenetic tree
const SCRIPT_LINE_MAX = 1000;                      // The maximum number of script lines to show

const ELE_REDUCE_FAC = 0.7;                        // Factor reduction when too many elements
const TABLE_ROW_MAX = 10000;                       // The maximum number of rows displayed on a table

const DEN_X = 200;                                 // Resolution of density p;ots
const SR_MIN = 5;                                  // Determines if DEN_X must be increased
const COMPARTMENT_CURVE = 0.25;                    // Determines the curvature put around compartments

const NORMAL = 99999989;                           // Denotes dot a source or sink
const SOURCE = 99999990;                           // Denotes source 
const SINK = 99999991;                             // Denotes sink
const OUT = 99999992;                              // Denotes individual out of system
const EV_ENTER = 99999993;                         // Denotes and enter event
const EV_LEAVE = 99999994;                         // Denotes a leave event
const EV_TRANS = 99999995;                         // Denotes a transition event
const EV_MOVE = 99999996;                          // Denotes a move event
const MIDPOINT = 99999997;                         // Denotes a midpoint in a transition line
const OUTSIDE_INF = 99999998;                      // Denotes an outside infection
const TRANS_OVER_RANGE = 0.5;                      // The width around which a transition is sensitive to mouse over
const SIM_VALUE_DASH = 2;                          // Dash value used when displaying sim values 
const SIM_VALUE_COL = BLACK;                       // The colour of the sim value line
const SIM_VALUE_THICK = 1;                         // The thickness of the sim value line
const TRANS_EXP_DASH = 2;                          // Dash for expected transitions
const TRANS_EXP_THICK = 2;                         // Thickness for expected transitions

const SPLIT_SIZE = 0.5;                            // Size of mouse movement needed to split transition
const GENRATION_COL = BLUE;                        // Used for generation plot
const GENRATION_CI_COL = LBLUE;                    // Used for generation plot
const GENRATION_THICK = 1;                         // The thickness of the sim value line
const GENRATION_DASH = 3;                          // Dash value used when displaying gen CIs 
const TRANS_POINT_R = 0.35;                        // The size of transitions points
const COLOUR_KEY_DIV = 100;                        // The numb er
const COMP_FILTER_MAX = 1000;                      // The maximum number of compartments for a filter
const TRANS_FILTER_MAX = 1000;                     // The maximum number of transitions for a filter

const endl = '\n';                                 // End line character
const missing_str = ".";                           // Denotes missing in data file
const not_alive_str = "!";                         // Denotes individul not alive
const set_str = "Please set";                      // String prompt to set
const select_str = "Please Select";                // String prompt to select         
const select_drop_str = "Select";                  // String prompt to set dropdown
const unset_type = "unset_type";                   // Denotes unset type
const no_elements = "No elements";                 // Message if no elements in array
const data_invalid = "Data source is invalid";     // Invalid data source message
const dist_matrix_name = "D";                      // The name of a distance matrix
const iden_matrix_name = "δ";                      // The name of identity matrix
const iden_matrix_name2 = "\\delta";               // The name of identity matrix
const RN_name = "RN";                              // The basic reproduction number function
const RNE_name = "RNE";                            // The effective reproduction number function
const RNC_name = "RNC";                            // The computational reproduction number function
const GT_name = "GT";                              // The basic generation time
const GTE_name = "GTE";                            // The effective generation time function
const GTC_name = "GTC";                            // The computational generation time function

// List of colours used for compartments
const collist = [LGREEN,LRED,LBLUE,LPURPLE,LORANGE,LYELLOW,LBROWN,LGREY,GREEN,RED,BLUE,PURPLE,ORANGE,YELLOW,BROWN,GREY,DGREEN,DRED,DBLUE,DPURPLE,DORANGE,DYELLOW,DBROWN,DGREY,DDGREEN,DDRED,DDBLUE,DDPURPLE,DDORANGE,DDYELLOW,DDBROWN,BLACK];

// List of colours used for lines
const auto_color = [BLUE,RED,GREEN,DRED,PURPLE,ORANGE,BROWN,GREY,GREEN,DRED,DBLUE,DPURPLE,DORANGE,DBROWN,DGREY,LGREEN,LRED,LBLUE,LPURPLE,LORANGE,LBROWN,LGREY,DDGREEN,DDRED,DDBLUE,DDPURPLE,DDORANGE,DDBROWN,BLACK];

// Colours for OK buttons
const OK_butcol = "#444499", OK_butcol2 = "#333399";
const OK_butcol3 = "#6666cc", OK_butcol4 = "#6666cc";


const functi = ["exp","cos","sin","log","pow","max","min","abs","sqrt","step","sig","thresh","ubound"];
const functi_dx = ["def","def","def","def","def","def","def","def",2.8,2.8,"def",3.4,3.7];

const opbut = ["_","^","+","-","\u00d7","\u2215","(",")","Σ","∫","'"];

const spline_radio_pos = [{value:"Linear"},{value:"Square"},{value:"Cubic +ve"},{value:"Cubic"}];

const numbut = ["0","1","2","3","4","5","6","7","8","9","."];

const greek = ["\u03B1","\u03B2","\u03B3","\u03B4","\u03B5","\u03B6","\u03B7","\u03B8","\u03B9",
"\u03Ba","\u03Bb","\u03Bc","\u03Bd","\u03Be","\u03Bf","\u03C0","\u03C1","\u03C3","\u03C4","\u03C5","\u03C6","\u03C7","\u03C8","\u03C9"];

const greek_latex = [["alpha","α"],["beta","β"],["gamma","γ"],["Gamma","Γ"],["delta","δ"],["Delta","Δ"],["epsilon","ε"],["zeta","ζ"],["eta","η"],["Eta","Η"],["theta","θ"],["Theta","Θ"],["iota","ι"],["kappa","κ"],["lambda","λ"],["Lambda","Λ"],["mu","μ"],["nu","ν"],["xi","ξ"],["Xi","Ξ"],["omicron","ο"],["pi","π"],["Pi","Π"],["rho","ρ"],["sigma","σ"],["tau","τ"],["upsilon","υ"],["phi","φ"],["Phi","Φ"],["chi","χ"],["psi","ψ"],["Psi","Ψ"],["omega","ω"],["Omega","Ω"],["sum","Σ"],["int","∫"]];
//["Sigma","Σ"],

const greek_capital = ["A","B","Γ","Δ","E","Z","H","Θ","I","K","Λ","M","N","Ξ","O","Π","P","Σ","T","Y","Φ","X","Ψ","Ω"];

const alphabet = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z"];

const like_name = ["L^markov","L^non-markov","L^ie","L^dist","L^obs","L^genetic-proc","L^genetic-obs","L^init","Prior"];

const trans_tree_name = ["N^origin","N^infected","N^mut-tree","N^mut-origin","N^unobs","t^root"];

// Numeric constants

const UNSET = 99999999;                            // Deotes unset
const LARGE = 10000000;                            // Denote a large quantity
const VLARGE = 100000000000;                       // Denote a very large quantity
const TINY = 0.000000001;                          // Denotes a tiny quantity
const VTINY = 0.00000000000001;                    // A very tiny number
const ALMOST_ONE = 0.9999999;                      // Denotes almost one
const ELEMENT_MAX = 1000;                          // The maximum number of elements which can be displayed
const PARAM_LIST_MAX = 1000;                       // The maximum number of parameters on list
const HASH_MAX = 10000;                            // The value used for the hash tables
const H_BIN = 10;                                  // Used for distributions in cumulative prob

const COR_MAX = 0.9;                               // Maximum correlation for individual effects

const SIM_NUM_DEFAULT = 1;                         // The default simulation number
const PPC_NUM_DEFAULT = 200;                       // The default number of ppc simulation
const ANNEAL_DEFAULT = "none";                     // Default annealing type
const ALG_DEFAULT = "DA-MCMC";                     // Default inference algorithm
const SEED_DEFAULT = 0;                            // The default seed
const ANNEAL_POWER_DEFAULT = "4";                  // Default annealing power
const ANNEAL_RATE_DEFAULT = "0.01";                // Rate at which annealing is done
const COORD_DEFAULT = "cartesian";                 // Default coordinate system
const INTERNAL_ERROR = " This is an internal error to BICI and not something you have done wrong! Please send the BICI script to Chris for diagnosis.";

const PARAM_OUTPUT_MAX_DEFAULT = 1000;             // The default maximum number of tensor elements to be output
const INDMAX_DEFAULT = 20000;                      // The default maximum number of individuals
const BURNIN_FRAC_DEFAULT = 20;                    // The default percentage burnin
const MCMC_SAMPLE_DEFAULT = 5000;                  // The default number of MCMC samples
const MCMC_OP_PARAM_DEFAULT = 1000;                // The default number of output parameters
const MCMC_OP_STATE_DEFAULT = 200;                 // The default number of output states  
const MCMC_CHAIN_DEFAULT = 3;                      // The default number of MCMC chains
const MCMC_CHAIN_PER_CORE_DEFAULT = 1;             // The default chains per core
const PAS_PART_DEFAULT = 10;                       // The default number of PAS particles
const PAS_GEN_UPDATE_DEFAULT = 1;                  // The default updates per generation percent  
const PAS_PART_PER_CORE_DEFAULT = 1;               // The default particles per core

const ABC_SAMPLE_DEFAULT = 1000;                    // The default number of ABC samples
const ABC_ACFRAC_DEFAULT = 0.1;                    // The default acceptance fraction for ABC
const ABCSMC_SAMPLE_DEFAULT = 1000;                 // The default number of ABCSMC samples
const ABCSMC_ACFRAC_DEFAULT = 0.5;                 // The default acceptance fraction for ABCSMC
const ABCSMC_KERNEL_DEFAULT = 0.5;                 // The default kernel size for ABCSMC 
const ABCSMC_GEN_DEFAULT = 5;                      // The number of generations for ABCSMC
const GRAPH_EXTRA = 1;                             // The extra space given on the right hand side of line plot
const COMP_NOISY_MAX = 10;                         // When simulating noisy compartmental observations gives max number
const si_anno = 1.4;
const size_annotation_default = 10;                // Default size of annotation text
const annotation_col_default = BLACK;              // Default colour for labels
const import_frac_pro = 0.5;

// Different line thicknesses
const NOLINE = 0, THINLINE = 0.5, NORMLINE = 1, MEDIUMLINE = 1.5, THICKLINE = 2, MTHICKLINE = 3, VTHICKLINE = 4;

const SD_MIN = TINY;                         // The minimum value for sd
const SD_MAX = LARGE;                        // The maximum value for sd
const CV_MIN = 0.01;                         // The minimum value for cv
const CV_MAX = 10.0;                         // The maximum value for cv
const MEAN_MIN = VTINY;                      // The minimum for dist mean
const MEAN_MAX = LARGE;                      // The maximim for dist mean
const NORM_MEAN_MIN = -LARGE;                // The minimum for normal mean
const NORM_MEAN_MAX = LARGE;                 // The maximum for normal mean
const SCALE_MIN = TINY;                      // The minimum for dist scale
const SCALE_MAX = LARGE;                     // The minimum for dist scale
const SHAPE_MIN = 0.1;                       // The minimum for dist shape
const SHAPE_MAX = 100.0;                     // The maximum for dist shape
const ALPBETA_MIN = 0.01;                    // The minumum for alpha or beta
const ALPBETA_MAX = 100;                     // The minumum for alpha or beta
const LAM_MIN = 0;                           // The minimum value for poisson
const LAM_MAX = LARGE;                       // The maximum value for poisson
const EXP_MEAN_MIN = TINY;                   // Minimum mean for exponential
const TIME_MIN = TINY;                       // The minimum time
const RATE_MIN = TINY;                       // The minimum rate
const P_MIN = TINY;                          // The minimum probability
const P_MAX = 1-TINY;                        // The maximum probability
 
// Sets diffent distribution for displaying text
const NORM_TE=0, LOGNORM_TE=1, WEIBULL_TE=2, GAMMA_TE=3, BETA_TE=4, NEGBINO_TE=5, BERN_TE=6, EXP_RATE_TE=7, EXP_MEAN_TE=8, POIS_TE=9, PERIOD_TE=10;

// Different quantities
const SD_QU=0, CV_QU=1, MEAN_QU=2, NORM_MEAN_QU=3, SHAPE_QU=4, SCALE_QU=5, ALPHA_QU=6, BETA_QU=7, P_QU=8, BERNP_QU=9, RATE_QU=10, EXP_MEAN_QU=11, POIS_QU=12, TIME_QU=13;

