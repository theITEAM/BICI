// Lists all the variables used in the graphical interface

var	width, height, widthold, heightold;                   // Define the size of the window

var nchrun=0, nchain=1;                                   // The number of MCMC chains used 

var runtype = "";                                         // The type of run (either "sim" for simulation or "inf" for inference)

var datanote, datanoteon, descnote;                       // Provide description of the data and model

var maincv, graphcv, cv, resultcv;                        // Used for drawing on canvases in HTML5

var intervalid;                                           // Id for Interval which automatically draws trace plots

var convid;                                               // Id for Timeout which triggers converegence test 

var corid;                                                // Id for TImeout for correlation calculation

var kernid;                                               // Id for TImeout for kernal estimation

var indplotst=[];                                         // Stores information about individual (do it doesn't have to be recalculated)

var obsloaded = 0;                                        // Determines if observations have been loaded

var modelstart = 0, modelsetup = 0;                       // Derermines if a new model has been started or setup

var examppic=[];                                          // Stores png images for examples on home page

var xmlhttp;                                              // Used to load models

var shx, shy, shfac;                                      // Describes the shifting and rescaling of model representation

var ncla=0, cla=[];                                       // Stores all information about the classifications in the model

var age=[];                                               // The list of age transitions

var time=[];                                              // The list of time transitions

var param=[], paramsim=[], paramnew=[], paramsimnew=[];   // Stores all the parameters or priors in the model

var paramagesmooth=[];                                    // Information about parameter age smoothing

var paramtimesmooth=[];                                   // Information about parameter time smoothing

var exporton = 0;                                         // Turms on when export submenu opens
var saveon = 0;                                           // Turms on when save submenu opens

var nsampevmax = 500;                                     // The maximum number of event samples stored 
var nsampmax = 10000;                                     // The maximum number of parameter samples stored 

var GRmax = 1.1;                                          // The maximum value for the Gelman rubin statistic
var ESSmin = 400;                                         // The minimum value for the effective sample size
var itermax = 10000;                                      // The number of iteration to run before termination
var termtype=0;                                           // The type of termination (0=None,1=Use ESS and GR)

var fiformat;                                             // Stores the file format

var fileToLoad, newfile;                                  // Name of loaded fileCreatedDate
var fileToLoadlist=[];                                    // Stores loaded files                        

var IDcol;                                                // Stores the column used to represent individual ID

var BFtext;                                               // Text output from Bayes factor calculation

var helptype = -1;                                        // Determines what help window is open

var exampst, examploaded;                                 // Temorarily stores name of example 

var simst, simsta;                                        // Indicate start of simulation

var tsimmin, tsimmax;                                     // The time range for simulation

var simty, simnumber = 10;                                // Determine the number of simulations to be performed

var indmaxnumber;                                         // The maximum allowable number of individuals

var errbar = 0, errbarscale = 1, errbarfix = 2;           // Determine error bars for population data 

var pageyst=[];                                           // Stores the scroll bar position
        
var radstyle = 0, radstyle2 = 1;                          // Determine plotting: colour/greyscale

var radlab = 2, radnamlab = 1;                            // Determine plotting: labels/not, name/not
 
var distCI=1, distmean = 1;                               // Determines if the credible interval or mean are calcualted                        

var advop;                                                // Detemines if the advanced options are shown

var subna;                                                // Stores the subpage name

var ytot;                                                 // The total height of a page

var derive=[];                                            // List of derived quantities
var dergensel;                                            // List of selected derived values

var adddepfl;                                             // Is set if adding a distribution in the model

var exe=[], child=[];                                     // Stores which chains are running    
                                   
var cornx, corny;                                         // The position of the canvas on the page

var gdropinfo=[], gdropsel=-1, gdropfrac, gdropnum;       // Used for dropdown menus
var gdropdy=20, gdropmyst, gdropfracst, gdropbut, gdropselop;

var probdetection;                                        // The probability of detection

var corlist=[], cormat=[], corx, cory;                    // Used for correlations

var ncpu;                                                 // The number of CPUs

var canmx, canmy;                                         // The mouse position when transformed onto the canvas

var nob, obx=[], oby=[], obty=[];                         // Objects (groups of buttons plotted together)
var obval=[], obval2=[], obval3=[], obval4=[]; 

var tablehist=[];                                         // Store changes to a table (to allow undo option)

var addingdata = 0;                                       // Determines is data is being added

var datashow = -1;                                        // The type of data being displayed

var datatemp, dataselected;                               // A specific data element is selected

var transty, transcl, transni, transnf;                   // Stores information about a selected transition

var modelty;                                              // The type of the model being drawn on page

var tpostmin, tpostmax;                                   // The posterior time range

var tlinetmin, tlinetmax;                                 // The time range for individual time lines

var inddata;                                              // Stores information about all the individuals in the data

var indsim;                                               // Stores information about all the individuals in the initial state in the simulation

var data=[];                                              // All the data used for inference

var simdata=[];                                           // Generated "data" used for simulation

var indshow=[];                                           // Determines which classification timelines are shown
 
var datatype;                                             // Stores the type of data being uploaded

var nbut, buttext=[], butx=[], buty=[];                   // Stores information about buttons
var butdx=[], butdy=[], butac=[], buttype=[];
var butover=[], butval=[], butval2=[];

var ncanbut, canbuttext=[], canbutx=[], canbuty=[];       // Stores information about canvas buttons
var canbutdx=[], canbutdy=[], canbutac=[], canbuttype=[];
var canbutover=[], canbutval=[], canbutval2=[];

var nvarlist, varlist=[], vardeplist=[];                  // The variables and their dependencies taken from an equation

var popflag;                                              // Set if the equation contains a population termination

var cladepflag=[];                                        // Set if the equation is dependent on a given classification

var ccanx, cany, anvasdx, canvasdy;                       // The position and size of the canvas

var nid, id=[];                                           // Store the ids of all the individuals

var simstarttext;                                         // Store how a simulation is started (i.e. if it is being rerun or not)

var page = HOMEPAGE, pagesub=[], pagesubsub=[];           // Stores the page/subpage being viewed

var over = -1;                                            // The button the mouse is over

var canover = -1;                                         // The canvas button the mouse is over

var fitype;                                               // The type of file loaded

var mx, my;                                               // The mouse position

var timedown, timeup;                                     // The timing the mouse was pushed down and went up

var arrow = 0;                                            // Determines the pointer (arrow or hand)

var nout, out=[];                                         // Stores a list of submenus

var overx=-1, overy, overdx, overdy;                      // The position of the overlay (used for help)

var nlines, lines=[], linesy=[], lineheight = 25;         // Used when drawing paragraphs

var mxst, myst, mxlast, mylast, drag, dragged, dragval;   // Used when dragging objects
var dragval2, dragtext, dragx, dragy;
var dragdx, dragdy, dragdh, dragww;

var clgl;                                                 // A selected classification

var transtempcl, transtempk, transoverj;                  // Stores details of a transition

var errmsg="", errormsg="", errmsgline=[];                // Error messages

var new_win = -1;                                         // Stores ID for newly openned window
 
var popdx, popdy;                                         // Used when openning a new window                        

var inpon = 0;                                            // Set if input is being given in bubble

var addtransmode, transtempi, addcompmode, addsourcemode; // Set when adding transition to model

var uoflag;                                               // Set if unobserved individuals in the inference

var induoflag=0;                                          // Set if there are unobserved individuals in the system

var running = 0;                                          // Set if inference is running

var percentrun, percentload, loading = 0;                 // Set if loading and what percent

var varCImin=[], varCImax=[];                             // The calculated credible interval

var varESS=[];                                            // The effective sample size

var varGR=[];                                             // The Gelmann Rubin statistic

var varmean=[];                                           // The mean variable value

var infres={result:-1}, ppcres;                           // Store results of inference/simulation
var simres={result:-1}, datares, res, rest;

var leftover=[];                                          // Used to buffer output from C++ code

var infpagename=[], simpagename=[];                       // The names of the inference and simulation pages

var tableyfr, tableyfrac;                                 // Used for the y scroll bar

var ytotright, ybeg, ymax;                                // Used for the y scroll bar on the right menu

var sliy1, sliy2;                                         // Range for y scroll bar

var vcalc, vcalced=[], vconvtest, timestart;              // Keeps track of statistics that have been calculated

var tempCI=[], CImin, CImax, mean, sd;                    // Stores statistics data

var col, col2;                                            // Stores colours used for plotting the table

var tags=[], tagst, typest, jsto;                         // Used when importing  

var searchterm="", searchres=[], searchresnum;            // Used when searching in a table

var replaceterm="", nreplace, ndel;                       // Used when replacing and deleting

var warning=[];                                           // Stores warnings generated when inference is started

var warnpic;                                              // The warning png image

var foc;                                                  // Stores focus for a text area

var selectbub, selectbubst, selectbubval, selectbubval2;  // Bubble for a selected element
var selectbubtext, selectbubx, selectbuby;
var selectbubdx, selectbubdy;

var bubx, buby, bubw, bubh, bubdx, bubdy;                 // Bubble size and position
var bubx1, bubx2, bubx3, buby1, buby2, buby3;

var cursx, cursy;                                         // The position when drawing in bubble

var ord = [];                                             // Used for ordering transition times

var calclist=[];                                          // Gives a list of all the event sequences need to be analysed to geenrate transition plot

var points=[];                                            // Stores bullet points

var tableheight, simheight;                               // Store the positions of sizes of objects on the page
var corheight, corwidth; 
var indtablewidth, indtableheight;
var siminitwidth, exampheight, priorheight, priorwidth;
var rightheight, tlinexmax;
var tablewidth, tottablewidth, tottableheight;
var obswidth, obsheight;
var tablewidthdesc, tableheightdesc;
var canvaswidth, canvasheight;
var rowmax, tableyfrac, tablexfrac;
var tablewidthbig, tableheightbig;
var derwidth, derheight;
var tablewidthbigger, tableheightdata; 
var graphdy, graphdx, graphdy;
var modeldx, modeldy, modeldytrans, modelplaydx;
var initmodeldx, rightval, wmax;

var clagendata=[], gentype = -1;                          // Varibles used when generating data
var detectprob, datagenttype, datagenttype2;
var tgenmin, tgenmax, dtgen, tgenuserdef;
var tlist=[], poplist=[], errbarlist=[];
var drawgennum = -1, datagen=[], Segen, Spgen;
var selclass=[], sellist=[], siminit="";
var transselcl, transseltr;
var whichind="all", obspd="all";

var ncol, ncoldef=0, ncoldefmax, colname=[];              // Information about tables
var colx=[], colw=[];
var nrow, row=[], rowwidth=[], rowcopy=[];

var chmin, chmax, smin, smax, plotinitfl = 0;             // Used when plotting
var sstart, plotmi, plotma, calcto;
var ncomp, compmult=[], compval=[], cpop=[];   
var NOTALIVE;
var simpopinitset=0, viewpopinitset=0, simpopinitlist=[];
var simindinit=[], showsimindinit = 0, popmax;
var claselect=[]; 
var depsel=[], histoplot=[], depop=[], depopsel;
var nticksize, ticksize=[]; 
var polypoint=[]; for(i = 0; i < 10; i++) polypoint[i] = [];
var gx, gy;
var axxmin, axxmax, axymin, axymax;
var axtickx, axticky;
var ntickx, tickx=[], nticky, ticky=[];
var popfilter=[], playtime, poptmin, poptmax, playstartt;
var playing=0, playtimeinit, poppage, poppagesub, poppagesubsub;
var JX = 200, Jbin=[], Jbinpriorfac=[];
var KX = 200, KY = KX, Kbin=[];
var scatsamp=[];
var POPX;
var nquant, quant=[], quanteqn=[], quantcol=[];
var quantmean=[], quantCImin=[], quantCImax=[], quantname=[], quantmaxi=[], quantmaxitot;
var popinit=[], popch=[], npop;
var lablist=[], labcollist=[];
var transpos=[], transsel=[], transmat=[], transindfilt, transindcheck=[], transcfilt=[];
var cusamp=[], culim=[];
var curvevar=[], curvech=[];
var xaxisauto = 1, xaxisfixmin, xaxisfixmax, yaxisauto = 1, yaxisfixmin, yaxisfixmax;
var ctrlkey;
var kde, kdest, kdemin, kdemax, smoothon;
var chooseaxis;
var transcfilt=[]; 
