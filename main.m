% -----------------------------------------------------------------
%       INTEGRATED CHOICE AND LATENT VARIABLE MODEL (ICLV)
%           Written by: Neeraj Saxena, 24th Jun 2017 
%            emailid: n.saxena@student.unsw.edu.au
%                 !! SHRI GANESHAAYE NAMAH !!
% -----------------------------------------------------------------

% INTRODUCTION: 
%
% -----------------------------------------------------------------
%                            CHANGE LOG
% None
% -----------------------------------------------------------------

% clear all

% DATA

% Number of people (decision-makers) in dataset.. CHANGE THIS FIELD !!!!
NP=100;        

% Number of choice tasks each one faces.. CHANGE THIS FIELD !!!!
NC=3;

% Number of choice situations in dataset. This is the number faced by all
% the people combined.. 
NCS=NP*NC;  

% Load and/or create XMAT, a matrix that contains the data.
%
% XMAT must contain one row of data for each alternative in each choice situation for each person.
% The rows are grouped by person, and by choice situations faced by each person.
% The number of rows in XMAT must be NROWS, specified above.
% The columns in XMAT are variable that describe the alternative.
% 
% The *first* column of XMAT identifies the person who faced this alternative. 
% The people must be numbered sequentially from 1 to NP, in ascending order.
% All alternatives for a given person must be grouped together.
% The *second* column of XMAT identifies the choice situation. The choice
% situations must be numbered sequentially from 1 to NCS.
% All alternatives for a given choice situation must be grouped together.
% The *third* column of XMAT identifies the chosen alternatives (1 for
% chosen, 0 for not). One and only one alternative must be chosen for each
% choice situation.
% The remaining columns of XMAT can be any variables.

% XMAT=csvread('Choicedata.csv');  %The variables are described below

% To help you keep up with the variables, list the variables in XMAT here.
% Start each line with % so that matlab sees that it is a comment rather than a command.
% NOTES for XMAT for Stop-&-go ICLV run:
% This dataset is for people's choice among two routes in stated-preference
% experiments. Each person faced with 3 experiments. Each
% experiment contained 2 alternatives representing two different routes, 
% each defined in terms of 4 attributes (TT, TTS, SnGo and VRC). The person stated which
% of the three routes he/she would take for their travel.
% The variables in XMAT are:
% 1. Person number (1-NP)            MUST BE THIS. DO NOT CHANGE.
% 2. Choice situation number (1-NCS) MUST BE THIS. DO NOT CHANGE.
% 3. Chosen alternative (1/0)        MUST BE THIS. DO NOT CHANGE.
% 4. Total travel time on a route (minutes)
% 5. Time spent in stop-&-go traffic (minutes)
% 6. Number of stop-&-go experienced (number)
% 7. Vehicle running cost for the trip (AU $)
% 8. Indicator variables 1..n
% 9. Socio-demographic variables 1..n

% NOTE: Refer to the input file layout before populating the following fields
%----------- SET THE ATTRIBUTES (XS) FOR THE STRUCTURAL MODEL ----------
IDS=[6;7;9;10];
NS=size(IDS,1);     
% 
%----------- SET THE INDICATORS (I) FOR MEASUREMENT MODEL ----------
IDIND=[11];
% IDIND=[14];
NIND=size(IDIND,1);     % Number of Indicators.. Measurement part.. 
%
%----------- SET THE ATTRIBUTES (XC) FOR THE CHOICE MODEL ----------
IDC=[5;8];
NCH=size(IDC,1);     
% 
% Rating for each Indicator 
RAT=5;

% OPTIMISATION 
% Maximum number of iterations for the optimization routine.
% The code will abort after ITERMAX iterations, even if convergence has
% not been achieved. The default is 400, which is used when MAXITERS=[];
MAXITERS=[];

% Convergence criterion based on the maximum change in parameters that is considered
% to represent convergence. If all the parameters change by less than PARAMTOL 
% from one iteration to the next, then the code considers convergence to have been
% achieved. The default is 0.000001, which is used when PARAMTOL=[];
PARAMTOL=0.000001;
% PARAMTOL=0.0001;

% Convergence criterion based on change in the log-likelihood that is
% considered to represent convergence. If the log-likelihood value changes
% less than LLTOL from one iteration to the next, then the optimization routine
% considers convergence to have been achieved. The default is 0.000001,
% which is used when LLTOL=[];
LLTOL=[];
% LLTOL=0.0001;

%Do not change the next line. It runs the model.
doitr

