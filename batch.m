% -----------------------------------------------------------------
%       INTEGRATED CHOICE AND LATENT VARIABLE MODEL (ICLV)
%           Written by: Neeraj Saxena, 22nd Jul 2017 
%            emailid: n.saxena@student.unsw.edu.au
%                 !! SHRI GANESHAAYE NAMAH !!
% -----------------------------------------------------------------

% INTRODUCTION

% This script executes the model in a batch and stores the results from
% each run.
% This file should be executed from now on and not main.m

% NOTE: The vector of Xs and I given in "InputX.csv" remains the same
% across samples. In other words, there is NO need to run GenerateX.m for
% each sample.
% Secondly, the file main.m from Datagen folder should be run
% Thirdly, the file should be read 
% Fourthly, main.m for model estimation should be run
% -----------------------------------------------------------------

% -----------------------------------------------------------------
%                            CHANGE LOG
% None
% -----------------------------------------------------------------

% OUTPUT FILE
% Put the name you want for your output file (including full path if not the current 
% working directory) after words "delete" and "diary".
% The 'diary off' and 'delete filename' commands close and delete the previous version 
% of the file created during your current matlab session or in any previous sessions. 
% If you want to append the new output to the old output, then 
% put % in front of the 'diary off' and 'delete filename' commands (or erase them).

tic
diary off
delete modelrun.out
diary modelrun.out

% Prepare the explanatory variables first. Its a one time thing
   cd('./DataGen')
   GenerateX
   cd('../')
   
   clear all 
   
   
% Number of samples to be generated 
sizesample=50;

for sample=1:sizesample

    % Set seed for the random number generator.
    SEED = 25834+randi(11000);
    % Setting the state of a random number generator. For more info read the
    % weblink: http://au.mathworks.com/matlabcentral/answers/17701-rand-state-11
    randn('state',SEED)  %For draws from normal
    rand('state',SEED)   %For draws from uniform
    
    % Declare GLOBAL variables
    % GLOBAL variables are all in caps
    % DO NOT CHANGE ANY OF THESE 'global' STATEMENTS
    global NP NC NCS NROWS NALTMAX NCSMAX
    global LV IDLV NLV
    global IDIND NIND RAT
    global XMAT
    global choice Xchoice I XS XC 
    global XD UpInd LowInd numele nfl nchc
    global NV NAMES
    global IDF NF IDF1 IDF2 IDF3 IDS NS NCH IDC
    global DRAWTYPE NDRAWS AMPL MASTERDRAWS MASTERDR CHOICEDR
    global NMEM NTAKES
    
    % TITLE
    % Put a title for the run in the quotes below, to be printed at the top of the output file.
    disp '**** Integrated Choice & latent Variable Model (ICLV) ****'
    disp(['-------> For Sample ID: ' num2str(sample) ' ']);
    
    % Generate the data file first... 
    cd('./DataGen')
    main
    
    % Once the data is generated, read the datafile "InputC.csv"
    XMAT=csvread('InputC.csv');  %The variables are described below

    % Call the script main.m to estimate the model
    cd('../')
    main
    
    % Once the run succeeds for a particular run, save the results into an array
    sampleparamhat(sample,:)=horzcat(sample,paramhat);
    samplestderr(sample,:)=horzcat(sample,transpose(stderr));
    
    % Clear variables from the memory    
    clearvars -except sample sizesample SEED sampleparamhat samplestderr
    
end

% Write the arrays into a csv file
csvwrite('Paramhats.csv', sampleparamhat);
csvwrite('Stderrs.csv', samplestderr);

toc
% Once the entire run finishes. 
disp(' ')
disp ('**** Simulation Ended ****')
disp(['Number of samples estimated: ' num2str(sample) ' ']);
disp(['Execution time: ' num2str(toc./3600) ' hours.']);

