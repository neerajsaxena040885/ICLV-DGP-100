% This code prepares the X vector to be fed into the Structural Part

% POINTS TO REMEMBER:
% 1. The mean of continuous variables should be within the range 2-7
% 2. This can be done by conducting arithmetic operations on the variable
% 2. Generate range of X values by using the formaula:
%    X = (max-min).*rand(size,1) + min;

% --------------------- SET VALUES SECTION ------------------------
% Number of respondents 
NP=100;

% Number of alternatives in each choice task ... 
NALT=2;

% Number of choice tasks in the experiment ... t
NC=3;

% Number of rows in the PANEL dataset .. 
NPANEL=NP*NALT*NC;

% Number of indicators... I
NIND=1;

% Seed value
SEED = 28593 + randi(11000);
rand('state',SEED)  %For draws from uniform
randn('state',SEED)  %For draws from normal

% ------------ PREPARING CONTINUOUS VARIABLES ---------------

% Time spent in Stop-&-go
meantts=1;
stdtts=1.5;

% Number of Stop-&-go experienced
meansng=1;
stdsng=2;

% Travel time
meantt=1;
stdtt=2.75;

% Vehicle Running Cost
meanvrc=1.5;
stdvrc=1.5;

% Generate vector of Xtts
Xtts=meantts + stdtts.*rand(NPANEL,1);
% Generate vector of SnG
Xsng=meansng + stdsng.*rand(NPANEL,1);
% Generate vector of Xtt
Xtt=meantt + stdtt.*rand(NPANEL,1);
% Generate vector of Xvrc
Xvrc=meanvrc + stdvrc.*rand(NPANEL,1);

% Check the means.. They must all be comparable
Xttsbar=mean(Xtts);
Xsngbar=mean(Xsng);
Xttbar=mean(Xtt);
Xvrcbar=mean(Xvrc);

% -------------------------------------------------------------------

% --------- PREPARING CASE SPECIFIC VARIABLES (CATEGORICAL) ---------

% Preparing age vector
cumage=[0.4];       % BINARY VARIABLE AT THE MOMENT
% Preparing driving experience vector
cumdr=[0.48];       % BINARY VARIABLE AT THE MOMENT

% Generate uniform draws
temp01=rand(NP,1);
temp02=rand(NP,1);

for n=1:NP
    Xage1(n,:)=(temp01(n,:)<=0.4);   % Majority selection is the base
    Xdr1(n,:)=(temp02(n,:)<=0.48);   % Majority selection is the base
end

% Check for means of Categorical valriables
Xagebar=mean(Xage1);
Xdrbar=mean(Xdr1);

Xage=kron(Xage1,ones(NALT*NC,1));  % Size of the vector is NPANEL X 1
Xdr=kron(Xdr1,ones(NALT*NC,1));    % Size of the vector is NPANEL X 1

% -------------------------------------------------------------------

% Write the variables to a file
% Formaty is: Alternate specific followed by case specific
inputdataX = horzcat(Xtt,Xtts,Xsng,Xvrc,Xage,Xdr);
csvwrite('inputX.csv', inputdataX);
