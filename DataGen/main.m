% -----------------------------------------------------------------
%         INTEGRATED CHOICE AND LATENT VARIABLE MODEL
%           Written by: Neeraj Saxena, 11TH aUG 2017 
%            emailid: n.saxena@student.unsw.edu.au
%                 !! SHRI GANESHAAYE NAMAH !!
% -----------------------------------------------------------------
%
% INTRODUCTION
% ICLV models are generally used in cases where the number of indicators is
% large making the resulting variance-covariance matrix as computationally
% expensive. Thus, the structural part is introduced which clubs indicators
% together to form latent constructs leading to a smaller VC matrix.
%   Since, I have got few indicators, it is better that i ignore the
% structural component of the ICLV model.

% Thus, the new model will have two components now
% 1. Measurement model --> Where I relate indicators in terms of
% socio-demographic information
% 2. Choice model --> Using indicators and attributes to determine the
% choice

% We will evaluate the following probability
% Probability of selecting a route conditional on the observed indicators

% This script prepares data for the model

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

% Number of latent variables... LV
LV=1;

% Read the X vector from the input file.
XDGP=csvread('inputX.csv');


% Make Standard Normal draws ...

% For the Structural Part.. It should be the same for an individual
z0=normrnd(0,1,LV,NP);  % Size of the vector is LV X NP
zs=kron(z0,ones(1,NALT*NC));  % Size of the vector is LV X NP*NC*NALT
zs=transpose(zs);   % Size of the vector is NP*NC*NALT X LV


% For the Measurement Part.. It should be different on each route for every individual
% z1=normrnd(0,1,1,NP);       % Size of the vector is 1 X NP
% zl=kron(z1,[1 0 1 0 1 0]);  % Size of the vector is 1 X NP*NC*NALT
% z2=normrnd(0,1,1,NP);       % Size of the vector is 1 X NP
% zr=kron(z2,[0 1 0 1 0 1]);  % Size of the vector is 1 X NP*NC*NALT
% zm=zl+zr;
zm=normrnd(0,1,1,NP*NALT*NC);   % Size of the vector is 1 X NP*NC*NALT
zm=transpose(zm);               % Size of the vector is NP*NC*NALT X 1


% For the Choice Part.. It should be the the same for an individual and route
% z3=normrnd(0,1,1,NP);       % Size of the vector is 1 X NP
% zc=kron(z3,[0 1 0 1 0 1]);  % Size of the vector is 1 X NP*NC*NALT
z3=normrnd(0,1,1,NP*NC);       % Size of the vector is 1 X NP*NC
zc=kron(z3,[0 1]);  % Size of the vector is 1 X NP*NC*NALT
zc=transpose(zc);           % Size of the vector is NP*NC*NALT X 1


% -------------------------------------------------------------------------
%                          FIX PARAMETER VALUES
% -------------------------------------------------------------------------

% Structual Part
a=[0.55 0.2 0.5 0.65];

% Measurement Part
bc=[-0.95;-1.10];
d=[1.05];

mu1=[-0.5 -0.6 -0.65];
mu2=[-0.4 -0.55 -0.7];
MCmu=vertcat(mu1,mu2);

% Choice Part
c=[-0.9 -0.7];
cl=[-0.8];


% -------------------------------------------------------------------------
%                     COMPUTE STRUCTURAL EQUATION
% -------------------------------------------------------------------------

XS=XDGP(:,[2;3;5;6]); % Vector of structural covariates

% Estimate the Latent Variable .. FLV
FLV = XS*a' + zs;   % Size of the vector is NP*NC*NALT X 1



% -------------------------------------------------------------------------
%                    COMPUTE MEASUREMENT EQUATION
% -------------------------------------------------------------------------

% Exponentiating the MCmu matrix
MCmu=exp(MCmu);
% Doing a cumulative of elements 
MCmu=cumsum(MCmu,2);
% Adding zero as the first element
MCmu=horzcat(zeros(2,1),MCmu);

% Estimate the Indicator latent variable .. Istar
Istar = repmat(bc,NP*NC,1) + d.*FLV + zm;       % Size of the vector is NP*NC*NALT X 1

% Generate Indicator Ratings
for n=1:NPANEL
    if rem(n,2)==0
       i=2;
    else
       i=1; 
    end    
    if Istar(n,1) <= MCmu(i,1)
        I(n,1)= 1;
    end
    if (Istar(n,1) > MCmu(i,1)) && (Istar(n,1) <= MCmu(i,2))
        I(n,1)= 2;
    end
    if (Istar(n,1) > MCmu(i,2)) && (Istar(n,1) <= MCmu(i,3))
        I(n,1)= 3;
    end
    if (Istar(n,1) > MCmu(i,3)) && (Istar(n,1) <= MCmu(i,4))
        I(n,1)= 4;
    end
    if Istar(n,1) > MCmu(i,4)
        I(n,1)= 5;
    end
end

% Check the data by using these commands
I1=length(find(I(:,1)==1));
I2=length(find(I(:,1)==2));
I3=length(find(I(:,1)==3));
I4=length(find(I(:,1)==4));
I5=length(find(I(:,1)==5));



% -------------------------------------------------------------------------
%                         COMPUTE CHOICE PART
% -------------------------------------------------------------------------

XC=XDGP(:,[1;4]); % Vector of choice variables

% Evaluate the observed utility
Util = XC*c' + cl.*FLV + zc;    % Size of the vector is NP*NC*NALT X 1

% Code the Utilities into observed binary outcome.
% Within each choice task, higher utility gets 1 and other as 0
k=1;
for n=1:NALT*NC:NPANEL
    for j=1:NC
        if Util(k,:)>Util(k+1,:)
            Choice(k,:)=1;      % Size is NPANEL X 1
            Choice(k+1,:)=0;
        else
            Choice(k,:)=0;
            Choice(k+1,:)=1;
        end
        k=k+2;
    end
end

% Preparing a column of PersonIDs
for n=1:NP
    for j=1:NALT*NC
      PID(NALT*NC*(n-1)+j,1)=n;
    end
end

% Preparing a column of Choice TaskIDs
for n=1:NP 
      SID(NALT*NC*(n-1)+1,1)=3*(n-1)+1;
      SID(NALT*NC*(n-1)+2,1)=3*(n-1)+1;
      SID(NALT*NC*(n-1)+3,1)=3*(n-1)+2;
      SID(NALT*NC*(n-1)+4,1)=3*(n-1)+2;
      SID(NALT*NC*(n-1)+5,1)=3*(n-1)+3;
      SID(NALT*NC*(n-1)+6,1)=3*(n-1)+3;
end

% Preparing a column of Alternative IDs
for n=1:NALT:NP*NALT*NC 
      ALTID(n,1)=1;
      ALTID(n+1,1)=2;
end

inputdataC = horzcat(PID,SID,ALTID,Choice,XDGP,I);
csvwrite('inputC.csv', inputdataC);
