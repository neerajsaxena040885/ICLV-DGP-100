% Calculates log-likelihood function value for the ICLV Model
%
% This code is input to Matlab's fminunc command
%
% Input param is a row vector of parameters, dimension "check doitr line 57"
%
% Output ll is the scalar value of the negative of the simulated log-likelihood 
% at the input parameters
% 
% PURPOSE of this code is to split the param array into beta values, which then go into other functions loglikmeasure and llgrad2
%

function [ll, g] =loglik(param)     % The contents of param array are --> {Fixed and Sigma}
                                    
global NIND RAT NP NS NALTMAX XS NCSMAX UpInd LowInd nfl nchc XC NCH choice

% Splitting the parameter array into ...

% Structural coefficients 
struct=param(1,1:NS);       % Size of vector is 1 X NS

% Measurement coefficients (constant)
if NIND>0
    Bc=param(1,NS+1:NS+NALTMAX);         
    Bc=transpose(Bc);              % size of Bc is NALTMAX X 1
else
    Bc=[];
end

% Measurement Factor Loading
if NIND>0
    Fac=param(1,NS+NALTMAX+1:NS+NALTMAX+nfl);    % size of Fac is 1 X nfl
else
    Fac=[];
end

% Measurement Thresholds
for n=1:NALTMAX
    eval(['wk' num2str(n) '= param(1,NS+NALTMAX+nfl+1+(n-1)*((RAT-1)-1):NS+NALTMAX+nfl+n*((RAT-1)-1));']);
    eval(['wk' num2str(n) '= exp(wk' num2str(n) ');']);     % All threshold points to be exponentiated
    eval(['wk' num2str(n) '= cumsum(wk' num2str(n) ');']);      % Taking cumulative of the exponentiated elements 
    eval(['wk' num2str(n) '= [-Inf 0 wk' num2str(n) ' Inf];']);      % Adding the first threshold point as ZERO. This is being done because of a revised utility specification for the ordered model.
end

% Concatenating wk arrays vertically
for n=1:NALTMAX
    if n == 1
        eval(['threshvec = wk' num2str(n) ';']);
    else
        eval(['threshvec = vertcat(threshvec,wk' num2str(n) ');']);   % Size of threshvec is NALTMAX X (RAT+1)
    end
end

% Choice coefficients
chcoeff=param(1,NS+NALTMAX+nfl+NALTMAX*((RAT-1)-1)+1:NS+NALTMAX+nfl+NALTMAX*((RAT-1)-1)+nchc);



% ---------------------------------------------------------------------------
%                       COMPUTING THE JOINT PROBABILITY
% ---------------------------------------------------------------------------

% Initialisation
PM=zeros(nchoosek(NALTMAX*NCSMAX,2),NP);    % Measurement Probability Matrix. Size is (NALTMAX*NCSMAX)C2 X NP
PC=zeros(1,NP);     % Choice Probability Matrix. Size is 1 X NP
PI=zeros(NALTMAX*NCSMAX*(NALTMAX-1)*NCSMAX,NP);    % Measurement-Choice Probability Matrix. Size is (NALTMAX*NCSMAX*(NALTMAX-1)*NCSMAX) X NP

% Calculate the joint probability for each individual
for n=1:NP
    % Set the starting and ending row indices
    startrowind=(n-1)*NALTMAX*NCSMAX+1;
    endrowind=n*NALTMAX*NCSMAX;
    
    
    % ---------------------------------------------------------------------------
    %                            PREPARING THE VECTORS
    % ---------------------------------------------------------------------------
    
    % Prepare the Latent Factor (Frustration) --> STRUCTURAL PART
    frustLV=XS(startrowind:endrowind,:)*transpose(struct);       % Size of the vector is NALTMAX*NCSMAX X 1
    
    % Form the measurement equations (Frustration Rating)   --> MEASUREMENT PART
    ystar=repmat(Bc,NCSMAX,1) + Fac.*frustLV;                    % Size of the vector is NALTMAX*NCSMAX X 1
    
    % Evaluate the Observed Utility of all the alternatives     --> CHOICE PART
    U=XC(startrowind:endrowind,:)*transpose(chcoeff(1,1:NCH)) + chcoeff(1,NCH+1:end).*frustLV;    % Size of U NALTMAX*NCSMAX X 1
    
    % Prepare the lambda coefficient... Refer to Bhat's paper
    lamda=repmat(chcoeff(1,NCH+1:end),NALTMAX*NCSMAX,1);        % Size of vector is NALTMAX*NCSMAX X 1
    
    % Vertically concatenate the measurement and choice vectors
    yvector=vertcat(ystar,U);    % Size of yvector is 2*NALTMAX*NCSMAX X 1   
    
    % Prepare the matrix of Mstar
    for k=1:NCSMAX
        if (choice((n-1)*NALTMAX*NCSMAX + (k-1)*NALTMAX + 1,:)==1) && (choice((n-1)*NALTMAX*NCSMAX + k*NALTMAX,:)==0)
            tempm(k,:)=[-1 1];       % The chosen alternative must be subtracted
        else
            tempm(k,:)=[1 -1];
        end
    end
    
    Mstar(1,:)=horzcat(tempm(1,:),zeros(1,4));
    Mstar(2,:)=horzcat(zeros(1,2),tempm(2,:),zeros(1,2));
    Mstar(3,:)=horzcat(zeros(1,4),tempm(3,:));                  % Size of Mstar is NCSMAX X NALTMAX*NCSMAX
    
    % Prepare the M matrix
    temp011=horzcat(eye(NALTMAX*NCSMAX),zeros(NALTMAX*NCSMAX,NALTMAX*NCSMAX));
    temp012=horzcat(zeros(NCSMAX,NALTMAX*NCSMAX),Mstar);
    M=vertcat(temp011,temp012);     % Size of M is NCSMAX+NALTMAX*NCSMAX X 2*NALTMAX*NCSMAX
    
    % Prepare the undifferenced error matrix
    sigma=zeros(NALTMAX*NCSMAX,NALTMAX*NCSMAX);
    sigma(2,2)=1;
    sigma(4,4)=1;
    sigma(6,6)=1;
    
    % Prepare the overall covariance matrix
    % Measurement Part
    vcm11 = eye(NALTMAX*NCSMAX) + repmat(Fac,NALTMAX*NCSMAX,1)*ones(nfl,nfl)*repmat(Fac,1,NALTMAX*NCSMAX);  % Size is (NALTMAX*NCSMAX) X (NALTMAX*NCSMAX)
    % Measurement-Choice Parts
    vcm12 = repmat(Fac,NALTMAX*NCSMAX,1)*ones(nfl,nfl)*lamda';  % Size is (NALTMAX*NCSMAX) X (NALTMAX*NCSMAX)
    vcm21 = lamda*ones(nfl,nfl)*repmat(Fac,1,NALTMAX*NCSMAX);   % Size is (NALTMAX*NCSMAX) X (NALTMAX*NCSMAX)
    % Choice part
    vcm22 = sigma + lamda*ones(nfl,nfl)*transpose(lamda);     % Size of vector is NALTMAX*NCSMAX X NALTMAX*NCSMAX

    vcmatrix=[vcm11,vcm12;vcm21,vcm22];         % Size of vcmatrix is 2*NALTMAX*NCSMAX X 2*NALTMAX*NCSMAX
    
    % Evaluate the overall mean vector of 9 observations
    B = M*yvector;  % Size of the vector is NCSMAX+NALTMAX*NCSMAX X 1
    
    % Evaluate the covariance matrix of 9 observations
    delta = M*vcmatrix*M';  % Size of vcmatrix is NCSMAX+NALTMAX*NCSMAX X NCSMAX+NALTMAX*NCSMAX
    
    
    
    % ---------------------------------------------------------------------------
    %                             MEASUREMENT PART
    % ---------------------------------------------------------------------------
    
    % Populate the lower threshold vector
    LT=LowInd(startrowind:endrowind,1);     % Size of the vector is NALTMAX*NCSMAX X 1
    LT=vertcat(LT,zeros(NCSMAX,1));         % Increasing the size of the matrix to NCSMAX+NALTMAX*NCSMAX X 1
    % Populate the upper threshold vector
    UT=UpInd(startrowind:endrowind,1);      % Size of the vector is NALTMAX*NCSMAX X 1
    UT=vertcat(UT,zeros(NCSMAX,1));         % Increasing the size of the matrix to NCSMAX+NALTMAX*NCSMAX X 1
    
    % We now evaluate pairs of joint probilities using the COMPOSITE MARGINAL LIKELIHOOD (CML) technique
    % Making 6C2 = 15 bivariate pairs of the joint 6 dimensional rating intergral
    
    % ----- COMPOSITE MARGINAL LIKELIHOOD (CML) (BIVARIATE) -----
    
    pairnumber=1;       % Count the number of bivatriate pairs
    
    for i=1:(NALTMAX*NCSMAX)-1
        for j=i+1:NALTMAX*NCSMAX
            
            % Prepare an R matrix which will pick the necessary pairs accordingly
            R=zeros(2,NCSMAX+NALTMAX*NCSMAX);  % The reason it is 2 is because it is a bi-variate CML
            R(1,i)=1;   % Setting ith column of the 1st row as 1
            R(2,j)=1;   % Setting jth column of the 2nd row as 1
            
            if rem(i,2)==0
               iind=2; 
            else
               iind=1; 
            end   
            
            if rem(j,2)==0
                jind=2;
            else
                jind=1;
            end    
            
            ystarsel=R*B;   % The propensity equations of the selected pair
            ystarsel=transpose(ystarsel);   % Size of the vector is 1 X 2
            ystarsel=repmat(ystarsel,2,1);      % Size of the vector is 2 X 2
            LTsel=transpose(R*LT);         % Lower threshold of the selected pair. Size of the vector is 1 X 2
            threshLTsel(1,1)=threshvec(iind,LTsel(1,1));
            threshLTsel(1,2)=threshvec(jind,LTsel(1,2));
            UTsel=transpose(R*UT);         % Upper threshold of the selected pair. Size of the vector is 1 X 2
            threshUTsel(1,1)=threshvec(iind,UTsel(1,1));
            threshUTsel(1,2)=threshvec(jind,UTsel(1,2));
            threshvecsel=vertcat(threshLTsel,threshUTsel);      % Size of the vector is 2 X 2
            
            E=R*delta*transpose(R);      % Size of E is 2 X 2. This is the VCM of the selected pair
            
            % Taking square root of the diagonal elements of vcmatrix
            SQ11=sqrt(E(1,1));
            SQ22=sqrt(E(2,2));
            SQ11m=repmat(SQ11,2,1);
            SQ22m=repmat(SQ22,2,1);
            SQ=horzcat(SQ11m,SQ22m);  % Size of the matrix is 2 X 2
            % Calculating coeff of correlation from this vcmatrix
            r=E(2,1)/(SQ11*SQ22);      % It is a scalar
            
            % Normalise the Observed Propensity
            A=(threshvecsel - ystarsel)./SQ;    % Size of the matrix is 2 X 2
            
            % Calculate the bivariate probability distribution of the selected pair
            pp=bvn(-Inf, A(2,1), -Inf, A(2,2), r) - bvn(-Inf, A(1,1), -Inf, A(2,2), r) - bvn(-Inf, A(2,1), -Inf, A(1,2), r) + bvn(-Inf, A(1,1), -Inf, A(1,2), r);
            
            % Store bivariate probability into PM Matrix. Size of the matrix is
            % (NALTMAX*NCSMAX)C2 X NP
            PM(pairnumber,n)=pp;
            pairnumber=pairnumber+1;
        end
    end
    
    
    
    % ---------------------------------------------------------------------------
    %                               CHOICE PART
    % ---------------------------------------------------------------------------  

    % Extract the covariance matrix of choice part from delta
    G=delta(NALTMAX*NCSMAX+1:end,NALTMAX*NCSMAX+1:end);
    
    % Calculate the Standardised Utility
    Us=-1.*B(NALTMAX*NCSMAX+1:end,1)./diag(sqrt(G));    % Size of vector is NCSMAX X 1
    Ustvn=transpose(Us);            % Size of vector is 1 X NCSMAX
    
    % Find the correlation matrix of G. This will be used in the trivariate normal function
    CorG=corrcov(G);    % Size of vector is NCSMAX X NCSMAX
    
    % Prepare the correlation coefficient vector for tvn matlab function.. Mind the Order!!
    CrCoef(1,1)=CorG(2,1);      % Corr coeff for Row 2 column 1
    CrCoef(1,2)=CorG(3,1);      % Corr coeff for Row 3 column 1
    CrCoef(1,3)=CorG(3,2);      % Corr coeff for Row 3 column 2
    
    % Evaluate the joint probability of observing all three route choices simultaneously
    PC(1,n)=tvn([-Inf, -Inf, -Inf], Ustvn, CrCoef);   % Size of the vector is 1 X NP
    
    
    % ---------------------------------------------------------------------------
    %                           MEASUREMENT CHOICE PART
    % ---------------------------------------------------------------------------
    
    MCpairnumber=1;
    
    for i=1:(NALTMAX*NCSMAX)
       for j=1:(NALTMAX-1)*NCSMAX 
                      
           % Prepare an RMC matrix which will pick the necessary pairs accordingly
            RMC=zeros(2,NCSMAX+NALTMAX*NCSMAX);  % The reason it is 2 is because it is a bi-variate CML
            RMC(1,i)=1;   % Setting ith column of the 1st row as 1
            RMC(2,NALTMAX*NCSMAX+j)=1;   % Setting jth column of the 2nd row as 1
            
            meanMC=RMC*B;   % The propensity equations of the selected pair
            meanMC=transpose(meanMC);   % Size of the vector is 1 X 2
            meanMC=repmat(meanMC,2,1);      % Size of the vector is 2 X 2
            if rem(i,NALTMAX)==0
               index=2; 
            else
               index=1; 
            end    
            MCLTsel=transpose(RMC*LT);         % Lower threshold of the selected pair. Size of the vector is 1 X 2
            MCthreshLTsel(1,1)=threshvec(index,MCLTsel(1,1));
            MCthreshLTsel(1,2)=-Inf;
            MCUTsel=transpose(RMC*UT);         % Upper threshold of the selected pair. Size of the vector is 1 X 2
            MCthreshUTsel(1,1)=threshvec(index,MCUTsel(1,1));
            MCthreshUTsel(1,2)=0;
            MCthreshvecsel=vertcat(MCthreshLTsel,MCthreshUTsel);      % Size of the vector is 2 X 2
                        
            MCE=RMC*delta*transpose(RMC);      % Size of E is 2 X 2. This is the VCM of the selected pair
            
            % Taking square root of the diagonal elements of vcmatrix
            MCSQ11=sqrt(MCE(1,1));
            MCSQ22=sqrt(MCE(2,2));
            MCSQ11m=repmat(MCSQ11,2,1);
            MCSQ22m=repmat(MCSQ22,2,1);
            MCSQ=horzcat(MCSQ11m,MCSQ22m);  % Size of the matrix is 2 X 2
            % Calculating coeff of correlation from this vcmatrix
            mcr=MCE(2,1)/(MCSQ11*MCSQ22);      % It is a scalar
            
            % Normalise the Observed Propensity
            MCA=(MCthreshvecsel - meanMC)./MCSQ;    % Size of the matrix is 2 X 2
            
            % Calculate the bivariate probability distribution of the selected pair
            mcpp=bvn(-Inf, MCA(2,1), -Inf, MCA(2,2), mcr) - bvn(-Inf, MCA(1,1), -Inf, MCA(2,2), mcr) - bvn(-Inf, MCA(2,1), -Inf, MCA(1,2), mcr) + bvn(-Inf, MCA(1,1), -Inf, MCA(1,2), mcr);
            
            % Store bivariate probability into PM Matrix. Size of the matrix is
            % (NALTMAX*NCSMAX)C2 X NP
            PI(MCpairnumber,n)=mcpp;
            MCpairnumber=MCpairnumber+1;
       end 
    end
end


% Taking Ln of PM and PC vectors
PM=log(PM);     % Size of the matrix is (NALTMAX*NCSMAX)C2 X NP
PM=sum(PM,1);   % Size of the matrix is 1 X NP
PC=log(PC);     % Size of the matrix is 1 X NP
PI=log(PI);     % Size of the matrix is (NALTMAX*NCSMAX)*(NALTMAX-1)*NCSMAX X NP
PI=sum(PI,1);   % Size of the matrix is 1 X NP


% --------------- EVALUATE THE JOINT PROBABILITY ----------------

% Jointly observing three set of measurement indicators and three choices from an individual
ptot=PM+PC+PI;
ptot=transpose(ptot);       % Size of the vector is NP X 1
ll=-sum(ptot,1);            % Size of the vector is 1 X 1. Extra -ve to make it into convex function for fminunc in Matlab

% % Check if the probability is too small which might give negative infinity
% % upon taking log
% ptot(ptot < 10^-20)=10^-20;   
% ptot(1,isnan(ptot))=1; %Change missing values to 1, as a precaution.
% ll=-sum(log(ptot),1);  % Sum of log-probabilities of all individuals... Extra -ve to make it into convex function for fminunc in Matlab

% Since log has already been taken, the ln(probability) just needs to be
% ADDED. Extra -ve to make it into convex function for fminunc in Matlab
% ll=-sum(ptot,2);    % ll is now a scalar

g=zeros(NS+NALTMAX+nfl+NALTMAX*((RAT-1)-1),NP);  % This line is just to give some value to g. Ideally this line should be inside llgrad2
g=-sum(g,2);        % Sum of gradient of all individuals... Extra -ve to make it into convex function for fminunc in Matlab



