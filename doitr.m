% Create Global variables to use in estimation
cp=XMAT(:,1); % person
% Create Choice column (4th column)
choice=XMAT(:,4);
% choice=XMAT(:,5);
% Create Structure columns
XS=XMAT(:,IDS(:,1));    % Size of vector is NROWS X NS
% Create I columns
I=XMAT(:,IDIND(:,1));   % Size of vector is NROWS X NIND
% Create Structure columns
XC=XMAT(:,IDC(:,1));    % Size of vector is NROWS X NCH

nn=zeros(NCS,1);
for n=1:NCS;
    nn(n,1)=sum(XMAT(:,2) == n,1);  %This condition counts the number of lines that have the same choice taskid
end;
NALTMAX=max(nn);    % This gives the number of alternatives 

% Total number of alternatives faced by all people in all choice situations combined.
% This is the number of rows of data in XMAT below.
NROWS=NCS*NALTMAX;

% Preparing indices for lower and upper threshold values
LowInd=I;
UpInd=I+ones(NROWS,1);

nn=zeros(NP,1);
for n=1:NP;
   k=(XMAT(:,1)==n);
   k=XMAT(k,2);
   nn(n,1)=1+k(end,1)-k(1,1);
end;
NCSMAX=max(nn);  %Maximum number of choice situations faced by any person

% --------------- NAME THE PARAMETERS TO BE ESTIMATED ---------------------
% Names of parameters for the study are listed in this array
% Mu_ij means the jth threshold value for indicator i
% NAMEST = {'A_TTS' 'A_SnGo' 'A_Male' 'A_AgeLT40' 'A_ExpLT8'};
NAMEST = {'A_TTS' 'A_SnGo' 'A_AgeLT40' 'A_ExpLT8'};
NAMEME = {'Bc_1' 'Bc_2' 'LV_Load' 'Mu_12' 'Mu_13' 'Mu_14' 'Mu_22' 'Mu_23' 'Mu_24'};
NAMECH = {'C_TT' 'C_VRC' 'C_LV'};

%--------------------------------------------------------------------------

% Set of starting values

% Structural coefficients
a=[0.1 0.1 -0.1 -0.1];

% Measurement coefficients (constant)
bc1=[0.5];
bc2=[0.5];
bc=horzcat(bc1,bc2);
% Latent construct Factor Loading
fl=[0.1];
nfl=size(fl,2);
% Threshold coefficients.. Number of mu arrays should be equal to NIND
mu1=[-0.1 -0.1 -0.1];
mu2=[-0.1 -0.1 -0.1];
mu=horzcat(mu1,mu2);

% Choice coefficients
% chc=[-0.1 -0.1 -0.1];
chc=[-0.1 -0.1 -0.1];
nchc=size(chc,2);

% Concatenate the coefficients horizontally. 
% NOTE: Right now RAT1 = RAT2... If that doesn't hold then taking Max of
% the two quanitities would be better.
param=horzcat(a,bc,fl,mu,chc);

% Adding starting values of 61 participant analysis 
% param=[1.0691 0.1905 -0.7905 0.4433 -2.5504 -2.1995 1.0513 0.4105 0.2207 0.4643 0.2413 0.0523 0.3768 -0.2140 -0.2313 -0.4290]; %   third last element

% BOUNDS.... Should be a vector of size 1 X [Number of Xs + Number of Ys + sigma]. It
% comprises of two components, Lower and Upper bounds
lb=[-100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100 -100];
ub=[100 100 100 100 100 100 100 100 100 100 100 100 100 100 100 100];

disp('Start estimation');
disp('The negative of the log-likelihood is minimized,');
disp('which is the same as maximizing the log-likelihood.');
tic;


% ADDED NEW CODE TO PRINT BETAS AT EVERY ITERATION
% Use OutputFcn to display results from every iteration of fmincon


% Assign one set of starting values for likelihood estimation 

options=optimset('LargeScale','off','Display','iter','OutputFcn',@outfun,'GradObj','off',...
    'MaxFunEvals',10000,'MaxIter',MAXITERS,'TolX',PARAMTOL,'TolFun',LLTOL,'DerivativeCheck','off');
[paramhat,fval,exitflag,output,lamda,grad,hessian]=fmincon(@loglik,param,[],[],[],[],lb,ub,@confun,options);    %@ is a handle that passes one function into another function.. Here loglik is a function defined in loglik.m.. Read more on http://au.mathworks.com/help/matlab/matlab_prog/symbol-reference.html

disp(' ');
disp(['Estimation took ' num2str(toc./60) ' minutes.']);
disp(' ');
if exitflag == 1
    disp('Convergence achieved.');
elseif exitflag == 2
    disp('Convergence achieved by criterion based on change in parameters.');
    if size(PARAMTOL,1)>0
        disp(['Parameters changed less than PARAMTOL= ' num2str(PARAMTOL)]);
    else
        disp('Parameters changed less than PARAMTOL=0.000001, set by default.');
    end
    disp('You might want to check whether this is actually convergence.');
    disp('The gradient vector is');
    grad
elseif exitflag == 3
    disp('Convergence achieved by criterion based on change in log-likelihood value.');
    if size(PARAMTOL,1)>0
        disp(['Log-likelihood value changed less than LLTOL= ' num2str(LLTOL)]);
    else
        disp('Log-likelihood changed less than LLTOL=0.000001, set by default.');
    end
    disp('You might want to check whether this is actually convergence.');
    disp('The gradient vector is');
    grad
else
    disp('Convergence not achieved.');
    disp('The current value of the parameters and hessian');
    disp('can be accesses as variables paramhat and hessian.');
    disp('Results are not printed because no convergence.');
    return
end

disp(['Value of the log-likelihood function at convergence: ' num2str(-fval)]);

% Adding the logic to evaluate AIC and BIC values wrt the final
% log-likelihood value and the number of parameters
numParam=length(paramhat);  % Determining the number of estimated parameters
[aic,bic] = aicbic(-1*fval,numParam,NROWS);   % -fval is the final log-likelihood value.... BIC calculation requires the sample size too... Lower the AIC/BIC value, better is the model
disp(' ');
disp('                      AIC                      BIC');
disp('              ------------------   -----------------------');
fprintf('                  %10.4f               %10.4f\n', aic,bic);
disp(' ');

%Calculate standard errors of parameters
disp(' ');
disp('Taking inverse of hessian for standard errors.');
disp(' ');
ihess=inv(hessian);
% ihess1 = inv(NP*hessian);
% grad1 = NP*grad;
% Jacobian = (grad1*grad1');
% covbhh = ihess*Jacobian*ihess;
% st_final = sqrt(diag(covbhh));
stderr=sqrt(diag(ihess));
disp(['The value of grad*inv(hessian)*grad is: ' num2str(grad'*ihess*grad)]);

% Finding the estimated parameters for each part

% Structural model
if nfl>0
    sthat=paramhat(1,1:NS);   % paramhat is a row vector
    stsd=stderr(1:NS,1);      % stderr is a column vector
end

% Measurement coefficients & thresholds
if NIND>0
    mehat=paramhat(1,NS+1:NS+NALTMAX+nfl+NALTMAX*((RAT-1)-1));   % paramhat is a row vector
    mesd=stderr(NS+1:NS+NALTMAX+nfl+NALTMAX*((RAT-1)-1),1);      % stderr is a column vector
end

% Choice model
chhat=paramhat(1,NS+NALTMAX+nfl+NALTMAX*((RAT-1)-1)+1:NS+NALTMAX+nfl+NALTMAX*((RAT-1)-1)+nchc);   % paramhat is a row vector
chsd=stderr(NS+NALTMAX+nfl+NALTMAX*((RAT-1)-1)+1:NS+NALTMAX+nfl+NALTMAX*((RAT-1)-1)+nchc,1);      % stderr is a column vector

% Taking transpose of standard error vectors
stsd=transpose(stsd);
mesd=transpose(mesd);
chsd=transpose(chsd);

disp(' ');
disp('-------------------------------------------------------------------');
disp('     RESULTS FROM THE INTEGRATED CHOICE & LATENT VARIABLE MODEL ');
disp('-------------------------------------------------------------------');

disp(' ');
disp(' ');
disp('STRUCTURAL MODEL');
disp(' ');
disp('                    Coeffs      ');
disp('              ------------------ ');
disp('                Est         SE       Z-Stat ');
for r=1:NS;
    fprintf('%-10s %10.4f %10.4f %10.4f\n', NAMEST{1,r}, [sthat(1,r) stsd(1,r) sthat(1,r)/stsd(1,r)]);
end

disp(' ');
disp(' ');

disp(' ');
disp(' ');
disp('MEASUREMENT MODEL');
disp(' ');
disp('                    Coeffs      ');
disp('              ------------------ ');
disp('                Est         SE       Z-Stat ');
for r=1:NALTMAX+nfl+NALTMAX*((RAT-1)-1);
    fprintf('%-10s %10.4f %10.4f %10.4f\n', NAMEME{1,r}, [mehat(1,r) mesd(1,r) mehat(1,r)/mesd(1,r)]);
end
disp(' ');
disp(' ');
disp('** The first threshold point is normalised to ZERO because a constant term is getting estimated. ');

disp(' ');
disp(' ');

disp(' ');
disp(' ');
disp('CHOICE MODEL');
disp(' ');
disp('                    Coeffs      ');
disp('              ------------------ ');
disp('                Est         SE       Z-Stat ');
for r=1:nchc;
    fprintf('%-10s %10.4f %10.4f %10.4f\n', NAMECH{1,r}, [chhat(1,r) chsd(1,r) chhat(1,r)/chsd(1,r)]);
end

disp(' ');
disp(' ');
disp('ESTIMATED VARIANCE COVARIANCE MATRIX.');
disp(' ');
disp('Measurement Model');
disp(' ');

factor=paramhat(1,NS+NALTMAX+1:NS+NALTMAX+nfl);
vcmatrix = eye(NALTMAX) + repmat(factor,NALTMAX,1)*ones(nfl,nfl)*repmat(factor,1,NALTMAX);  % Size is 2 X 2

vcmatrix

disp(' ');
disp(' ');
disp('ESTIMATED PARAMETERS AND FULL COVARIANCE MATRIX.');
disp('The estimated values of the parameters are:');
paramhat
disp('The covariance matrix for these parameters is:');
ihess

disp(' ');
disp('You can access the estimated parameters as variable paramhat,');
disp('the gradient of the negative of the log-likelihood function as variable grad,');
disp('the hessian of the negative of the log-likelihood function as variable hessian,');
disp('and the inverse of the hessian as variable ihess.');
disp('The hessian is calculated by the BFGS updating procedure.');
