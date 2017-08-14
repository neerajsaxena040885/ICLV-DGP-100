 function stop = outfun(x,optimValues,state)
     stop = false;
 
     % NOTE: The optimValues structure contains results from the iteration.
     % The fields in this structure are as follows:
     % 1. iteration --> Iteration number
     % 2. funccount --> ???
     % 3. fval --> The log-likelihood value
     % 4. constrviolation --> Any constraints violated
     % 5. stepsize --> The step length made to the next value
     % 6. gradient --> The numerical gradient vector of parameters
     % 7. firstorderopt --> The first order optimality condition
     % 8. trustregionradius --> Radius of the Trust Region
     % 9. cgiterations --> ????
     
     switch state
         case 'iter'
             % Concatenate current point and objective function
             % value with history. x must be a row vector.
             dlmwrite('historyX.csv',horzcat(optimValues.iteration,x),'-append')
             dlmwrite('historyF.csv',horzcat(optimValues.iteration,optimValues.fval),'-append')
             dlmwrite('historyFFO.csv',horzcat(optimValues.iteration,optimValues.firstorderopt),'-append')
             dlmwrite('historySS.csv',horzcat(optimValues.iteration,optimValues.stepsize),'-append')
             dlmwrite('historyGr.csv',horzcat(optimValues.iteration,transpose(optimValues.gradient)),'-append')
         otherwise
     end
 end