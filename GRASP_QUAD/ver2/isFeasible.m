function feasibleSol = isFeasible(itrGRASP, xtemp, fvaltemp, exitflag, P, S) %return 1 if feasible, 0 if not feasible
    feasibleSol = 0;
    if(exitflag>0) %viable linear solution!
       % TE = getTE(xtemp, P, S);
        %if(TE <= 0.1 && TE >= 0.05)
         %   AS = getActiveShare(xtemp, P, S);
          %  if(AS <= 0.4)
                %fprintf('ITER: %d\t FEASIBLE SOL FOUND!', itrGRASP);
               % fprintf('\t0.05 <= TE = %f <= 0.1\tAS = %f <= 0.4\t(fval = %f)\n', TE, AS, fvaltemp);
                feasibleSol = 1;
         %   else
          %      fprintf('ITER: %d\t INFEASIBLE POINT: AS = %f\n', itrGRASP, AS);
        %    end
       % else
         %  fprintf('ITER: %d\t INFEASIBLE POINT: TE = %f\n', itrGRASP, TE);
        %end
    else
       fprintf('ITER: %d\t INFEASIBLE POINT: DO NOT OBEY LINEAR CONSTRAINTS\n', itrGRASP);
    end

end