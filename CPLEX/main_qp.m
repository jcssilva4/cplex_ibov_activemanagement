%{
solves unconstrained mv portfolio problem
constraints:
0 <= wi <= 1
sum(wi) = 1;
%}
function [x,fval,exitflag] = main_miqcp(P, lambda, COVIJ) 
    %get general parameters
    assets = P.initialAssetList;
    nVars = 2*length(assets);
    %get objfunction (H, f)
    H = zeros(nVars);
    for i = 1:length(P.initialAssetList)
        for j = 1:length(P.initialAssetList)
            H(i,j) = COVIJ(i,j);
        end
    end
     %using di as a variable
    %since cplex qcqp solver deals with a vectorized obj function
    %get i
    f = zeros(1, nVars); %initializie linear coeff vec: 
    for i=1:length(P.initialAssetList)
        f(i) = -(lambda*assets{1,i}.alphaScore);
    end
    f = f'; %avoids: Dimensions of matrices being concatenated are not consistent.
    %get linear constraints (Aineq, bineq, Aeq, beq)
    Aeq = zeros(1+length(P.initialAssetList), nVars);
    beq = zeros(1,1+length(P.initialAssetList));
    %sum(wi) = 1
    Aeq(1,length(P.initialAssetList)+1:nVars) = 1;
    beq(1) = 1;
    % -di + wi = wbi
    for i = 1:length(P.initialAssetList)
        Aeq(i+1,i+1) = -1; %-di
        Aeq(i+1,i+length(P.initialAssetList)) = 1;%+wi
        beq(i+1) = P.initialAssetList{1,i}.benchWeight;
    end
    %[Aineq, bineq] = get_qplincon(P);
    beq = beq'; %avoids: Dimensions of matrices being concatenated are not consistent.
    %get boundaries
    lb = zeros(1,nVars); lb(1:length(P.initialAssetList)) = -0.05;
    ub = ones(1,nVars); ub(1:length(P.initialAssetList)) = 0.05;
    %set options
    options = cplexoptimset;
    options.Display = 'off';
    %solve
    %QP
    [x, fval, exitflag, output]= cplexqp(H,f,[],[],Aeq,beq,[],[],[],options);
    %maxCOV = max(max(H))
    fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
end

%{
%original formulation
%{
solves unconstrained mv portfolio problem
constraints:
0 <= wi <= 1
sum(wi) = 1;
%}
function [x,fval,exitflag] = main_miqcp(P, lambda, COVIJ) 
    %get general parameters
    assets = P.initialAssetList;
    nVars = length(assets);
    %get objfunction (H, f)
    H = COVIJ;
     %using di as a variable
    %since cplex qcqp solver deals with a vectorized obj function
    %get i
    f = zeros(1, nVars); %initializie linear coeff vec: 
    for i=1:nVars 
        f(i) = -(lambda*assets{1,i}.alphaScore);
    end
    f = f'; %avoids: Dimensions of matrices being concatenated are not consistent.
    %get linear constraints (Aineq, bineq, Aeq, beq)
    %0 <= (wi) <= 1
    Aeq = ones(1, nVars); 
    sumbench = 1; %sum of benchmark weights
    beq = 1 - sumbench; %beq is equal to 0
    [Aineq, bineq] = get_qplincon(P);
     bineq = bineq'; %avoids: Dimensions of matrices being concatenated are not consistent.
    %get boundaries
    %lb = -0.05*(ones(1,nVars));
    %ub = 0.05*(ones(1,nVars));
    %set options
    options = cplexoptimset;
    options.Display = 'off';
    %solve
    %QP
    [x, fval, exitflag, output]= cplexqp(H,f,Aineq,bineq,Aeq,beq,[],[],[],options);
    %maxCOV = max(max(H))
    fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
end
%}