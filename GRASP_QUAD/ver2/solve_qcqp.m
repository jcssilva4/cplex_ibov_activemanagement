function [x, fval, exitflag] = solve_qcqp(Portfolios, lambda, S, COVIJ, numOfPastPeriods, t)
 %{
       min      0.5*x'*H*x + f*x or f*x
       st.      Aineq*x      <= bineq
                Aeq*x         = beq
                l*x + x'*Q*x <= r
                lb <= x <= ub
 %}  
%FIRST: you need to get an index mapping to search the K stocks in P
    %now the problem is reduced from |P.initialList|x|P.initialList| to |K|x|K|
    %get general parameters
    %we are using di as variable
    nVars = 2*length(S); %two type of vars di(1:length(assets)) and wi(length(assets)+1:nVars)
    P = Portfolios{1, t};
    %get objfunction (H, f)
    COVmat = (getCovMat(P, S, COVIJ));%reduced COVIJ
    H = zeros(nVars);
    for i = 1:length(S)
        for j = 1:length(S)
            H(i,j) = COVmat(i,j);
        end
    end
    f = getObjF_qp(P, lambda, S, H);
    f = f';
    %get linear constraints (Aineq, bineq, Aeq, beq)
    %sum(wi) = 1
    Aeq = zeros(1+length(S), nVars);
    beq = zeros(1,1+length(S));
    Aeq(1,length(S)+1:nVars) = 1;
    beq(1) = 1;
    % -di + wi = wbi
    for i = 1:length(S)
        Aeq(i+1,i+1) = -1; %-di
        Aeq(i+1,i+length(S)) = 1;%+wi
        beq(i+1) = P.initialAssetList{1,S(i)}.benchWeight;
    end
    %[Aineq, bineq] = get_qplincon(P);
    beq = beq'; %avoids: Dimensions of matrices being concatenated are not consistent.
    %get boundaries
    lb = 0.001*zeros(1,nVars); lb(1:length(S)) = -0.05;
    ub = ones(1,nVars); ub(1:length(S)) = 0.05;
    [Aineq, bineq] = getLinCon_qcqp(Portfolios, S, t, numOfPastPeriods);
    bineq = bineq';
    %get quadratic constraints (l, Q, r)
    [~, Q, r] = getNonlinCon_qcqp(Portfolios, COVmat, S, t-numOfPastPeriods, t-1);
    l = 0; %you need to initialize l, otherwise cplex will understand that l, Q and r are not consistent
    %set options
    options = cplexoptimset;
    options.Display = 'off';
    %solve QP
    [x, fval, exitflag, ~]= cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub,[],options);
    
%{
ORIGINAL FORMULATION
    %FIRST: you need to get an index mapping to search the K stocks in P
    %now the problem is reduced from |P.initialList|x|P.initialList| to |K|x|K|
    %get general parameters
    nVars = length(S);
    %get objfunction (H, f)
    H = getCovMat(P, S, COVIJ);%reduced COVIJ
    f = getObjF_qcqp(P, lambda, S, H, nVars);
    f = f';
    %get linear constraints (Aineq, bineq, Aeq, beq)
    %sum of ws = 1
    Aeq = ones(1, nVars); 
    beq = 1;
    [Aineq, bineq] = getLinCon_qcqp(P, S);
    bineq = bineq';
    %get quadratic constraints (l, Q, r)
    [l, Q, r] = getNonlinCon_qcqp(P, lambda, S, H, f);
    %get boundaries
    %lb: wi >= -0.05 + wb (r5.1)  (but wi>=0.001) because it is included in the portfolio (GRASP SOLUTION)
    lb = zeros(1,nVars);
    ub = ones(1,length(lb));
    inclusion_lb = 0.001;
    for i = 1:nVars
        candidatelb = -0.05 + P.initialAssetList{1,S(i)}.benchWeight;  
        if(candidatelb < inclusion_lb), candidatelb = inclusion_lb;  end %if candidatelb is less than inc_lb, then we use inc_lb!!
        lb(i) = candidatelb;
    end
    %ub: wi <= 0.05 + wb (r5.2)  
    for i = 1:nVars, ub(i) = 0.05 + P.initialAssetList{1,S(i)}.benchWeight; end
    %set options
    options = cplexoptimset;
    options.Display = 'off';
    %options.MaxNodes = 2;
    %options.MaxTime =  10; %in seconds, we need this to apply constraint 10
    %solve
    x0 = (getx0(P))'; %get benchmark as a initial point   
    %QCQP
    [x, fval, exitflag, output]= cplexqcp(H,f,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,[],options);
    %maxCOV = max(max(H))
    
    %fprintf ('\nSolution status = %s \n', output.cplexstatusstring)
  %}  
end