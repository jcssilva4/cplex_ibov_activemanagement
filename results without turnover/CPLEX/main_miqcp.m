function [x,fval,exitflag] = main_miqcp(P, lambda)
 %{
    [x,fval,exitflag,output]=cplexmiqcp(H, f, Aineq, bineq, Aeq, beq,...
        l, Q, r, sostype, sosind, soswt, lb, ub, ctype, x0, options)

       min      0.5*x'*H*x + f*x or f*x
       st.      Aineq*x      <= bineq
                Aeq*x         = beq
                l*x + x'*Q*x <= r
                lb <= x <= ub
 %}
%{
%WITH BINVARS
    %get general parameters
    nAssets = length(P.initialAssetList);
    nVars = 2*nAssets;
    %get objfunction (H, f)
    [H, f] = getObjFunc(P, lambda);
    %get linear constraints (Aineq, bineq, Aeq, beq)
    %sum of ws = 1
    Aeq = zeros(1, nVars); Aeq(1,1:nAssets) = 1;
    beq = 1;
    [Aineq, bineq] = getLinCon(P);
    bineq = bineq';
    %get quadratic constraints (l, Q, r)
    [l, Q, r] = getNonlinCon(P, lambda, H, f);
    %get boundaries
    lb = zeros(1,nVars);
    ub = ones(1,length(lb));
    %get BISCN indicator
    ctype = 'S'; %continuous vars
    for i = 1 : nAssets-1
        ctype = strcat(ctype,'S');
    end
    for i = nAssets + 1 : nVars %bin vars
        ctype = strcat(ctype,'B');
    end
    %set options
    options = cplexoptimset;
    options.Display = 'on';
    options.MaxNodes = 2;
    %options.MaxTime =  10; %in seconds, we need this to apply constraint 10
    %solve
    %{
    [x,fval,exitflag,output]=cplexmiqcp(H, f, Aineq, bineq, Aeq, beq,...
        l, Q, r, sostype, sosind, soswt, lb, ub, ctype, x0, options)
    %}
    %x = 0; fval = 0; exitflag =0;
    x0 = (getx0(P))'; %get benchmark as a initial point
    %{
    x0 = getx0(P);
    MaxIters = 10;
    tempAineq = 0; %auxiliar A
    tempbineq = 0; %auxiliar b
    TE = zeros(1, MaxIters); %tracking error
    for (iter = 1 : MaxIters)
        fprintf("\niter: %d", iter)
        TE(iter) = getTE(x0, P);
        if(iter > 1)    x0 = x; end %start from the last obtained solution
        %{
        [Acon10,bcon10] = getcon10(x0, P); %get constraint 10
        tempAineq = vertcat(Aineq, Acon10); %joint all constraints with constraint 10
        tempbineq = vertcat(bineq, bcon10); %joint all constraints with constraint 10
        size(tempAineq)
        size(tempbineq)
        [x,fval,exitflag,output] = cplexmiqcp(H, f, tempAineq, tempbineq, Aeq, beq, ...
            [],[], [], [], [], [], lb, ub, ctype, x0, options);
        %}
          [x,fval,exitflag,output] = cplexmiqcp(H, f, Aineq, bineq, Aeq, beq, ...
            [],[], [], [], [], [], lb, ub, ctype, [], options);
    end
    for(i = 1:MaxIters)
        fprintf("\niter: %d\tTE: %f\n", i, TE(i));    
    end
    %}
    [x,fval,exitflag,output] = cplexmiqcp(H, f, Aineq, bineq, Aeq, beq, ...
            [],[], [], [], [], [], lb, ub, ctype, [], options);
    %maxCOV = max(max(H))
    TE = getTE(x, P)
        
    fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
  %}  
    %get general parameters
    nAssets = length(P.initialAssetList);
    nVars = nAssets;
    %get objfunction (H, f)
    [H, f] = getObjFunc(P, lambda);
    f = f';
    %get linear constraints (Aineq, bineq, Aeq, beq)
    %sum of ws = 1
    Aeq = zeros(1, nVars); Aeq(1,1:nAssets) = 1;
    beq = 1;
    [Aineq, bineq] = getContLinCon(P);
    bineq = bineq';
    %get quadratic constraints (l, Q, r)
    [l, Q, r] = getNonlinCon(P, lambda, H, f);
    %get boundaries
    lb = zeros(1,nVars);
    ub = ones(1,length(lb));
    %set options
    options = cplexoptimset;
    options.Display = 'on';
    %options.MaxNodes = 2;
    %options.MaxTime =  10; %in seconds, we need this to apply constraint 10
    %solve
    x0 = (getx0(P))'; %get benchmark as a initial point   
    %{
    %MIQCQP
    [x,fval,exitflag,output] = cplexmiqcp(H, f, Aineq, bineq, Aeq, beq, ...
            l,Q, r, [], [], [], lb, ub, ctype, [], options);
    %}
    %QCQP
    [x, fval, exitflag, output]= cplexqcp(H,f,Aineq,bineq,Aeq,beq,l,Q,r,lb,ub,[],options);
    %maxCOV = max(max(H))
    TE = getTE(x, P)
    fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
end