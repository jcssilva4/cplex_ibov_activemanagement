%function [w,Sfinal, cpx] = CPLEX_classical(P, lambda) 
function [w,Sfinal] = CPLEX_classical(P, lambda) 
    %original
    %get general parameters
    assets = P.initialAssetList;
    nAssets = length(P.initialAssetList);
    nVars = 2*nAssets;
    %get objfunction (H, f)
    H = zeros(nVars);  %initialize coeff matrix:
    COV = getCovMat(P);
    for i = 1:nAssets
        for j = 1:nAssets
            H(i,j) = COV(i,j);
        end
    end
    f = zeros(1, nVars); %initialzie linear coeff vec: 
    for i=1:nAssets  
        f(i) = -(lambda*assets{1,i}.alphaScore); 
    end
    f = f';
    %get linear constraints (Aineq, bineq, Aeq, beq)
    %sum of ws = 1
    Aeq = zeros(1, nVars); Aeq(1,1:nAssets) = 1;
    beq = 1;
    [Aineq, bineq] = getLinCon_Classical(P);
    bineq = bineq';
    %get boundaries
    lb = zeros(1,nVars);
    ub = (ones(1,length(lb)));  
    %get BISCN indicator
    ctype = 'C'; %continuous vars
    for i = 1 : nAssets-1
        ctype = strcat(ctype,'C');
    end
    for i = nAssets + 1 : nVars %bin vars
        ctype = strcat(ctype,'B');
    end
    %set options
    options = cplexoptimset('cplex');
    options.emphasis.mip = 2;
    %options.mip.pool.absgap = 0;
    %options.mip.pool.intensity = 4;
    options.optimalitytarget = 1;
    %options.mip.strategy.dive = 3;
    %options.mip.limits.populate = 2100000000;
    %options.mip.strategy.heuristicfreq = -1;
    %options.mip.strategy.lbheur = 1;
    %options.mip.strategy.nodeselect = 2;
    options.mip.tolerances.absmipgap = 0;
    options.display = 'on';
    %options.mip.strategy.search = 1;
    %options.MaxTime =  10; %in seconds, we need this to apply constraint 10
    [x,fval,exitflag,output] = cplexmiqp(H, f, Aineq, bineq, Aeq, beq, ...
        [], [], [], lb, ub, ctype, [], options)
%}
%{
%optimize via cplex class
    %get general parameters
    assets = P.initialAssetList;
    nAssets = length(P.initialAssetList);
    nVars = 2*nAssets;
    %get objfunction (H, f)
    H = zeros(nVars);  %initialize coeff matrix:
    COV = getCovMat(P);
    for i = 1:nAssets
        for j = 1:nAssets
            H(i,j) = COV(i,j);
        end
    end
    f = zeros(1, nVars); %initialzie linear coeff vec: 
    for i=1:nAssets  
        f(i) = -(lambda*assets{1,i}.alphaScore); 
    end
    f = f';
    %get linear constraints (Aineq, bineq, Aeq, beq)
    %sum of ws = 1
    Aeq = zeros(1, nVars); Aeq(1,1:nAssets) = 1;
    beq = 1;
    [A, lhs, rhs] = getLinCon_Classical(P);
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
    %options = cplexoptimset('cplex');
    %options.emphasis.mip = 2;
    %options.mip.pool.absgap = 0;
    %options.mip.pool.intensity = 4;
    %options.optimalitytarget = 1;
    %options.mip.strategy.dive = 3;
    %options.mip.strategy.heuristicfreq = -1;
    %options.mip.strategy.lbheur = 1;
    %options.mip.strategy.nodeselect = 2;
    %options.mip.tolerances.absmipgap = 1e-12;
    %options.mip.pool.intensity = 4;
    %options.display = 'on';
    %options.mip.strategy.search = 1;
    %options.MaxTime =  10; %in seconds, we need this to apply constraint 10
    
    cpx = Cplex();
    cpx.Model.Q = H;
    cpx.Model.obj = f;
    cpx.Model.lb = lb';
    cpx.Model.ub = ub';
    cpx.Model.A = A;
    cpx.Model.lhs = lhs';
    cpx.Model.rhs = rhs';
    cpx.Model.sense = 'minimize';
    cpx.Model.ctype = ctype;
    
    cpx.Param.timelimit.Cur = 60;
    cpx.Param.optimalitytarget.Cur = 1;
    %cpx.Param.mip.pool.intensity.Cur = 4;
    cpx.Param.mip.limits.populate.Cur = 210000;
    %cpx.Param.WorkMem.Cur = 3024; %3GB
    %cpx.Param.mip.pool.absgap.Cur = 0;
    cpx.Param.mip.tolerances.absmipgap.Cur = 0;
    %cpx.Param.absmipgap.Cur = 1e-09;
    %cpx.populate();
    solstruct = cpx.solve();
    x = cpx.Solution.x(1:nVars);
    %}
     %%get results
    Sfinal = []; %assets included
    for i = nAssets+1:nVars
        if(x(i)), Sfinal = [Sfinal i-nAssets]; end
    end
    w = zeros(1, length(Sfinal));
    for i = 1 : nAssets
        if(ismember(i,Sfinal))
            w(i) = x(i);
        end
    end
%}
end