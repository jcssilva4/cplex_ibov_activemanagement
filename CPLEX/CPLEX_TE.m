function [w,Sfinal] = CPLEX_TE(P, lambda)
    %optimize via cplex class
    %get general parameters
    assets = P.initialAssetList;
    nAssets = length(P.initialAssetList);
    nVars = 3*nAssets; %di, wi and sigma 
    %[di-> 1:nAssets, wi-> nAssets+1:2*nAssets, sigma_i-> 2*nAssets+1:3*nAssets]
    %get objfunction (H, f)
    H = zeros(nVars);  %initialize coeff matrix:
    COV = getCovMat(P);
    for i = 1:nAssets
        for j = 1:nAssets
            H(i,j) = COV(i,j);
        end
    end
    f = zeros(1, nVars); %initialize linear coeff vec: 
    %get alpha-> vector of expected excess returns
    expectedBenchReturn = 0;
    for i=1:nAssets  
        expectedBenchReturn = (assets{1,i}.benchWeight*assets{1,i}.alphaScore); 
    end
    %get fobj part2
    for i=1:nAssets  
        f(i) = -(lambda*((assets{1,i}.alphaScore)-expectedBenchReturn)); 
    end
    f = f';
    %get linear constraints (Aineq, bineq, Aeq, beq)
    %sum of ws = 1
    Aeq = zeros(1+length(P.initialAssetList), nVars);
    beq = zeros(1,1+length(P.initialAssetList));
    %Aeq(1,1:nAssets) = 1; %sumdi = 0
    Aeq(1,nAssets+1:2*nAssets) = 1; %sumwi = 1
    beq(1) = 1; %sumwi = 1
    % -di + wi = wbi
    for i = 1:nAssets
        Aeq(i+1,i) = -1; %-di
        Aeq(i+1,i+nAssets) = 1;%+wi
        beq(i+1) = P.initialAssetList{1,i}.benchWeight;
    end
    beq = beq'; %avoids: Dimensions of matrices being concatenated are not consistent.
    [A, b] = getLinCon_TE(P);
    b = b';%avoids: Dimensions of matrices being concatenated are not consistent.
    %get boundaries
    lb = zeros(1,nVars); lb(1:nAssets) = -0.1; 
    ub = ones(1,length(lb)); ub(1:nAssets) = 0.1;  %ub(nAssets+1:2*nAssets) = 0.15;
    %get BISCN indicator
    ctype = 'C'; %continuous vars
    for i = 1 : (2*nAssets)-1
        ctype = strcat(ctype,'C');
    end
    for i = (2*nAssets) + 1 : nVars %bin vars
        ctype = strcat(ctype,'B');
    end
    %set options
    options = cplexoptimset('cplex');
    %options.emphasis.mip = 2;
    %options.mip.pool.absgap = 0;
    %options.mip.pool.intensity = 4;
    options.optimalitytarget = 1;
    %options.mip.strategy.dive = 3;
    %options.mip.strategy.heuristicfreq = -1;
    %options.mip.strategy.lbheur = 1;
    %options.mip.strategy.nodeselect = 2;
    options.mip.tolerances.absmipgap = 0;
    %options.mip.pool.intensity = 4;
    options.display = 'on';
    %options.mip.strategy.search = 1;
    %options.MaxTime =  10; %in seconds, we need this to apply constraint 10
    [x,fval,exitflag,output] = cplexmiqp(H, f, A, b, Aeq, beq, ...
        [], [], [], lb, ub, ctype, [], options)
    %{
    %%show results
    fprintf('\n\nfval: %f', fval);
    fprintf('\noutput: '); output
    fprintf('\n\t\tdi\t\t\t\twi\t\t\tsigma')
    for i = 1:nAssets
        fprintf('\n\t%f\t\t%f\t\t%d', x(i), x(i+nAssets), x(i+(2*nAssets)));
    end
    fprintf('\n\nsum(di)\t\tsum(wi)\t\t(number of stocks)');
    fprintf('\n\n%d\t\t%d\t\t%d', sum(x(i:nAssets)), ...
        sum(x(nAssets+1:2*nAssets)),sum(x((2*nAssets)+1:nVars)));
        %}
        
    %%get results
    Sfinal = []; %assets included
    for i = (2*nAssets)+1:nVars
        if(x(i)), Sfinal = [Sfinal i-(2*nAssets)]; end
    end
    w = zeros(1, length(Sfinal));
    for i = 1:nAssets
        %x(i+nAssets)
        if(ismember(i, Sfinal) >= 0.001)
            w(i) = x(i+nAssets);
        end
    end
end