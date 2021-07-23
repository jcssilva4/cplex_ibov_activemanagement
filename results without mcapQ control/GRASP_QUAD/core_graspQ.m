%% this version buys and sells stocks 
function [d,w, fval, Sfinal, COVIJ, covij] = core_graspQ(P, lambda, iterMAX, P0)
%{
    Some problems:
    -CPU TIME per iter is relatively high. Create a vector Sc
     for each sector and mcapQ for random preselection? select 2*|secSET|+
     2*|mcapSET| stocks and then use a normal RCL vector, doing this you'll
     guarantee that S contains all sectors and mcapQ. You need to select 
     two stocks from each sector and mcapq in order to guarantee that MiLS will not
     draw a sector or mcapq. (DONE!)
      -We need to get better mean fval
%}
    [SInitial, ScInitial, witPRE] = getGQInitialSolution(P0, P); %get t-1 portfolio and witPRE
    Sfinal = []; %final list of portfolios
    x = [];
    fval = 1000;
    meanT = 0;% mean time per iteration
    totalT = 0;%total elapsed time
    %configure GRASP input Parameters
    itrMAX_GRASP = iterMAX;
    SIZERCLBUY = round(0.5*length(P.initialAssetList)); %fixed approximately to the half of the |P.initialAssetList|
    SIZERCLSELL = round(0.2*length(P0.initialAssetList)); %fixed approximately to the half of the |P.initialAssetList|
    %configure Minor Local Search Parameters
    minLocalSearch = 0.7; %probability of Min Local Search
    %markowitz MV problem inputs
    COVIJ = getCovMat(P);
    %Start GRASP-QUAD
    itrGRASP = 1;
    T0 = cputime; %instant where GRASP was initiated (Used to calculate totalT)     
    acumfval = 0; %used to calculate mean_fval
    %[d, ~, ~] = solve_qp(P, lambda, 1:length(P.initialAssetList), COVIJ);
    while (itrGRASP <= itrMAX_GRASP)
        a = rand(); b = rand();
        K = round(60+(a*10)-(b*10)); %size of the solution (number of stocks to be selected)% remember that: 50<=K<=70 (r9.1)
        t0 = cputime; %tic
        S = SInitial; %initialize(t-1 portfolio)the list(tuple) of the selected stocks
        Sc = ScInitial; % S complement (list of all stocks - si
        %stock selection (solution construction) phase
        %K selection: let the BEST K WINS!!!!!!! The best K have stocks in
        %every sector and mcapQ!!
        while(K ~= length(S))
            %fprintf('\n\nK -> %d, sizeP0 -> %d  >>> ', K, length(S));
            if(K > length(S)) %buy stocks with the highest benchmark weight
                fprintf('buy stocks...');
                sizeRCL = SIZERCLBUY;
                [S, Sc] = select_newstocks(P, lambda, COVIJ, sizeRCL, K, S, ScInitial);
            else
                if(K < length(S)) %sell stocks with the lowest turnover contribution
                    fprintf('sell stocks...');
                    sizeRCL = SIZERCLSELL;
                    S = remove_stocks(P, sizeRCL, K, SInitial, witPRE); %still not programmed
                end
                %AND IF K == length(SInitial) Mils will do all the work swapping stocks.. GRASP
                %don't need to buy or sell
            end
            getNumStockperSec(P, S)
        end
        mcapdiversified = 1; %assume a full diversified portfolio
        numStockperQ = getnumStocksperQ(P,S); %check mcapq diversification
        if(ismember(0, numStockperQ))
            mcapdiversified = 0;
        end
        %solve quadratic program
        S = sort(S,'ascend');
        %fprintf('\nsolving QCQP:\t(CARD CONSTRAINT SETTLED!)\n');
        [xtemp, fvaltemp, exitflag] = solve_qp(P, lambda, S, COVIJ);
        feasible = isFeasible(itrGRASP, xtemp, fvaltemp, exitflag, P, S);
        if(feasible)
            %Minor Local Search
            if(rand() <= minLocalSearch)
                %fprintf('applying MiLS...\n')
                [improved, St, xt, fvalt] = MiLS(itrGRASP, round(0.1*K),...
                    P, lambda, S, Sc, xtemp, fvaltemp, COVIJ);
                if(improved) %check if fval was improved
                    fprintf('\nMiLS improved the solution!');
                    S = St;
                    xtemp = xt; %xt is the best trial solution obtained with MiLS
                    fvaltemp = fvalt;
                end
                acumfval = acumfval + fvaltemp;
            elseif(itrGRASP == itrMAX_GRASP) %do MaLS
            end
            if(fvaltemp < fval) %compare fval temp with fval (minimize fval)
                if(mcapdiversified) %check mcapQ diversification
                    fval = fvaltemp;
                    x = xtemp;
                    Sfinal = S;
                else
                    fprintf('\nmcapQ diversification not satisfied...')
                end
            end
        end
        %iteration done!
        itrGRASP = itrGRASP + 1; 
        fprintf('\nfval -> %f (best sol -> fval = %f)       mean fval = %f', fvaltemp, fval, acumfval/(itrGRASP-1));
        toc = cputime - t0; 
        fprintf('\niteration %d elapsed time: %fs\tMeanTPI: %fs\n\n', itrGRASP-1, toc, (cputime - T0)/(itrGRASP-1));
    end
    fprintf('\ntotal elapsed time: %fs\tMeanTPI: %fs\n', cputime - T0,(cputime - T0)/(itrGRASP-1));
    dreduced = x(1:length(Sfinal));
    wreduced = x(length(Sfinal)+1:2*length(Sfinal));
    d = zeros(1,length(P.initialAssetList));
    w = zeros(1,length(P.initialAssetList));
     for i = 1:length(P.initialAssetList)
        if(ismember(i,Sfinal)) %if stock i is included in the Portfolio S
            d(i) = dreduced(find(ismember(Sfinal,i)));
            w(i) = wreduced(find(ismember(Sfinal,i)));
        end
     end
    covij = getCovMat(P, S, COVIJ);
end