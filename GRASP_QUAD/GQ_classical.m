%% GRASP QUAD for the classical model
function [d,w, fval, Sfinal, COVIJ, covij] = GQ_classical(P, lambda, iterMAX)
%{
    Some problems:
    -CPU TIME per iter is relatively high. Create a vector Sc
     for each sector and mcapQ for random preselection? select 2*|secSET|+
     2*|mcapSET| stocks and then use a normal RCL vector, doing this you'll
     guarantee that S contains all sectors and mcapQ. You need to select 
     two stocks from each sector and mcapq in order to guarantee that MiLS will not
     draw a sector or mcapq
      -We need to get better mean fval
%}
    x = [];
    fval = 1000;
    Sfinal = [];
    meanT = 0;% mean time per iteration
    totalT = 0;%total elapsed time
    %configure GRASP input Parameters
    itrMAX_GRASP = iterMAX;
    a = rand(); b = rand();
    SIZERCL = round(0.5*length(P.initialAssetList)); %fixed approximately to the half of the |P.initialAssetList|
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
        K = round(30+(a*10)-(b*10)); %size of the solution (number of stocks to be selected)% remember that: 50<=K<=70 (r9.1)
        t0 = cputime; %tic
        sizeRCL = SIZERCL;
        S = []; %initialize(empty list)the list(tuple) of the selected stocks
        Sc = 1:length(P.initialAssetList); % S complement (list of all stocks - si
        RCL = []; %Restricted Candidate List (stocks are randomly selected from this list)  
        %stock selection (solution construction) phase
        %K selection: let the BEST K WINS!!!!!!! The best K have stocks in
        %every sector and mcapQ!!
        while (length(S) < K) %K is the number of stocks to be selected (we need to vary in [50,70]) 
            if(length(S) < 1) 
                    %fprintf('\nstart stock selection with: sizeRCL = %d   sec_mcapq_threshold = %d  K = %d',sizeRCL, sec_mcapq_Threshold,K); 
            end
            %evaluate greedy values Fs for each candidate stock s (a candidate stock is not included in S)
            if(length(S) >= length(P.sectorSet)) %USE NORMAL RCL
                [Fs_pos,Fs_neg] = getGreedyVals_benchweight(P, lambda, Sc, COVIJ);%get sorted positive(Fs_pos) and negative(Fs_neg) indexes of stocks relative to the greedyvals
                if(sizeRCL <= length(Fs_pos))
                    RCL = Fs_pos(1:sizeRCL); %fcl restricted candidate list
                else
                    RCL = Fs_pos;
                    RCL = [RCL Fs_neg(1:length(Fs_neg))];
                end
                %select a random stock from RCL
                s = randsample(length(RCL),1);
                S = [S RCL(s)]; %S U s
                Sc = setdiff(Sc, RCL(s)); %Sc - s
                %fprintf('\n(%d)|S| + (%d)|Sc| -> %d   #sizeRLC -> %d',length(S), length(Sc), length(S) + length(Sc), sizeRCL);
            else %first: use a Sc that contain every sector and mcapq to ensure diversification
                %select a random stock from RCL
                [SDiversf, Sc] = getSDiversf(P, Sc,lambda, COVIJ);
                S = SDiversf; %S U s
            end
        end
        %solve quadratic program
        S = sort(S,'ascend');
        %fprintf('\nsolving QCQP:\t(CARD CONSTRAINT SETTLED!)\n');
        [xtemp, fvaltemp, exitflag] = solve_classical(P, lambda, S, COVIJ);
        if(isFeasible(itrGRASP, xtemp, fvaltemp, exitflag, P, S))
            %Minor Local Search
            if(rand() <= minLocalSearch)
                %fprintf('applying MiLS...\n')
                [improved, St, xt, fvalt] = MiLS(itrGRASP, round(0.1*K),...
                    P, lambda, S, Sc, xtemp, fvaltemp, COVIJ);
                if(improved) %check if fval was improved
                    %fprintf('MiLS improved the solution!\n');
                    S = St;
                    xtemp = xt; %xt is the best trial solution obtained with MiLS
                    fvaltemp = fvalt;
                end
                acumfval = acumfval + fvaltemp;
            elseif(itrGRASP == itrMAX_GRASP)
            end
            if(fvaltemp < fval) %compare fval temp with fval (minimize fval)
                fval = fvaltemp;
                x = xtemp;
                Sfinal = S;
            end
        end
        %iteration done!
        itrGRASP = itrGRASP + 1; 
        %fprintf('best sol -> fval = %f       mean fval = %f\n', fval, acumfval/(itrGRASP-1));
        toc = cputime - t0; 
        %fprintf('iteration elapsed time: %fs\tMeanTPI: %fs\n\n', toc, (cputime - T0)/(itrGRASP-1));
    end
    %fprintf('\ntotal elapsed time: %fs\tMeanTPI: %fs\n', cputime - T0,(cputime - T0)/(itrGRASP-1));
    wreduced = x(1:length(Sfinal));
    %wreduced = x(length(S)+1:2*length(Sfinal));
    d = zeros(1,length(P.initialAssetList));
    w = zeros(1,length(P.initialAssetList));
     for i = 1:length(P.initialAssetList)
        if(ismember(i,Sfinal)) %if stock i is included in the Portfolio S
            %d(i) = dreduced(find(ismember(Sfinal,i)));
            w(i) = wreduced(find(ismember(Sfinal,i)));
        end
     end
    covij = getCovMat(P, S, COVIJ);
end