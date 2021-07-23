%% 
function [x, fval, Sfinal, COVIJ, covij] = core_graspQ(P, lambda)
%{
    Some problems:
    -CPU TIME per iter is relatively high. Create a vector Sc
     for each sector and mcapQ for random preselection? select |secSET|+
     |mcapSET| stocks and then use a normal RCL vector, doing this you'll
     guarantee that S contains all sectors and mcapQ. You need to select 
     two stocks from each sector and mcapq in order to guarantee that MiLS will not
     draw a sector or mcapq
      -We need to get better mean fval
%}
%function [A,b,x, fval] = core_graspQ(P, lambda)
    x = [];
    fval = 1000;
    Sfinal = [];
    meanT = 0;% mean time per iteration
    totalT = 0;%total elapsed time
    %configure GRASP input Parameters
    itrMAX_GRASP = 10;
    a = rand(); b = rand();
    K = round(60+(a*10)-(b*10)); %size of the solution (number of stocks to be selected)% remember that: 50<=K<=70 (r9.1)
    SIZERCL = round(length(P.initialAssetList)/2); %fixed approximately to the half of the |P.initialAssetList|
    %configure Minor Local Search Parameters
    minLocalSearch = 0.7; %probability of Min Local Search
    %markowitz MV problem inputs
    COVIJ = getCovMat(P);
    %Start GRASP-QUAD
    itrGRASP = 1;
    T0 = cputime; %instant where GRASP was initiated (Used to calculate totalT)     
    acumfval = 0; %used to calculate mean_fval
    while (itrGRASP <= itrMAX_GRASP)
        t0 = cputime; %tic
        proceedGRASP = 0; 
        sizeRCL = SIZERCL;
        S = []; %initialize(empty list)the list(tuple) of the selected stocks
        Sc = 1:length(P.initialAssetList); % S complement (list of all stocks - size N)
        RCL = []; %Restricted Candidate List (stocks are randomly selected from this list)  
        sec_mcapq_Threshold = round(0.7*K);
        setGreedySelection = 0; %use this after a random selection
        %stock selection (solution construction) phase
        %K selection: let the BEST K WINS!!!!!!! The best K have stocks in
        %every sector and mcapQ!!
        while (length(S) < K) %K is the number of stocks to be selected (we need to vary in [50,70]) 
            if(length(S) < 1) 
                    %fprintf('\nstart stock selection with: sizeRCL = %d   sec_mcapq_threshold = %d  K = %d',sizeRCL, sec_mcapq_Threshold,K); 
            end
            %evaluate greedy values Fs for each candidate stock s (a candidate stock is not included in S)
            [Fs_pos,Fs_neg] = getGreedyVals(P, lambda, Sc, COVIJ);%get sorted positive(Fs_pos) and negative(Fs_neg) indexes of stocks relative to the greedyvals
            if(sizeRCL <= length(Fs_pos))
                RCL = Fs_pos(1:sizeRCL); %fcl restricted candidate list
            else
                RCL = Fs_pos;
                RCL = [RCL Fs_neg(1:length(Fs_neg))];
            end
            %check S sector set cardinality and mcapq set cardinality (check when about 90% of the sol. is constructed)
            if(length(S) > sec_mcapq_Threshold) 
                if(~proceedGRASP) %when proceedGRASP = 1, you jump to the else and continue GRASP normally
                    %fprintf('\nchecking S sector set cardinality and mcapq set cardinality\n ');
                    [S, Sc, proceedGRASP] = adjustSecMcapQ(P, S, Sc, RCL); %check
                    if( ~proceedGRASP)%if ~mcap_ok nor sector_ok then do a random selection
                        sizeRCL = length(P.initialAssetList);
                        %reset solution and do random selection
                     %   fprintf('\nswitching to random selection')
                        setGreedySelection = 1; %use this to compensate loss of greediness
                        S = []; %initialize(empty list)the list(tuple) of the selected stocks
                        Sc = 1:length(P.initialAssetList); % S complement (list of all stocks - size N)
                        RCL = []; %Restricted Candidate List (stocks are randomly selected from this list) 
                        a = rand(); b = rand();
                        %get another K candidate
                        K = round(60+(a*10)-(b*10)); % remember that: 50<=K<=70 (r9.1)
                    end
                else
                    %select a random stock from RCL
                    s = randsample(length(RCL),1);
                    S = [S RCL(s)]; %S U s
                    Sc = setdiff(Sc, RCL(s)); %Sc - s
                    if(setGreedySelection)
                      %  fprintf('\ncompensate randomness...switching to greedy selection');
                        sizeRCL= round(0.2*K);
                        setGreedySelection = 0;
                    end
                end
            else
                %select a random stock from RCL
                s = randsample(length(RCL),1);
                S = [S RCL(s)]; %S U s
                Sc = setdiff(Sc, RCL(s)); %Sc - s
            end
            %fprintf('\n(%d)|S| + (%d)|Sc| -> %d   #sizeRLC -> %d',length(S), length(Sc), length(S) + length(Sc), sizeRCL);
        end
        %solve quadratic program
        S = sort(S,'ascend');
        %fprintf('\nsolving QCQP:\t(CARD CONSTRAINT SETTLED!)\n');
        [xtemp, fvaltemp, exitflag] = solve_qcqp(P, lambda, S, COVIJ);
        if(isFeasible(itrGRASP, xtemp, fvaltemp, exitflag, P, S))
            %Minor Local Search
            if(rand() <= minLocalSearch)
                %fprintf('applying MiLS...\n')
                [improved, St, xt, fvalt] = MiLS(itrGRASP, round(0.1*K),...
                    P, lambda, S, Sc, xtemp, fvaltemp, COVIJ);
                if(improved) %check if fval was improved
                    fprintf('MiLS improved the solution!\n');
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
        fprintf('best sol -> fval = %f       mean fval = %f\n', fval, acumfval/(itrGRASP-1));
        toc = cputime - t0; 
        fprintf('iteration elapsed time: %fs\tMeanTPI: %fs\n\n', toc, (cputime - T0)/(itrGRASP-1));
    end
    fprintf('\ntotal elapsed time: %fs\tMeanTPI: %fs\n', cputime - T0,(cputime - T0)/(itrGRASP-1));
    xreduced = x;
    x = zeros(1,length(P.initialAssetList));
     for i = 1:length(P.initialAssetList)
        if(ismember(i,Sfinal)) %if stock i is included in the Portfolio S
            x(i) = xreduced(find(ismember(Sfinal,i)));
        end
     end
    covij = getCovMat(P, S, COVIJ);
end