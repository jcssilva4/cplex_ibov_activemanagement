function [improved, SMiLS, xMiLS, fvalMiLS] = MiLS(itrGRASP, UB_minLocalSearch,...
            P, lambda, S, Sc, x, fval, COVIJ)
    
    SMiLS = S;
    xMiLS = x;
    fvalMiLS = fval; 
    improved = 0;
    i = 1;
    while (i < UB_minLocalSearch)
        s1 = randsample(length(S),1); 
        s2 = randsample(length(Sc),1);
        original_s1 = S(s1); %use this to retrieve the original S(s1)
        original_s2 = Sc(s2); %use this to retrieve the original Sc(s2)
        %exchange s1 by s2
        S(s1) = original_s2; 
        Sc(s2) = original_s1;
        S = sort(S,'ascend');
        Sc = sort(Sc,'ascend');
        [xt, fvalt, exitflag] = solve_qp(P, lambda, S, COVIJ);
        if(isFeasible(itrGRASP, xt, fvalt, exitflag, P, S))
            if(fvalt < fvalMiLS)
                improved = 1;  
                fvalMiLS = fvalt;
                xMiLS = xt;
                SMiLS = S;
            end
        end
        i = i + 1;
        %back to the original solution
        s1 = find(ismember(S,original_s2)); % we need to use this because we sort S, and indexes may change...
        S(s1) = original_s1;
        s2 = find(ismember(Sc,original_s1)); % we need to use this because we sort Sc, and indexes may change...
        Sc(s2) = original_s2;
        S = sort(S,'ascend');
        Sc = sort(Sc,'ascend');
    end
end

function [sector_satisfied, mcapq_satisfied, sectorSet, mcapqSet] = checkSecMcapQ(P, S, sec_ok, mcap_ok)
    sector_satisfied = 0;
    mcapq_satisfied = 0;
    sectorSet = [];
    mcapqSet = [];
    if(~sec_ok)
        %get sector set
        sectorSet = P.initialAssetList{1,S(1)}.sector; %convert this array of chars to string
        for i = 2:length(S)
           %compare strings!
           currentSector = P.initialAssetList{1,S(i)}.sector;
           if(~ismember(currentSector,sectorSet))
               sectorSet = [sectorSet currentSector]; %add new sector
           end
        end
        if(~(length(sectorSet) < length(P.sectorSet)))
            sector_satisfied = 1; %S contains all sectors
        end
    else
        sector_satisfied = 1; %S already contains all sectors
    end
    if(~mcap_ok)
        %get mcapq set
        mcapqSet = P.initialAssetList{1,S(1)}.mCapQ; %convert this array of chars to string
        for i = 2:length(S)
           %compare strings!
           currentmcapq = P.initialAssetList{1,S(i)}.mCapQ;
           if(~ismember(currentmcapq,mcapqSet))
               mcapqSet = [mcapqSet currentmcapq]; %add new sector
           end
        end
        if(~(length(mcapqSet) < 5))
            mcapq_satisfied = 1; %S contains all mcapq
        end
    else
        mcapq_satisfied = 1; %S already contains all mcapq 
        mcapqSet = 1;
    end
end