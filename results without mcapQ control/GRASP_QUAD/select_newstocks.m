function [S, Sc] = select_newstocks(P, lambda, COVIJ, sizeRCL, K, S0, Sc0)
    S = S0;
    Sc = Sc0;
    Scdiv = Sc;
    RCL = []; %Restricted Candidate List (stocks are randomly selected from this list)  
    %PRIORIZATION: you need at least two assets in each sector
    secOrder = []; %num of purchase orders for each sector
    numStockperSec = getNumStockperSec(P, S);
    auxSec =1;
    for sec = 1:length(numStockperSec)
        if(numStockperSec(sec) < 2)
            %fprintf('\npriority: buy %d assets from %s',secOrder(auxSec,2),P.sectorSet{1,sec});
            [Fs_pos,Fs_neg] = getBuyGreedyVals(P, lambda, P.idxSec{1,sec}, COVIJ,1);%get sorted positive(Fs_pos) and negative(Fs_neg) indexes of stocks relative to the greedyvals
            sizeRCL = 0.5*P.idxSec{1,sec};
            RCL = [];
            if(sizeRCL <= length(Fs_pos))
                RCL = Fs_pos(1:sizeRCL); %fcl restricted candidate list
            else
                RCL = Fs_pos;
                RCL = [RCL Fs_neg(1:length(Fs_neg))];
            end
            for i = 1:2-numStockperSec(sec)
                s = randsample(length(RCL),1);
                S = [S RCL(s)]; %S U s
                Scdiv = setdiff(Scdiv, RCL(s)); %Sc - s1
                RCL = setdiff(RCL, RCL(s));
            end
        end
    end
    Sc = Scdiv;
    while (length(S) < K) %K is the number of stocks to be selected (we need to vary in [50,70]) 
        %evaluate greedy values Fs for each candidate stock s (a candidate stock is not included in S)
        %if(length(S) >= length(P.sectorSet)) %USE NORMAL RCL
            [Fs_pos,Fs_neg] = getBuyGreedyVals(P, lambda, Sc, COVIJ);%get sorted positive(Fs_pos) and negative(Fs_neg) indexes of stocks relative to the greedyvals
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
            %{
            %already diversified
        else %first: use a Sc that contain every sector and mcapq to ensure diversification
            %select a random stock from RCL
            [SDiversf, Sc] = getSDiversf(P, Sc,lambda, COVIJ);
            S = SDiversf; %S U s
        end
            %}
    end
end