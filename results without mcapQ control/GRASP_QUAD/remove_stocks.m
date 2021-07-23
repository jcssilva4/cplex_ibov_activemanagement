function S = remove_stocks(P, sizeRCL, K, S0, witPRE)
    S = S0;
    Ssell = []; %withdrawn stock list
    RCL = []; %Restricted Candidate List (stocks are randomly selected from this list) 
    criticalsectors = 0; %vector containing sectors with 2 assets
    operations = 0;
    numStockperSec = [];
    while (length(S) > K) %K is the number of stocks to be selected (we need to vary in [50,70]) 
        %evaluate greedy values Fs for each candidate stock s (a candidate stock is not included in S)
        [Fs_pos,Fs_neg] = getSellGreedyVals(P, S, witPRE);%get sorted positive(Fs_pos) and negative(Fs_neg) indexes of stocks relative to the greedyvals
        if(sizeRCL <= length(Fs_pos))
            RCL = Fs_pos(1:sizeRCL); %rcl restricted candidate list
        else
            RCL = Fs_pos;
            RCL = [RCL Fs_neg(1:length(Fs_neg))];
        end
        %you need at least two assets in each sector
        numStockperSec = getNumStockperSec(P, S); %use this when removing stocks to guarantee that every sector has at least one associated stock
        if(~isempty(find(ismember(numStockperSec,2), 1)))
            criticalsectors = find(ismember(numStockperSec,2));
            for sec = 1:length(criticalsectors)
                RCL = setdiff(RCL, P.idxSec{1,criticalsectors(sec)});
            end
        end
        %select a random stock from RCL
        s = randsample(length(RCL),1);
        S = setdiff(S, RCL(s)); %Sc - s
        Ssell = [Ssell RCL(s)];
        RCL = setdiff(RCL, RCL(s)); %Sc - s 
        %fprintf('\n(%d)|S| + (%d)|Ssell| -> %d   #sizeRLC -> %d',length(S), length(Ssell), length(S) + length(Ssell), sizeRCL);
    end
    %if there's a sector with <2 number of stocks
    if(~isempty(find(ismember(numStockperSec,1), 1)))
        criticalsectors = find(ismember(numStockperSec,1));
        for sec = 1:length(criticalsectors)
            %remove 1 asset from S for each sector with only 1 asset (to induce an asset purchase for each sector)
            s = randsample(length(RCL),1); %select a random stock from RCL
            S = setdiff(S, RCL(s)); %Sc - s 
            Ssell = [Ssell RCL(s)];
            RCL = setdiff(RCL, RCL(s)); %Sc - s 
        end
    else
        if(~isempty(find(ismember(numStockperSec,0), 1)))
            criticalsectors = find(ismember(numStockperSec,0));
            for sec = 1:length(criticalsectors)
                %fprintf('\nremove 2 assets from S for each sector with 0 assets (to induce 2 asset purchases for each sector)');
                s = randsample(length(RCL),1); %select a random stock from RCL
                RCLs = RCL(s);
                S = setdiff(S, RCL(s)); %Sc - s 
                Ssell = [Ssell RCL(s)];
                RCL = setdiff(RCL, RCL(s));
                s = randsample(length(RCL),1); %select a random stock from RCL
                RCLs = RCL(s);
                S = setdiff(S, RCL(s)); %Sc - s 
                Ssell = [Ssell RCL(s)];
            end
        end
    end
end

