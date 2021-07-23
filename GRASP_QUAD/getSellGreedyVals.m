function [Fs_pos,Fs_neg] = getSellGreedyVals(P, S, witPRE)
    Fs = zeros(1, length(S));
    Fs_pos = []; %positive greedy vals indexes
    Fs_neg = []; %negative greedy vals indexes
    npos = 0; %num of positive
    nneg = 0; %num of negative
    %fprintf('second option')
    for i = 1 : length(S)
        S;
        witPRE(:,1);
        witIdx = find(ismember(witPRE(:,1),S(i)));
        Fs(i) = witPRE(witIdx,2); %get turnover contribution ;
        if(Fs(i) < 0) %if Fs is negative
            %Fs_neg = [Fs_neg Fs];   %then Fs_neg U Fs 
            nneg = nneg + 1;
        else
            %Fs_pos = [Fs_pos Fs];   %then Fs_pos U Fs 
            npos = npos + 1;
        end
    end
    %sort Fs_pos and Fs_neg by ascending order
    [witVals, FsIdx] = sort(Fs, 'ascend');
    
    %to avoid index repetitions, like S = [10 2 344 20 10 11 12...] we need
    %to retrive indexes associated with Sc
    if(npos>0)  
        for i = nneg+1: nneg+npos %positive values comes after negative values in Fs: Fs = [-2 -1 2 1 3]
            witIdx = find(ismember(witPRE(:,2),witVals(i)));
            SIdx = find(ismember(S,witPRE(witIdx,1)));
            Fs_pos = [Fs_pos S(SIdx)];  %retrive S associated indexes
        end
    end 
    if(nneg>0)
        for i = 1 : nneg 
            Fs_neg = [Fs_neg witPRE(FsIdx(i),1)];    
        end
    end 
end

function getTurnover
    
end