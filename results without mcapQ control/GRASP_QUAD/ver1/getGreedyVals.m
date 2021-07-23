function [Fs_pos,Fs_neg] = getGreedyVals(P, lambda, Sc, COV)
    
    assets = P.initialAssetList;
    Fs = zeros(1, length(Sc));
    Fs_pos = []; %positive greedy vals indexes
    Fs_neg = []; %negative greedy vals indexes
    npos = 0; %num of positive
    nneg = 0; %num of negative
    for i = 1 : length(Sc)
            COVrowi = COV(Sc(i),:);
            A = sum(COVrowi)/length(Sc);
            B = lambda*assets{1,Sc(i)}.alphaScore;
            Fs(i) = (1+A)/(1+B);
            if(Fs(i) < 0) %if Fs is negative
                %Fs_neg = [Fs_neg Fs];   %then Fs_neg U Fs 
                nneg = nneg + 1;
            else
                %Fs_pos = [Fs_pos Fs];   %then Fs_pos U Fs 
                npos = npos + 1;
            end
    end
    %sort Fs_pos and Fs_neg by ascending order
    [trash, FsIdx] = sort(Fs, 'ascend');
    %to avoid index repetitions, like S = [10 2 344 20 10 11 12...] we need
    %to retrive indexes associated with Sc
    if(npos>0)  
        for i = nneg+1: nneg+npos %positive values comes after negative values in Fs: Fs = [-2 -1 2 1 3]
            Fs_pos = [Fs_pos Sc(FsIdx(i))];  %retrive Sc associated indexes
        end
    end 
    if(nneg>0)
        for i = 1 : nneg 
            Fs_neg = [Fs_neg Sc(FsIdx(i))];    
        end
    end 
end