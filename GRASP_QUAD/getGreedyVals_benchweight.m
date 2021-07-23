function [Fs_pos,Fs_neg] = getGreedyVals_benchweight(P, lambda, Sc, COV,d)
    assets = P.initialAssetList;
    Fs = zeros(1, length(Sc));
    Fs_pos = []; %positive greedy vals indexes
    Fs_neg = []; %negative greedy vals indexes
    npos = 0; %num of positive
    nneg = 0; %num of negative
    if(nargin < 5)
        for i = 1 : length(Sc)
            COVrowi = COV(Sc(i),:);
            A = sum(COVrowi)/length(Sc);
            B = lambda*assets{1,Sc(i)}.r*assets{1,Sc(i)}.alphaScore;
            %Fs(i) = (1+B)/(1+A);
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
    else%USE topK_WEIGHTBENCH method
        %fprintf('second option')
        A = 0;
        for i = 1 : length(Sc)
            A = -P.initialAssetList{1,Sc(i)}.benchWeight; 
            %B = 1+ lambda*d(Sc(i))*assets{1,Sc(i)}.alphaScore;
            Fs(i) = A;%((1+A)/length(Sc))/B ;
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
    end
    
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