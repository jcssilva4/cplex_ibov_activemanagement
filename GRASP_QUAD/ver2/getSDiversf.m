function [SDiversf, Scdiv] = getSDiversf(P, Sc, lambda, COVIJ)
    SDiversf = [];
    Scdiv = Sc;
    %first...get all sectors!
    for s = 1:length(P.sectorSet) 
        [Fs_pos,Fs_neg] = getGreedyVals(P, lambda, P.idxSec{1,s}, COVIJ,1);%get sorted positive(Fs_pos) and negative(Fs_neg) indexes of stocks relative to the greedyvals
        sizeRCL = 0.5*P.idxSec{1,s};
        RCL = [];
        if(sizeRCL <= length(Fs_pos))
            RCL = Fs_pos(1:sizeRCL); %fcl restricted candidate list
        else
            RCL = Fs_pos;
            RCL = [RCL Fs_neg(1:length(Fs_neg))];
        end
        
        s = randsample(length(RCL),1);
        SDiversf = [SDiversf RCL(s)]; %S U s
        Scdiv = setdiff(Scdiv, RCL(s)); %Sc - s1
        RCL = setdiff(RCL, RCL(s));
        s = randsample(length(RCL),1);
        SDiversf = [SDiversf RCL(s)]; %S U s
        Scdiv = setdiff(Scdiv, RCL(s)); %Sc - s1
        RCL = setdiff(RCL, RCL(s));
    end
    %second...get all quintiles!
    for q = 1:5 
        [Fs_pos,Fs_neg] = getGreedyVals(P, lambda, P.idxQ{1,q}, COVIJ,1);%get sorted positive(Fs_pos) and negative(Fs_neg) indexes of stocks relative to the greedyvals
        sizeRCL = 0.5*P.idxQ{1,q};
        RCL = [];
        if(sizeRCL <= length(Fs_pos))
            RCL = Fs_pos(1:sizeRCL); %fcl restricted candidate list
        else
            RCL = Fs_pos;
            RCL = [RCL Fs_neg(1:length(Fs_neg))];
        end
        s = randsample(length(RCL),1);
        if(~ismember(RCL(s),SDiversf)) %check if this stock is already contained in SDiversf
            SDiversf = [SDiversf RCL(s)]; %S U s
            Scdiv = setdiff(Scdiv, RCL(s)); %Sc - s1
            RCL = setdiff(RCL, RCL(s));
        end  
        s = randsample(length(RCL),1);
        if(~ismember(RCL(s),SDiversf)) %check if this stock is already contained in SDiversf
            SDiversf = [SDiversf RCL(s)]; %S U s
            Scdiv = setdiff(Scdiv, RCL(s)); %Sc - s1
            RCL = setdiff(RCL, RCL(s)); 
        end
    end
end
