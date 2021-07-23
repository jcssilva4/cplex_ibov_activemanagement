function [SInitial, ScInitial, witPRE] = getGQInitialSolution(P0, P1)
    SInitial = []; %initial asset list (maybe t-1 portfolio)
    ScInitial = 1:length(P1.initialAssetList); %initial sp500 list
    %get persistent portfolios contained in P0 and in S&P500
    assetsP0 = P0.initialAssetList;
    SP500assets = P1.initialAssetList;
    witPRE2 = 0; %witPRE denominator: witPRE pt 2 for persistent stocks and withdrawn sotcks
    witPRE = [];
    %getSP500 sedol list
    SP_SEDOLS = strings(1,length(SP500assets));
    for s1 = 1:length(SP500assets)
        SP_SEDOLS(s1) = SP500assets{1,s1}.SEDOLidentifier;
    end
    for s0 = 1:length(assetsP0)%loop over all P0 stocks
        value = 0;
        %fprintf('current s0.SEDOL: %s\n',assetsP0{1,s0}.SEDOLidentifier) %this
        value = assetsP0{1,s0}.weight*(1+assetsP0{1,s0}.r);
        if(ismember(assetsP0{1,s0}.SEDOLidentifier,SP_SEDOLS)) %if persistent
            SPidx = find(ismember(SP_SEDOLS,assetsP0{1,s0}.SEDOLidentifier)); %get SP500 idx associated with s0
            %fprintf('%s (w_P0 = %f) is also included in P1-> idx: %d\n',assetsP0{1,s0}.SEDOLidentifier,assetsP0{1,s0}.weight, SPidx);
            SInitial = [SInitial SPidx]; %add it to t portfolio
            ScInitial = setdiff(ScInitial, SPidx); 
            witPRE = [witPRE; [SPidx,value]]; %store index and (w_it-1*r_t-1)
        else %this stock is not included in SP500_t (withdraw)
            
        end
        witPRE2 = witPRE2 + value;
    end
    %get witPREi
    for i = 1:length(SInitial)
        witPRE(i,2) =  witPRE(i,2)/witPRE2;
    end
end