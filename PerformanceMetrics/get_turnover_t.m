function [TOopt_t] = get_turnover_t(P0, P1)
    %{
    THE TURNOVER IS DIVIDED IN 3 PARTS:
    part 1: withdrawn stocks:
        turnover = turnover + witPRE_i, where i belongs to P0
    part 2: persistent stocks:
        turnover = turnover + abs(wit - witPRE_i), where i belongs to both P0 AND P1
    part 3: new stocks:
        turnover = turnover + wit, where i belongs to P1  
    %}
    assetsP0 = P0.initialAssetList;
    assetsP1 = P1.initialAssetList;
    TOopt_t = 0;
    witPRE1_excl = []; %witPRE numerator: witPRE pt 1 for withdrawn stocks (withdrawn: s only belongs to P0)
    witPRE1_sed = []; %store persistent stocks sedols
    witPRE1_prs = []; %witPRE numerator: witPRE pt 1 for persistent stocks (persistent: s belongs to both P0 and P1)
    witPRE2 = 0; %witPRE denominator: witPRE pt 2 for persistent stocks and withdrawn sotcks
    %get witPRE for persistent and withdrawn stocks
    %fprintf('\n'); %this
    %getP1 sedol list
    P1_SEDOLS = strings(1,length(assetsP1));
    for s1 = 1:length(assetsP1)
        P1_SEDOLS(s1) = assetsP1{1,s1}.SEDOLidentifier;
    end
    for s0 = 1:length(assetsP0)%loop over all P0 stocks
        %fprintf('current s0.SEDOL: %s\n',assetsP0{1,s0}.SEDOLidentifier) %this
        value = 0;
        if(ismember(assetsP0{1,s0}.SEDOLidentifier,P1_SEDOLS)) %if persistent
            %fprintf('%s (w_P0 = %f) is also included in P1-> idx: %d\n',assetsP0{1,s0}.SEDOLidentifier,assetsP0{1,s0}.weight, s1);
            witPRE1_sed = [witPRE1_sed assetsP0{1,s0}.SEDOLidentifier]; %build list of persistent indexes
            value = assetsP0{1,s0}.weight*(1+assetsP0{1,s0}.r);
            witPRE1_prs = [witPRE1_prs value];
        else %this stock is not included in P1 (withdraw)
            value = assetsP0{1,s0}.weight*(1+assetsP0{1,s0}.r);
            witPRE1_excl = [witPRE1_excl value];
        end
        witPRE2 = witPRE2 + value;
    end
    %witPRE1_prs = [witPRE_sed; witPRE1_prs]; %concatenate  SEDOL      A     BRYGH1  BRYGH1
    %now you have:                                         witPRE1  0.001    0.04   0.007 
    %now you can calculate parts 1 and 2 of the turnover
    %compute turnover part 1
    TOopt1_t = 0;
    for s0 = 1:length(witPRE1_excl) %compute part 1
        TOopt1_t = TOopt1_t + (witPRE1_excl(s0)/witPRE2);
    end
    %fprintf('\npart 1: %f', TOopt1_t); %this
    %compute part 2
    TOopt2_t = 0;
    if(~isempty(witPRE1_prs))
        for s0 = 1:length(witPRE1_prs(1,:)) %for every P0 SEDOL included in P1
            P1_idx = find(ismember(P1_SEDOLS,witPRE1_sed(s0))); %included associated assetsP1 index of s0
            witPRE = witPRE1_prs(s0)/witPRE2;
            TOopt2_t = TOopt2_t + abs(assetsP1{1,P1_idx}.weight - witPRE); 
        end
    else
        witPRE1_sed = 'ALLNEW';
    end
    %fprintf('\npart 2: %f', TOopt2_t); %this
    %compute turnover part 3
    TOopt3_t = 0;
    for s1 = 1:length(assetsP1)%loop over all P1 stocks
        if(~ismember(assetsP1{1,s1}.SEDOLidentifier, witPRE1_sed)) %check if s1 is a new stock
            TOopt3_t = TOopt3_t + assetsP1{1,s1}.weight;
        end
    end
    %fprintf('\npart 3: %f', TOopt3_t); %this
    TOopt_t = TOopt1_t + TOopt2_t + TOopt3_t;
    %fprintf('\nturnover =  %f\n', TOopt_t); %this
end
        


