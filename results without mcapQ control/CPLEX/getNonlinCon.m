function [l, Q, r] = getNonlinCon(P,lambda, H, f)
 nAssets = length(P.initialAssetList);
    assets = P.initialAssetList;
    % YOU NEED TO SUBTRACT ALPHA*lambda terms from f: f = f + alpha
    for (i = 1:nAssets)
        f(i,1) = f(i,1) + (lambda*assets{1,P.i_asset(i)}.alphaScore);
    end    
    %get wbi*COVIJ*wbj and add it to the rhs...remember that when we expand
    % di*dj and di^2 the scalar wbi*COVIJ*wbj remains...
    brisk = 0; %benchmark risk
    cov = 0;
    bi = 0;
    bj = 0;
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t))
            cov = 2*P.covij(t); %sum covij and covji
        else
            cov = P.covij(t);
        end
        bi = assets{1,P.i_asset(t)}.benchWeight;
        bj = assets{1,P.j_asset(t)}.benchWeight;
        brisk= brisk + bi*cov*bj;
    end
    r1 = brisk - (0.05^2);
    r2 = 0.1^2 - brisk;
    Q = H; %this
    l = f; %this
    r = r2; %THIS
%{
    nNonlinConstr = 2;
    nAssets = length(P.initialAssetList);
    assets = P.initialAssetList;
    Q = cell (1,nNonlinConstr);
    % YOU NEED TO SUBTRACT ALPHA*lambda terms from f: f = f + alpha
    for (i = 1:nAssets)
        f(i,1) = f(i,1) + (lambda*assets{1,P.i_asset(i)}.alphaScore);
    end
    %constraint(11.1: We need to transform >= 0.05 in <= 0, multiplying by -1  )
    Q{1,1} = -H;
    l1 = -f; 
    %constraint(11.2: We need to transform <= 0.1 in <= 0)
    Q{1,2} = H;
    l2 = f;
    l = horzcat(l1,l2);
    %get wbi*COVIJ*wbj and add it to the rhs...remember that when we expand
    % di*dj and di^2 the scalar wbi*COVIJ*wbj remains...
    brisk = 0; %benchmark risk
    cov = 0;
    bi = 0;
    bj = 0;
    for t=1:length(P.covij)  %loop over covij vector
        if(P.i_asset(t) ~= P.j_asset(t))
            cov = 2*P.covij(t); %sum covij and covji
        else
            cov = P.covij(t);
        end
        bi = assets{1,P.i_asset(t)}.benchWeight;
        bj = assets{1,P.j_asset(t)}.benchWeight;
        brisk= brisk + bi*cov*bj;
    end
    r1 = brisk - (0.05^2);
    r2 = 0.1^2 - brisk;
    r = [r1 r2]; %rhs, you need to add the scalar part
%}
end