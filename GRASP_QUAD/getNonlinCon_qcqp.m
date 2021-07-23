function [l, Q, r] = getNonlinCon_qcqp(P,lambda, S, H, f)
    nVars = length(S);
    assets = P.initialAssetList;
    % YOU NEED TO SUBTRACT ALPHA*lambda terms from f: f = f + alpha
    for i = 1:nVars
        f(i,1) = f(i,1) + (lambda*assets{1,S(i)}.alphaScore);
    end    
    %get wbi*COVIJ*wbj and add it to the rhs...remember that when we expand
    % di*dj and di^2 the scalar wbi*COVIJ*wbj remains...
    brisk = 0; %benchmark risk
    bi = 0;
    bj = 0;
    cov = 0;
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
    Q = H; %this
    l = f; %this
    r = 0.1^2 - brisk; %THIS
end