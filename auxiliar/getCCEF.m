function [ropt, stdopt, AR, TE, MRpareto, VARpareto,...
    ARpareto , TEpareto] = getCCEF(P, iterMAX)
    Ssols = [];
    W = [];
    D = [];
    LAMBDA = [];
    TE = [];
    AR = [];
    lambda = 0;
    frontierPoints = 20;
    COV = getCovMat(P);
    ropt = [];
    stdopt = [];
    while(lambda <= 10)
        fprintf('lambda-> %f\n\n', lambda)
        [d,w,~,~,~,~] = core_graspQ(P, lambda, iterMAX); 
        D = [D; d];
        W = [W; w];
        LAMBDA = [LAMBDA lambda];
        lambda = lambda + 10/frontierPoints;
    end
    for i = 1 : frontierPoints
        d = D(i,:);
        w = W(i,:);
        sum(w)
        ropt = [ropt getOverallR(w, P)];
        stdopt = [stdopt getStdDev(w, P)];
        TE = [TE getTEFinal(d, P)];
        AR = [AR getActiveR(d, P)];
    end
    [VARpareto, MRpareto, ~] = getNDPoints(stdopt, ropt);
    [TEpareto, ARpareto, ~] = getNDPoints(TE, AR);
end

function [TEpareto, ARpareto, nonDominatedPoints] = getNDPoints(TE, AR) 
    %get non-dominated solutions
    nsols = length(TE) ;
    dominated = 0;
    countcrit = 0; %count losses in a objective
    nonDominated = []; %list containing non dominated sols!
    TEpareto = [];
    ARpareto = [];
    nonDominatedPoints = [];
    for i = 1:nsols
        for j = 1:nsols
            if(~dominated) 
                %j dominates i?
                if(TE(i) - TE(j) > 0)
                    countcrit = countcrit + 1; fprintf('%d dominates %d in respect to TE\n' , j, i);
                end  %TE(min): if TE(j) < TE(i), i is dominated in relation to this objective
                if(AR(i) - AR(j) < 0)
                    countcrit = countcrit + 1; fprintf('%d dominates %d in respect to AR\n' , j, i);
                end   %TE(min): if AR(j) > AR(i), i is dominated in relation to this objective
                if(countcrit > 1) 
                    fprintf('%d dominates %d\n' , j, i); dominated = 1;
                end
            end
            countcrit = 0;
        end
        if(~dominated), nonDominated = [nonDominated i]; end
        dominated = 0;
    end
    for i = 1:length(nonDominated)
        TEpareto(i) = TE(nonDominated(i));
        ARpareto(i) = AR(nonDominated(i));
    end
end