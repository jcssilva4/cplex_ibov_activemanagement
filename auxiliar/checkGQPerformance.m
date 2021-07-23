function [] = checkGQPerformance(itr)
    if(length(itr)==3)
        plot_UTEEF_CCCEF(itr(1), itr(2), itr(3))
    end
    
end

function [] = plot_UTEEF_CCCEF(itr1, itr2, itr3)
    %%%%%%%%%get CARDINALITY CONSTRAINED Efic frontier%%%%%%%
    iter1 = itr1;
    [ropt_grasp1, stdopt_grasp1, AR_grasp1, TE_grasp1, MRpareto_grasp1,...
        VARpareto_grasp1, ARpareto_grasp1, TEpareto_grasp1] = getCCEF(P, iter1); 
    iter2 = itr2;
    [ropt_grasp2, stdopt_grasp2, AR_grasp2, TE_grasp2, MRpareto_grasp2,...
        VARpareto_grasp2, ARpareto_grasp2, TEpareto_grasp2] = getCCEF(P, iter2);
    iter3 = itr3;
    [ropt_grasp3, stdopt_grasp3, AR_grasp3, TE_grasp3, MRpareto_grasp3,...
        VARpareto_grasp3, ARpareto_grasp3, TEpareto_grasp3] = getCCEF(P, iter3); 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  

    %%%%%%%%%%get UNCONSTRAINED Efic Frontier%%%%%%%%%%%%%%%%%
    [ropt, stdopt,  AR, TE, MRpareto, VARpareto, ...
        ARpareto, TEpareto, LAMBDApareto_MV, LAMBDApareto_Active] = getUEF(P);
    % Find max TE grasp
    maxTEgrasp = -1;
    candidates = [max(TE_grasp1) max(TE_grasp2) max(TE_grasp3)];
    for i = 1:3
        if(maxTEgrasp < candidates(i))
            maxTEgrasp = candidates(i);
        end
    end
    newTEVals = TE(1):maxTEgrasp/100:maxTEgrasp; %min of ef frontier and max of grasp TE solution( intervals of maxTEgrasp/100)
    newAR_graspvals1 = interp1(TE_grasp1, AR_grasp1, newTEVals);
    newAR_graspvals2 = interp1(TE_grasp2, AR_grasp2, newTEVals);
    newAR_graspvals3 = interp1(TE_grasp3, AR_grasp3, newTEVals);
    newAR = interp1(TE, AR, newTEVals); %UEF interpolated
    % Compare the original values to these interpolated values
    figure
    plot(newTEVals, newAR,'b');
    hold on
    plot(newTEVals, newAR_graspvals1, 'k--o');
    hold on
    plot(newTEVals, newAR_graspvals2, 'r--o');
    hold on
    plot(newTEVals, newAR_graspvals3, 'c--o');
    hold off

    legend('Unconstrained TE Efficient Frontier',...
        sprintf('GRASP-QUAD CCEF - %diters',itr1),...
        sprintf('GRASP-QUAD CCEF - %diters',itr2),...
        sprintf('GRASP-QUAD CCEF - %diters',itr3))
end