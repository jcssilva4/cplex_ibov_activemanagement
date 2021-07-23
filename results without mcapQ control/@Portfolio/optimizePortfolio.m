function postAssetList = optimizePortfolio(P, rawCovMatFile)
    postAssetList = [];
    if(strcmp(P.solver.solverName, 'GA'))
        P.solver.fitness = P.solver.singleopt_fitness(P.initiallAssetList,...
            rawCovMatFile, P.solver.CovMatIdxMemorized);
        if(~P.solver.CovMatIdxMemorized)
            P.solver.CovMatIdxMemorized = 1; %after first use
        end
    end
end