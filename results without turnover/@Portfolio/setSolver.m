function solver = setSolver(solverName)
    if(strcmp(solverName, 'GA'))
        solver = GA_Solver();
    end
end