classdef Data
    properties
        tsSet %S&P500 timeseries set of Data
        returns %historical returns
        initialPortfolio %vector of weights
    end
    methods
        function obj = Data()
        end
    end
    methods(Static)
        tsSet = readSP500file(filepath); %read S&P500 file
        returns = readResultsTempFile(filepath); %read results template file
        initialPortfolio = readInitPortFile(filepath); %read initialPortfolio file
    end
end