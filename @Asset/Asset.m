classdef Asset
    properties
        sedol;%unique stock identifier
        sector;%information about the industrial sector each asset belongs to
        alphaScore;%expected return score, higher values higher future returns expected
        name; %asset's name
        benchWeight; %proportion of each asset in the benchmark
        weight; %proportion of each asset in the portfolio
        r; %asset return
        sectorwrite; %excel only accepts a vector of chars, not a string
        SEDOLidentifier;
    end
    methods
        function obj = Asset(sed, n, sec, bw, alpha, ret, w)
            obj.sedol = sed;
            obj.SEDOLidentifier = convertCharsToStrings(sed);
            obj.sectorwrite = sec;
            obj.sector = convertCharsToStrings(sec);
            obj.alphaScore = alpha;
            obj.name = n;
            obj.benchWeight = bw; 
            obj.r = ret;
            obj.weight = w;
            %{
            fprintf('%s\t%s\t%f\t%s\t%f\n', obj.sedol, obj.sector,...
                obj.alphaScore, obj.name, obj.benchWeight); %print data
            %}
        end
    end
    methods(Static)
       [list, eIdx, refDate] = getAll(tsSet, bIdx) %returns a list of assets and an 
       %index indicating where the date changed
    end
end