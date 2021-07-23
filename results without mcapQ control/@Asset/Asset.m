classdef Asset
    properties
        sedol;%unique stock identifier
        sector;%information about the industrial sector each asset belongs to
        beta;%Measure of an asset’s return sensitivity relative the equity market return as a whole
        alphaScore;%expected return score, higher values higher future returns expected
        name; %asset's name
        benchWeight; %proportion of each asset in the benchmark
        mCapQ;%MarketCapQuintile->measure of asset's capitalisation size relative to the market (distinguish large companies from small ones)
        weight; %proportion of each asset in the portfolio
        r; %asset return
        sectorwrite; %excel only accepts a vector of chars, not a string
        SEDOLidentifier;
    end
    methods
        function obj = Asset(sed, sec, b, alpha, n, bw, mcq, ret, w)
            obj.sedol = sed;
            obj.SEDOLidentifier = convertCharsToStrings(sed);
            obj.sectorwrite = sec;
            obj.sector = convertCharsToStrings(sec);
            obj.beta = b;
            obj.alphaScore = alpha;
            obj.name = n;
            obj.benchWeight = bw;
            obj.mCapQ = mcq;
            obj.r = ret;
            obj.weight = w;
            %{
            fprintf('%s\t%s\t%f\t%f\t%s\t%f\t%f\n', obj.sedol, obj.sector,...
                obj.beta, obj.alphaScore, obj.name, obj.benchWeight,...
                obj.mCapQ); %print data
            %}
        end
    end
    methods(Static)
       [list, eIdx] = getAll(tsSet, returns, bIdx) %returns a list of assets and an 
       %index indicating where the date changed
    end
end