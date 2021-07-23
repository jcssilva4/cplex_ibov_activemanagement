function tsSet = readSP500file(filepath)
    fprintf('\nreading S&P500 file...');
    fileID = fopen(filepath, 'r');%open timeseries file
    tsSet = textscan(fileID, '%s %s %s %s %s %s %s %s', 'delimiter', '\t');
    fclose(fileID);
    fprintf('ok\n');
end