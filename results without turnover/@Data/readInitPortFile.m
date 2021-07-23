function initialPortfolio = readInitPortFile(filepath)
    fprintf('\nreading initialPorfolio file...');
    fileID = fopen(filepath, 'r');%open timeseries file
    raw = textscan(fileID, '%f');
    initialPortfolio = (raw{1,1})';
    fclose(fileID);
    fprintf('ok\n');
end