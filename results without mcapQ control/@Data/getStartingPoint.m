function w = getStartingPoint(filepath)
    fprintf('\nreading starting point(initial portfolio) file...');
    fileID = fopen(filepath, 'r');%open timeseries file
    w = textscan(fileID, '%f');
    fclose(fileID);
    fprintf('ok\n');
end