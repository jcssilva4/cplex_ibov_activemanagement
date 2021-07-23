function returns = readResultsTempFile(filepath)
    fprintf('\nreading Results Template file...');
    rfid = fopen(filepath, 'rt'); %open returns file
    rawResTemp = textscan(rfid, '%s'); %create a cell containing raw data
    returns = zeros(1,length(rawResTemp)-1); %create a vector that will store returns
    rawLine = cell(1,8); %temp cell that contains a line of rawResTemps
    auxcounter = 0; %aids construction of the initial portfolio weight vector
    for i=1:(length(rawResTemp{1,1})-1)
        rawLine = textscan(rawResTemp{1,1}{i+1}, '%s %s %s %s %s %s %s',... 
        'delimiter', ';');
        returns(i) = str2double(rawLine{1,4}); %get returns (4th rawLine column)
    end
    fprintf('ok\n');
    %REMEMBER YOU STILL NEED TO USE rfid TO GET THE INITIAL PORTFOLIO
    fclose(rfid);
end