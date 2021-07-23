function [date] = formatDate(currYear, currMonth, currDay)
    yy = sprintf('%d',currYear);
    if(currMonth >= 10)
        mm = sprintf('%d',currMonth);
    else
        mm = sprintf('0%d',currMonth);
    end
    if(currDay >= 10)
        dd = sprintf('%d',currDay);
    else
        dd = sprintf('0%d',currDay);
    end
    date = sprintf('%s-%s-%s',yy,mm,dd);
end