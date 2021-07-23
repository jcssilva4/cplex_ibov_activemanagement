function [y, m, d] = getNextDate(currYear, currMonth, currDay, dDay)
    %set current calendar coords
    [i, j, c] = setCcoords(currYear, currMonth, currDay); 
    %get next date
    for daycounter = 1:dDay
         if (j == 7) %then go to the next week
                j = 1;
                i = i + 1;
            else
                j = j + 1; %then go to the next day
        end
        if (c(i,j) == 0)%then go to the next month
            if(currMonth == 12) %then go to the next year
                currMonth = 1;
                currYear = currYear + 1;
            else
                currMonth = currMonth + 1;
            end
        %get new calendar starting coords (day 1)
        [i, j, c] = setCcoords(currYear, currMonth, 1); 
        end
    end
    currDay = c(i,j);
    y = currYear;
    m = currMonth;
    d = currDay;
end