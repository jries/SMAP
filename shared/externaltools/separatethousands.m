function numstring = separatethousands(number, separator, decimalspaces)
    % ----------------------
    % Author: Eero Kuusi
    % Created 15th July 2010
    % ----------------------
    %
    % SEPARATETHOUSANDS
    %   Separates thousands (and millions etc) with the given separator.
    %
    %   USE:
    %       separatethousands(number, separator, decimalspaces)
    %   INPUTS:
    %       number          (double)    = number to convert
    %       separator       (string)    = separator sign
    %      (decimalspaces   (int)       = number of decimal spaces)
    %
    %
    %   EXAMPLES:
    %
    %       separatethousands(1212.112, ' ', 2)
    %           ans = '1 212.11'
    %
    %       separatethousands(12121212.12, ',', 0)
    %           ans = '12,121,212'
    %
    %
    
    % Default decimalspaces is zero
    if exist('decimalspaces','var') == 0
        decimalspaces = 0;
    end
    
    if isempty(number)
        numstring='';
        return
    end
    if isinf(number)
        numstring='inf';
        return
    end
    if isnan(number)
        numstring='NaN';
        return
    end
    % Check if negative
    if number < 0
        modified_number(1) = -number;
        minus = 1;
    else
        modified_number(1) = number;
        minus = 0;
    end
    
    % Pull out decimals
    if decimalspaces > 0
        decimals = modified_number-fix(modified_number);
        modified_number = fix(modified_number);
    end
    
    % Get each block of numbers
    count = 1;
    while modified_number(count) >= 1000
        modified_number(count+1) = fix(modified_number(count)/1000);
        modified_number(count) = modified_number(count)...
            -modified_number(count+1)*1000;
        count = count+1;
    end
    
    % Sign if negative
    if minus
        numstring = sprintf('-%1.0f',modified_number(end));
    else
        numstring = sprintf('%1.0f',modified_number(end));
    end
    
    % Add blocks
    for i = 1:count-1
        numstring = sprintf('%s%s%03.0f',...
            numstring,separator,modified_number(end-i));
    end
    
    % Add decimals
    if decimalspaces > 0
        decimalformat = sprintf('%%0.%1.0ff',decimalspaces);
        decimalstring = sprintf(decimalformat,decimals);
        decimalstring = decimalstring(2:end); % Get rid of the leading zero
        numstring = [numstring decimalstring];
    end
   
    % End of function
    
end