function out = isEven(in)
%ISEVEN checks whether a number is even
%
% SYNOPSIS out = isEven(in)
%
% INPUT    in :  input (array) of numbers to be tested. 
% OUTPUT   out:  array of size(in) with 
%                   1 for even integers and zero
%                   0 for odd integers
%                 NaN for non-integers
%                out is a logical array as long as the input is all integers.
%
% c: jonas 5/05
% Last modified 11/24/2009 - Jonas Dorn
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

out = mod(in+1, 2);
% Set NaN for non-integer data, because they are neither odd or even
out((out ~= 0) & (out ~= 1)) = NaN;

% since doubles cannot be used for logical indexing, we should convert to
% logicals if possible. 
if all(isfinite(out(:)))
    out = logical(out);
end