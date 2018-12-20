function out = repeatEntries(val,kTimes)
%REPEATENTRIES fills a matrix with k repeats the rows of the input matrix
%
% SYNOPSIS out = repeatEntries(val,kTimes)
%
% INPUT    val    : matrix (or vectors) containing the rows to repeat (works for strings, too)
%          kTimes : number of repeats of each row (scalar or vector of size(vlaues,1))
%
% OUTPUT   out    : matrix of size [sum(kTimes) size(values,2)] containing
%                   repeated entries specified with k
%
% EXAMPLES     repeatEntries([1;2;3;4],[2;3;1;1]) returns [1;1;2;2;2;3;4]
%
%              repeatEntries([1;2;3;4],2) returns [1;1;2;2;3;3;4;4]
%
% c: jonas, 2/04
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% note: in case we need to speed this up: adapt the code below
% nn = cellfun(@numel,points);
% a = find(nn);
% index = zeros(sum(nn),1);
% index([1;cumsum(nn(a(1:end-1)))+1])=1;
% 
% % get the indices
% ii = a(cumsum(index));

%===========
% test input
%===========

% nargin
if nargin ~= 2 || isempty(val) || isempty(kTimes)
    error('two non-empty input arguments are needed!')
end

% size
valSize = size(val);
if length(valSize)>2
    error('only 2D arrays supported for val')
end



% decide whether we have scalar k
numK = length(kTimes);
if numK == 1
    scalarK = 1;
elseif numK ~= valSize(1)
    error('vector k must have the same length as the number of rows in val or be a scalar')
else
    % check again whether we could use scalar k
    if all(kTimes(1) == kTimes)
        scalarK = 1;
        kTimes = kTimes(1);
    else
        scalarK = 0;
    end
end

% do not care about size of k: we want to make a col vector out of it - and
% this vector should only contain nonzero positive integers
kTimes = round(kTimes(:));
% if there are any negative values or zeros, remove the entry
if scalarK && kTimes < 1
    out = [];
    return
end
if ~scalarK
    badK = kTimes < 1;
    kTimes(badK) = [];
    val(badK,:) = [];
    % update valSize
    valSize = size(val);
    if any(valSize==0)
        out = [];
        return
    end
end
%kTimes = max(kTimes,ones(size(kTimes)));


%============
% fill in out
%============

% first the elegant case: scalar k
if scalarK

    % build repeat index matrix idxMat
    idxMat = meshgrid( 1:valSize(1), 1:kTimes(1) );
    idxMat = idxMat(:); % returns [1;1...2;2;... etc]

    out = val(idxMat,:);

    % second: the loop
else

    % init out, init counter
    if iscell(val)
        out = cell(sum(kTimes) , valSize(2));
    else
    out = zeros( sum(kTimes), valSize(2) );
    end
    endct = 0;

    if valSize(2) == 1

        % vector: fill directly

        % loop and fill
        for i = 1:valSize(1)
            startct = endct + 1;
            endct   = endct + kTimes(i);
            out(startct:endct,:) = val(i);
        end % for i=1:valSize(1)

    else

        % matrix: fill via index list

        idxMat = zeros(sum(kTimes),1);

        for i = 1:valSize(1)
            startct = endct + 1;
            endct   = endct + kTimes(i);
            idxMat(startct:endct) = i;
        end % for i=1:valSize(1)
        out = val(idxMat,:);

    end

    % check for strings and transform if necessary
    if ischar(val)
        out = char(out);
    end

end % if doScalar


