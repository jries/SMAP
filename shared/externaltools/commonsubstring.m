function [substr,loc1,loc2] = commonsubstring(str1,str2,substringlength)
% commonsubstring - finds the longest common substring(s) betwen two strings
% usage: [substr,loc1,loc2] = commonsubstring(str1,str2)
% usage: [substr,loc1,loc2] = commonsubstring(str1,str2,substringlength)
%
% Find the longest common substrings between a pair of strings,
% or if the substring length is indicated, finds all common
% substrings with that specified length.
%
% Arguments: (input)
%  str1, str2 - character strings that will be searched for
%       common substrings.
%
%  substringlength - (OPTIONAL) - scalar positive integer 
%       if provided. This specifies the length of the common
%       substrings to be generated. If not provided or empty
%       the maximum length common substrings are returned.
%
%       Note, if str1 and str2 contain many thousands of
%       elements, then searching for a very long substring
%       of a specific length may be a time and memory consuming
%       operation.
%
% Arguments: (output)
%  substr - character string that contains the longest common
%       substring. If more than one string has that maximum
%       length, then substr will be a matrix of strings.
%       If multiple solutions were found, then each row of
%       that character matrix will be a solution.
%
%  loc1, loc2 - the respective starting locations of these
%       substrings in each of the original strings provided.
%       Since a given substring may be found multiple times,
%       these will be cell arrays, possibly containing vectors
%       of indices
%
% Example:
% % Generate a pair of long random letter sequences, then
% % determine the longest common substring between them
%  bases = 'acgt';
%  str1 = bases(ceil(rand(1,100000)*4));
%  str2 = bases(ceil(rand(1,100000)*4));
%
%  tic,[substr,ind1,ind2] = commonsubstring(str1,str2);toc
% % Elapsed time is 16.650532 seconds.
%
% % There were two substrings of the maximum length
% % (16) characters.
%  substr
% % substr =
% % gctttagggcgtacgc
% % cttcggataccttgtt
%
% ind1
% % ind1 = 
% %     [22189]
% %     [74425]
%
% ind2
% % ind2 = 
% %     [64948]
% %     [32833]
%
% Example:
% % Find all common substrings of a given fixed length.
% str1 = char('a' + round(rand(1,100)*1.5))
% % str1 =
% % bbbabbbbbabbbbbbabbbbbabbababbbabbabbbbbbabbbbbbbbbbbbbabbaaabbbbbaabbbbbbbbbbabbbbbaaabbabbaabbbbbb
% str2 = 'aaabbabbb';
%
% [substr,ind1,ind2] = commonsubstring(str1,str2,3)
% % substr =
% % aaa
% % aab
% % abb
% % bab
% % bba
% % bbb
% 
% % ind1 = 
% %     [1x2  double]
% %     [1x4  double]
% %     [1x15 double]
% %     [1x11 double]
% %     [1x15 double]
% %     [1x19 double]
%
% % ind2 = 
% %     [         1]
% %     [         2]
% %     [1x2 double]
% %     [         5]
% %     [         4]
% %     [         7]
%
% See also: regexp, strfind, findstr
%
% Author: John D'Errico
% e-mail: woodchips@rochester.rr.com
% Release date: 4/28/2010

if (nargin < 2) || (nargin > 3)
  error('COMMONSUBSTRING:improperarguments', ...
    'Either 2 or 3 arguments must be supplied')
end

if (nargin < 3)
  substringlength = [];
elseif ~isnumeric(substringlength) || (numel(substringlength) > 1) || ...
    isinf(substringlength) || isnan(substringlength) || ...
    (rem(substringlength,1) ~= 0) || (substringlength < 1)
  
  error('COMMONSUBSTRING:improperarguments', ...
    'substringlength must be positive, scalar, finite integer if provided')
end

if ~ischar(str1) || ~ischar(str2)
  error('COMMONSUBSTRING:improperarguments', ...
    'Both arguments must be strings')
elseif isempty(str1) || isempty(str2)
  % If either of these strings are empty, the
  % common substring must be also
  substr = '';
  loc1 = {};
  loc2 = {};
  return
end

% unroll the strings into row vectors
str1 = str1(:).';
str2 = str2(:).';
nstr1 = numel(str1);
nstr2 = numel(str2);
nstr = min(numel(str1),numel(str2));

% was a specific substring length indicated?
if ~isempty(substringlength)
  % yes
  if substringlength > nstr
    substr = '';
    loc1 = {};
    loc2 = {};
  end
  
  % get the unique substrings of that length 
  if nstr1 <= nstr2
    % str1 is the shorter of the two, or they
    % have the same lengths. In either event,
    % take the list of possible substrings
    % from str1.
    substr = substrings(str1,substringlength,1);
  else
    % str2 is the shorter of the two
    substr = substrings(str2,substringlength,1);
  end
  
  % get the locations of each of these strings
  cellsubstr = cellstr(substr);
  loc1 = regexp(str1,cellsubstr);
  loc2 = regexp(str2,cellsubstr);
  
  % delete those substrings thet were not found
  % in both original strings.
  k = cellfun(@isempty,loc1) | cellfun(@isempty,loc2);
  substr(k,:) = [];
  loc1(k) = [];
  loc2(k) = [];
  
  % we can exit now
  return
end

% extract the single character string elements of str1
substr = unique(str1');
cellsubstr = cellstr(substr);

% Location of these single character substrings in str2
% regext will give this quite efficiently
ind2 = regexp(str2,cellsubstr);

% drop those characters that are not located in str2
k2 = cellfun(@isempty,ind2);
substr(k2) = [];
ind2(k2) = [];
cellsubstr(k2) = [];

% Location of these single character substrings in str1
ind1 = regexp(str1,cellsubstr);

% were there any common single digit strings?
% if not, then we are done already
if isempty(substr)
  substr = '';
  loc1 = {};
  loc2 = {};
  return
end

% now, see if we can extend any of these substrings
% by just a single character at a time
flag = true;
substrlength = 1;
while flag
  nsubstr = size(substr,1);
  ext = cell(nsubstr,3);
  L = 0;
  for i = 1:nsubstr
    substri = substr(i,:);
    
    % Where did we find this string in each of
    % str1 and str2?
    loc1 = ind1{i};
    loc2 = ind2{i};
    
    % we can drop any occurrences that lie
    % at the very beginning of str1 or str2
    loc1(loc1 == 1) = [];
    loc2(loc2 == 1) = [];
    
    % find the very next characters that lie just
    % prior to this substring in str1 and str2.
    % pc refers to the previous character, to be
    % then prepended.
    pc1 = str1(loc1 - 1);
    pc2 = str2(loc2 - 1);
    
    % are there any common characters? If
    % so, then we can extend the corresponding
    % substrings
    pccommon = intersect(pc1,pc2);
    ncommon = numel(pccommon);
    if ncommon > 0
      extind1 = cell(ncommon,1);
      extind2 = cell(ncommon,1);
      L = L + 1;
      for j = 1:ncommon
        extind1{j} = loc1(pc1 == pccommon(j)) - 1;
        extind2{j} = loc2(pc2 == pccommon(j)) - 1;
      end
      
      ext{L,1} = [pccommon',repmat(substri,ncommon,1)];
      ext{L,2} = extind1;
      ext{L,3} = extind2;
    end % if ncommon > 0
  end % for i = 1:nsubstr
  
  
  % were we able to prepend any new characters to extend
  % this list of substrings?
  if (L >= 1) && (substrlength < nstr)
    % yes,so continue prepending
    substrlength = substrlength + 1;
    
    substr = cell2mat(ext(1:L,1));
    ind1 = cat(1,ext{1:L,2});
    ind2 = cat(1,ext{1:L,3});
  else
    % no, so terminate with the previous results
    loc1 = ind1;
    loc2 = ind2;
    flag = false;
  end
  
end % while flag


