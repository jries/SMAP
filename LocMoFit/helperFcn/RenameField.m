function S = RenameField(S, Old, New)
% RenameField - Rename a field of a struct
% T = RenameField(S, Old, New)
% INPUT:
%  INPUT:
%    S:    Struct or struct array.
%    Old:  String or cell string, name of the fields to be renamed. If Old is
%          not existing in S, the output T equals the input S.
%    New:  String or cell string, new field name, which must be a valid Matlab
%          symbol: up to 63 characters, first character is a letter, the others
%          are alpha-numeric or the underscore.
%  OUTPUT:
%    T:    Struct S with renamed fields.
%
% EXAMPLES:
%   S.A = 1; S.B = 2;
%   T = RenameField(S, 'B', 'C');  %  >>  T.A = 1, T.C = 2
%
% NOTE: Hardcore programmers can omit the validity checks of the new name.
%   Then all names with up to 63 characters are allowed. Although this does
%   not crash Matlab, the effects can be rather strange: you can rename a
%   field to '*', ' ' and even ''. Such fields can be accessed by dynamic field
%   names: S.('') works!
%   If the checking is disabled, the field names are not necessarily unqiue
%   Then the dynamic field name access picks the first occurence of a name. But
%   e.g. RMFIELD will stop with an error. But other functions might crash.
%
% NOTE: This function was created after some discussions in Loren's blog:
%   http://blogs.mathworks.com/
%          loren/2010/05/13/rename-a-field-in-a-structure-array
%
% COMPILATION: See RenameField.c
% Run uTest_RenameField to check validity and speed of the Mex function.
%
% Tested: Matlab 6.5, 7.7, 7.8, WinXP
% Author: Jan Simon, Heidelberg, (C) 2006-2011 matlab.THISYEAR(a)nMINUSsimon.de
%
% See also CELL2STRUCT, STRUCT, GENVARNAME, RMFIELD.

% $JRev: R-f V:005 Sum:rsC5ZfSTA7dM Date:11-Feb-2011 00:17:01 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $File: Tools\GLStruct\RenameField.m $
% History:
% 001: 19-Aug-2010 00:11, Created after discussion in Loren's blog.

% Initialize: ==================================================================
% Do the work: =================================================================

% This is an implementation as M-code. Prefer the mex file, which is 50% (S has
% 1 field only) to 95% (S has 1000 fields) faster.

% Comment this out, if you want to use the M-version:
error(['JSimon:', mfilename, ':NoMex'], 'Cannot find compiled Mex file!');

% Under Matlab 6.5 CELL2STRUCT accpets names with more than 63 characters. But
% the later recognition fails!

if isempty(S) && isa(S, 'double')  % Accept [] as empty struct without fields
   return;
end

Data  = struct2cell(S);
Field = fieldnames(S);
if ischar(Old)
   Field(strcmp(Field, Old)) = {New};
elseif iscellstr(Old)
   for iField = 1:numel(Old)
      Field(strcmp(Field, Old{iField})) = New(iField);
   end
else
   error(['JSimon:', mfilename, ':BadInputType'], ...
      'Fields must be a string or cell string!');
end

S = cell2struct(Data, Field);

return;
