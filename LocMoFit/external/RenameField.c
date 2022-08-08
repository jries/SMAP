// RenameField.c
// RenameField: Rename a fieldname of a struct
// T = RenameField(S, Old, New)
// INPUT:
//   S: Struct or struct array.
//   Old, New: CHAR vectors or cell string. In the M-Version of this function
//      strings are allowed also. The names must be valid Matlab symbols:
//      <= 63 characters, first character is a letter, others characters are
//      alpha-numeric or the underscore.
//      If a name in Old exist in S, it is renamed to the corresponding name in
//      New. Not existing names are ignored.
// OUTPUT:
//   T: Struct S with renamed fields.
//
// EXAMPLES:
//   S.A = 1; S.B = 2;
//   T = RenameField(S, 'B', 'new');  %  >>  T.A = 1, T.new = 2
//   U = RenameField(S, 'Q', 'C');    %  >>  U.A = 1, U.B = 2
//
// * The C-Mex-function is much faster, but does not accept strings, only CHAR
//   vectors and cell string.
// * Hardcore programmers can omit the validity checks of the new name in the
//   C-Mex function: All names up to 63 ASCII characters are allowed, even '*',
//   ' ' and the empty char vector ''. This does not crash Matlab in my
//   experiments. In R2018b such fields can be accessed with dynamic fieldnames:
//   S.(''), S.('*') works!
// * The check for multiple field names can be omitted also. Dynamic field names
//   pick the first occurence in this case, but e.g. RMFIELD stops with an
//   error. Other functions might crash.
// * This function was created after a discussion in Loren's blog:
//   http://blogs.mathworks.com/
//          loren/2010/05/13/rename-a-field-in-a-structure-array
//
// NOTE: This function was created after some discussions in Loren's blog:
//   http://blogs.mathworks.com/
//        loren/2010/05/13/rename-a-field-in-a-structure-array/
//
// COMPILATION:
//   mex -O RenameField.c
// Omit validity check of names (funny, but cruel):
//   mex -DNO_NAME_CHECK -O RenameField.c
// If the undocumented SharedDataCopy does not work anymore:
//   mex -O RenameField.c -DNO_SHAREDDATACOPY
// Consider C99 comments on Linux:
//   mex -O CFLAGS="\$CFLAGS -std=c99" RenameField.c
// Pre-compiled Mex: http://www.n-simon.de/mex
// Run the unit-test uTest_RenameField after compiling.
//
// Tested: Matlab 2009a, 2015b(32/64), 2016b, 2018b, Win7/10
//         Compiler: WinSDK7.1, MSVC 2008/2010/2017
// Assumed Compatibility: higher Matlab versions, Mac, Linux
// Author: Jan Simon, Heidelberg, (C) 2006-2022 matlab.2010(a)n(MINUS)simon.de
//
// See also CELL2STRUCT, STRUCT, GENVARNAME, RMFIELD.

/*
% $JRev: R-v V:021 Sum:4lJJm6RZYqY3 Date:09-Jun-2022 01:08:18 $
% $License: BSD (use/copy/change/redistribute on own risk, mention the author) $
% $UnitTest: uTest_RenameField $
% $File: Tools\Mex\Source\RenameField.c $
% History:
% 001: 16-Aug-2010 21:43, Created after Loren's blog.
% 009: 12-Feb-2011 00:45, Cell string input.
% 018: 14-Jan-2017 14:41, STRCMP=> safer STRNCMP. Nicer error messages.
%      mxCreateSharedDataCopy=> faster mxCreateReference (Thanks to James Tursa)
%      Not tested with R6.5 anymore.
% 020: 31-Mar-2019, Back to mxCreateSharedDataCopy.
%      mxCreateReference was not faster to my surprise.
*/

// =============================================================================
#include "mex.h"
#include <string.h>

// Error messages do not contain the function name in Matlab 6.5! This is not
// necessary in Matlab 7, but it does not bother:
#define ERR_ID   "JSimon:RenameField:"
#define ERR_HEAD "*** RenameField[mex]: "
#define ERROR_2(id,msg)   mexErrMsgIdAndTxt(ERR_ID id, ERR_HEAD msg);
#define ERROR_3(id,msg,p) mexErrMsgIdAndTxt(ERR_ID id, ERR_HEAD msg, p);

// There is an undocumented method to create a shared data copy. This is much
// faster, if the replied object is not changed, because it does not duplicate
// the contents of the array in the memory.
#if !defined(NO_SHAREDDATACOPY)
   mxArray *mxCreateSharedDataCopy(const mxArray *mx);
#  define COPY_ARRAY mxCreateSharedDataCopy
#else
#  define COPY_ARRAY mxDuplicateArray     // slower, but documented
#endif

// Allow the disabling of checking for bad names:
#if !defined(NO_NAME_CHECK)
#define CHECK_SYMBOL 1
#else
#define CHECK_SYMBOL 0
#endif

// Prototypes:
mxLogical CheckSymbol(const char *S);
void GetString(const mxArray *S, char *C);
mxArray *StringMethod(const mxArray *S,
                    const mxArray *OldM, const mxArray *NewM);
mxArray *CellMethod(const mxArray *S,
                    const mxArray *OldC, const mxArray *NewC);
int InsertNewName(const char **FieldList, mwSize nField, char *Old, char *New);

// Main function ===============================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  const mxArray *S;
  
  // Check number and type of arguments:
  if (nrhs != 3) {
     ERROR_2("BadNInput", "3 inputs required.");
  }
  if (nlhs > 1) {
     ERROR_2("BadNOutput", "1 output allowed.");
  }
  
  // Type of input arguments: Struct or empty matrix
  S = prhs[0];
  if (!mxIsStruct(S)) {
     // Allow empty matrix as input - nothing to rename:
     if (mxIsDouble(S) && mxGetNumberOfElements(S) == 0) {
        plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
        return;
     }
     ERROR_2("BadTypeInput1", "1st input must be a struct.");
  }
  
  // Replace field names:
  if (mxIsChar(prhs[1]) && mxIsChar(prhs[2])) {
     plhs[0] = StringMethod(S, prhs[1], prhs[2]);
     
  } else if (mxIsCell(prhs[1]) && mxIsCell(prhs[2])) {
     if (mxGetNumberOfElements(prhs[1]) != mxGetNumberOfElements(prhs[2])) {
         ERROR_2("BadInputType",
                 "[Old] and [New] must have the same numer of elements.");
     }

     plhs[0] = CellMethod(S, prhs[1], prhs[2]);
     
  } else {  // Bad input:
     ERROR_2("BadInputType",
             "[Old] and [New] must be strings or cell strings.");
  }
  
  return;
}

// =============================================================================
mxLogical CheckSymbol(const char *S)
{
  // Reply TRUE if the string is a valid symbol.
  // The 1st character must be a letter, the following alphanumerical or the
  // undescore. The string length must not exceed 63 characters.
  
  // This function is published as "isValidSymbol" in the FEX.
   
  static char Table[256], Filled = 0;
  unsigned char i, maxLen = mxMAXNAM - 1;
  
  // Create a lookup table for all byte values:
  if (!Filled) {
     memset(Table, (char) 0, 256);
     
     // Letters:
     for (i = 'A'; i <= 'Z'; Table[i++] = 1) ;   // empty loop
     for (i = 'a'; i <= 'z'; Table[i++] = 1) ;   // empty loop
     
     // Numbers and underscore:
     for (i = '0'; i <= '9'; Table[i++] = 2) ;   // empty loop
     Table['_'] = 2;
     
     Filled = 1;
  }
     
  // First letter must be a letter:
  if (Table[*S] != 1) {
     return (false);
  }

  // The other characters must be alpha-numerical or the underscore:
  for (i = 1; i < maxLen; i++) {   // i < 63
     if (Table[S[i]] == 0) {
        return ((mxLogical) (S[i] == '\0'));
     }
  }
  
  return (false);
}


// =============================================================================
void GetString(const mxArray *S, char *C)
{
  // Copy contents of Matlab String to C-String.
  
  // Check for unitialized cell element:
  if (S == NULL) {
     ERROR_2("ElementIsNULL", "Uninitialized cell element.");
  }
  
  if (mxGetString(S, C, mxMAXNAM)) {
     if (mxGetNumberOfElements(S) > mxMAXNAM) {
        ERROR_2("NameTooLong", "Field name has too many characters.");
     } else if (!mxIsChar(S)) {
        ERROR_2("NoString", "Field name must be a string.");
     } else {
        ERROR_2("NoMemory", "Cannot convert field name to C-string.");
     }
  }
  
  return;
}

// =============================================================================
mxArray *StringMethod(const mxArray *S,
                      const mxArray *OldM, const mxArray *NewM)
{
  // [OldM] and [NewM] are CHAR vectors.
  
  char       Old[mxMAXNAM], New[mxMAXNAM];   // mxMAXNAM = 64, with terminator
  const char **FieldList;
  mwSize     FieldNumber, FieldNumberNew, iField, nField, nElem, iElem;
  mxArray    *R;

  // Dimension of the struct:
  nField = mxGetNumberOfFields(S);
  nElem  = mxGetNumberOfElements(S);
  
  // No fields => nothing to rename:
  if (nField == 0) {
     return COPY_ARRAY(S);
  }
  
  // Get list of pointers to field names:
  FieldList = (const char **) mxMalloc(nField * sizeof(char *));
  if (FieldList == NULL) {
     ERROR_2("NoMemory", "Cannot get memory for FieldList.");
  }
  for (iField = 0; iField < nField; iField++) {
     FieldList[iField] = mxGetFieldNameByNumber(S, iField);
  }
  
  // Obtain old and new field name (Matlab string -> C string):
  GetString(OldM, Old);
  GetString(NewM, New);
  
  // Check the validity new name:
  // The length of the name is still limited to 63 characters, but you can
  // rename a field to '%', ' ' or even '', if you omit the check. Then these
  // fields can be access by dynamic field names only: S.('') works!
# if CHECK_SYMBOL
  if (!CheckSymbol(New)) {
     ERROR_2("InvalidSymbol", "[New] is no valid Matlab symbol.");
  }
# endif

  // Nothing to do if Old is not a field of S:
  if ((FieldNumber = mxGetFieldNumber(S, Old)) == -1) {
     return COPY_ARRAY(S);
  }
  
  // Check if new name is existing already, but accept New==Old:
  FieldNumberNew = mxGetFieldNumber(S, New);
  if (FieldNumberNew != -1) {
     if (FieldNumberNew != FieldNumber) {
        ERROR_3("NewNameExisting", "S.('%s') is existing already.", New);
     }
     
     // New equals Old - reply the unchanged input struct:
     return COPY_ARRAY(S);
  }
     
  // Get list of field names:
  FieldList = (const char **) mxMalloc(nField * sizeof(char *));
  if (FieldList == NULL) {
     ERROR_2("NoMemory", "Cannot get memory for FieldList.");
  }
     
  // Copy pointers to field names and replace the found one:
  for (iField = 0; iField < FieldNumber; iField++) {
     FieldList[iField] = mxGetFieldNameByNumber(S, iField);
  }
  
  FieldList[FieldNumber] = New;
  
  for (iField = FieldNumber + 1; iField < nField; iField++) {
     FieldList[iField] = mxGetFieldNameByNumber(S, iField);
  }
  
  // Create the output struct:
  R = mxCreateStructArray(mxGetNumberOfDimensions(S),
                          mxGetDimensions(S), nField, FieldList);
                          
  // Copy fields for each element of the struct array S:
  for (iElem = 0; iElem < nElem; iElem++) {
     for (iField = 0; iField < nField; iField++) {
        mxSetFieldByNumber(R, iElem, iField,
                           COPY_ARRAY(mxGetFieldByNumber(S, iElem, iField)));
     }
  }
  
  // Cleanup:
  mxFree(FieldList);

  return R;
}

// =============================================================================
mxArray *CellMethod(const mxArray *S,
                    const mxArray *OldC, const mxArray *NewC)
{
  // [OldC] and [NewC] are cell strings of the same size.

  char       Old[mxMAXNAM],   // mxMAXNAM = 64, with terminator
             *NewNameChar, *aNew;
  const char **FieldList;
  mwSize     nC, iC, iField, nField, nElem, iElem;
  mxArray    *R;
  
  // Dimension of the struct:
  nField = mxGetNumberOfFields(S);
  nElem  = mxGetNumberOfElements(S);
  
  // Number of strings:
  nC = mxGetNumberOfElements(OldC);
  
  // No fields => nothing to rename:
  if (nField == 0 || nC == 0) {
     return COPY_ARRAY(S);
  }
  
  // Get list of pointers to field names of the struct:
  FieldList = (const char **) mxMalloc(nField * sizeof(char *));
  if (FieldList == NULL) {
     ERROR_2("NoMemory", "Cannot get memory for FieldList.");
  }
  for (iField = 0; iField < nField; iField++) {
     FieldList[iField] = mxGetFieldNameByNumber(S, iField);
  }
  
  // Create array to carry all new field names:
  NewNameChar = (char *) mxMalloc(nC * mxMAXNAM * sizeof(char));
  if (NewNameChar == NULL) {
     ERROR_2("NoMemory", "Cannot get memory for NewNameChar.");
  }
  
  // Loop over elements of the cell strings [OldC] and [NewC]:
  aNew = NewNameChar;
  for (iC = 0; iC < nC; iC++) {
     // Obtain old and new field name (Matlab string -> C string):
     GetString(mxGetCell(OldC, iC), Old);
     GetString(mxGetCell(NewC, iC), aNew);
     
     // Check the validity new name:
#    if CHECK_SYMBOL
        if (!CheckSymbol(aNew)) {
           ERROR_2("InvalidSymbol",
                   "New field name is no valid Matlab symbol.");
        }
#    endif
     
     // Replace matching field name:
     if (InsertNewName(FieldList, nField, Old, aNew)) {
        aNew += mxMAXNAM;
     }
  }  // end: for(iC)
  
  // No new fields:
  if (aNew == NewNameChar) {
     return COPY_ARRAY(S);
  }
  
  // Create the output struct:
  R = mxCreateStructArray(mxGetNumberOfDimensions(S),
                          mxGetDimensions(S), nField, FieldList);
                          
  // Copy fields for each element of the struct array S:
  for (iElem = 0; iElem < nElem; iElem++) {
     for (iField = 0; iField < nField; iField++) {
        mxSetFieldByNumber(R, iElem, iField,
                           COPY_ARRAY(mxGetFieldByNumber(S, iElem, iField)));
     }
  }
  
  // Cleanup:
  mxFree(FieldList);
  mxFree(NewNameChar);
  
  return R;
}

// =============================================================================
int InsertNewName(const char **FieldList, mwSize nField, char *Old, char *New)
{
  // Search for the Old name in the list of field names and insert the pointer
  // to the New name in case of a match. Check if the new name equals any
  // already existing name.
  // NOTE: mxMAXNAM is 64, which is 63 charcatcers and a string terminator \0.
  mwSize i;
  int    Reply = 0;
  
  for (i = 0; i < nField; i++) {
     // Replace matching field name:
     if (strncmp(FieldList[i], Old, mxMAXNAM) == 0) {
        FieldList[i] = New;
        Reply = 1;
        break;                // BREAK: for (i)
     }
     
     // Reject existing field name:
     if (strncmp(FieldList[i], New, mxMAXNAM) == 0) {
        ERROR_3("MultipleName", "Field [%s] is existing already.", New);
     }
  }
  
  // Check later field names for already existing name:
  while (++i < nField) {
     if (strncmp(FieldList[i], New, mxMAXNAM) == 0) {
        ERROR_3("MultipleName",
                "Field [%s] is existing already.", New);
     }
  }
    
  return Reply;
}
