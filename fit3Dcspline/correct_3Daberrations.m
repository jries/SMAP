%  Copyright (c)2017 Ries Lab, European Molecular Biology Laboratory,
%  Heidelberg.
%  
%  
%  
%  Additional permission under GNU GPL version 3 section 7
%  
%  If you modify this Program, or any covered work, by
%  linking or combining it with libraries required for interaction
%  with analysis programs such as Igor Pro or Matlab,
%  the licensors of this Program grant you additional permission
%  to convey the resulting work.
%%
function zout=correct_3Daberrations(zcorr,zin,objectivepos)
% zcorr structure returned by calibrate3Daberrations
% zin: uncorrected z (all units: nm)
% objectivepos: position of the focal plane above the coverslip in
% nanometers
% zout: correctec z-positions
zout=zin+zcorr(ones(size(zin))*objectivepos,zin);