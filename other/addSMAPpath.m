function addSMAPpath

ph=pwd;
ind=strfind(ph,'git');

psmap=[ph(1:ind) 'it' filesep 'SMAP'];
addpath(psmap)

addpath(genpath([psmap filesep 'shared']))
addpath([psmap filesep 'plugins' filesep 'shared'])