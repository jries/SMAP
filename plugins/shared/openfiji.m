function ijm=openfiji(obj)
disp('Due to a Fiji update, the ImageJ-MATLAB plugin does not work properly any more. I am working on fixing it.');

ijclass='IJM';
ijclass='MIJ';

ijm=obj.getPar('IJM');
if isempty(ijm) %open fiji
    
    fijipath=obj.getGlobalSetting('fijipath');  
    if ~exist(fijipath,'dir')
        mij=[];
        errordlg('cannot find Fiji, please select Fiji directory in menu SMAP/Preferences...')
        return
    end
    
    dir=pwd;
    obj.setPar('status','open Fiji');
    if ~isdeployed
        addpath(fijipath)
    end
    try
    ImageJ
    ijm=evalin('base','IJM');
%     Miji();
%     mij=MIJ;
    if ~isdeployed
        cd(dir);
    end
    obj.setPar('IJM',ijm);
    catch err
        warning('ImageJ script from ImageJ-MATLAB could not be opened')
        err
    end
    obj.setPar('IJ',ij.IJ);
%     obj.setPar('MIJ',mij);
end
end