function ijm=openfiji(obj)
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
    ImageJ
    ijm=evalin('base','IJM');
%     Miji();
%     mij=MIJ;
    if ~isdeployed
        cd(dir);
    end
    obj.setPar('IJM',ijm);
    obj.setPar('IJ',ij.IJ);
%     obj.setPar('MIJ',mij);
end
end