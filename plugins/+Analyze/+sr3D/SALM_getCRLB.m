classdef SALM_getCRLB<interfaces.DialogProcessor
%     calclates CRLB for photon counts, necessary for zSALM error
%     calculation
    properties
        
    end
    methods
        function obj=SALM_getCRLB(varargin)   
            obj@interfaces.DialogProcessor(varargin{:}) ;
        end
        function out=run(obj,p)
            out=[];
            % get fields
            fsa=p.assignfield1.selection;
            fua=p.assignfield2.selection;
            fsabg=p.bgfield1.selection;
            fsubg=p.bgfield2.selection;
            
            if isfield(obj.locData.loc,'znm_original')
                zfield='znm_original';
            elseif isfield(obj.locData.loc,'znm_a')
                zfield='znm_a';
            else 
                zfield='znm';
            end
            
            Ns=obj.locData.loc.(fsa);
            Nu=obj.locData.loc.(fua);
            bgs=obj.locData.loc.(fsabg);
            bgu=obj.locData.loc.(fsubg);
            zas=obj.locData.loc.(zfield);
            % rescale back to original units without refractive index
            % factor to get right CRLB
            if obj.locData.files.file.savefit.fitparameters.MLE_GPU_Yiming.userefractive_index_mismatch;
                rif=obj.locData.files.file.savefit.fitparameters.MLE_GPU_Yiming.refractive_index_mismatch;
            else 
                rif=1;
            end
            zas=zas/rif; 
        
            psf_sa=splinePSF;
            psf_sa.loadmodel(p.cal_3Dfile,2);
            
            psf_ua=splinePSF;
            psf_ua.loadmodel(p.cal_3Dfile,1);          
            
            roi_ua=obj.locData.files.file.savefit.fitparameters.RoiCutterWF.loc_ROIsize;
            crlbu=(psf_ua.crlb(Nu,bgu,zas,roi_ua));
            crlbs=(psf_sa.crlb(Ns,bgs,zas,roi_ua));
            Nserr=sqrt(crlbs(:,3));
            Nuerr=sqrt(crlbu(:,3));
            obj.locData.setloc([fsa 'err'], single(Nserr));
            obj.locData.setloc([fua 'err'], single(Nuerr));
            obj.locData.setloc('znm_a', single(zas*p.RIFHighNA));
            obj.locData.setloc('locprecznm_a',obj.locData.loc.locprecznm*p.RIFHighNA/rif);
            obj.locData.regroup;       
        end
        
        function initGuiFinal(obj)  
            fields={'assignfield1','assignfield2','bgfield1','bgfield2'};
            fieldnames={'psf_ns','psf_nu','psf_bgs','psf_bgu'};
            for k=1:length(fields)
                stri1=obj.guihandles.(fields{k}).String;
                indg=find(contains(stri1,fieldnames{k}));
                if ~isempty(indg)
                    obj.guihandles.(fields{k}).Value=indg(1);
                end
            end   
        end
        
        function pard=guidef(obj)
            pard=guidef(obj);
        end
    end
end



function loadcall_callback(a,b,obj)
p=obj.getAllParameters;
if isempty(p.cal_3Dfile)
    path=fileparts(obj.getPar('lastSMLFile'));
%     path=obj.getGlobalSetting('DataDirectory');
%     fh=obj.getPar('loc_fileinfo');
%     if ~isempty(fh) && ~isempty(fh.imagefile)
%         path=fileparts(fh.imagefile);
%     end  
    p.cal_3Dfile=[path filesep '*3dcal.mat'];
end
[f,p]=uigetfile(p.cal_3Dfile);
if f
    l=load([p f]);
    if ~isfield(l,'outforfit') && ~isfield(l,'SXY') && ~isfield(l,'cspline')
        msgbox('no 3D data recognized. Select other file.');
    end
    obj.setGuiParameters(struct('cal_3Dfile',[p f]));
    obj.setPar('cal_3Dfile',[p f]);
    
end
end

function pard=guidef(obj)


pard.t1.object=struct('String','Intensity SA','Style','text');
pard.t1.position=[1,1];
pard.t2.object=struct('String','Intensity UA','Style','text');
pard.t2.position=[2,1];

pard.assignfield1.object=struct('Style','popupmenu','String','');
pard.assignfield1.position=[1,2];
pard.assignfield2.object=struct('Style','popupmenu','String','');
pard.assignfield2.position=[2,2];

pard.t3.object=struct('String','bg SA','Style','text');
pard.t3.position=[1,3.5];
pard.t4.object=struct('String','bg UA','Style','text');
pard.t4.position=[2,3.5];

pard.bgfield1.object=struct('Style','popupmenu','String','');
pard.bgfield1.position=[1,4];
pard.bgfield2.object=struct('Style','popupmenu','String','');
pard.bgfield2.position=[2,4];

pard.rift.object=struct('String','refractive index mismatch factor new high NA','Style','text');
pard.rift.position=[3,1];
pard.rift.Width=2.5;
pard.RIFHighNA.object=struct('Style','edit','String','0.6');
pard.RIFHighNA.position=[3,3.5];
pard.RIFHighNA.Width=0.5;

pard.loadcal.object=struct('Style','pushbutton','String','Load 3D cal','Callback',{{@loadcall_callback,obj}});
pard.loadcal.position=[4,1];
pard.loadcal.Width=.75;
pard.cal_3Dfile.object=struct('Style','edit','String','');
pard.cal_3Dfile.position=[4,1.75];
pard.cal_3Dfile.Width=3.25;
pard.cal_3Dfile.TooltipString=sprintf('3D calibration file for astigmtic 3D. \n Generate from bead stacks with plugin: Analyze/sr3D/CalibrateAstig');


 pard.syncParameters={{'locFields','assignfield1',{'String'}},{'locFields','assignfield2',{'String'}},...
     {'locFields','bgfield1',{'String'}},{'locFields','bgfield2',{'String'}},{'cal_3Dfile','cal_3Dfile',{'String'}}};
            
%             obj.addSynchronization('locFields',[],[],@obj.updateLocFields)
%             obj.addSynchronization('filelist_short',obj.guihandles.dataselect,'String')
pard.plugininfo.type='ProcessorPlugin';

pard.plugininfo.description='converts reltive intensities of super- and undercritical channel into z position.';
end

function err=zerrSALM(fitp,Ns,Nu,bgs,bgu)
indb=Ns<1 | Nu<10;
% Ns(Ns<1)=1;
% Nu(Nu<1)=1;
c=max(fitp.c,0);
r=Ns./Nu;
dru=-Ns./Nu.^2;
drs=1./Ns;
dR2=drs.^2.*errN2(Ns,bgs)+dru.^2.*errN2(Nu,bgu);
dz2=dR2./(-fitp.b.*(r-c)).^2;
dz2(r<c)=inf;
dz2(indb)=inf;
% dz2=(a./r).^2.*dR2;
% dz2=dR2/b^2./r.^2;
err=sqrt(dz2);
end
function err=zerrSALMspline(obj,dz_dr,Ns,Nu,bgs,bgu,zas)
% XXX replace by CRLB from splinePSF
% psf_sa=splinePSF;
% psf_sa.loadmodel(cal_3Dfile,2);
% 
% psf_ua=splinePSF;
% psf_ua.loadmodel(cal_3Dfile,1);

% roi_ua=obj.locData.files.file.savefit.fitparameters.RoiCutterWF.loc_ROIsize;
% errNu2=(psf_ua.crlb(Nu,bgu,zas,roi_ua));

indb=Ns<1 | Nu<10;

% c=max(fitp.c,0);
r=Ns./Nu;
dru=-Ns./Nu.^2;
drs=1./Ns;
dR2=drs.^2.*errN2(Ns,bgs)+dru.^2.*errN2(Nu,bgu);

dz2=dz_dr(r).^2.*dR2;
% dz2=dR2./(-fitp.b.*(r-c)).^2;
% dz2(r<c)=inf;
dz2(indb)=inf;
% dz2=(a./r).^2.*dR2;
% dz2=dR2/b^2./r.^2;
err=sqrt(dz2);
end

function dN2=errN2(N,bg)
bg(bg<0)=0;
s_a=1; %sigmapsf/pixelsize
tau=2*pi*(bg)*(s_a^2+1/12)./N;
dN2=N.*(1+4*tau+sqrt(tau./(14*(1+2*tau))));
end