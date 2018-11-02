function info=minfoparse(minfo)

sstrEV={'ROI"','Evolve-ReadoutRate":','Evolve-Exposure":','Evolve-Gain":','Evolve-Offset":',...
    'Evolve-MultiplierGain":','Evolve-ChipName":','Evolve-Actual Gain e/ADU":','Evolve-Port":','RO'};
sstrAndor={'ROI"','"Andor-ReadoutMode":','"Andor-Exposure":','Andor-Pre-Amp-Gain":','Andor-Offset":',...
     '"Andor-Gain":','"Andor-Camera":','Evolve-Actual Gain e/ADU":','"Andor-Output_Amplifier":'};
sstrED={'ROI"','Evolve Delta-ReadoutRate":','Evolve Delta-Exposure":','Evolve Delta-Gain":','Evolve Delta-Offset":',...
    'Evolve Delta-MultiplierGain":','Evolve Delta-ChipName":','Evolve Delta-Actual Gain e/ADU":','Evolve Delta-Port":','RO'};

 
if ~isempty(strfindfast(minfo,'Delta'))
    disp('Evolve')
    sstr=sstrED;
    cam=1;
elseif ~isempty(strfindfast(minfo,'Evolve'))
    disp('Evolve')
    sstr=sstrEV;
    cam=1;
elseif ~isempty(strfindfast(minfo,'Andor'))
    disp('Andor')
    sstr=sstrAndor;
    cam=1;
else
    disp('cannot read metadata');
    info.readoutrate=1;
    info.exposure=1;
    info.gain=1;
    info.offset=100;
    info.emgain=1; 
    info.actualconv=0;
    cam=0;
    info.port=' ';
    info.timediff=1;
end
 
if cam
for k=1:length(sstr)
    ir=strfindfast(minfo,sstr{k});
    xo{k}=getval(minfo,ir+length(sstr{k})+1);
end
    info.readoutrate=xo{1};
    info.exposure=str2double(xo{2});
    info.emgain=str2double(xo{5});
   
    info.offset=str2double(xo{4});
    
    info.chip=xo{6};
    info.actualconv=str2double(xo{7});
    info.port=xo{8};
    
    if ~isempty(strfindfast(minfo,'Evolve'))
         info.gain=str2double(xo{3});
    elseif ~isempty(strfindfast(minfo,'Andor'))
         info.gain=str2double(xo{3}(5:end));
    end
    
    %rois=xo{9};
    rois=[1 1 1 1];
%     info.roi=str2num(rois(2:end-2));
    info.roi=rois;
    %XXXXXXXXXXXXXXXXXXX 
    ir=strfindfast(minfo,'FrameKey-499-0-0":');
    timefind=499;
    if isempty(ir)
        info.timediff=1;
    else
    smi=minfo(ir+35:ir+550);
    ir=strfind(smi,'ElapsedTime-ms"');
    smi2=smi(ir+16:ir+25);
    timee=str2double(smi2);
%     timee=sscanf(smi,'Time-ms": %f')
    info.timediff=timee/timefind;
    end

%     xo{7}


if isnan(info.exposure)
    info.readoutrate=1;
    info.exposure=1;
    info.gain=1;
    info.offset=1000;
    info.emgain=300;  
    
    info.chip='not determined';
    disp('minfoparse:meta info could not be parsed, camera could not be identified. CameraCalibration.xls')
%     errordlg('no meta info present')

%      info.actualconv=0;
end
%time difference
% ir=strfind(minfo,'ElapsedTime-ms');
% expt=zeros(length(ir),1);
% for k=1:length(ir)
%     t=minfo(ir(k)+16:ir(k)+23);
%     expt(k)=str2double(t);
%     
%     
% end
% info.timediff=max(expt)/length(expt);
% if isempty(info.timediff)
% info.timediff=1;
% end

if isnan(info.offset)
    info.offset=100; %Andor: check!
end
end
    
function txt=getval(minfo,index)


while minfo(index)~='"'
    index=index+1;
end
i1=index+1;
index=index+1;
while minfo(index)~='"'
    index=index+1;
end
txt=minfo(i1:index-1);

function ind=strfindfast(bigs,finds)
ind1=1;
dind=1000000;
lenf=length(finds);
lenbig=length(bigs);
ind2=ind1+dind;
ind=[];
while  ind2<=lenbig+dind
    
    indf=strfind(bigs(ind1:min(ind2,lenbig)),finds);
    if ~isempty(indf)
        ind=ind1+indf(1)-1;
        break;
    end
    ind1=ind2-lenf;
    ind2=ind1+dind;
end
