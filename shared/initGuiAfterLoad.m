function initGuiAfterLoad(obj, resetview)
if isempty(obj.locData.loc)||isempty(obj.locData.loc.frame)
    return
end
if nargin<2
    loc=obj.locData.getloc({'xnm','ynm'},'layer',find(obj.getPar('sr_layerson')),'Position','fov');
    resetview=isempty(loc.xnm);

end
% disp('init gui after load')
obj.status('grouping');drawnow;
obj.locData.regroup;
% obj.locData.filter;
% fmax=max(obj.locData.loc.frame);
%this calls histogram maker etc
% obj.setPar('frame_max',fmax);
% obj.setPar('frame_min',1);

fl={obj.locData.files.file(:).name};
for k=1:length(fl)
    flh=strrep(fl{k},'\',filesep);
    flh=strrep(flh,'/',filesep);
    [~,fls{k}]=fileparts(flh);
end

obj.setPar('filelist_long',fl,'String');
    obj.setPar('cam_pixelsize_um',obj.locData.files.file(1).info.cam_pixelsize_um)
    obj.setPar('cam_pixelsize_nm',obj.locData.files.file(1).info.cam_pixelsize_um*1000)
% obj.locData.filter;
if resetview

    file=obj.locData.files(1).file;
    info=file.info;
    roi=info.roi;
    sr_pos(1)=(roi(1)+roi(3)/2)*info.cam_pixelsize_um(1)*1000;
    sr_pos(2)=(roi(2)+roi(4)/2)*info.cam_pixelsize_um(end)*1000;  
    sr_size=roi(3:4).*info.cam_pixelsize_um*1000/2;
    obj.setPar('sr_pos',sr_pos);
    obj.setPar('sr_size',sr_size);
    pixrec=sr_size./obj.getPar('sr_imagesize')*2;
    if resetview
        obj.setPar('sr_pixrec',pixrec(1));
    end
% obj.setPar('mainfile',fl{1});
end


locfields=fieldnames(obj.locData.loc);
% if ~isempty(obj.locData.loc.(locfields{1}))
obj.setPar('locFields',locfields,'String');
% end

obj.setPar('filelist_short',fls,'String');
fsx=[{'layer','all'} fls];
obj.setPar('filelist_short_ext',fsx,'String');
            
obj.status('filter');drawnow;
obj.locData.filter;
obj.setPar('currentfileinfo',obj.locData.files.file(1).info) %triggers format to update stuff. and redraw.

%set logLikelyhood? or in channel?
%         sf={fn,s{data.Indices(1),2},s{data.Indices(1),6},s{data.Indices(1),7},true,true};
%         obj.setPar('selectedField',sf,'layer',obj.layer)
        
        
%obj.locData.SE.addFile(obj.locData.files.file(k).name,obj.locData.files.file(k).number,obj.locData.files.file(k).info)

if ~isempty(strfind(obj.getPar('mainfile'),'_sml'))
   tg=obj.getPar('mainGui').guihandles.maintab;
   tg.SelectedTab=tg.Children(3);
end

se=obj.locData.SE;
% fnames={se.files(:).name};
if isempty(se.files) || length(obj.locData.files.file)>length(se.files) %later: test if sites were added already, otherwise run to updata
    for k=1:length(obj.locData.files.file)
%         if ~contains(fnames,obj.locData.files.file(k).name)
            se.addFile(obj.locData.files.file(k).name,obj.locData.files.file(k).number,obj.locData.files.file(k).info)
%         end
    end
end

obj.status('loading done');drawnow;
end

