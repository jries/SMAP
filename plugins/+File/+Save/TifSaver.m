classdef TifSaver<interfaces.DialogProcessor
    methods
        function obj=TifSaver(varargin)  
                obj@interfaces.DialogProcessor(varargin{:}) ;
                obj.inputParameters={'filelist_long','mainfile','mainGui','numberOfLayers','sr_layerson','sr_pixrec','layers','sr_image','sr_pos','group_dt','group_dx'};
        end
        
        function out=save(obj,p)
            p2=obj.getGuiParameters;
            obj.status('save tif file')
            fn=p.filelist_long.selection;
            [path,file]=fileparts(fn);
            of=[path filesep file p2.img_ext.selection];
            ind=2;
            while exist(of,'file')
                of=[path filesep file '_' num2str(ind) p2.img_ext.selection];
                ind=ind+1;
            end
            
            
            serchstr={['*' p2.img_ext.selection];['*' strjoin(p2.img_ext.String,';*')]};
            [f,path]=uiputfile(serchstr,'select output file for image', of);
            if f
                img=obj.getPar('sr_image');
                res=ones(2,1)/p.sr_pixrec*2.5e7/1e6;%rescale to macroscopic values, otherwise size in PPT is zero  
                p.scalebarnm=img.scalebarnm;
                txt=filterpar(obj,p);
                description=sprintf(txt);
                [~,~,ext]=fileparts(f);
%                 switch ext
%                     case '.tif'
                        imwrite(img.image,[path f],'Description',description);
%                     case '.png'
%                          imwrite(img.image,[path f],'Description',description,'XResolution',res(1),'YResolution',res(2),'ResolutionUnit','inch');
%                         imwrite(img.image,[path f],'Description',description,'XResolution',res(1),'YResolution',res(2),'ResolutionUnit','inch');
%                 end
            end
            obj.status('save done')
          
        end
        function pard=guidef(obj)
           pard.plugininfo.type='SaverPlugin';
           
            pard.img_ext.object=struct('Style','popupmenu','String',{{'.tif','.png'}});
            pard.img_ext.position=[1,1];
            pard.img_ext.Width=2;
        end
        function run(obj,p)
            obj.save(p)
        end        

    end
end

function txt=filterpar(obj,p)
txt='SMAP \n';
txt=[txt 'pixelsize(nm) \t' num2str(p.sr_pixrec) '\n'];
txt=[txt 'scalebar (nm) \t' num2str(p.scalebarnm) '\n'];
txt=[txt 'position (nm) \t' num2str(p.sr_pos (1:2)) '\n'];
txt=[txt 'group_dx (nm) \t' num2str(p.group_dx) '\n'];
txt=[txt 'group_dt (frames) \t' num2str(p.group_dt) '\n'];
for k=1:length(p.sr_layerson)
    if p.sr_layerson(k)
        lp=['layer' num2str(k)];
        txt=[txt lp ':\n'];
        txt=[txt p.([lp '_']).ch_filelist.selection '\n'];
        
        filn=p.([lp '_']).ch_filelist.Value;
        txt=[txt p.filelist_long.String{filn} '\n'];
        try
            fn=obj.locData.files.file(filn).info.filename;
            fn=strrep(fn,'\','/');
            txt=[txt fn '\n'];
        catch
        end
        txt=[txt 'rendermode: ' num2str(p.([lp '_']).rendermode.selection) '\n'];
        txt=[txt 'channels: ' num2str(p.([lp '_']).channels) '\n'];
        txt=[txt 'grouping: ' num2str(p.([lp '_']).groupcheck) '\n'];
        txt=[txt 'lut: ' p.([lp '_']).lut.selection '\n'];
        txt=[txt 'quantile/Imax: ' num2str(p.([lp '_']).imax_min) '\n'];
        txt=[txt 'color range: \t' num2str(p.([lp '_']).colorfield_min) ' : \t' num2str(p.([lp '_']).colorfield_max) '\n'];
        txt=[txt 'remove outside c-range: ' num2str(p.([lp '_']).remout) '\n'];
        if strcmp(p.([lp '_']).renderfield.selection,'field')
            txt=[txt 'render field: ' (p.([lp '_']).render_colormode.selection) '\n'];
        end
        if strcmp(p.([lp '_']).renderfield.selection,'z')
            txt=[txt 'render field: znm'  '\n'];
        end
        txt=[txt 'shift x,y,z (nm): ' num2str(p.([lp '_']).shiftxy_min) ', ' num2str(p.([lp '_']).shiftxy_max)  ', ' num2str(p.([lp '_']).shiftxy_z) '\n'];
        s=obj.getPar(['layer' num2str(k) '_filtertable']);
        ss=size(s);
        txt=[txt 'filters:'  '\n'];
        for l=1:ss(1)
            if s{l,7} &&~strcmp(s{l,1},'xnm')&&~strcmp(s{l,1},'ynm')&&~strcmp(s{l,1},'filenumber')
                line=s(l,1:6);
                th=[line{1} ':\t' num2str(line{2}) ' : \t' num2str(line{6}) '\n'];
                txt=[txt th];
            end
            
        end
        
        txt=[txt 'min gauss nm: ' num2str(p.([lp '_']).mingaussnm) '\n'];
        txt=[txt 'min gauss pix: ' num2str(p.([lp '_']).mingausspix) '\n'];
        txt=[txt 'gauss factor: ' num2str(p.([lp '_']).gaussfac) '\n'];
        txt=[txt 'gamma: ' num2str(p.([lp '_']).gamma) '\n'];
%         disp(sprintf(txt));
    end
    
end
end