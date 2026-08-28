classdef export_tiffs<interfaces.DialogProcessor&interfaces.SEProcessor
%     export superresolution reconstructions of selected ROIs
    methods
        function obj=export_tiffs(varargin)        
                obj@interfaces.DialogProcessor(varargin{:});
            obj.inputParameters={};
        end
        
        function out=run(obj,p)  
            out=[];
            sites=obj.SE.sites;
            keep=obj.getPar('se_keeptempimages');
            obj.setPar('se_keeptempimages',true);
            if p.export_selected
                sev=obj.getPar('se_viewer');
                pv=sev.getAllParameters;
                selectedsites=pv.sitelist.Value;
            else
              
                selectedsites=1:length(sites);
            end
            mainfile=obj.getPar('mainfile');
            path=fileparts(mainfile);
            prefix='img_.tif';
            [f,path]=uiputfile([path filesep prefix]);
            
            if ~f
                return
            end
            [~,f]=fileparts(f);

        
     
            for k=selectedsites
                site=sites(k);
                imold=site.image.image;
                site.image=[];
                site.image=obj.SE.plotsite(site,-1);
                site.image.image=imold;
                options.color=true;
                options.comp='lzw';
                pixelsize=site.image.parameters.sr_pixrec;
                xResolution = 1 / pixelsize;
                tags.ResolutionUnit = Tiff.ResolutionUnit.Centimeter;
                tags.XResolution = xResolution * 10000000;
                tags.YResolution = xResolution * 10000000;
                

                filen=[f strrep(site.name,'.','_')];
                fhere= [filen '.tif'];
                if p.export_scalebar
                    imsel=site.image.image;
                else
                    imsel=site.image.composite;
                end
                imout=uint8(imsel*255);
                
                saveastiff(imout,[path fhere],options,tags);
                if length(site.image.layers)>1
                    for ll=1:length(site.image.layers)
                        
                        iml=site.image.layers(ll).images.renderimages.image;
                        filell= [filen '_' int2str(ll) '.tif'];
                        imoutll=uint8(iml*255);
                        saveastiff(imoutll,[path filell],options);
                    end
                end
                site.image.layers=[];site.image.composite=[];
            end
            
            if p.export_cells
                cells=obj.SE.cells;
                for k=1:length(cells)
                    cell=cells(k);
                    imold=cell.image.image;
                    cell.image=[];
                    cell.image=obj.SE.plotsite(cell,-1);
                    cell.image.image=imold;
                    options.color=true;
                    options.comp='lzw';
                    filen=[f 'cell_C' num2str(cell.ID)  '_F' num2str(cell.info.filenumber)];
                    fhere= [filen '.tif'];
                    imout=uint8(cell.image.image*255);

                    saveastiff(imout,[path fhere],options);
                end
            end
            
            if p.export_files
                files=obj.SE.files;
                for k=1:length(files)
                    file=files(k);
       
                    options.color=true;
                    options.comp='lzw';
                    filen=[f  '_F' num2str(file.ID)];
                    fhere= [filen '.tif'];
                    imout=uint8(file.image.image*255);

                    saveastiff(imout,[path fhere],options);
                end
             end
            
            out=0;
            obj.setPar('se_keeptempimages',keep);
        end
        function pard=guidef(obj)
            pard=guidef;
        end
    end
end




function pard=guidef

pard.export_selected.object=struct('String','export only selected sites','Style','checkbox');
pard.export_selected.position=[1,1];
pard.export_selected.Width=2;

pard.export_cells.object=struct('String','export cell images as well','Style','checkbox');
pard.export_cells.position=[2,1];
pard.export_cells.Width=2;

pard.export_files.object=struct('String','export file images as well','Style','checkbox');
pard.export_files.position=[3,1];
pard.export_files.Width=2;

pard.export_scalebar.object=struct('String','scale bar and LUT','Style','checkbox','Value',1);
pard.export_scalebar.position=[1,3];
pard.export_scalebar.Width=2;



pard.plugininfo.type='ROI_Analyze';
pard.plugininfo.description='export superresolution reconstructions of selected ROIs';

end