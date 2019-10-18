classdef CombineFiles<interfaces.DialogProcessor
    %combines several loaded files into one file
    methods
        function obj=CombineFiles(varargin)      
            obj@interfaces.DialogProcessor(varargin{:}) ;  
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','CombineFiles');
            notify(obj.P,'backup4undo');
            %adjust channel
            
            dat=p.table.Data;
            use=[dat{:,3}];
            ftarget=find(use,1,'first');
            
            for k=1:length(use)
                if use(k)
                    ind=obj.locData.loc.filenumber==k;
                    if p.addchannel
                        obj.locData.loc.channel(ind)=obj.locData.loc.channel(ind)+dat{k,2};
                    else
                        obj.locData.loc.channel(ind)=dat{k,2};
                    end
                    obj.locData.loc.filenumber(ind)=ftarget;
                    obj.locData.files.file(ftarget).tif=[obj.locData.files.file(ftarget).tif obj.locData.files.file(k).tif];
                    if isfield(obj.locData.files.file(ftarget),'raw')
                        obj.locData.files.file(ftarget).raw=[obj.locData.files.file(ftarget).raw obj.locData.files.file(k).raw];
                    end
                end
            end
            

            %rename file
            obj.locData.files.file(ftarget).name=strrep(obj.locData.files.file(ftarget).name,dat{ftarget,1},p.newfilename);
            %adjust file info (e.g. number of frames, ROI etc), combine tif
            
%              obj.locData.files.file(ftarget).numberOfTif=length(obj.locData.files.file(ftarget).tif);
%             obj.locData.files.file(ftarget).raw=[obj.locData.files.file(ftarget).raw obj.locData.files.file(fnt).raw];
            %etc, if different pixelsize: warn
            %remove second file.
            delete=use;delete(ftarget)=false;
            obj.locData.files.file(delete)=[];
            obj.locData.files.filenumberEnd=obj.locData.files.filenumberEnd-sum(delete);
            %clear SE
            delfiles=sort(find(delete));
            for k=length(delfiles):-1:1
                obj.locData.SE.removeFile(delfiles(k));
            end
            %update filelist
            obj.locData.regroup;
            initGuiAfterLoad(obj)
             
        end
        function initGui(obj)
            fl=obj.getPar('filelist_short');
            numfiles=length(fl.String);
            data=cell(numfiles,3);
            data(:,1)=fl.String;
            data(:,2)=num2cell((1:numfiles)');
            data(:,3)=num2cell(true(numfiles,1));
            pos=obj.guihandles.table.Position;
            ht=uitable(obj.handle,'Position',pos,'Data',data);
            ht.ColumnEditable=[false ,true,true];
%             ht.ColumnFormat{1}=fl.String';
            ht.ColumnName={'File','Channel','include'};
            wt=pos(3);
            ht.ColumnWidth=num2cell(round(wt*[.7 , .1 ,.1]));
            obj.addSynchronization('filelist_short',[],[],{@obj.updatefiles})
            obj.guihandles.table=ht;
            %synch filelist short to table
        end
        function pard=guidef(obj)
            pard=guidef(obj);
        end
        function updatefiles(obj)
        end
    end
end

function newfilenameb_callback(a,b,obj)
p=obj.getGuiParameters;
dat=p.table.Data;

f1=dat{1,1};

fo=['C' num2str(sum([dat{:,3}])) '_' f1];
po.newfilename=fo;
obj.setGuiParameters(po)

end

function pard=guidef(obj)
pard.table.object=struct('String','Channel','Style','text');
pard.table.position=[5,1];
pard.table.Width=4;
pard.table.Height=5;

pard.newfilenameb.object=struct('String','New filename:','Style','pushbutton','Callback',{{@newfilenameb_callback,obj}});
pard.newfilenameb.position=[7,1];

pard.newfilename.object=struct('String','newfile_sml','Style','edit');
pard.newfilename.position=[7,2];
pard.newfilename.Width=3;

pard.addchannel.object=struct('String','add Channel number to existing','Style','checkbox');
pard.addchannel.position=[6,3];
pard.addchannel.Width=2;
pard.plugininfo.description='combines several loaded files into one file';
pard.plugininfo.type='ProcessorPlugin';
end