classdef CombineFiles<interfaces.DialogProcessor
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
            
%             fn1=p.file1.Value;
%             fn2=p.file2.Value;
%             ftarget=min(fn2,fn1);
%             fnt=max(fn2,fn1);
%             
%             ind1=obj.locData.loc.filenumber==fn1;
%             ind2=obj.locData.loc.filenumber==fn2;
%             if p.addchannel1
%                 obj.locData.loc.channel(ind1)=obj.locData.loc.channel(ind1)+p.channel1;
%             else
%                 obj.locData.loc.channel(ind1)=p.channel1;
%             end
%             if p.addchannel2
%                 obj.locData.loc.channel(ind2)=obj.locData.loc.channel(ind2)+p.channel2;
%             else
%                 obj.locData.loc.channel(ind2)=p.channel2;
%             end
%             
%             % adjust filenumber
%             obj.locData.loc.filenumber(ind2)=ftarget;
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
% f2=dat{2,1};
% re='[^_]*';
% [i1,i2]=regexp(f1,re);
% for k=1:length(i1)
%     if i2(k)-i1(k)>2
%     f2=strrep(f2,f1(i1(k):i2(k)),'');
%     end
% end
% 
% fo=['C' num2str(dat{1,2}) '_' f2 'C' num2str(dat{2,2}) '_' f1];
% fo=regexprep(fo,'[_]*','_');
fo=['C' num2str(sum([dat{:,3}])) '_' f1];
po.newfilename=fo;
obj.setGuiParameters(po)

end

function pard=guidef(obj)
pard.table.object=struct('String','Channel','Style','text');
pard.table.position=[5,1];
pard.table.Width=4;
pard.table.Height=5;
% pard.t1.object=struct('String','Channel','Style','text');
% pard.t1.position=[1,4];
% pard.t1.Width=0.5;
% 
% pard.t3.object=struct('String','add ','Style','text');
% pard.t3.position=[1,4.5];
% pard.t3.Width=0.5;
% 
% pard.file1.object=struct('Style','popupmenu','String','File');
% pard.file1.position=[2,1];
% pard.file1.object.TooltipString='choose localization file data set';
% pard.file1.Width=3;
% 
% pard.file2.object=struct('Style','popupmenu','String','File');
% pard.file2.position=[3,1];
% pard.file2.object.TooltipString='choose localization file data set';
% pard.file2.Width=3;
% 
% 
% pard.channel1.object=struct('String','1','Style','edit');
% pard.channel1.position=[2,4];
% pard.channel1.Width=0.5;
% 
% pard.channel2.object=struct('String','2','Style','edit');
% pard.channel2.position=[3,4];
% pard.channel2.Width=0.5;
% 
% pard.addchannel1.object=struct('String','','Style','checkbox');
% pard.addchannel1.position=[2,4.5];
% pard.addchannel1.Width=0.5;
% 
% pard.addchannel2.object=struct('String','','Style','checkbox');
% pard.addchannel2.position=[3,4.5];
% pard.addchannel2.Width=0.5;

pard.newfilenameb.object=struct('String','New filename:','Style','pushbutton','Callback',{{@newfilenameb_callback,obj}});
pard.newfilenameb.position=[7,1];

pard.newfilename.object=struct('String','newfile_sml','Style','edit');
pard.newfilename.position=[7,2];
pard.newfilename.Width=3;

pard.addchannel.object=struct('String','add Channel number to existing','Style','checkbox');
pard.addchannel.position=[6,3];
pard.addchannel.Width=2;

% pard.syncParameters={{'filelist_short','file1',{'String'}},{'filelist_short','file2',{'String'}}};

pard.plugininfo.type='ProcessorPlugin';
end