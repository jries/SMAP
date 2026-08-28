classdef AlternatingExcitationCorrect<interfaces.DialogProcessor
%     Renames localization fields
    methods
        function obj=AlternatingExcitationCorrect(varargin)     
            obj@interfaces.DialogProcessor(varargin{:});
        end
        
        function out=run(obj,p)
            out=[];
            obj.setPar('undoModule','AlternatingExcitationCorrect');
            notify(obj.P,'backup4undo');
            frames=obj.locData.loc.frame;
            ch=obj.locData.loc.channel;
            if p.removech
                inda=find(ch==p.cha);
                indb=find(ch==p.chb);
                fevena=mod(frames(inda),2)==0; fodda=mod(frames(inda),2)==1;
                fevenb=mod(frames(indb),2)==0; foddb=mod(frames(indb),2)==1;
                if sum(fevena)/length(inda)+sum(foddb)/length(indb)<sum(fodda)/length(inda)+sum(fevenb)/length(indb)
                    inx=inda; inda=indb; indb=indx;
                    fevena=mod(frames(inda),2)==0; fodda=mod(frames(inda),2)==1;
                    fevenb=mod(frames(indb),2)==0; foddb=mod(frames(indb),2)==1;
                end
                indremove=sort(vertcat(inda(fodda),indb(fevenb)));
                obj.locData.removelocs(indremove);
            end
            obj.locData.loc.frame=ceil(obj.locData.loc.frame/2);
            obj.locData.regroup;

            % p.layers

%             Remove wrong channel checkbox: Find for each channel: bright and dark frames. Remove
% Dual: frame -> ceil(frame/2).
% For consistency: start always with channel 1 and combine with next ch2?

   
        end
        function pard=guidef(obj)
            pard=guidef;
        end

        % function updateGui(obj,event,data)
        %     if ~isempty(obj.locData.loc)
        %     fn=fieldnames(obj.locData.loc);
        %     obj.guihandles.fieldselect.String=fn;
        %     end
        % end       
    end
end




function pard=guidef
pard.removech.object=struct('String','remove wrong ch:','Style','checkbox');
pard.removech.position=[1,1];
pard.removech.Width=1.5;

pard.channelsused.object=struct('String','Channels used:','Style','text');
pard.channelsused.position=[1,2.5];
pard.cha.object=struct('String','1','Style','edit');
pard.cha.position=[1,3.5];
pard.cha.Width=0.5;

pard.chb.object=struct('String','2','Style','edit');
pard.chb.position=[1,4];
pard.chb.Width=0.5;

pard.plugininfo.type='ProcessorPlugin';
pard.plugininfo.description='Alternating excitation to standard representation';
end