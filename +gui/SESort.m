classdef SESort< interfaces.SEProcessor
    properties
        menue
    end
    methods
        function obj=SESort(varargin)
            obj@interfaces.SEProcessor;
            if nargin>0
                obj.handle=varargin{1};  
            end        
        end
        function pard=guidef(obj)
            pard=guidef;
        end
        function initGui(obj)
            initGui@interfaces.SEProcessor(obj);
            men=menutext;
            obj.menue=men;
            for k=1:4
                obj.guihandles.(['sortprop' num2str(k)]).String={men(:).name};
                obj.guihandles.(['sortspec' num2str(k)]).Visible='off';
                obj.guihandles.(['sortbutton' num2str(k)]).Visible='off';
%                 obj.guihandles.(['sortedit' num2str(k)]).Visible='off';
                
                
                obj.guihandles.(['sortprop' num2str(k)]).Callback={@prop_callback,obj,k};
                obj.guihandles.(['sortspec' num2str(k)]).Callback={@spec_callback,obj,k};
                obj.guihandles.(['sortbutton' num2str(k)]).Callback={@button_callback,obj,k};
                obj.guihandles.(['sortedit' num2str(k)]).Callback={@edit_callback,obj,k};
                obj.guihandles.sort.Callback={@sort_callback,obj};
            end
%             set(obj.guihandles.showSE,'Callback',{@make_siteexplorer,obj})
%             set(obj.guihandles.redrawall,'Callback',{@redrawall_callback,obj})
%             set(obj.guihandles.clearall,'Callback',{@clearall_callback,obj})
%             addlistener(obj.SE.locData,'loaded',@obj.loaded_notify);
        end

%         function loaded_notify(obj,lb,eventdata)
%             obj.updateParameters;
%         end
    end
end

function prop_callback(object,action,obj,num)
nums=num2str(num);
specstr=obj.menue(object.Value).spec;
if ~isempty(specstr)
    obj.guihandles.(['sortspec' nums]).String=specstr;
    obj.guihandles.(['sortspec' nums]).Visible='on';
    obj.guihandles.(['sortspec' nums]).Value=1;
    spec_callback(obj.guihandles.(['sortspec' nums]),0,obj,num)
else
    obj.guihandles.(['sortspec' nums]).Visible='off';
    obj.guihandles.(['sortedit' nums]).String=obj.menue(object.Value).field;
end

    
end

function spec_callback(object,data,obj,num)
nums=num2str(num);
specstr=object.String{object.Value};
obj.guihandles.(['sortbutton' nums]).Visible='off';
switch specstr
    case 'Edit'
        obj.guihandles.(['sortbutton' nums]).Visible='on';
        obj.guihandles.(['sortedit' nums]).Visible='on';
    case 'File'
        obj.guihandles.(['sortedit' nums]).String='info.filenumber';
    case 'Cell'
        obj.guihandles.(['sortedit' nums]).String='info.cell';
    case 'Site'
        obj.guihandles.(['sortedit' nums]).String='ID';
    otherwise
         pr=obj.guihandles.(['sortprop' nums]);
        field=obj.menue(pr.Value).field;

        str2=[field '.' (specstr)];% '.value'];
        obj.guihandles.(['sortedit' nums]).String=str2;
%         obj.guihandles.(['sortedit' nums]).Visible='off';        
end
end

function button_callback(a,b,obj,num)
nums=num2str(num);
spec=obj.guihandles.(['sortprop' nums]);
selection=spec.String{spec.Value};
site=obj.SE.currentsite;
siteh=site;
switch selection
    case 'Annotation'
        field='annotation';
    case 'Evaluation'
        field='evaluation';
    case 'Other'
        warning('off','MATLAB:structOnObject')
        
        s.x=struct(site);
        warning('on','MATLAB:structOnObject')
        siteh=s;
        field='x';
    otherwise
        field='';
end

 if isempty(field)||isstruct(siteh.(field))
    str=browsefields(siteh,field);

    if ~isempty(str)
        if strcmpi(selection,'Other')
            str=str(3:end);
        end
%         obj.guihandles.(['sortedit' nums]).String=str;
    end
 else
     str='field';
 end
 if ~isempty(str)
    val=eval(['site.' str]);
    if isa(val,'double')&&length(val)>1
        numind=length(val(:));
        for k=1:numind
        numstr{k}=num2str(k);
        end
        answ=listdlg('ListString',numstr,'Name','Index','PromptString','select index','SelectionMode','single');
        if ~isempty(answ)
            str=[str '(' num2str(answ) ')'];
        end     
    end
    obj.guihandles.(['sortedit' nums]).String=str;
 end
end

function str=browsefields(prop,field)
if isempty(field)
    fn=fieldnames(prop);
    
elseif isstruct(prop.(field))
    fn=fieldnames(prop.(field));
else
    str=field;
    return;
end
    answ=listdlg('ListString',fn,'SelectionMode','single');
    if isempty(answ)
        str='';
    else
    field2=fn{answ};

    str=[field '.' browsefields(prop.(field),field2)];
    end

% if isstruc(prop)
% fn=fieldnames(prop)
% selectfield(fieldnames(obj.
end

function edit_callback(a,b,obj,num)
end

function sort_callback(a,b,obj)
par=obj.getAllParameters;
sites=obj.SE.sites;
sortmatrix=zeros(length(sites),4);
for k=1:4
    if isfield(par,['sortedit' num2str(k)])
        field{k}=par.(['sortedit' num2str(k)]);
    end
end
for k=1:length(sites)
    for s=1:length(field)
%         if isfield(sites(k),field{s})
%          try
         if ~isempty(field{s})
            evalstring=['sites(k).' field{s}];
            val=eval(evalstring);
            if isstruct(val)
                val=val.value;
            end
            sortmatrix(k,s)=val(1);
         end
%          end
    end
end

for s=1:4
    if strcmpi(par.(['direction' num2str(s)]).selection,'descend') || strcmpi(par.(['direction' num2str(s)]).selection,'decend')
        sortmatrix(:,s)=-sortmatrix(:,s);
    end
end
[~,sortind]=sortrows(sortmatrix);
obj.SE.sites=obj.SE.sites(sortind);
obj.SE.setIndList;
obj.SE.processors.preview.updateSitelist;

end

function pard=guidef


pard.title1.object=struct('String','1st sort','Style','text');
pard.title1.position=[1,1];

pard.sortprop1.object=struct('String','sortpar','Style','popupmenu');
pard.sortprop1.position=[1,2];

pard.sortspec1.object=struct('String','sortspec','Style','popupmenu');
pard.sortspec1.position=[1,3];

pard.sortbutton1.object=struct('String','select','Style','pushbutton');
pard.sortbutton1.position=[1,4];

pard.sortedit1.object=struct('String','','Style','edit');
pard.sortedit1.position=[2,2];
pard.sortedit1.Width=3;
pard.direction1.object=struct('String','ascend|descend','Style','popupmenu');
pard.direction1.position=[2,1.05];
pard.direction1.Width=.95;

pard.title2.object=struct('String','2nd sort','Style','text');
pard.title2.position=[3,1];

pard.sortprop2.object=struct('String','sortpar','Style','popupmenu');
pard.sortprop2.position=[3,2];

pard.sortspec2.object=struct('String','sortspec','Style','popupmenu');
pard.sortspec2.position=[3,3];

pard.sortbutton2.object=struct('String','select','Style','pushbutton');
pard.sortbutton2.position=[3,4];

pard.sortedit2.object=struct('String','','Style','edit');
pard.sortedit2.position=[4,2];
pard.sortedit2.Width=3;
pard.direction2.object=struct('String','ascend|descend','Style','popupmenu');
pard.direction2.position=[4,1.05];
pard.direction2.Width=.95;

pard.title3.object=struct('String','3rd sort','Style','text');
pard.title3.position=[5,1];

pard.sortprop3.object=struct('String','sortpar','Style','popupmenu');
pard.sortprop3.position=[5,2];

pard.sortspec3.object=struct('String','sortspec','Style','popupmenu');
pard.sortspec3.position=[5,3];

pard.sortbutton3.object=struct('String','select','Style','pushbutton');
pard.sortbutton3.position=[5,4];

pard.sortedit3.object=struct('String','','Style','edit');
pard.sortedit3.position=[6,2];
pard.sortedit3.Width=3;

pard.direction3.object=struct('String','ascend|descend','Style','popupmenu');
pard.direction3.position=[6,1.05];
pard.direction3.Width=.95;

pard.title4.object=struct('String','4th sort','Style','text');
pard.title4.position=[7,1];

pard.sortprop4.object=struct('String','sortpar','Style','popupmenu');
pard.sortprop4.position=[7,2];

pard.sortspec4.object=struct('String','sortspec','Style','popupmenu');
pard.sortspec4.position=[7,3];

pard.sortbutton4.object=struct('String','select','Style','pushbutton');
pard.sortbutton4.position=[7,4];

pard.sortedit4.object=struct('String','','Style','edit');
pard.sortedit4.position=[8,2];
pard.sortedit4.Width=3;

pard.direction4.object=struct('String','ascend|descend','Style','popupmenu');
pard.direction4.position=[8,1.05];
pard.direction4.Width=.95;

pard.sort.object=struct('String','Sort','Style','pushbutton');
pard.sort.position=[10,3];
pard.sort.Height=2;


end


function men=menutext
men(1).name='None';
men(1).field='';
men(1).spec={};
men(1).button='';
men(1).edit='';

nmen=8;

for k=2:nmen
    men(k)=men(1);
end
men(2).name='Hierarchy';
% men(2).field='info.filenumber';
men(2).spec={'File','Cell','Site'};
men(3).name='Statistics';
men(3).field='evaluation.generalStatistics';
men(3).spec={'Nphot','locplayers','PSFlayers','Edit'};
% men(3).field='info.cell';
% men(4).name='Site';
% men(4).field='ID';
 men(4).name='List';
men(4).field='annotation';
men(4).spec={'list1','list2','list3','list4'};

men(5).name='Annotation';
men(5).field='annotation';
men(5).spec={'Line1','Line2','Edit'};

men(6).name='Evaluation';
men(6).spec={'Edit'};

men(7).name='Other';
men(7).spec={'Edit'};
end