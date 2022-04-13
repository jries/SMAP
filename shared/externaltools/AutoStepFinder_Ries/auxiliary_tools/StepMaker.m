
function varargout = StepMaker(varargin);
% 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StepMaker_OpeningFcn, ...
                   'gui_OutputFcn',  @StepMaker_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
               
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


function Stepmaker(handles)   
%GUI parameters
initval.flatstep        =get(handles.flatstep,'Value');              %Flat distr steps
initval.minstep         =str2double(get(handles.minstep,'String'));
initval.maxstep         =str2double(get(handles.maxstep,'String'));
initval.gausstep        =get(handles.gausstep,'Value');              %Gaus distr steps
initval.meanstep        =str2double(get(handles.meanstep,'String'));
initval.sigmastep       =str2double(get(handles.sigmastep,'String'));
initval.expstep         =get(handles.expstep,'Value');               %Exp distr steps
initval.decaystep       =str2double(get(handles.decaystep,'String'));
initval.flatdwell       =get(handles.flatdwell,'Value');             %Flat dwell dwell
initval.mindwell        =str2double(get(handles.mindwell,'String'));
initval.gausdwell       =get(handles.gausdwell,'Value');             %Gausd distr dwell
initval.maxdwell        =str2double(get(handles.maxdwell,'String'));
initval.meandwell       =str2double(get(handles.meandwell,'String'));
initval.expdwell        =get(handles.expdwell,'Value');              %Flat distr dwell
initval.sigmadwell      =str2double(get(handles.sigmadwell,'String'));
initval.decaydwell      =str2double(get(handles.decaydwell,'String'));
initval.stepsnumber     =str2double(get(handles.stepsnumber,'String'));
initval.noisesteps      =str2double(get(handles.noisesteps,'String'));
initval.repeats         =str2double(get(handles.repeats,'String'));
initval.addbase         =str2double(get(handles.addbase,'String'));
initval.repeatsteps     =str2double(get(handles.traces,'String'));
stepmaker_mainloop(initval)
 
     function stepmaker_mainloop(initval)

%JWJK_B:-------------------------------------------------------------------
%Title:Simulate Step Traces
%Summary: %this function  creates a simple stepped trace with noise. 
%Approach: Steps  are taken as random samples from a 
%stepsize - and dwell time distribution. 
%This distribution can be refined in the 'Set_type_of_signal' function

%Dwell times can be set as 'dependent' on the step size or 'independent'.
%In the latter case they are taken as random samples from a flat or
%exponential dwell time distribution.

%By varying the settings, a wide range of traces types can be generated

%Below, example settings are given for case studies

%Input: none
%Output: step trace saved as text column in root path

%Jacob Kerssemakers, Cees Dekker Lab, 2018-20
%JWJK_B:-------------------------------------------------------------------

% if nargin<1
%     %definition of steps per type
%     initval.flatstep=1;       %type flat
%     initval.minstep=8;         
%     initval.maxstep=8;         
%     %---------------------------------
%     initval.gausstep=0;       %type gauss
%     initval.meanstep=10;      
%     initval.sigmastep=20;       
%     %---------------------------------
%     initval.expstep=0;         %type exponential
%     initval.decaystep=50;      
%     %---------------------------------
%     
%     %definition of dwells per type
%     initval.flatdwell=1;      %type Flat 
%     initval.mindwell=1; 
%     initval.maxdwell=100;       
%     %---------------------------------
%     initval.gausdwell=0;      %type Gauss 
%     initval.meandwell=100;      
%     initval.sigmadwell=100;      
%     %----------------------------------
%     initval.expdwell=0;      %type exponential
%     initval.decaydwell=100;      
%     %---------------------------------
% 
%     %trace properties
%     initval.stepsnumber=6;     %number of steps
%     initval.noisesteps=1;      %noise level   
%     initval.addbase=-700;        %add a L-points baseline before (-L) or after (+L)
%     initval.repeatsteps=3;    %repeat the number of 'blocks' per trace
%     initval.repeats=5;        %repeat the generation of traces
% end

%% distribution generation

        %Build dwell distribution    
        [dwell_distribution,binaxis_dwells,initval]=Set_dwell_distributions(initval,50)          
        summed_dwell_distribution=cumsum(dwell_distribution);

        %%Build step distribution
        [step_distribution,binaxis_steps,initval]=Set_step_distributions(initval,50);         
        summed_stepdistribution=cumsum(step_distribution); 

%% Build lists picking a random sample from the step and dwell distributions        
for rp=1:initval.repeats
    Concat_Curve=[];
    for bl=1:initval.repeatsteps 
        % Build step size list  
        StepSizes=zeros(initval.stepsnumber,1);  
        for jj=1:initval.stepsnumber+1  
                step_sample=Pick_from_distribution(binaxis_steps,summed_stepdistribution);
                StepSizes(jj)=step_sample;
        end 
       % Build dwell-time list 
       Dwelltimes=zeros(initval.stepsnumber,1); 
        for jj=1:initval.stepsnumber+1
            dwell_sample=Pick_from_distribution(binaxis_dwells,summed_dwell_distribution);
            DwellTimes(jj)=round(dwell_sample);
        end


        %% Build curve     
        Curve=[];
        Level=0;
        for jj=1:initval.stepsnumber+1
            Dwell=DwellTimes(jj);
            Step=StepSizes(jj);        
            NewSection=Level+Step*ones(Dwell,1);        
            Curve=[Curve ; NewSection];
            Level=Curve(end);
        end  

        % Add start or tail
        if initval.addbase<0,Curve=[zeros(abs(initval.addbase),1)+Curve(1) ; Curve];end
        if initval.addbase>0,Curve=[Curve ; zeros(abs(initval.addbase),1)+Curve(end)];end       

        Concat_Curve=[Concat_Curve; Curve];
     end



    %% Finalize:
    %Add time axis, noise
    Ld=length(Concat_Curve);
    TimAx=(1:Ld)';
    Noise=initval.noisesteps*randn(Ld,1);
    Trace=Concat_Curve+Noise; 
    data=[TimAx Trace];







    %% plot
    close(findobj('type','figure','name','Distributions and Trace'));        %close S-curve plots --> for batch mode
    figure('Name','Distributions and Trace','NumberTitle','off','units', 'normalized', 'position', [0.05 0.35 0.5 0.5]);
    %figure(169);
    subplot(2,2,1);
        bar(binaxis_steps,step_distribution, 'w','Linewidth',1);
        title('Step Sizes');
        xlabel('step size, a.u.');
        ylabel('occurence, a.u.');
        axis tight;
    subplot(2,2,2); 
        bar(binaxis_dwells,dwell_distribution, 'w','Linewidth',1);
        title('Dwell Times');
        xlabel('dwell time, pts');
        ylabel('occurence, a.u.');
        axis tight;
    subplot(2,1,2);
        plot(data(:,1),data(:,2), 'b-','Linewidth',1);
        title('Trace');
    pause(0.1);



    %% saveu
    save_it=1;   
    writepath=[pwd, '\test_traces'];
    if ~isdir(writepath), mkdir(writepath);,end
    SaveName=strcat('testdata',...
        '_step_distribution',initval.TypeOfSteps,...
        '_dwell_distribution',initval.TypeOfSteps,...
        '_steps',num2str(initval.stepsnumber),...
        '_blocks',num2str(bl),...
        '_repeat',num2str(bl),num2str(rp),...
        '.txt');
    if save_it  
        dlmwrite(strcat(writepath,'/',SaveName), data);
        %saveas(gcf,strcat(writepath,'/',SaveName(1:end-4),'.jpg'));
    end
end

 function Sample=Pick_from_distribution(BinAxis,SumCurve)
            %pick from distribution
            pickval=rand(1);  %pick a random percentile
            sel_st=find(SumCurve>0);
            NonZeroSumCurve=SumCurve(sel_st);
            NonZeroBinAxis=BinAxis(sel_st);                        
            %nearest occupied bin            
            [~,idx]=min(abs(NonZeroSumCurve-pickval));
            %will return first index of plateau which is correct
            Sample=NonZeroBinAxis(idx);  


   function [StepSizeCurve,binaxis_steps,initval]=Set_step_distributions(initval,binz)
    %set up global range of stepsizes
    if initval.flatstep==1,initval.TypeOfSteps='Flat';end
    if initval.gausstep==1,initval.TypeOfSteps='Gaussian';end
    if initval.expstep==1,initval.TypeOfSteps='Exponential';end  
    switch initval.TypeOfSteps
       case 'Flat'
           %each step in the specified range has an equal chance of
           %occurring
           binaxis_steps=unique(linspace(initval.minstep,initval.maxstep,binz));
           StepSizeCurve=(1+0*(binaxis_steps));         
        case 'Gaussian'
            minstep=initval.meanstep-2*initval.sigmastep;
            maxstep=initval.meanstep+2*initval.sigmastep;
            binaxis_steps=unique(linspace(minstep,maxstep,binz));
            StepSizeCurve=exp(-( (binaxis_steps-initval.meanstep)/(initval.sigmastep)).^2);
        case 'Exponential'  
            binaxis_steps=unique(linspace(0,2*initval.decaystep,binz));
            StepSizeCurve=exp(-binaxis_steps/(initval.decaystep)); 
    end
    StepSizeCurve=StepSizeCurve/sum(StepSizeCurve);  %normalized dstribution

    function [DwellTimeCurve,binaxis_dwells,initval]=Set_dwell_distributions(initval,binz)
    %set up global range of stepsizes
    if initval.flatdwell==1,initval.TypeOfDwells='Flat';end
    if initval.gausdwell==1,initval.TypeOfDwells='Gaussian';end
    if initval.expdwell==1, initval.TypeOfDwells='Exponential';end

    switch initval.TypeOfDwells
       case 'Flat'
           %each step in the specified range has an equal chance of
           %occurring
           mindwell=max([1 initval.mindwell]);
           maxdwell=max([mindwell initval.maxdwell]);
           binaxis_dwells=unique(linspace(mindwell,maxdwell,binz));
           DwellTimeCurve=(1+0*(binaxis_dwells));         
        case 'Gaussian'  %should be>0
            mindwell=max([1 initval.meandwell-2*initval.sigmadwell]);;
            maxdwell=max([mindwell initval.meandwell+2*initval.sigmadwell]);
            binaxis_dwells=unique(linspace(mindwell,maxdwell,binz));
            DwellTimeCurve=exp(-( (binaxis_dwells-initval.meandwell)/(initval.sigmadwell)).^2);
        case 'Exponential'  
            decay=max([1 initval.decaydwell]);
            binaxis_dwells=unique(linspace(1,2*decay,binz));
            DwellTimeCurve=exp(-binaxis_dwells/(decay)); 
    end
    DwellTimeCurve=DwellTimeCurve/sum(DwellTimeCurve);  %normalized dstribution



    
% --- Executes just before StepMaker is made visible.
function StepMaker_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StepMaker (see VARARGIN)

% Choose default command line output for StepMaker
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);





% --- Outputs from this function are returned to the command line.
function varargout = StepMaker_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function stepsnumber_Callback(hObject, ~, ~)
checkmax_stepsnumber=get(hObject,'String');
checkNan_stepsnumber=isnan(str2double(checkmax_stepsnumber));
checkmax_stepsnumber=str2double(checkmax_stepsnumber);

     if checkNan_stepsnumber==1
         msgbox('The number of steps is NaN.','ERROR', 'error')
         set(hObject,'String',20);
     return;
     end
if checkmax_stepsnumber < 1
         msgbox('The number of steps is smaller than 1. The input value has been set to 1','ERROR', 'error')
         set(hObject,'String',1);
     return;     
end 
 

function noisesteps_Callback(hObject, ~, ~)
checkmax_noisesteps=get(hObject,'String');
checkmax_noisesteps=isnan(str2double(checkmax_noisesteps));
     if checkmax_noisesteps==1
         msgbox('The noise setting is NaN.','ERROR', 'error')
         set(hObject,'String',3);
     return;
     end

   
function minstep_Callback(hObject, ~, ~)
checkmax_minstep=get(hObject,'String');
checkmax_minstep=isnan(str2double(checkmax_minstep));
     if checkmax_minstep==1
         msgbox('The Min stepsize setting is NaN.','ERROR', 'error')
         set(hObject,'String',-10);
     return;
     end


function GenerateData_Callback(~, handles, ~)
% hObject    handle to GenerateData (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(StepMaker);
Stepmaker(handles)

function decaystep_Callback(hObject, ~, ~)
checkmax_decaystep=get(hObject,'String');
checkmax_decaystep=isnan(str2double(checkmax_decaystep));
     if checkmax_decaystep==1
         msgbox('The decay setting is NaN.','ERROR', 'error')
         set(hObject,'String',100);
     return;
     end



function mindwell_Callback(hObject, ~, ~)
checkmax_mindwell=get(hObject,'String');
checkmax_mindwell=isnan(str2double(checkmax_mindwell));
     if checkmax_mindwell==1
         msgbox('The Min dwell time setting is NaN.','ERROR', 'error')
         set(hObject,'String',50);
     return;
     end
   



function meandwell_Callback(hObject, ~, ~)
checkmax_meandwell=get(hObject,'String');
checkmax_meandwell=isnan(str2double(checkmax_meandwell));
     if checkmax_meandwell==1
         msgbox('The mean dwell time setting is NaN.','ERROR', 'error')
         set(hObject,'String',75);
     return;
     end
if checkmax_meandwell < 2
         msgbox('The number of steps is smaller than 1. The input value has been set to 1','ERROR', 'error')
         set(hObject,'String',2);
     return;     
end 



function decaydwell_Callback(hObject, ~, ~)
checkmax_decaydwell=get(hObject,'String');
checkmax_decaydwell=isnan(str2double(checkmax_decaydwell));
     if checkmax_decaydwell==1
         msgbox('The decay setting is NaN.','ERROR', 'error')
         set(hObject,'String',100);
     return;
     end




function maxdwell_Callback(hObject, ~, ~)
checkmax_maxdwell=get(hObject,'String');
checkmax_maxdwell=isnan(str2double(checkmax_maxdwell));
     if checkmax_maxdwell==1
         msgbox('The max dwell time setting is NaN.','ERROR', 'error')
         set(hObject,'String',100);
     return;
     end
     if checkmax_maxdwell < 2
         msgbox('The max dwelltime is smaller than 2. The input value has been set to 2','ERROR', 'error')
         set(hObject,'String',2);
     return;     
end 





function sigmadwell_Callback(hObject, ~, ~)
checkmax_sigmadwell=get(hObject,'String');
checkmax_sigmadwell=isnan(str2double(checkmax_sigmadwell));
     if checkmax_sigmadwell==1
         msgbox('The sigma dwell time setting is NaN.','ERROR', 'error')
         set(hObject,'String',25);
     return;
     end



function addbase_Callback(hObject, ~, ~)
checkmax_addbase=get(hObject,'String');
checkmax_addbase=isnan(str2double(checkmax_addbase));
     if checkmax_addbase==1
         msgbox('The add baseline setting is NaN.','ERROR', 'error')
         set(hObject,'String',0);
     return;
     end

function meanstep_Callback(hObject, ~, ~)
checkmax_meanstep=get(hObject,'String');
checkmax_meanstep=isnan(str2double(checkmax_meanstep));
     if checkmax_meanstep==1
         msgbox('The mean stepsize setting is NaN.','ERROR', 'error')
         set(hObject,'String',10);
     return;
     end

function sigmastep_Callback(hObject, ~, ~)
checkmax_sigmastep=get(hObject,'String');
checkmax_sigmastep=isnan(str2double(checkmax_sigmastep));
     if checkmax_sigmastep==1
         msgbox('The sigma stepsize setting is NaN.','ERROR', 'error')
         set(hObject,'String',10);
     return;
     end



function repeats_Callback(hObject, ~, ~)
checkmax_repeats=get(hObject,'String');
checkmax_repeats=isnan(str2double(checkmax_repeats));
     if checkmax_repeats==1
         msgbox('The number of repeats is NaN.','ERROR', 'error')
         set(hObject,'String',1);
     return;
     end
checkmax_repeats1=str2double(get(hObject,'String'));     
     if checkmax_repeats1 < 1
         msgbox('The number of repeats is smaller than 1. The input value has been set to 1','ERROR', 'error')
         set(hObject,'String',1);
     return;     
end 
 


function traces_Callback(hObject, ~, ~)
checkmax_traces=get(hObject,'String');
checkmax_traces2=isnan(str2double(checkmax_traces));
     if checkmax_traces2==1
         msgbox('The number of traces is NaN.','ERROR', 'error')
         set(hObject,'String',1);
     return;
     end
checkmax_traces1=str2double(get(hObject,'String'));
if checkmax_traces1 < 1
         msgbox('The number of traces is smaller than 1. The input value has been set to 1','ERROR', 'error')
         set(hObject,'String',1);
     return;     
end 
 

function maxstep_Callback(hObject, ~, ~)
checkmax_maxstep=get(hObject,'String');
checkmax_maxstep=isnan(str2double(checkmax_maxstep));
     if checkmax_maxstep==1
         msgbox('The max stepsize setting is NaN.','ERROR', 'error')
         set(hObject,'String',1);
     return;
     end


% --- Executes on button press in flatstep.
function flatstep_Callback(~, ~, handles)
% hObject    handle to flatstep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.maxstep, 'Enable','On');
set(handles.minstep, 'Enable','On');
set(handles.meanstep, 'Enable','Off');
set(handles.sigmastep, 'Enable','Off');
set(handles.decaystep, 'Enable','Off');



% --- Executes on button press in gausstep.
function gausstep_Callback(~, ~, handles)
set(handles.maxstep, 'Enable','Off');
set(handles.minstep, 'Enable','Off');
set(handles.meanstep, 'Enable','On');
set(handles.sigmastep, 'Enable','On');
set(handles.decaystep, 'Enable','Off');



% --- Executes on button press in expstep.
function expstep_Callback(~, ~, handles)
set(handles.maxstep, 'Enable','Off');
set(handles.minstep, 'Enable','Off');
set(handles.meanstep, 'Enable','Off');
set(handles.sigmastep, 'Enable','Off');
set(handles.decaystep, 'Enable','On');


% --- Executes on button press in flatdwell.
function flatdwell_Callback(~, ~, handles)
set(handles.maxdwell, 'Enable','On');
set(handles.mindwell, 'Enable','On');
set(handles.meandwell, 'Enable','Off');
set(handles.sigmadwell, 'Enable','Off');
set(handles.decaydwell, 'Enable','Off');


% --- Executes on button press in gausdwell.
function gausdwell_Callback(~, ~, handles)
set(handles.maxdwell, 'Enable','Off');
set(handles.mindwell, 'Enable','Off');
set(handles.meandwell, 'Enable','On');
set(handles.sigmadwell, 'Enable','On');
set(handles.decaydwell, 'Enable','Off');


% --- Executes on button press in expdwell.
function expdwell_Callback(~, ~, handles)
set(handles.maxdwell, 'Enable','Off');
set(handles.mindwell, 'Enable','Off');
set(handles.meandwell, 'Enable','Off');
set(handles.sigmadwell, 'Enable','Off');
set(handles.decaydwell, 'Enable','On');
