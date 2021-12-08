function varargout = StepMerger(varargin)
% STEPMERGER MATLAB code for StepMerger.fig
%      STEPMERGER, by itself, creates a new STEPMERGER or raises the existing
%      singleton*.
%
%      H = STEPMERGER returns the handle to a new STEPMERGER or the handle to
%      the existing singleton*.
%
%      STEPMERGER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in STEPMERGER.M with the given input arguments.
%
%      STEPMERGER('Property','Value',...) creates a new STEPMERGER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before StepMerger_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to StepMerger_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help StepMerger

% Last Modified by GUIDE v2.5 10-Feb-2021 18:31:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @StepMerger_OpeningFcn, ...
                   'gui_OutputFcn',  @StepMerger_OutputFcn, ...
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

function Step_merger_main(handles)
    
%% settings for despiking
    init.path               =get(handles.directory, 'String');
    init.despike            =get(handles.despikeon,'Value');            %Despiking on
    init.spikemaxwidth      =str2double(get(handles.width,'String'));  
    init.updownmargin       =str2double(get(handles.margin,'String'));  %fraction that steps can be different
    init.spikeup            =get(handles.dirup,'Value');
    init.spikedown          =get(handles.dirdown,'Value');
    init.spikeboth          =get(handles.dirboth,'Value');
%% settings for slope merging
    init.slopemerge         =get(handles.mergeon,'Value');                  %Merging on
    init.wmin               =str2double(get(handles.widthmerge,'String'));  %Max width
%% Bootstrap settings
    init.booton             =get(handles.erroreston,'Value'); %bootstrapping on
    Step_merger_main_mainloop(init)
    
function Step_merger_main_mainloop(init)
%This function combines options for merging small step trains, or
init.codepath=pwd;

if init.spikeup,  init.spikesign=1;else
if init.spikedown,   init.spikesign=-1;else
if init.spikeboth, init.spikesign=0;end; end; end

if init.booton
    init.bootstraprepeats=1000; 
else
    init.bootstraprepeats=0;
end



%run_modus='demo'; 'quickload'; 
%run_modus='quickload';
run_modus='default';
switch run_modus
    case 'demo' %DEMO
        DeSpiker;
        SlopeMerger;
         FileName='test_fits.txt';
         PathName=[pwd, '/'];
    case 'default'
        %% load
        cd(init.path);
        [FileName,PathName] = uigetfile('*_fits.txt', 'Open the fit of interest');
        data=double(dlmread(strcat(PathName,'/',FileName),',',2,0));
        cd(init.codepath);
        TT=data(:,1);  
        XX=data(:,2);
        FitX=data(:,3);
        oriFitX=FitX;
        %% run
        if init.despike, 
            [FitX,step_props]=DeSpiker(init.spikemaxwidth, init.updownmargin, init.spikesign,XX,FitX);
        end  
        if init.slopemerge, 
            [FitX,step_props]=SlopeMerger(init.wmin,XX,FitX);
        end  
end   

%% Re-process errors
%needs adding of residual error and bootstrap here 
step_props=AddStep_Errors(XX,step_props);

 if init.bootstraprepeats>0
    step_idxes=step_props(:,1);
    [err_st, err_t, ~,~ ]=bootstrap_get_errors(XX, step_idxes,init.bootstraprepeats);
    step_props(:,9)=err_st;
    step_props(:,10)=err_t;
end
 


%% save
if ~strcmp(run_modus, 'demo')
    if init.despike==1,
        SaveName=[FileName(1:end-9), '_despiked'];
    end
    if init.slopemerge==1,
        SaveName=[FileName(1:end-9), '_merged'];
    end
        
    Save_AS_formats(TT,XX,FitX,step_props,SaveName, PathName);
end
%% plot
close(findobj('type','figure','name','Merged Trace'));        %close S-curve plots --> for batch mode
figure('Name','Merged Trace','NumberTitle','off','units', 'normalized', 'position', [0.02 0.5 0.6 0.3]);

plot(XX,'b-','LineWidth', 2); hold on;
plot(oriFitX, 'r-','LineWidth', 2);
plot(FitX, 'y-', 'LineWidth', 2);
legend('data','original fit','processed fit');
 

function [FitX_merged,StepFitProps_out]=DeSpiker(spikemaxwidth, updownmargin, spikesign,XX,FitX)
%% Remove spikes
%Aim: spikes (short up-and-down) or blinks (reverse) may be of less
%interest to a user. This tool finds and removes these. 

%Description: 
%1) DeSpiker: no input, runs a demo trace on trains
%4) DeSpiker(spikemaxwidth, updownmargin, spikesign): load data and process according to modus.
    %with this option, output is saved with extension '_merged'
%3) DeSpiker(spikemaxwidth, updownmargin, spikesign,XX,FitX)

%Input:
    %XX is a single data column, FitX is its fit.
    %spikemaxwidth(5), : max width for spike to consider
    %updownmargin(0.30): relative magnitude margin per up and down step 
    %spikesign(0): 
%         1: remove positive spikes, 
%        -1: remove negative spikes or 'blinks'
%         0 remove all signs


%Output: new fit and step properties

%Extras: includes autorun demo option of a  noisy step trace

%Reference: Auxiliairy tool of AutoStepfinder: 
%A fast and automated step detection method for single-molecule analysis.
%Luuk Loeff*, Jacob Kerssemakers*, Chirlmin Joo & Cees Dekker.
% * Equal contribution
% Last update by Jacob Kerssemakers, October 2020
%--------------------------------------------------------------------------

%% Get trace
if (nargin==3)  %Load standard Autostepfinder output
    [FileName,PathName] = uigetfile('*_fits.txt', 'Open the fit of interest');
    data=double(dlmread(strcat(PathName,'/',FileName),',',2,1));
    XX=data(:,1);
    FitX=data(:,2);
end  
 %% Demo: make a simple trace with a number of  spikes
if nargin==0
    spikemaxwidth=7;
    updownmargin=0.4, 
    spikesign=0;
    [XX,FitX]=MakeDemoSpikeTrace;   
end
% ---------------------------------------------
 
 Lx=length(XX);
 TT=(1:Lx)';
[step_props,~,~]=Get_StepsFromFit(TT,FitX);              

%% remove spikes; 
%run twice to clean ''comb sections''
 FitX_merged=FitX;
 for i=1:2
    [~,isolated,FitX_merged]=remove_spikes(step_props(:,1), step_props(:,5), spikemaxwidth, updownmargin, spikesign,FitX_merged);
    [step_props,~,~]=Get_StepsFromFit(TT,FitX_merged);
 end
[StepFitProps_out,~,~]=Get_StepsFromFit(TT,FitX_merged);

if nargin <5
    figure(67);
    plot(XX,'k-','LineWidth', 2); hold on;
    plot(FitX, 'c-','LineWidth', 1.5);
      plot(FitX_merged, 'r-', 'LineWidth', 2);
     legend('data','original fit','final fit');
     title('DeSpiker');
     dum=1;
end

if nargin==3
     Axz=(1:length(XX))';
     data=[Axz XX FitX];
     outname1=strcat(PathName,'/',FileName(1:end-4),'_despiked.txt');
     dlmwrite(outname1,data);
     outname2=strcat(PathName,'/',FileName(1:end-4),'_despiked_props.txt');
     dlmwrite(outname2,StepFitProps_out);
end

 function  [spikes,isolated,FitX]=remove_spikes(index, steps, spikemaxwidth, updownmargin, spikesign,FitX);
%find all spikes; output is sets of indices:
%'spikes'  pairs of indices
%'isolated'

%indices and steps
OLi=index(1:end-3);     OLs=steps(1:end-3);
CLi=index(2:end-2);     CLs=steps(2:end-2);
CRi=index(3:end-1);     CRs=steps(3:end-1);
ORi=index(4:end);       ORs=steps(4:end);
%differences
dt_L=CLi-OLi;           ds_L=CLs-OLs;   %outer left, center left etc
dt_C=CRi-CLi;           ds_C=CLs-CRs;
dt_R=ORi-CRi;           ds_L=ORs-CRs; 

narrow=(dt_C<=spikemaxwidth);
opposed=(sign(CLs)+sign(CRs)==0);
nullingratio=abs((CRs+CLs))./(0.5*(abs(CRs)+abs(CLs)));
nulling=(nullingratio<updownmargin);
if spikesign~=0
    propersign=(sign(CLs)==spikesign);
else
    propersign=1+0*sign(CLs);
end

sel_L=find((narrow&opposed&nulling&propersign)==1)+1;   
spike_index_L=index(sel_L);
spike_index_R=index(sel_L+1);
spikes=sort(unique([spike_index_L spike_index_R]));
isolated=index(~ismember(index,spikes));

%remove the spikes
for spi=1:length(spike_index_L)    
    lft_i=spike_index_L(spi);
    rgt_i=spike_index_R(spi);
    padval=FitX(lft_i-1);
    FitX(lft_i:rgt_i)=padval;
end

function [XX,FitX]=MakeDemoSpikeTrace
%Demo Trace
Lxi=2000;
noize=10;                   %noise level
XX=[]; FitX=[];
stepsizes=[-10 10 30 50 5];
Nspikes=30; 
spike_ampli=100;
isolated=[1 [10:50:Lxi]];  %just a step train
indices=isolated;
Li=length(indices);
for ii=1:Li
    dice=ceil(length(stepsizes)*(rand(1,1)));
    levs(ii)=stepsizes(dice);  %add random step levels
end
%build initial trace
for ii=1:Li-1
    L_sect=indices(ii+1)-indices(ii)+1;
    XX=[XX ; noize*randn(L_sect, 1)+levs(ii)];
    FitX=[FitX ; zeros(L_sect, 1)+levs(ii)];
end 
%add spikes

ups=round(linspace(105,Lxi-105,Nspikes));
downs=ups+2;
for spi=1:length(ups)
    peak=spike_ampli; %*randn(1,1);
    peaksign=2*round(rand(1,1))-1;
    XX(ups(spi):downs(spi))=XX(ups(spi):downs(spi))+peaksign*peak;
    FitX(ups(spi):downs(spi))=FitX(ups(spi):downs(spi))+peaksign*peak;
end

         
function [StepsX,levelX, histX]=Get_StepsFromFit(T,FitX)
%This function builds tables of steps or levels properties from a step fit
%Values are based on Averages of plateaus
lx=length(FitX);
difX=FitX(2:lx)-FitX(1:lx-1);
sel=find(difX~=0);  %note: index points to last point before step
lsel=length(sel);

%dwell time after
dwellT=T(sel(2:lsel))-T(sel(1:lsel-1)); 
dwellTafter=[dwellT' T(lx)-T(sel(lsel))]';
dwellTbefore=[T(sel(1))-T(1) dwellT']'; 
StepsX=[sel T(sel) FitX(sel) FitX(sel+1) difX(sel) dwellTbefore dwellTafter]; 
%list of stepsizes: [index time levelbefore levelafter step dwelltime before dwelltimeafter]

Levstartidx=[1; sel+1];  %first index of level
Levstopidx=[sel; lx];    %last index of level
LevstartT=T(Levstartidx);
LevstopT=T(Levstopidx);
LevLevel=[FitX(sel); FitX(lx)];
LevDwell=Levstopidx-Levstartidx+1;
LevStepBefore=[0; difX(sel)];
LevStepAfter=[difX(sel); 0];
levelX=[Levstartidx Levstopidx LevstartT LevstopT LevLevel LevDwell LevStepBefore LevStepAfter];
%list of levels: [startindex starttime stopindex stoptime level dwell stepbefore stepafter]   

%histogram
stepsizes=StepsX(:,4)-StepsX(:,3);
mx=max(stepsizes); mn=min(stepsizes);
stp=(mx-mn)/20;
hx=(mn:stp:mx);
if ~isempty(hx)
    histX=hist(stepsizes,hx)';
    histX=[hx' histX];
else
    histX=[1 1];
end
        
    
function FitX=Get_FitFromSteps(X,indexlist,modus)       
% This function builds plateau data
%list of levels: [startindex  stopindex starttime stoptime level dwell stepbefore stepafter]
lx=length(X);
lsel=length(indexlist); %note: index points to last point before step
%Build a 'FitX' based on the median levels (ot on the averages)
idxes=[0 ; indexlist ; lx];
FitX=0*X;
for ii=1:lsel+1
    ixlo=idxes(ii)+1;  %first index of plateau
    ixhi=idxes(ii+1);  %last index
    switch modus
        case 'mean', FitX(ixlo:ixhi)=nanmean(X(ixlo:ixhi));
        case 'median', FitX(ixlo:ixhi)=nanmedian(X(ixlo:ixhi));
    end
end 
        
   function [FitX,step_props]=SlopeMerger(wmin,XX,FitX)
%% Remove slope series by merging
%Aim 1: The Autostepfinder tends to fit trains of small steps to smooth
%slopes of non-instant steps, for example if the data is low-pass filtered.
%This tools finds such step trains and merges
%them by retaining the middle step index. After that, fit is re-done

%Description: 
%1) SlopeMerger: no input, runs a demo trace
%2) SlopeMerger(wmin): load data and process.
    %with this option, output is saved with extension '_merged'
%3) SlopeMerger(wmin,XX,FitX)  

%Output:
 %XX is a single data column, FitX is its fit. 
  %wmin(3) is the maximum distance aqual-signed steps should be spaced to 
  %be considered for merging   
    
%Output: new fit and step properties

%Extras: includes autorun demo option of a  noisy step trace

%Reference: Auxiliairy tool of AutoStepfinder: 
%A fast and automated step detection method for single-molecule analysis.
%Luuk Loeff*, Jacob Kerssemakers*, Chirlmin Joo & Cees Dekker.
% * Equal contribution
% Last update by Jacob Kerssemakers, October 2020
%--------------------------------------------------------------------------

%% Get trace
if (nargin==1)  %Load standard Autostepfinder output
    [FileName,PathName] = uigetfile('*_fits.txt', 'Open the fit of interest');
    data=double(dlmread(strcat(PathName,'/',FileName),',',2,1));
    XX=data(:,1);
    FitX=data(:,2);
end  
 %% Demo: make a simple trace with a number of step trains  or spikes
demo=0;
if nargin==0, demo=1;end
if demo 
    wmin=3;
    [XX,FitX]=MakeSlopedStepsTrace;   
end
% ---------------------------------------------

 Lx=length(XX);
 TT=(1:Lx)';
[step_props,~,~]=Get_StepsFromFit(TT,FitX);  
oriFitX=FitX;

%% remove slopes
for rp=1:2
    [~,~,~,together]=find_slopesteps(step_props(:,1), step_props(:,5), wmin);  
    FitX=Get_FitFromSteps(XX,together,'median');
    [step_props,~,~]=Get_StepsFromFit(TT,FitX);
end

if nargin <3
    figure(63);
    close(gcf);
    figure(63);
    plot(XX,'k-','LineWidth', 2); hold on;
    plot(oriFitX, 'c-','LineWidth', 1.5);
      plot(FitX, 'r-', 'LineWidth', 2);
     legend('data','original fit','slope-merged fit');
     title('SlopeMerger');
     dum=1;
end

if nargin==1
     Axz=(1:length(XX))';
     data=[Axz XX FitX];
     outname1=strcat(PathName,'/',FileName(1:end-4),'_merged.txt');
     dlmwrite(outname1,data);
     outname2=strcat(PathName,'/',FileName(1:end-4),'_merged_props.txt');
     dlmwrite(outname2,step_props);
end

function  [steptrains,isolated,merged,together]=find_slopesteps(index, steps, wmin)
%find all connected trains of steps; output is a structure 'steptrains' and
%sets of indices: 
%'isolated', 
%'merged'  set of one index of thelargest step per step train
%together: 'isolated' and 'merged' combined
%----------------------------------------------

ls=length(steps);
%get indices close to at least one another AND same sign
dt_fw=diff(index);
ds_fw=diff(sign(steps));
sel=find((dt_fw<=wmin)&ds_fw==0);                   %nearby with same sign
if ~isempty(sel)
narrows=unique(sort([index(sel) ; index(sel+1)]));  %these are potential trains
isolated=index(find(~ismember(index,narrows)));     %these are isolated
merged=[];
%proceed with points nearby another
%initialize first set
first_ix2add=narrows(1);                    %first index of set              %
sign2add=sign(steps(index==first_ix2add));  %sign of set
steptrains(1).ix=first_ix2add;  
steptrains(1).sign=sign2add;
set_index=1; 
keep=narrows(find(first_ix2add~=narrows));
keep_steps=steps(find(ismember(index,keep)));    
while length(keep)>0
    set_grows=1;
while set_grows   
    current_ones=steptrains(set_index).ix;
    Lset=length( current_ones);
    new_ones=[];
    for ii=1:Lset  %for all set members
        this_i=steptrains(set_index).ix(ii);
        this_trainsign=steptrains(set_index).sign;
        dt=abs(keep-this_i);                    %time difference
        ds=sign(keep_steps)-this_trainsign;     %sign difference
        sel=find((dt<=wmin)&(ds==0));           %includes same-sign
        %%find one new point per existing to the train ; update remainder:
        if ~isempty(sel)  
            one2add=keep(sel(1));
            new_ones=[new_ones one2add];   
            keep=keep(find(keep~=one2add));  %remove from remainder
            keep_steps=steps(find(ismember(index,keep)));            
        end
    end
    %add the set of new points to the train
    if ~isempty(new_ones)
        new_ones=unique(new_ones);
        steptrains(set_index).ix=unique(sort([current_ones new_ones]));
    else  %set complete; get step closest to middle and go to next set
        this_train=steptrains(set_index).ix';
        these_steps=abs(steps(ismember(this_train,index)));
        %1) find largest step       
        %[~,idx]=max(these_steps);
        %2) closest to middle
        [~,idx]=min(abs(this_train-mean(this_train)));
        
        
        merged_train=this_train(idx);
        merged=[merged ; merged_train];
        %prepare for next
        set_grows=0;
        set_index=set_index+1; 
        %initialize the next step train:
        if ~isempty(keep)
            first_ix2add=keep(1);                 %get one new point
            sign2add=sign(steps(index==first_ix2add));  %sign of set           
            keep=keep(find(keep~=first_ix2add));  %remove from remainder
            keep_steps=steps(find(ismember(index,keep)));
            steptrains(set_index).ix=first_ix2add; %first index of next set
            steptrains(set_index).sign=sign2add;   %sign of set
        end
    end
end
end
together=sort(unique([isolated ; merged]));
else
    steptrains=[];
    merged=[];
    isolated=index;
    together=index;
end

function [XX,FitX]=MakeSlopedStepsTrace
%For demo purposes
Lxi=1000;
noize=1;                   %noise level
XX=[]; FitX=[];
stepsizes=[-5:2:5];
isolated=[1 [10:50:Lxi]];  %just a step train
indices=isolated;
Li=length(indices);
for ii=1:Li
    dice=ceil(length(stepsizes)*(rand(1,1)));
    levs(ii)=stepsizes(dice);  %add random step levels
end
%build initial trace
for ii=1:Li-1
    L_sect=indices(ii+1)-indices(ii)+1;
    XX=[XX ; noize*randn(L_sect, 1)+levs(ii)];
    FitX=[FitX ; zeros(L_sect, 1)+levs(ii)];
end 

%smooth both
XX=smooth(XX,10);
FitX=round(smooth(FitX,10));

function [idx, avl, avr,rankit,errorcurve]=Splitfast(Segment)              %
%this function also adresses a one-dim array 'Segment'
%and determines the best step-fit there
%To save time, functions like 'mean' are avoided
    w=length(Segment);   
    Chisq=(1:w-1)*0;  
    if w>3
		Chisq=(1:w-1)*0;                           
        AvL=Segment(1);    AvR=sum(Segment(2:w))/(w-1); AvAll=sum(Segment)/w;  
        for t=2:w-2
            AvL=(AvL*(t-1)+Segment(t))/t;     AvR=(AvR*(w-t+1)-Segment(t))/(w-t);
            DStepL=AvL-AvAll;           DStepR=AvR-AvAll;
            DeltaChisq=((DStepL.^2)*t+(DStepR.^2)*(w-t));            
            Chisq(t)=-DeltaChisq;       
        end
         [~,idx]=min(Chisq(2:w-2)); idx=idx+1;
         avl=mean(Segment(1:idx));                    avr=mean(Segment(idx+1:w));
         %rankit=(avr-avl)^2/(1/(idx-1)+1/(w-idx-1));  %quantity expresing predicted accuracy step (squared)
         rankit=(avr-avl)^2*w;
    else                                            %low-limit cases
         rankit=0;               
         switch w
             case 3 
                a=(Segment(2)+Segment(3))/2-Segment(1);            b=Segment(3)-(Segment(1)+Segment(2))/2;
                cL=[Segment(1) , (Segment(1)+Segment(2))/2];      cR=[(Segment(2)+Segment(3))/2 , Segment(3)]; 
                [~,idx]=max([a b]);    avl=cL(idx);  avr=cR(idx);
            case 2
                idx=1;  avl=Segment(1); avr=Segment(2); 
        end
    end
    errorcurve=Chisq/(w-1);
    
    function [error_st_boot, error_t_boot, t_refit,error_t_refit]=bootstrap_get_errors(data, indices,bootstraprepeats)
% bootstrap 
%Aim: bootstrap to get location error (plus a time estimate)

%use: [error_t_boot,error_t_refit,error_st_boot]=get_errors_by_bootstrap(data, indices)

%Input: 
%data: original data, sigle column, 
%indices: locations of stepfits, in pts
%

%Output:
% error_t_boot: 95% confidence range of error in time (pts) by bootstrapping
% t_refit: location obtained by re-fitting
% error_t_refit: difference of input step location and location obtained by re-fitting
% error_st_boot:  95% confidence range of error in sep size by bootstrapping

%References: 
%following method:
%[1] Received 18 Jul 2015 | Accepted 16 Nov 2015 | Published 17 Dec 2015 
% ATP hydrolysis assists phosphate release and promotes reaction ordering in F1-ATPase 
% Chun-Biu Li1, Hiroshi Ueno2, Rikiya Watanabe2,3,4, Hiroyuki Noji2,4 & Tamiki Komatsuzaki1


%% bootstrap main cycle
%loop all segments
N_steps=length(indices-2);
error_t_refit=zeros(N_steps,1);
error_t_boot=zeros(N_steps,1);
error_st_boot=zeros(N_steps,1);

t_refit=zeros(N_steps,1);

indices_ext=[0; indices; length(data)];

for ii=1:length(indices)
    %% 1) our usual step fit to get step location (plus extra export of error curve)
    Segment=data(indices_ext(ii)+1:indices_ext(ii+2));
    idx_old=indices(ii);  %input location
    [idx, ~, ~,~, error_curve]=Splitfast(Segment) ;
    
    %note that we also obtain a new estimate for the best fit in this
    %segment:
    t_refit(ii)=indices_ext(ii)+idx;
    error_t_refit(ii,1)=abs(t_refit(ii)-idx_old);
    
    %% 2 bootstrapping
    % repeat plateau fits left and right by bootstrapping 
    %(location is kept constant)
    tic
     %repeat many times:
    leftpart=Segment(1:idx);        L_left=length(leftpart);
    rightpart=Segment(idx+1:end);   L_right=length(rightpart);

    if ((L_left>20)&(L_right>20))    
        [bootstat_av_left,bootsam_left] = bootstrp(bootstraprepeats,@mean,leftpart);  %resample left, get many left averages
        [bootstat_av_right,bootsam_right] = bootstrp(bootstraprepeats,@mean,rightpart); %resample right, get many right averages
        stepsize_boot=bootstat_av_right-bootstat_av_left;

        %% 3 Get the value Chi-square(idx) for all resamplings.
        %To gain time, work matrix-wise :
        %expand the plateau results in a block, such that every column is a new step-fit
        %build with the resampled left and right averages
        newstep_indices=[bootsam_left; idx+bootsam_right];  %matrix of indices
        NewSegments=Segment(newstep_indices);
        NewFits=[repmat(bootstat_av_left,1,L_left) repmat(bootstat_av_right,1,L_right)]';

        all_chi_squares=(mean((NewSegments-NewFits).^2)).^0.5;
        bootstrap_error_of_minimum_value=1.96*std(all_chi_squares);



        %% 4 get estimate of 95% confidence range of value of Errcurv at minimum
        %to increase precision, we interpolate 10-fold near the minimum
        LL=length(error_curve);
        nearmin_lox=max([1 idx-50]);
        nearmin_hix=min([LL idx+50]);
        axz=nearmin_lox:nearmin_hix;
        axz_ip=nearmin_lox:0.1:nearmin_hix;
        Errcurv_blowup=interp1(axz,error_curve(axz),axz_ip);

        sel=find(abs(Errcurv_blowup-error_curve(idx))<=bootstrap_error_of_minimum_value);
        lox_ip=min(sel);  
        hix_ip=max(sel);

        lox=axz_ip(lox_ip); 
        hix=axz_ip(hix_ip);

        
        error_t_boot(ii,1)=(hix-lox)/2;
        error_st_boot(ii,1)=1.96*std(stepsize_boot)/(2^0.5);
    else %too short plateaus, do not bootstrap
        error_t_boot(ii,1)=NaN;
        error_st_boot(ii,1)=NaN;
    end
end

function Save_AS_formats(TT,XX,FitX,step_props,SaveName, PathName);
%save data identical to AutoStepfinder format
    [~, cols]=size(step_props);
     FinalSteps=step_props;
     curpth=pwd;
     cd(PathName);
     time_resolution=median(diff(TT));
        
    %Step properties
      IndexStep                 =  FinalSteps(:,1);             %Index where step occured
      TimeStep                  =  FinalSteps(:,2)...           %Time when step occured
                                   *time_resolution; 
      LevelBefore               =  FinalSteps(:,3);             %Level before step
      LevelAfter                =  FinalSteps(:,4);             %Level after step
      StepSize                  =  LevelAfter - LevelBefore;            %Size of step
      DwellTimeStepBefore       =  FinalSteps(:,6)...           %Dwelltime step before
                                   *time_resolution;       
      DwellTimeStepAfter        =  FinalSteps(:,7)...           %Dwelltime step after
                                   *time_resolution;  
                               
      if cols==7;                       
      properties_table          = table(IndexStep,TimeStep,...          %Save variables in table
                                  LevelBefore,LevelAfter,StepSize,...
                                  DwellTimeStepBefore,DwellTimeStepAfter);
      end                        
                               
      if cols==8;                  
      StepError                 =  FinalSteps(:,8);             %Error of each step      
      properties_table          = table(IndexStep,TimeStep,...          %Save variables in table
                                  LevelBefore,LevelAfter,StepSize,...
                                  DwellTimeStepBefore,DwellTimeStepAfter,StepError);
      end
      if cols>8     %bootstrap step errrors 
          StepError                 =  FinalSteps(:,8);             %Error of each step  
          StepError_boot       =  FinalSteps(:,9);             %Error of each step
          DwellError_boot      =  FinalSteps(:,10)*time_resolution;      %Dwelltime step after
          properties_table          = table(IndexStep,TimeStep,...          %Save variables in table
                                  LevelBefore,LevelAfter,StepSize,...
                                  DwellTimeStepBefore,DwellTimeStepAfter,StepError,...
                                  StepError_boot, DwellError_boot);  
      end
                       
      fits_table                = table(TT, XX, FitX);          %Save variables in table                   

      %% Save
      cd(PathName);       
      writetable(fits_table, [SaveName,'_fits.txt']);                   %Save table containing fits                       
      writetable(properties_table, [SaveName,'_properties.txt']);       %Save table containing properties                                
      cd(curpth);

function StepsX=AddStep_Errors(XX,StepsX)   
%This function calculates step errors associated with the steps.
%from the standard deviation of the adjacent plateaus

[ls,col]=size(StepsX); 
i1=0;
for ii=1:ls
    i2=StepsX(ii,1);
    if ii<ls
        i3=StepsX(ii+1,1);
    else 
        i3=length(XX);
    end    
    Nbefore=i2-i1;    
    Nafter=i3-i2;
      
    rmsbefore=std(XX(i1+1:i2));
    rmsafter=std(XX(i2+1:i3)) ;
    StepsX(ii,col+1)=2*((rmsbefore^2/Nbefore+rmsafter^2/Nafter)^0.5)/2^0.5; %plus minus 95%
    i1=i2;
end


% --- Executes just before StepMerger is made visible.
function StepMerger_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to StepMerger (see VARARGIN)

% Choose default command line output for StepMerger
handles.output = hObject;
pwd=cd;
set(handles.dirboth,'enable','On')
set(handles.dirup,'enable','On')
set(handles.dirdown,'enable','On')
set(handles.width,'enable','On')
set(handles.margin,'enable','On')
set(handles.widthmerge,'enable','Off')
set(handles.errorestoff,'Value',1)
set(handles.erroreston,'Value',0)
set(handles.directory,'String',pwd)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes StepMerger wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = StepMerger_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function widthmerge_Callback(hObject, eventdata, handles)
% hObject    handle to widthmerge (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of widthmerge as text
%        str2double(get(hObject,'String')) returns contents of widthmerge as a double


% --- Executes on button press in go.
function go_Callback(~, ~, handles)
% hObject    handle to go (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

Step_merger_main(handles)


function directory_Callback(hObject, eventdata, handles)
chckfldr=get(hObject,'String');
chckfldr= exist(chckfldr);
if chckfldr ~= 7,  msgbox('The provided directory is not valid.','ERROR', 'error')
    pwd=cd;
    set(hObject, 'String', cd);
    return; end



function margin_Callback(hObject, eventdata, handles)
% hObject    handle to margin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of margin as text
%        str2double(get(hObject,'String')) returns contents of margin as a double




function width_Callback(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width as text
%        str2double(get(hObject,'String')) returns contents of width as a double


% --- Executes on button press in despikeon.
function despikeon_Callback(~, ~, handles)
set(handles.dirboth,'enable','On')
set(handles.dirup,'enable','On')
set(handles.dirdown,'enable','On')
set(handles.width,'enable','On')
set(handles.margin,'enable','On')
set(handles.widthmerge,'enable','Off')

% --- Executes on button press in mergeon.
function mergeon_Callback(~, ~, handles)
set(handles.dirboth,'enable','Off')
set(handles.dirup,'enable','Off')
set(handles.dirdown,'enable','Off')
set(handles.width,'enable','Off')
set(handles.margin,'enable','Off')
set(handles.widthmerge,'enable','On')
