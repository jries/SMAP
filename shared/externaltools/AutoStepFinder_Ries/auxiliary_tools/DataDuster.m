%% AutoStepfinder: A fast and automated step detection method for single-molecule analysis.
%Luuk Loeff*, Jacob Kerssemakers*, Chirlmin Joo & Cees Dekker.
% * Equal contribution
% DataDuster algorithm: remove inf, nan and other artefacts from data
% Last update: November 2018
%% Concise de overview: for details and explanation refer to main text.
%Line 13-67:    Initialization of GUI   
%Line 68-136:	Load data
%Line 137-152:  Check input files for errors
%Line 153-166:  Run mode DataDuster
%Line 167-215:  Identify indices of improperly formatted datapoint and clean up
%Line 216-357:  GUI 
%% Initialization of GUI  
function varargout = DataDuster(varargin)
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @DataDuster_OpeningFcn, ...
                   'gui_OutputFcn',  @DataDuster_OutputFcn, ...
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

% --- Executes just before DataDuster is made visible.
function DataDuster_OpeningFcn(hObject, ~, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to DataDuster (see VARARGIN)
set(handles.neighbor,'Value',1);
set(handles.single,'Value',1);
set(handles.columns,'enable','Off');
set(handles.checkbox1,'Value',1);
set(handles.single,'Value',1);
initval.columsanalysis=str2double(get(handles.columns,...               %Threshold for second round of fitting
                                 'string'));


% Choose default command line output for DataDuster
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function path_CreateFcn(hObject, ~, ~)

    data_directory = pwd;                           %insert your default directory here
    set(hObject,'String', num2str(data_directory));

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function datacleaner_func(handles)
%% Read and clean data
clc
%1) single or multiple files from GUI
disp('Loading..');
initval.datapath=get(handles.path,...               
 'string');
initval.single=get(handles.single,...            
 'value'); 
if initval.single == 1
    initval.hand_load=1;
else
    initval.hand_load=2;
    initval.datapath=uigetdir(initval.datapath);
end
initval.nextfile=1;   
CurrentFolder=initval.datapath; 

%2) Work single or multiple files
while initval.nextfile>0      
        switch initval.hand_load
            case 1                                              %single file
            try
            cd(initval.datapath);
            catch
            errorpath=['Data path is not valid.'];              %Error for datapath
            msgbox(errorpath,'ERROR', 'error')
            display('Please provide an existing datapath')
            return
            end
            [FileName,PathName] = uigetfile('*.*','Select the signal file'); 
            cd(CurrentFolder);
            source=strcat(PathName,FileName);
            try
            data=double(dlmread(source));
            catch
            errorfile=[num2str(FileName),' is not formatted properly.']; %Error that occurs when data is not properly formatted
            msgbox(errorfile,'ERROR', 'error')
            display('Please provide an .txt file that is free of characters')
            return
            end
            SaveName=FileName(1:length(FileName)-4);
            initval.nextfile=0;      
            case 2                                      %Batch style processing
            fileindex=initval.nextfile;
            try
            cd(initval.datapath);
            catch
            return
            end
            PathName=initval.datapath;
            AllFileNames=dir('*.txt');
            cd(CurrentFolder);
            if isempty(AllFileNames)  
                msgbox('The provided input folder is empty.','ERROR', 'error')
            return; end;
            AllFiles=length(AllFileNames);
            FileName=AllFileNames(fileindex).name;
           try
            data=double(dlmread(strcat(initval.datapath,'\',FileName)));
            catch
            errorfile=[num2str(FileName),' is not formatted properly.'];
            msgbox(errorfile,'ERROR', 'error')
            display('Please provide an .txt file that is free of characters')
             return
           end
            if initval.nextfile==AllFiles
               initval.nextfile=0;
                else
                 initval.nextfile=fileindex+1;          
            end
        end


        %% Check columns of input file: cleaning is performed per column if more than one is specified
        %missing=NaN
        %data=fillmissing(data,'constant',missing)
        [initval.rowssfile,initval.columnsfile]=size(data);       %get data file size
        pads=num2str(max([1 (ceil(log10(initval.columnsfile)))]));
        precisionstring=strcat('%0',pads,'.0f');  %proper range of zeros
        
        initval.columsanalysis=str2double(get(handles.columns,... %read input for number of columns
                                              'string')); 
        initval.allcols=get(handles.checkbox1,'value');            %read input for number of columns                                               
        if ~initval.allcols
               if initval.columsanalysis > initval.columnsfile
                msgbox('The input for the number of columns exceeds the columns in the file. Please specify a different column number.','ERROR', 'error')
               end 
               data=data(:,initval.columsanalysis) ;     
        end
        %% Which mode to run?
        initval.mean=get(handles.mean,'value');               %read input for number of mean                        
        if initval.mean == 1,   initval.replacemode=1;     lbl='mean';     end
        initval.median=get(handles.median,'value');           %read input for number of median                       
            if initval.median == 1, initval.replacemode=2; lbl='median';     end
        initval.neighbor=get(handles.neighbor,'value');       %read input for number of neighbor                       
        if initval.neighbor == 1,   initval.replacemode=3; lbl='neighbours';    end
        initval.remove=get(handles.remove,'value');           %read input for number of remove                       
        if initval.remove == 1,     initval.replacemode=4; lbl='removal';     end
        
        outdir=strcat(PathName, '\cleaned_data_using_',lbl,'\');
        if ~isdir(outdir), mkdir(outdir);end
        

        %% Clean-up data per column
        [~, ncols]=size(data);       
        for cc=1:ncols
            thiscoldata=data(:,cc);
            %first, check if this column a time axis-this will be excluded from further
            %cleaning
            dif=thiscoldata(2:end)-thiscoldata(1:end-1);
            iftime=sum(dif==median(dif))>0.9*length(dif);
            if ~iftime  %proceed cleaning
                thiscoloridata=thiscoldata;  %keep it around for comparison
                specificbadval=NaN;  %option for specific value search - silenced 
                bad_idx=find(isinf(thiscoldata)|isnan(thiscoldata)|(thiscoldata==specificbadval));
                good_idx=find(~isinf(thiscoldata)&~isnan(thiscoldata)&(thiscoldata~=specificbadval));
                [nfoundrows, nfoundcols]=size(bad_idx);
                nfoundrows_rem=[num2str(nfoundrows), ' value(s) have been removed from your data'];
                nfoundrows_repl=[num2str(nfoundrows), ' value(s) have been replaced within your data'];
                switch initval.replacemode
                    case 1
                        thiscoldata(bad_idx)=mean(thiscoldata(good_idx));
                        disp(nfoundrows_repl);
                    case 2
                        thiscoldata(bad_idx)=median(thiscoldata(good_idx)); 
                        disp(nfoundrows_repl);
                    case 3
                        LB=length(bad_idx);
                        for ii=1:LB
                            left_i=max(find(good_idx<bad_idx(ii)));  %first good one on left
                            right_i=min(find(good_idx>bad_idx(ii)));  %first good one on right
                            padval=mean(thiscoloridata([good_idx(left_i) good_idx(right_i)]));
                            thiscoldata(bad_idx(ii))=padval;
                        end
                        disp(nfoundrows_repl);
                    case 4
                        thiscoldata=thiscoldata(good_idx);
                        disp(nfoundrows_rem);
                end  
                if initval.allcols
                    filename_loop=['cleaned_',FileName(1:end-4),'_col_',num2str(cc,precisionstring),'.txt'];
                else
                    filename_loop=['cleaned_',FileName(1:end-4),'_col_',...
                        num2str(initval.columsanalysis,precisionstring),'.txt'];
                end
                save([outdir, filename_loop],'thiscoldata','-ascii');
            end
    end
end
msgbox('Your data has been cleaned','DataDuster', 'help')
display('Your data has been cleaned')
% --- Outputs from this function are returned to the command line.


function varargout = DataDuster_OutputFcn(~, ~, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in cleandata.
function cleandata_Callback(~, ~, handles)
% hObject    handle to cleandata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
datacleaner_func(handles);

% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, ~, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
columsanalysis=get(hObject,'Value');
if columsanalysis==0
set(handles.columns,'enable','On');
else
set(handles.columns,'enable','Off');
end

function columns_Callback(hObject, ~, handles)
% hObject    handle to columns (see GCBO)
initval.columsanalysis=str2double(get(handles.columns,...               %Threshold for second round of fitting
                                 'string'));
                             manualmode=get(hObject,'String');
checkcl=isnan(initval.columsanalysis);
if checkcl==1
         msgbox('The input for the number of columns is NaN.','ERROR', 'error')
     return;
end

% Hints: get(hObject,'String') returns contents of columns as text
%        str2double(get(hObject,'String')) returns contents of columns as a double

% --- Executes during object creation, after setting all properties.
function columns_CreateFcn(hObject, ~, ~)
% hObject    handle to columns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function columns_DeleteFcn(~, ~, ~)
% hObject    handle to columns (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function path_Callback(hObject, ~, ~)
% hObject    handle to path (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of path as text
%        str2double(get(hObject,'String')) returns contents of path as a double
chckfldr=get(hObject,'String');
chckfldr= exist(chckfldr);
if chckfldr ~= 7,  msgbox('The provided directory is not valid.','ERROR', 'error')
     return; end;

% --- Executes on button press in mean.
function mean_Callback(hObject, ~, ~)
% hObject    handle to mean (see GCBO)
initval.mean=get(hObject,'Value');
if initval.mean == 1
    initval.replacemode=1;
end

% --- Executes on button press in median.
function median_Callback(hObject, ~, ~)
% hObject    handle to median (see GCBO)
initval.median=get(hObject,'Value');
if initval.median == 1
    initval.replacemode=2;
end

% --- Executes on button press in neighbor.
function neighbor_Callback(hObject, ~, ~)
% hObject    handle to neighbor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
initval.neighbor=get(hObject,'Value');
if initval.neighbor == 1
     initval.replacemode=3;
end

% --- Executes on button press in remove.
function remove_Callback(hObject, ~, ~)
% hObject    handle to remove (see GCBO)
initval.remove=get(hObject,'Value');
if initval.remove == 1
      initval.replacemode=4;
end


% --- Executes during object creation, after setting all properties.
function cleandata_CreateFcn(~, ~, ~)
% hObject    handle to cleandata (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in single.
function single_Callback(hObject, ~, ~)
% hObject    handle to single (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
initval.single=get(hObject,'Value');
if initval.single == 1
    initval.hand_load=1;
end
% Hint: get(hObject,'Value') returns toggle state of single


% --- Executes on button press in batch.
function batch_Callback(~, ~, ~)
% hObject    handle to batch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% initval.batch=get(hObject,'Value');
% if initval.batch == 1
%     initval.hand_load=2;
% end
% Hint: get(hObject,'Value') returns toggle state of batch
