function varargout = mb_GUI(varargin)
dbstop if error
addpath(genpath('functions'))
warning('off','all')
% MB_GUI M-file for mb_GUI.fig
%% Primary Objective: This is the interface for users to process glaciological data collected by the US Geological Surveys Benchmark Glacier Program
%   Program has five primary functions
%   A) Initalize Guided User Interface
%       1) Start GUI
%       2) Load glacier site map, sites, and years for glacier selected by user

%   B) Load glacier site map, sites, and years for glacier selected by user

%   C) Update Functions
%       1) Fill data gaps in primary weather station data (idealy a weather
%          station in the glaciers basin
%       2) Recalibrate Mass Balance Model
%           a) Calibrate precipitation model first, since it requires no
%               a priori assumptions
%           b) Calibrate snow melt model with calibrated precipitation
%               model. Requires no a priori assumptions for ice melt.
%           c) Calibrate ice melt model, using calibrated precipitation and
%               snow melt models
%       3) Fill missing glaciological observations using calibrated mass
%           balance model
%   D) Calculated Glacier-wide Mass balance
%       1) Users selects sites to use, years to proces, lapse rate to
%          apply, integration method(index or gradient), integration surface
%          (Conventional or Reference-Surface), time-system to use
%          (stratigraphic, Fixed-date Hydrological year, or Fixed-date
%          stratigraphic), and geodetic calibration (none, universal, or
%          piece-wise; See section D part 3
%       2) Glacier-wdie seasonal and annual balances are calculated using
%          the inputs described in part one of section D
%       3) Cumulative annual mass balance time-series is geodetically
%       calibrated using either a universal (best-fit to all data)
%       correction to glaciological balances, or a piece-wise (absolute fit
%       per geodetic measurement interval)
%       weather station in the 
%   E) Plot Glacier-wide mass balance time-series
%       1) Generate two plots to be displayed in the GUI once parts D is
%          complete.
%           a) Cumulative Glacier wide annual time-series uncalibrated and
%              calibrate(if selected).
%           b) Seasonal and Anuual balances per year
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   A) Initalize Guided User Interface
%       1) Start GUI 
%
%      MB_GUI, by itself, creates a new MB_GUI or raises the existing
%      singleton*.
%
%      H = MB_GUI returns the handle to a new MB_GUI or the handle to
%      the existing singleton*.
%
%      MB_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MB_GUI.M with the given input arguments.
%
%      MB_GUI('Property','Value',...) creates a new MB_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mb_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mb_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Edit the above text to modify the response to help mb_GUI
% Last Modified by GUIDE v2.5 13-Jan-2019 14:56:34
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @mb_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @mb_GUI_OutputFcn, ...
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

%     2) Load glacier site map, sites, and years for glacier selected by user
function mb_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for mb_GUI
addpath functions/ %addpath to functions
handles.output = hObject;
mydir = pwd; % get working directory
% set(handles.MYDIR_EDIT,'String',mydir)
handles.params.mydir = mydir;
% Choose default command line output for setup_inputs
handles.output = hObject;

%%%%%% setup site map %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Gulkana by default %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
axes(handles.Site_Map);
image(imread(['data/Gulkana/Input/Input_Gulkana_Glaciological_Sites.jpg']));%import glacier map
axis off          
axis image 
        
Glaciological_data=readtable(['data/Gulkana/Input/Input_Gulkana_Glaciological_Data.csv'],'Format','%f%s%s%s%f%f%f%f%f');%import glaciological data
all_years=unique(Glaciological_data.Year);%get years of record from glaciological data
selectable_years=[num2str(all_years(1)),':',num2str(all_years(end))];%string of years for display in GUI
set(handles.YEAR_EDIT,'String',selectable_years)
handles.my_yr = selectable_years;
        
all_sites=unique(Glaciological_data.site_name);%get list of sites from glaciological data
if length(all_sites)>36
    fprintf(1,'ERROR: Need to add more buttons to match the number of sites')
else
    for site=1:36
        if site <= length(all_sites);
        else
            all_sites=[all_sites;'  '];
        end
    end
end
%set handles for site names for GUI display
for i=1:length(all_sites)
    set(getfield(handles, ['SITE',num2str(i)]),'string',all_sites(i,1));
end
    
% Update handles structure
guidata(hObject, handles);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% B) Load glacier site map, sites, and years for glacier selected by user
% --- Executes on button press in GLACIER_ID.
function GLACIER_ID_Callback(hObject, eventdata, handles)
% hObject    handle to GLACIER_ID (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
glaciers = get(handles.Glacier_Selection,'String');%names of selectable glaciers in the GUI
glacier=cell2mat(glaciers(get(handles.Glacier_Selection,'value')));%name of glacier selected by users
Glaciological_data=readtable(['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv'],'Format','%f%s%s%s%f%f%f%f%f');%import glaciological data for that glacier       
default_settings=readtable(['data/',glacier,'/Input/Input_',glacier,'_defaultsettings.csv']);
record_years=unique(Glaciological_data.Year);%years of glaciological record
years=[num2str(record_years(1)),':',num2str(record_years(end))];%years selected by user

set(handles.YEAR_EDIT,'String',years)
handles.my_yr = years;
all_sites=unique(Glaciological_data.site_name);%get list of sites from glaciological data
if length(all_sites)>36
    fprintf(1,'ERROR: Need to add more buttons to match the number of sites')
    [s,v] = listdlg('PromptString','Select Sites',...
                'SelectionMode','multiple',...
                'ListString',all_sites);
    all_sites=all_sites(s);
    for i=1:36
        set(getfield(handles, ['SITE',num2str(i)]),'string','');
        set(getfield(handles, ['SITE',num2str(i)]),'value',0);
    end

else
    %set handles for site names for GUI display
    for i=1:36
        if i<=length(all_sites)
            set(getfield(handles, ['SITE',num2str(i)]),'string',all_sites(i,1));
        else
            set(getfield(handles, ['SITE',num2str(i)]),'string','');
        end
    end
end
    
axes(handles.Site_Map);
image(imread(['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Sites.jpg']));%import glacier map
axis off          
axis image        
       
       
    
guidata(hObject,handles)
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C) Update Functions
function UPDATE_BUTTON_Callback(hObject, eventdata, handles)
% Get user input from GUI
current_directory = pwd;                                            %gets present working directory
glaciers = get(handles.Glacier_Selection,'String');
glacier=cell2mat(glaciers(get(handles.Glacier_Selection,'value')));  %get glacier entered
integration_method = get(handles.GRAD_POPUP,'Value');   %get gradients you are using
%Glaciological_data=readtable(['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv'],'Format','%f%s%s%s%f%f%f%f%f');
Glaciological_data=readtable(['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv']);


record_years=unique(Glaciological_data.Year);
years=[record_years(1),':',record_years(end)];

AAD = importdata(['data/',glacier,'/Input/Input_',glacier,'_Area_Altitude_Distribution.csv']); 

lapse_rates=get(handles.Lapse_Rate_Selection,'String');
lapse_rate=str2num(cell2mat(lapse_rates(get(handles.Lapse_Rate_Selection,'Value'))));

cd(current_directory)
addpath functions
%
switch get(handles.UPDATE_POPUP,'Value')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C) Update Functions
%%    1) Fill data gaps in primary weather station data (idealy a weather
%          station in the glaciers basin

    case 1 %correct temp and precip data
        disp('        %%%%%   Filling data gaps in weather data    %%%%%          ')
        [ready]=fillWxData(glacier); %fill any missing data within primary weather station data
        if ready==1
            set(handles.UPDATE_BUTTON,'string','You can now Calibrate the Mass Balance Model');
        else
            set(handles.UPDATE_BUTTON,'string','Weather data has NaNs');
        end
        cd(current_directory)
        guidata(hObject, handles);
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% C) Update Functions
%%     2) Recalibrate Mass Balance Model
%           a) Calibrate precipitation model first, since it requires no
%               a priori assumptions
%           b) Calibrate snow melt model with calibrated precipitation
%               model. Requires no a priori assumptions for ice melt.
%           c) Calibrate ice melt model, using calibrated precipitation and
%               snow melt models

    case 2 %Calibrate the MAss Balance Model%  
        if get(handles.PLOTABLATIONMODEL,'Value')==1 %if plot ablation model has been selected warn user that this will produce a lot of figures
           promptMessage = sprintf('You selected to plot the ablation model for more than one year. This probably blow up your machine. Do you still want to?');
            titleBarCaption = 'Plotting ablation model';
            choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know');
                    if strcmpi(choice, 'Yes')
                        plot_ablation_model = 1;
                        msgbox('Ok... Take cover!')
                    elseif strcmpi(choice, 'No')
                        plot_ablation_model=0;
                    end
        else
                plot_ablation_model=0;
        end
        years = eval(get(handles.YEAR_EDIT,'String'));%get years input from user
        lapse_rates=get(handles.Lapse_Rate_Selection,'String');%get list of lapse rates available for user
        lapse_rate=str2num(cell2mat(lapse_rates(get(handles.Lapse_Rate_Selection,'Value'))));%get lapse rate selected by user
        Weather_data=readtable(['data/',glacier,'/Intermediate/',glacier,'FilledWx.csv']);%import weather data for glacier selected
        Weather_data.Date=datenum(Weather_data.Date);%convert dates to datenumbers
        disp('        %%%%%   Calibrating Mass Balance Model    %%%%%          ')
        plot_calibration=1;
        [ready,~, precipitation_ratio_table,~,~,Degree_Day_Factor_table]= Calibrate_Precipitation_and_Ablation_models(glacier,years,Glaciological_data,Weather_data,AAD,lapse_rate,plot_ablation_model,1) ;
        Degree_Day_Factors_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Degree_Day_Factors.csv'];
        precipitation_ratios_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Precipitation_Ratios.csv'];%path for precipitation ratios
        writetable(precipitation_ratio_table,precipitation_ratios_path)
        writetable(Degree_Day_Factor_table,Degree_Day_Factors_path)
        if ready==1
            set(handles.UPDATE_BUTTON,'string','Can now estimate missing obs.');
        else            
            set(handles.UPDATE_BUTTON,'string','Uh oh. Your data is f!*$@%');
        end
        cd(current_directory)
        guidata(hObject, handles);
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%% C) Update Functions
%%        3) Fill missing glaciological observations using calibrated mass balance model for glacier selected
            %import calibrated DDFs              
    case 3 %estimate missing observations     
            Degree_Day_Factors_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Degree_Day_Factors.csv'];
            if ~exist(Degree_Day_Factors_path) %if file doesn't exist then we need approximate DDFs to start
                disp('ERROR: You need to Calibrate the mass balance model!')
            else
                meltrates=readtable(Degree_Day_Factors_path);%otherwise import previously calibrated DDFs
                ks=meltrates.ks(end);   %DDF for snow
                ki=meltrates.ki(end);   %DDF for ice
            end 
            %import calibrated precipitation ratios
            precipitation_ratios_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Precipitation_Ratios.csv'];%path for precipitation ratios
            if ~exist(precipitation_ratios_path)%if file doesn't exist then we assume prcipitation rates are one-to-one for initial value
                disp('ERROR: You need to Calibrate the mass balance model!')
            else
                precipitation_ratios_table=readtable((precipitation_ratios_path));
            end
            %if users selected to plot 'Find_Mass_Maximum_and_Minimum.m'
            %ask them if they really want that mank plots
        if get(handles.PLOTABLATIONMODEL,'Value')==1
           promptMessage = sprintf('You selected to plot the ablation model for more than one year. This probably blow up your machine. Do you still want to?');
            titleBarCaption = 'Plotting ablation model';
            choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know');
                    if strcmpi(choice, 'Yes')
                        plot_ablation_model = 1;
                        msgbox('Ok... Take cover!')
                    elseif strcmpi(choice, 'No')
                        plot_ablation_model=0;
                    end
           
        else
                plot_ablation_model=0;
        end
        if get(handles.PLOTINTEGRATION,'Value')==1
               promptMessage = sprintf('You selected to plot integrations for more than one year. This might take a while. Do you still want to plot integrations?');
                titleBarCaption = 'Plotting Integrations';
                choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know');
                        if strcmpi(choice, 'Yes')
                            plot_integration = 1;
                        elseif strcmpi(choice, 'No')
                            plot_integration=0;
                        end

        else
                plot_integration=0;
        end
        years = eval(get(handles.YEAR_EDIT,'String'));%get years selected by user
        lapse_rates=get(handles.Lapse_Rate_Selection,'String');%list lapse rates available for user to select
        lapse_rate=str2num(cell2mat(lapse_rates(get(handles.Lapse_Rate_Selection,'Value'))));%lapse rate selected by user
        AAD = importdata(['data/',glacier,'/Input/Input_',glacier,'_Area_Altitude_Distribution.csv']); %import AAD for selected glacier
        Weather_data=readtable(['data/',glacier,'/Intermediate/',glacier,'FilledWx.csv']);%import weather data for selected glacier
        Weather_data.Date=datenum(Weather_data.Date);%convert dates to date number
        %get sites selected by user
        if length(unique(Glaciological_data.site_name))<36
            all_sites={};
            for i=1:length(unique(Glaciological_data.site_name))
                if get(getfield(handles, ['SITE',num2str(i)]),'Value')==1
                    all_sites=[all_sites;get(getfield(handles, ['SITE',num2str(i)]),'string')];
                end
            end
        elseif length(unique(Glaciological_data.site_name))>36 && get(handles.ALLSITES,'Value')==0
            all_sites=unique(Glaciological_data.site_name);
            [s,v] = listdlg('PromptString','Select Sites',...
                'SelectionMode','multiple',...
                'ListString',all_sites);
            all_sites=all_sites(s);
        elseif get(handles.ALLSITES,'Value')==1
            all_sites=unique(Glaciological_data.site_name);
        end
                
        disp('                                                 ')
        disp(['%%% You have selected to create data for all years for the following sites %%%'])
        disp(all_sites')
        
        disp('        %%%%%   Calculating Missing Glaciological Observations    %%%%%          ')
        Glaciological_data=[Glaciological_data table(zeros(height(Glaciological_data),1),zeros(height(Glaciological_data),1),'VariableNames',{'bw_fill', 'ba_fill'})];
        if exist(['data/',glacier,'/Input/Input_',glacier,'_Normalization_Site.csv'])%check if input file for normalization site exists. If it does then we can first fill using a mean balance profile
            normal_site=readtable(['data/',glacier,'/Input/Input_',glacier,'_Normalization_Site.csv']);
            plot_profiles = 1;
            Profile_Filled_Glaciological_data = Profile_Fill_Missing_Glaciological_Observations(glacier,Glaciological_data,years,normal_site,all_sites,AAD,plot_profiles);% fill any missing glaciological balance        
            for i=1:height(Profile_Filled_Glaciological_data)
                if isnan(Glaciological_data.bw(Glaciological_data.Year==Profile_Filled_Glaciological_data.Year(i)&strcmp(Glaciological_data.site_name,Profile_Filled_Glaciological_data.site_name(i))))
                    Glaciological_data.bw(Glaciological_data.Year==Profile_Filled_Glaciological_data.Year(i)&strcmp(Glaciological_data.site_name,Profile_Filled_Glaciological_data.site_name(i)))=Profile_Filled_Glaciological_data.bw(i);
                    Glaciological_data.bw_fill(Glaciological_data.Year==Profile_Filled_Glaciological_data.Year(i)&strcmp(Glaciological_data.site_name,Profile_Filled_Glaciological_data.site_name(i)))=Profile_Filled_Glaciological_data.bw_fill(i);
                elseif isnan(Glaciological_data.ba(Glaciological_data.Year==Profile_Filled_Glaciological_data.Year(i)&strcmp(Glaciological_data.site_name,Profile_Filled_Glaciological_data.site_name(i))))
                    Glaciological_data.ba(Glaciological_data.Year==Profile_Filled_Glaciological_data.Year(i)&strcmp(Glaciological_data.site_name,Profile_Filled_Glaciological_data.site_name(i)))=Profile_Filled_Glaciological_data.ba(i);
                    Glaciological_data.ba_fill(Glaciological_data.Year==Profile_Filled_Glaciological_data.Year(i)&strcmp(Glaciological_data.site_name,Profile_Filled_Glaciological_data.site_name(i)))=Profile_Filled_Glaciological_data.ba_fill(i);
                end
            end
        end
            
            Degree_Day_Factors_path = ['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Degree_Day_Factors.csv'];
            if ~exist(Degree_Day_Factors_path) %if file doesn't exist then we need approximate DDFs to start
                disp('ERROR: You need to Calibrate the mass balance model!')
            else
                meltrates = readtable(Degree_Day_Factors_path);%otherwise import previously calibrated DDFs
                ks = meltrates.ks(end);   %DDF for snow
                ki = meltrates.ki(end);   %DDF for ice
            end 
            %import calibrated precipitation ratios
            precipitation_ratios_path = ['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Precipitation_Ratios.csv'];%path for precipitation ratios
            if ~exist(precipitation_ratios_path)%if file doesn't exist then we assume prcipitation rates are one-to-one for initial value
                disp('ERROR: You need to Calibrate the mass balance model!')
            else
                precipitation_ratios_table = readtable((precipitation_ratios_path));
            end
            %if users selected to plot 'Find_Mass_Maximum_and_Minimum.m'
            %ask them if they really want that mank plots
        if get(handles.PLOTABLATIONMODEL,'Value')==1
           promptMessage = sprintf('You selected to plot the ablation model for more than one year. This probably blow up your machine. Do you still want to?');
            titleBarCaption = 'Plotting ablation model';
            choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know');
                    if strcmpi(choice, 'Yes')
                        plot_ablation_model = 1;
                        msgbox('Ok... Take cover!')
                    elseif strcmpi(choice, 'No')
                        plot_ablation_model = 0;
                    end
           
        else
                plot_ablation_model = 0;
        end
        if get(handles.PLOTINTEGRATION,'Value')==1
               promptMessage = sprintf('You selected to plot integrations for more than one year. This might take a while. Do you still want to plot integrations?');
                titleBarCaption = 'Plotting Integrations';
                choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know');
                        if strcmpi(choice, 'Yes')
                            plot_integration = 1;
                        elseif strcmpi(choice, 'No')
                            plot_integration = 0;
                        end

        else
                plot_integration = 0;
        end
        
        lapse_rates = get(handles.Lapse_Rate_Selection,'String');%list lapse rates available for user to select
        lapse_rate = str2num(cell2mat(lapse_rates(get(handles.Lapse_Rate_Selection,'Value'))));%lapse rate selected by user
        Weather_data = readtable(['data/',glacier,'/Intermediate/',glacier,'FilledWx.csv']);%import weather data for selected glacier
        Weather_data.Date = datenum(Weather_data.Date);%convert dates to date number
        
        

        if get(handles.incorporate_TSLs,'Value')==1
            if ~exist(['data/',glacier,'/Input/Input_',glacier,'_TSL_Data.csv'])
                disp('!!! You selected to incorporate transient snow line data but not input file exists !!!')
            else
                tsl_data = readtable(['data/',glacier,'/Input/Input_',glacier,'_TSL_Data.csv']);
                site_name = cell(height(tsl_data),1);
                for i = 1:height(tsl_data)
                    site_name(i,1) = {['TSL_',num2str(tsl_data.Elevation_m(i))]};
                    fall_dates(i,:) = datetime(['09/30/',datestr(tsl_data.Date_mmddyyyy(i),'yyyy')],'format','MM/dd/yyyy');
                end
                year = str2num(datestr(tsl_data.Date_mmddyyyy,'yyyy'));
                spring_dates = tsl_data.Date_mmddyyyy;
                elevation = tsl_data.Elevation_m;
                ba = nan*ones(height(tsl_data),1);
                bw = zeros(height(tsl_data),1);
                tsl_table = table(year,site_name,spring_dates,fall_dates,elevation,bw,ba,nan*ones(height(tsl_data),1),nan*ones(height(tsl_data),1),zeros(height(tsl_data),1),zeros(height(tsl_data),1),'VariableNames',Glaciological_data.Properties.VariableNames);
                Glaciological_data = [Glaciological_data;tsl_table];
                Glaciological_data = sortrows(Glaciological_data,'Year');    
                all_sites = [all_sites(~strcmp(all_sites,'  '));{'TSL'}];
            end
        end
        
        Filled_Glaciological_data = Model_Missing_Glaciological_Observations(glacier,Glaciological_data,years,all_sites,Weather_data,AAD,ks,ki,lapse_rate,precipitation_ratios_table,NaN,NaN,plot_ablation_model);% fill any missing glaciological balance        
        writetable(Filled_Glaciological_data,['data/',glacier,'/Intermediate/Filled_',glacier,'_Glaciological_Data.csv']);
        set(handles.UPDATE_BUTTON,'string','Ready to GO!');
        
        disp('        %%%%%    Missing Glaciological Observations Filled   %%%%%          ')


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   D) Calculated Glacier-wide Mass balance
function GO_BUTTON_Callback(hObject, eventdata, handles)
%       1) Users selects sites to use, years to proces, lapse rate to
%          apply, integration method(index or gradient), integration surface
%          (Conventional or Reference-Surface), time-system to use
%          (stratigraphic, Fixed-date Hydrological year, or Fixed-date
%          stratigraphic), and geodetic calibration (none, universal, or
%          piece-wise; See section D part 3
%       2) Calcukate glacier wide seasonal and annual mass balance time
%          series based on user defined time-system, hypsometry, etc...
%       3) Geodetically calibrate glacier wide time series. If selected,
%          include previously published solutions as well



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   D) Calculated Glacier-wide Mass balance
%%      1) Get user inputs
%           Users selects sites to use, years to proces, lapse rate to
%          apply, integration method(index or gradient), integration surface
%          (Conventional or Reference-Surface), time-system to use
%          (stratigraphic, Fixed-date Hydrological year, or Fixed-date
%          stratigraphic), and geodetic calibration (none, universal, or
%          piece-wise; See section D part 3
warning('off','all')
current_directory = pwd;%get(handles.MYDIR_EDIT,'string');
glaciers = get(handles.Glacier_Selection,'String');
glacier=cell2mat(glaciers(get(handles.Glacier_Selection,'value')));% glacier selected by user
if isempty(glacier) %if user did not enter a glacier name
    x = inputdlg('Enter glacier name',...
             'Sample', [1 50]);
glacier = str2num(x{:}); 
end
disp(['%%%%%%%%%%%%   Processing ',glacier,' Glacier Data   %%%%%%%%%%%%%%%%'])
integration_method = get(handles.GRAD_POPUP,'Value');                                                       %Integration method selected by user (Index, linear, gradient, piecwise linear, etc...)
Glaciological_data=readtable(['data/',glacier,'/Intermediate/Filled_',glacier,'_Glaciological_Data.csv']);  %import glaciological data for selected glacier
years = eval(get(handles.YEAR_EDIT,'String'));                                                              %years selected by user
lapse_rates=get(handles.Lapse_Rate_Selection,'String');
lapse_rate=str2num(cell2mat(lapse_rates(get(handles.Lapse_Rate_Selection,'Value'))));                       %lapse rate selected by user
AAD = importdata(['data/',glacier,'/Input/Input_',glacier,'_Area_Altitude_Distribution.csv']);              %import AAD for selected glacier
Weather_data=readtable(['data/',glacier,'/Intermediate/',glacier,'FilledWx.csv']);                          %weather data for selected glacier
Weather_data.Date=datenum(Weather_data.Date);                                                               %convert date to date numbers
Degree_Day_Factors_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Degree_Day_Factors.csv'];
if ~exist(Degree_Day_Factors_path)          %if file doesn't exist then we need to calibrate model
    disp('ERROR: You need to Calibrate the mass balance model!')
else
    meltrates=readtable(Degree_Day_Factors_path);%otherwise import previously calibrated DDFs
    ks=meltrates.ks(end);   %DDF for snow
    ki=meltrates.ki(end);   %DDF for ice
end 
precipitation_ratios_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Precipitation_Ratios.csv'];%path for precipitation ratios
if ~exist(precipitation_ratios_path)%if file doesn't exist then we need to calibrate model
    disp('ERROR: You need to Calibrate the mass balance model!')
else
    precipitation_ratios_table=readtable((precipitation_ratios_path));
end
dbstop if error
col1=[1,.5,0];
if length(unique(Glaciological_data.site_name))<=36 && get(handles.ALLSITES,'Value')==0
    all_sites={};
    site_names=unique(Glaciological_data.site_name);
    site_names=site_names(~contains(site_names,'TSL'))
    for i=1:length(site_names)
        if get(getfield(handles, ['SITE',num2str(i)]),'Value')==1
            all_sites=[all_sites;get(getfield(handles, ['SITE',num2str(i)]),'string')];
        end
    end
elseif length(unique(Glaciological_data.site_name))>36 && get(handles.ALLSITES,'Value')==0
        all_sites=unique(Glaciological_data.site_name);
        all_sites=all_sites(~contains(all_sites,'TSL'))
        [s,v] = listdlg('PromptString','Select Sites',...
            'SelectionMode','multiple',...
            'ListString',all_sites);
        all_sites=all_sites(s);
        if isempty(all_sites)
            disp(' You need to select sites!!')
        end
elseif get(handles.ALLSITES,'Value')==1
    all_sites=unique(Glaciological_data.site_name);
    all_sites=all_sites(~contains(all_sites,'TSL'));
end

        if get(handles.incorporate_TSLs,'Value')==1
            if ~contains(unique(Glaciological_data.site_name),'TSL')
                disp('%%%%%%%   TSLs observations are not available   %%%%%%%%%')
            end
            all_sites=[all_sites(~strcmp(all_sites,'  '));{'TSL'}];
        end

disp('                                                 ')
disp(['%%%%%%%%%% You Have Selected To Use The Following %%%%%%%%%%%'])
disp(all_sites')

my_yr = eval(get(handles.YEAR_EDIT,'String'));
numYrs=length(my_yr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Defining Time-system used  %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch get(handles.BAL_POPUP,'Value') 
    case 1 %Stratigraphic time-system used at each site. Site mass min to site mass min. Site balance is models using P/PDD model if missing
        time_system=1;
        y_axis_label='Stratigraphic Year (m w.e.)';
        time_system_text='Stratigraphic Mass Balance';
        time_system_name='Stratigraphic';
    case 2 %Stratigraphic time-system used at each site. Site mass min to site mass min. Site balance is models using P/PDD model if missing
        time_system=2;  
        y_axis_label='Stratigraphic Measurement Interval (m w.e.)';
        time_system_text='Stratigraphic Mass Balance Measurement Interval';
        time_system_name='Stratigraphic_Measurement_Interval';
    case 3 %glacier-wide from October 1 thru Sept 30. AKA the Hydrologic year
        time_system=3;
        y_axis_label='Hyrologic Year (m w.e.)';
        time_system_text= 'Hyrologic Year Mass Balance';
        time_system_name='Hyrologic'
    case 4 %Combined date system uses the the P/PDD model to find the weighted mean mass minimum date to calculate the glacier-wide balnce for one date
        time_system=4;
        y_axis_label='Combined Fixed-Date/Stratigraphic (m w.e.)';
        time_system_text= 'Fixed-Date Mass bBalance';
        time_system_name='Fixed_Date';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define Integration Surface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch get(handles.SURF_POPUP,'Value') 
    case 1 %use convention hypsometry == conventional
        integration_surface=1;
    case 2 %Use oldest hypsometry == reference
        integration_surface=2;
    case 3 % Use hypsometry of single year defined by user
        integration_surface(1,1)=3;
        x = inputdlg('Enter year for specifc hypsometry','Sample', [1 50]);     
        integration_surface(1,2) = str2num(x{:});
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select if integration of point balance over the glacier hypsometry should
% be plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.PLOTINTEGRATION,'Value')==1
%     close all
   if length(my_yr)>1 %warn user plotting integrations for every year might take a while
       promptMessage = sprintf('You selected to plot integrations for more than one year. This might take a while. Do you still want to plot integrations?');
        titleBarCaption = 'Plotting Integrations';
        choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know');
                if strcmpi(choice, 'Yes')
                    plot_integration = 1;
                elseif strcmpi(choice, 'No')
                    plot_integration=0;
                end
   else      
        plot_integration=1;
   end
else
        plot_integration=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select if 'Find_Mass_Maximum_and_Minimum_Adjustments.m' should be plotted
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if get(handles.PLOTABLATIONMODEL,'Value')==1
       % close all
       if length(my_yr)>1 %warn user plotting integrations for every year might take a while
       promptMessage = sprintf('You selected to plot the ablation model for more than one year. This probably blow up your machine. Do you still want to?');
        titleBarCaption = 'Plotting ablation model';
        choice = questdlg(promptMessage,titleBarCaption,'Yes','No','Don''t know');
                if strcmpi(choice, 'Yes')
                    plot_ablation_model = 1;
                    msgbox('Ok... Take cover!')
                elseif strcmpi(choice, 'No')
                    plot_ablation_model=0;
                end
   else      
        plot_ablation_model=1;
   end
else
        plot_ablation_model=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include Previous Glacierwide Solutions?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Include_Previous_Glacierwide_Solutions=get(handles.Include_Previous_Glacierwide_Solutions,'Value');


cd(current_directory)
addpath functions 
if integration_surface(1,1)==1 % for screen output
    surf2='CONVENTIONAL'; 
    Integration_surface_string='Conventional';
    marker='.-';
    Previous_solutions_marker='.--';
%     color=[round(rand,1) round(rand,1) round(rand,1)];
    color=[1 0 0];
elseif integration_surface(1,1)==2
    surf2='REFERENCE-SURFACE';
    Integration_surface_string='Reference_Surface';
    marker='.-';
    Previous_solutions_marker='.--';
    color=[round(rand,1) round(rand,1) round(rand,1)];
elseif integration_surface(1,1)==3
    surf2='SPECIFIC-HYPSOMETRY';
    Integration_surface_string='Specific_Hypsometry';
    marker='.-';
    Previous_solutions_marker='.--';
    color=[round(rand,1) round(rand,1) round(rand,1)];
end
if integration_method==1
    integration_method_string='Index_Method';
elseif integration_method==2
    integration_method_string='Gradient_Method';
elseif integration_method==3
    integration_method_string='Piecewise_Gradient_Method';
elseif integration_method==4
    integration_method_string='Non_Parametric_smoother_Method';
end

if time_system==1
    time_system_string='Stratigraphic';
elseif time_system==2
    time_system_string='Stratigraphic_Measurement';
elseif time_system==3
    time_system_string='Fixed_Date_Hydrologic';
elseif time_system==4
    time_system_string='Fixed_Date_Stratigraphic';
end
nan_incomplete_glaciological_data=get(handles.nan_incomplete_glaciological_data,'Value');
Geodetic_Calibration_index=get(handles.Geodetic_Calibration_Selection,'Value');    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   D) Calculated Glacier-wide Mass balance
%%      2) Glacier-wdie seasonal and annual balances are calculated using the inputs described in part one of section D
%           Calcukate glacier wide seasonal and annual mass balance time
%          series based on user defined time-system, hypsometry, etc...      

[final_point_balances,final_geodetic_data,Final_Glacier_Wide_solutions] = USGS_BenchmarkGlacier_Analysis(glacier,all_sites,years,Glaciological_data,Weather_data,AAD,ks,ki,precipitation_ratios_table,lapse_rate,time_system,integration_surface,integration_method,nan_incomplete_glaciological_data,Geodetic_Calibration_index,plot_integration,plot_ablation_model,Include_Previous_Glacierwide_Solutions)
writetable(final_point_balances,['data/',glacier,'/Output/','Output_',glacier,'_Adjusted_Point_Balances.csv']);

% % if user selected to plot ELAs against Ba values
if get(handles.BAELA,'Value')==1
    ela_ba_lm=fitlm(Final_Glacier_Wide_solutions.ELA_m,Final_Glacier_Wide_solutions.Ba_mwe);
    figure(1);hold on
    scatter(Final_Glacier_Wide_solutions.ELA_m,Final_Glacier_Wide_solutions.Ba_mwe,'k','filled')
    lsline
    text(min(Glacier_Wide_solutions_table.ELA_m)+100,-2,['r^2 = ',num2str(round(ela_ba_lm.Rsquared.Ordinary,2))]) 
    title([glacier,' Glacier ELA B_a Regression'],'fontname','arial ','fontsize',14,'fontweight','bold')
    set(gca,'fontname','arial ','fontsize',14,'fontweight','bold','TickLength',[0.025 0.025],'linewidth',2)
    axis square
    box on
    set(gcf, 'PaperPositionMode', 'auto');
    print -depsc2 gates_epoch2.eps 
    hold off
end    

writetable(Final_Glacier_Wide_solutions(2:end,:),['data/',glacier,'/Output/Output_',glacier,'_Glacier_Wide_solutions_calibrated.csv']);

writetable(final_geodetic_data,['data/',glacier,'/Output/Output_',glacier,'_Geodetics_Data.csv']);
axes(handles.CUMU_AXES)

uncalibrated_error=(0.2^2*abs(ones(height(Final_Glacier_Wide_solutions),1))).^0.5;

Calibrated_Cumulative_Ba=cumsum(Final_Glacier_Wide_solutions.Ba_mwe);
errorbar(Final_Glacier_Wide_solutions.Year,cumsum(Final_Glacier_Wide_solutions.Ba_mwe),uncalibrated_error,marker,'color',color);hold on%color); hold on
 if Geodetic_Calibration_index~=1
     first_geodetic_year_of_reanalyzed_solutions=str2num(datestr(final_geodetic_data.Date(find(double(strcmp(final_geodetic_data.Source,'DEM'))==1 & str2num(datestr(final_geodetic_data.Date,'yyyy'))>=Final_Glacier_Wide_solutions.Year(1),1,'first'),:),'yyyy'));
     final_geodetic_data=final_geodetic_data(str2num(datestr(final_geodetic_data.Date,'yyyy'))>=first_geodetic_year_of_reanalyzed_solutions,:);
     final_geodetic_data.Mass_Change_mwe=final_geodetic_data.Mass_Change_mwe-final_geodetic_data.Mass_Change_mwe(1);
     pin_reanalyzed_solutions=sum(Final_Glacier_Wide_solutions.Ba_mwe(Final_Glacier_Wide_solutions.Year<=first_geodetic_year_of_reanalyzed_solutions))-final_geodetic_data.Mass_Change_mwe(str2num(datestr(final_geodetic_data.Date,'yyyy'))==first_geodetic_year_of_reanalyzed_solutions);
     uncalibrated_Glacier_Wide_solutions=Final_Glacier_Wide_solutions.Ba_mwe-Final_Glacier_Wide_solutions.Calibration;
     uncalibrated_Glacier_Wide_solutions(1,1)=0;
     uncalibrated_Cumulative_Ba=cumsum(uncalibrated_Glacier_Wide_solutions);
     calibrated_error=(0.2^2*abs(Final_Glacier_Wide_solutions.Year-first_geodetic_year_of_reanalyzed_solutions*ones(height(Final_Glacier_Wide_solutions),1))).^0.5;
     errorbar(Final_Glacier_Wide_solutions.Year,uncalibrated_Cumulative_Ba,calibrated_error,marker,'color',[round(rand,1) round(rand,1) round(rand,1)]); hold on
 end
sources=unique(final_geodetic_data.Source);
final_geodetic_data.Mass_Change_mwe=final_geodetic_data.Mass_Change_mwe(1)-final_geodetic_data.Mass_Change_mwe;
if Geodetic_Calibration_index~=1
    for source=1:length(sources)
        data=final_geodetic_data(strcmp(final_geodetic_data.Source,sources(source)),:);
        if strcmp(sources(source),'DEM')
            marker='ks';
            data.Mass_Change_mwe=data.Mass_Change_mwe+data.Mass_Difference;
            data.Mass_Change_mwe=data.Mass_Change_mwe+data.Mass_Change_mwe(1);
        elseif strcmp(sources(source),'Altimetry')
            marker='ko';
        end
        data_start_Ba=Calibrated_Cumulative_Ba(Final_Glacier_Wide_solutions.Year==str2num(datestr(data.Date(1,:),'yyyy')))-data.Mass_Change_mwe(1);
        errorbar(str2num(datestr(data.Date,'yyyy')),data.Mass_Change_mwe+data_start_Ba,data.Uncertainty,marker,'markersize',8,'linewidth',1);        
    end
end
max_Ba=max(Calibrated_Cumulative_Ba);
min_Ba=min(Calibrated_Cumulative_Ba);
ylim([min_Ba-5 max_Ba+5])
set(handles.CUMU_AXES,'YMinorTick','on')
xlabel('Year')
ylabel('Cumulative balance (m w.e.)')
grid on
legend(Integration_surface_string,'Location','SW')
hold off;

axes(handles.ANN_NET_AXES)
r= rand(1,6);
bar(Final_Glacier_Wide_solutions.Year(2:end),Final_Glacier_Wide_solutions.Bw_mwe(2:end),1,'FaceColor', [r(1),r(2),1])
hold on;
bar(Final_Glacier_Wide_solutions.Year(2:end),Final_Glacier_Wide_solutions.Bs_mwe(2:end),1,'FaceColor', [1,r(3),r(4)])
hold on;
bar(Final_Glacier_Wide_solutions.Year(2:end),Final_Glacier_Wide_solutions.Ba_mwe(2:end),1,'FaceColor', [r(5),1,r(6)]) 
hold off; 
legend('B_w','B_s','B_a','Location','EO')
ylabel(y_axis_label);
set(handles.ANN_NET_AXES,'YMinorTick','on')
xlabel('Year')
xlim([Final_Glacier_Wide_solutions.Year(1)-.5 Final_Glacier_Wide_solutions.Year(end)+.5])
grid on



%
function varargout = mb_GUI_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;

function GRAD_POPUP_Callback(hObject, eventdata, handles)
%
function GRAD_POPUP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in Glacier_Select.
function Glacier_Select_Callback(hObject, eventdata, handles)
% hObject    handle to Glacier_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Glacier_Select contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Glacier_Select


% --- Executes during object creation, after setting all properties.
function Glacier_Select_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Glacier_Select (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%
function YEAR_EDIT_Callback(hObject, eventdata, handles)
%
function YEAR_EDIT_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% 
% 
%function will write a file that says which sites will use to compute
%balances
%
%
function SITE1_Callback(hObject, eventdata, handles)
function SITE2_Callback(hObject, eventdata, handles)
function SITE3_Callback(hObject, eventdata, handles)
function SITE4_Callback(hObject, eventdata, handles)
function SITE5_Callback(hObject, eventdata, handles)
function SITE6_Callback(hObject, eventdata, handles)
    
function SITE7_Callback(hObject, eventdata, handles)
function SITE8_Callback(hObject, eventdata, handles)
function SITE9_Callback(hObject, eventdata, handles)
function SITE10_Callback(hObject, eventdata, handles)
function SITE11_Callback(hObject, eventdata, handles)
function SITE12_Callback(hObject, eventdata, handles)
function SITE13_Callback(hObject, eventdata, handles)
function SITE14_Callback(hObject, eventdata, handles)
function ALLSITES_Callback(hObject, eventdata, handles)
    
%
%
%
function SURF_POPUP_Callback(hObject, eventdata, handles)
%
function SURF_POPUP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% 
function BAL_POPUP_Callback(hObject, eventdata, handles)
%
function BAL_POPUP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%
%
function ERRORBAR_CHECK_Callback(hObject, eventdata, handles)
%
%

function UPDATE_BUTTON_CreateFcn(hObject, eventdata, handles)
%
function UPDATE_POPUP_Callback(hObject, eventdata, handles)
%
function UPDATE_POPUP_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in BAELA.
function BAELA_Callback(hObject, eventdata, handles)
% hObject    handle to BAELA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BAELA


% --- Executes on button press in PLOTABLATIONMODEL.
function PLOTABLATIONMODEL_Callback(hObject, eventdata, handles)
% hObject    handle to PLOTABLATIONMODEL (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PLOTABLATIONMODEL



% --- Executes on selection change in Glacier_Selection.
function Glacier_Selection_Callback(hObject, eventdata, handles)
% hObject    handle to Glacier_Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Glacier_Selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Glacier_Selection


% --- Executes during object creation, after setting all properties.
function Glacier_Selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Glacier_Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Geodetic_Calibration_Selection.
function Geodetic_Calibration_Selection_Callback(hObject, eventdata, handles)
% hObject    handle to Geodetic_Calibration_Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Geodetic_Calibration_Selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Geodetic_Calibration_Selection


% --- Executes during object creation, after setting all properties.
function Geodetic_Calibration_Selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Geodetic_Calibration_Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in Lapse_Rate_Selection.
function Lapse_Rate_Selection_Callback(hObject, eventdata, handles)
% hObject    handle to Lapse_Rate_Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns Lapse_Rate_Selection contents as cell array
%        contents{get(hObject,'Value')} returns selected item from Lapse_Rate_Selection


% --- Executes during object creation, after setting all properties.
function Lapse_Rate_Selection_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Lapse_Rate_Selection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Include_Previous_Glacierwide_Solutions.
function Include_Previous_Glacierwide_Solutions_Callback(hObject, eventdata, handles)
% hObject    handle to Include_Previous_Glacierwide_Solutions (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of Include_Previous_Glacierwide_Solutions


% --- Executes on button press in nan_incomplete_glaciological_data.
function nan_incomplete_glaciological_data_Callback(hObject, eventdata, handles)
% hObject    handle to nan_incomplete_glaciological_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of nan_incomplete_glaciological_data


% --- Executes on button press in incorporate_TSLs.
function incorporate_TSLs_Callback(hObject, eventdata, handles)
% hObject    handle to incorporate_TSLs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of incorporate_TSLs


% --- Executes on button press in PLOTINTEGRATION.
function PLOTINTEGRATION_Callback(hObject, eventdata, handles)
% hObject    handle to PLOTINTEGRATION (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of PLOTINTEGRATION


% --- Executes on button press in SITE15.
function SITE15_Callback(hObject, eventdata, handles)
% hObject    handle to SITE15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE15


% --- Executes on button press in SITE16.
function SITE16_Callback(hObject, eventdata, handles)
% hObject    handle to SITE16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE16


% --- Executes on button press in SITE17.
function SITE17_Callback(hObject, eventdata, handles)
% hObject    handle to SITE17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE17


% --- Executes on button press in SITE18.
function SITE18_Callback(hObject, eventdata, handles)
% hObject    handle to SITE18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE18


% --- Executes on button press in SITE19.
function SITE19_Callback(hObject, eventdata, handles)
% hObject    handle to SITE19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE19


% --- Executes on button press in SITE20.
function SITE20_Callback(hObject, eventdata, handles)
% hObject    handle to SITE20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE20


% --- Executes on button press in SITE21.
function SITE21_Callback(hObject, eventdata, handles)
% hObject    handle to SITE21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE21


% --- Executes on button press in SITE22.
function SITE22_Callback(hObject, eventdata, handles)
% hObject    handle to SITE22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE22


% --- Executes on button press in SITE23.
function SITE23_Callback(hObject, eventdata, handles)
% hObject    handle to SITE23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE23


% --- Executes on button press in SITE24.
function SITE24_Callback(hObject, eventdata, handles)
% hObject    handle to SITE24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE24


% --- Executes on button press in SITE25.
function SITE25_Callback(hObject, eventdata, handles)
% hObject    handle to SITE25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE25


% --- Executes on button press in SITE26.
function SITE26_Callback(hObject, eventdata, handles)
% hObject    handle to SITE26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE26


% --- Executes on button press in SITE27.
function SITE27_Callback(hObject, eventdata, handles)
% hObject    handle to SITE27 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE27


% --- Executes on button press in SITE28.
function SITE28_Callback(hObject, eventdata, handles)
% hObject    handle to SITE28 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE28


% --- Executes on button press in SITE29.
function SITE29_Callback(hObject, eventdata, handles)
% hObject    handle to SITE29 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE29


% --- Executes on button press in SITE30.
function SITE30_Callback(hObject, eventdata, handles)
% hObject    handle to SITE30 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE30


% --- Executes on button press in SITE31.
function SITE31_Callback(hObject, eventdata, handles)
% hObject    handle to SITE31 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE31


% --- Executes on button press in SITE32.
function SITE32_Callback(hObject, eventdata, handles)
% hObject    handle to SITE32 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE32


% --- Executes on button press in SITE33.
function SITE33_Callback(hObject, eventdata, handles)
% hObject    handle to SITE33 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE33


% --- Executes on button press in SITE34.
function SITE34_Callback(hObject, eventdata, handles)
% hObject    handle to SITE34 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE34


% --- Executes on button press in SITE35.
function SITE35_Callback(hObject, eventdata, handles)
% hObject    handle to SITE35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE35


% --- Executes on button press in SITE36.
function SITE36_Callback(hObject, eventdata, handles)
% hObject    handle to SITE36 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of SITE36
