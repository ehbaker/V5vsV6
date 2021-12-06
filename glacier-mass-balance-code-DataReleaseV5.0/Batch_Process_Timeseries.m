clear all
close all
addpath functions
addpath data

%glaciers={'Taku';'SouthCascade';'Sperry'};
% glaciers={'LemonCreek'};
glaciers={'Taku'};
% glaciers={'Gulkana';'LemonCreek';'SouthCascade';'Sperry'};
master_dir=pwd;
cd(master_dir)


snow_table=table([],[],[],[],[],'VariableNames',{'Glacier','n','k','R2','MAE'});
ice_table=table([],[],[],[],[],'VariableNames',{'Glacier','n','k','R2','MAE'});
for glacier=1:length(glaciers)
clearvars -except snow_table ice_table glacier glaciers    
%set glacier name    
glacier_index=glacier;
glacier=cell2mat(glaciers(glacier));
plot_ablation_model=0;
plot_integration=0;
data_for_calibration='All';
%import input tables
Glaciological_data=readtable(['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv']);
AAD = importdata(['data/',glacier,'/Input/Input_',glacier,'_Area_Altitude_Distribution.csv']);
geodetic_data=readtable(['data/',glacier,'/Input/Input_',glacier,'_Geodetics.csv']);

%import default settings
settings=readtable(['data/',glacier,'/Input/Input_',glacier,'_defaultsettings.csv']);
index_sites=split(cell2mat(settings.Setting(strcmp(settings.Option,'index_sites'))),';');
nan_incomplete_glaciological_data=str2num(cell2mat(settings.Setting(strcmp(settings.Option,'nan_incomplete_glaciological_data'))));
Include_Previous_Glacierwide_Solutions=str2num(cell2mat(settings.Setting(strcmp(settings.Option,'Include_Previous_Glacierwide_Solutions'))));
years=str2num(cell2mat(settings.Setting(strcmp(settings.Option,'years'))));
break_point=str2num(cell2mat(settings.Setting(strcmp(settings.Option,'calibration_breakpoints'))));
incorporate_TSLs=str2num(cell2mat(settings.Setting(strcmp(settings.Option,'incorporate_TSLs'))));
normalization_site=split(cell2mat(settings.Setting(strcmp(settings.Option,'normalization_site'))),';')'
normalization_site=table(normalization_site(1),str2num(cell2mat(normalization_site(2))),str2num(cell2mat(normalization_site(3))),....
    'VariableNames',{'Normalization_Site' 'Year_1' 'Year_2'});

%set desired parameters
Geodetic_Calibration_index=4;
integration_surface=1;
integration_method=3;
time_system=4;
lapse_rate=-6.5;

integration_list={'Index';'LinearProfile';'PiecewiseProfile'};
calibration_list={'None';'BestFit';'Piecewise';'BreakPoints'};
timesystem_list={'stratigraphic';'measurement';'hydroyear';'floatingdate_stratigraphic'};
surface_list={'conventional';'reference';'specificdate'};

if integration_method==1
    all_sites=index_sites;
else
    all_sites=unique(Glaciological_data.site_name);
    if incorporate_TSLs==1
        all_sites=all_sites(~contains(all_sites,'TSL'));
        all_sites =[all_sites(~strcmp(all_sites,'  '));{'TSL'}];
    end
end

%% update weather data
[ready]=fillWxData(glacier); %fill any missing data within primary weather station data
Weather_data=readtable(['data/',glacier,'/Intermediate/',glacier,'FilledWx.csv']);
Weather_data.Date=datenum(Weather_data.Date);
%% update mass balance model calibration
plot_calibration=0;
% [~,~, precipitation_ratios_table,~,~,Degree_Day_Factor_table]= Calibrate_Precipitation_and_Ablation_models(glacier,[],Glaciological_data,Weather_data,AAD,lapse_rate,plot_ablation_model, plot_calibration);
[~,accum_data, precipitation_ratios_table,snow,ice,Degree_Day_Factor_table]= Calibrate_Precipitation_and_Ablation_models(glacier,[],Glaciological_data,Weather_data,AAD,lapse_rate,plot_ablation_model, plot_calibration);

lm=fitlm(snow(:,1),snow(:,2));
snow
snow_table=[snow_table;table(glaciers(glacier_index),length(snow),round(Degree_Day_Factor_table.ks(end),4),round(lm.Rsquared.Ordinary,2),round(nanmean(abs(lm.Residuals.Raw)),2),'VariableNames',{'Glacier','n','k','R2','MAE'})];
lm=fitlm(ice(:,1),ice(:,2));
ice_table=[ice_table;table(glaciers(glacier_index),length(ice),round(Degree_Day_Factor_table.ki(end),4),round(lm.Rsquared.Ordinary,2),round(nanmean(abs(lm.Residuals.Raw)),2),'VariableNames',{'Glacier','n','k','R2','MAE'})];


Degree_Day_Factors_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Degree_Day_Factors.csv'];
precipitation_ratios_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Precipitation_Ratios.csv'];%path for precipitation ratios
writetable(precipitation_ratios_table,precipitation_ratios_path)
writetable(Degree_Day_Factor_table,Degree_Day_Factors_path)
if ~exist(Degree_Day_Factors_path) %if file doesn't exist then we need approximate DDFs to start
    disp('ERROR: You need to Calibrate the mass balance model!')
else
    meltrates=readtable(Degree_Day_Factors_path);%otherwise import previously calibrated DDFs
    ks=meltrates.ks(end);   %DDF for snow
    ki=meltrates.ki(end);   %DDF for ice
end 
%% fill data gaps in glaciological data
Glaciological_data=[Glaciological_data table(zeros(height(Glaciological_data),1),zeros(height(Glaciological_data),1),'VariableNames',{'bw_fill', 'ba_fill'})];
if ~strcmp(normalization_site.Normalization_Site,'NaN')
    normal_site=normalization_site;
    plot_profiles = 1;
    Profile_Filled_Glaciological_data = Profile_Fill_Missing_Glaciological_Observations(glacier,Glaciological_data,years,normal_site,[],AAD,plot_profiles);% fill any missing glaciological balance        
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
if incorporate_TSLs==1
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
%         all_sites = [all_sites(~strcmp(all_sites,'  '));{'TSL'}];
    end
end
Filled_Glaciological_data = Model_Missing_Glaciological_Observations(glacier,Glaciological_data,[],{},Weather_data,AAD,ks,ki,lapse_rate,precipitation_ratios_table,NaN,NaN,plot_ablation_model);% fill any missing glaciological balance        
writetable(Filled_Glaciological_data,['data/',glacier,'/Intermediate/Filled_',glacier,'_Glaciological_Data.csv']);        
Glaciological_data=readtable(['data/',glacier,'/Intermediate/Filled_',glacier,'_Glaciological_Data.csv']);

%% Calculate Glacier time series
[final_point_balances,final_geodetic_data,Final_Glacier_Wide_solutions] = USGS_BenchmarkGlacier_Analysis(glacier,all_sites,years,Glaciological_data,Weather_data,AAD,ks,ki,precipitation_ratios_table,lapse_rate,time_system,integration_surface,integration_method,nan_incomplete_glaciological_data,Geodetic_Calibration_index,plot_integration,plot_ablation_model,Include_Previous_Glacierwide_Solutions);

writetable(Final_Glacier_Wide_solutions(2:end,:),['data/',glacier,'/Output/Output_',glacier,'_Glacier_Wide_solutions_calibrated.csv']);
writetable(final_point_balances,['data/',glacier,'/Output/','Output_',glacier,'_Adjusted_Point_Balances.csv']);
writetable(final_geodetic_data,['data/',glacier,'/Output/Output_',glacier,'_Geodetics_Data.csv']);




end