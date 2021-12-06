clear all
close all
addpath functions
dbstop if error

 glacier='Sperry';
glacier='Gulkana';
 glacier='Wolverine';
% glacier='SouthCascade';
% glacier='LemonCreek';


if strcmp(glacier,'Wolverine')
    index_sites={'A';'AU';'B';'C'};
elseif strcmp(glacier,'Gulkana')
    index_sites={'A';'AU';'B';'C';'D'};
end
plot_ablation_model=0;
plot_integration=0;
data_for_calibration='All';
bas=[];
bws=[];
bss=[];
all_Ba_solutions=[];
Glaciological_data=readtable(['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv']);
AAD = importdata(['data/',glacier,'/Input/Input_',glacier,'_Area_Altitude_Distribution.csv']);
if strcmp(glacier,'SouthCascade')
    index_sites=unique(Glaciological_data.site_name);
    nan_incomplete_glaciological_data=1;
    Include_Previous_Glacierwide_Solutions=1;
    years=unique(Glaciological_data.Year);
    incorporate_TSLs=0;
elseif strcmp(glacier,'Sperry')
    index_sites=unique(Glaciological_data.site_name);
    Include_Previous_Glacierwide_Solutions=0;
    nan_incomplete_glaciological_data=0;
    years=unique(Glaciological_data.Year);
    incorporate_TSLs=0;
elseif strcmp(glacier,'LemonCreek')
    index_sites={'A';'B';'C';'D';'E';'F'};
    nan_incomplete_glaciological_data=0;
    Include_Previous_Glacierwide_Solutions=0;
    years=1998:2018;
    incorporate_TSLs=1;
elseif strcmp(glacier,'Gulkana')
    index_sites={'A';'AU';'B';'C';'D'};
    nan_incomplete_glaciological_data=0;
    Include_Previous_Glacierwide_Solutions=0;
    years=unique(Glaciological_data.Year);
    incorporate_TSLs=0;
elseif strcmp(glacier,'Wolverine')
    index_sites={'A';'AU';'B';'C'};
    nan_incomplete_glaciological_data=0;
    Include_Previous_Glacierwide_Solutions=0;
    years=unique(Glaciological_data.Year);
    incorporate_TSLs=0;
end

integration_list={'Index';'LinearProfile';'PiecewiseProfile'};
calibration_list={'None';'BestFit';'Piecewise';'BreakPoints'};
timesystem_list={'stratigraphic';'measurement';'hydroyear';'floatingdate_stratigraphic'};
surface_list={'conventional';'reference';'specificdate'};
figure();hold on

for lapse_rate=-6.5%:0.5:-4.5
    Weather_data=readtable(['data/',glacier,'/Intermediate/',glacier,'FilledWx.csv']);
    Weather_data.Date=datenum(Weather_data.Date);
    plot_calibration=0;
    [~,~, precipitation_ratio_table,~,~,Degree_Day_Factor_table]= Calibrate_Precipitation_and_Ablation_models(glacier,years,Glaciological_data,Weather_data,AAD,lapse_rate,plot_ablation_model, plot_calibration);
    Degree_Day_Factors_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Degree_Day_Factors.csv'];
    precipitation_ratios_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Precipitation_Ratios.csv'];%path for precipitation ratios
    writetable(precipitation_ratio_table,precipitation_ratios_path)
    writetable(Degree_Day_Factor_table,Degree_Day_Factors_path)
    if ~exist(Degree_Day_Factors_path) %if file doesn't exist then we need approximate DDFs to start
        disp('ERROR: You need to Calibrate the mass balance model!')
    else
        meltrates=readtable(Degree_Day_Factors_path);%otherwise import previously calibrated DDFs
        ks=meltrates.ks(end);   %DDF for snow
        ki=meltrates.ki(end);   %DDF for ice
    end 
    precipitation_ratios_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Precipitation_Ratios.csv'];%path for precipitation ratios
    if ~exist(precipitation_ratios_path)%if file doesn't exist then we assume prcipitation rates are one-to-one for initial value
        disp('ERROR: You need to Calibrate the mass balance model!')
    else
        precipitation_ratios_table=readtable((precipitation_ratios_path));
    end
    all_sites={};
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
                all_sites = [all_sites(~strcmp(all_sites,'  '));{'TSL'}];
            end
        end
    
    
    Filled_Glaciological_data = Model_Missing_Glaciological_Observations(glacier,Glaciological_data,years,all_sites,Weather_data,AAD,ks,ki,lapse_rate,precipitation_ratios_table,NaN,NaN,plot_ablation_model);% fill any missing glaciological balance        
    writetable(Filled_Glaciological_data,['data/',glacier,'/Intermediate/Filled_',glacier,'_Glaciological_Data.csv']);        
    Glaciological_data=readtable(['data/',glacier,'/Intermediate/Filled_',glacier,'_Glaciological_Data.csv']);
    break_point=nan;
 
            
    if strcmp(glacier,'SouthCascade')
        nan_incomplete_glaciological_data=1;
    elseif strcmp(glacier,'Sperry')
        index_sites=unique(Glaciological_data.site_name);
    elseif strcmp(glacier,'LemonCreek')
        index_sites=unique(Glaciological_data.site_name);
    elseif strcmp(glacier,'Gulkana')
        index_sites={'A';'AU';'B';'C';'D'};
    elseif strcmp(glacier,'Wolverine')
        index_sites={'A';'AU';'B';'C'};
    end
    %Set geodetic list reference to only include None, BestFit for Sperry
    
    geodetic_data=readtable(['data/',glacier,'/Input/Input_',glacier,'_Geodetics.csv']);
    if height(geodetic_data)>2
        geodetic_list_ref=4;
    else
        geodetic_list_ref=2;
    end
    for Geodetic_Calibration_index=1%:geodetic_list_ref
        for integration_surface=1%:2
            for integration_method=1:3
                if integration_method==1
                    all_sites=index_sites;
                else
                    all_sites=unique(Glaciological_data.site_name);
                    all_sites=all_sites(~contains(all_sites,'TSL'));
                    all_sites =[all_sites(~strcmp(all_sites,'  '));{'TSL'}];
                end
                for time_system=4
                      [final_point_balances,final_geodetic_data,Final_Glacier_Wide_solutions] = USGS_BenchmarkGlacier_Analysis(glacier,all_sites,years,Glaciological_data,Weather_data,AAD,ks,ki,precipitation_ratios_table,lapse_rate,time_system,integration_surface,integration_method,nan_incomplete_glaciological_data,Geodetic_Calibration_index,plot_integration,plot_ablation_model,Include_Previous_Glacierwide_Solutions);
%                         Final_Glacier_Wide_solutions.Bs_mwe=Final_Glacier_Wide_solutions.Ba_mwe-Final_Glacier_Wide_solutions.Bw_mwe;
                end
            

                if isempty(all_Ba_solutions)
                    all_Ba_solutions=table(Final_Glacier_Wide_solutions.Year,Final_Glacier_Wide_solutions.Ba_mwe,'VariableNames',{'Year' [cell2mat(integration_list(integration_method)),cell2mat(calibration_list(Geodetic_Calibration_index))]});
                else
                    all_Ba_solutions=[all_Ba_solutions table(Final_Glacier_Wide_solutions.Ba_mwe,'VariableNames',{[cell2mat(integration_list(integration_method)),cell2mat(calibration_list(Geodetic_Calibration_index))]})];
                end
            end
            
        end
        
    end
    writetable(all_Ba_solutions,['data/',glacier,'/Output/',glacier,'_all_solutions.csv']);
end
%%
for i=2:height(all_Ba_solutions)
    if sum(years==all_Ba_solutions.Year(i))
        stdevs(i-1,1)=iqr(table2array(all_Ba_solutions(i,2:end)));
    else
        stdevs(i-1,1)=nan;
    end
end
sensitivity=nanmean(stdevs)
    