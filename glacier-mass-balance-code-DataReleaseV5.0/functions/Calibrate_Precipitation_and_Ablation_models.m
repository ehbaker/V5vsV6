function [ready, all_precipitation_data, precipitation_ratio_table, all_snow_ablation_data, all_ice_ablation_data,Degree_Day_Factor_table]= Calibrate_Precipitation_and_Ablation_models(glacier,years,Glaciological_data,Weather_data,AAD,lapse_rate,plot_ablation_model, plot_calibration) 
%% Calibrate_Precipitation_and_Ablation_models.m
% This function calibrates positive degree day and precipiation coefficents 
% for a temperature index and precipiation modeling

%% Inputs:

%glacier - glacier name

%years - range of years of interest

%Glaciological_data - table of seasonal and annual glaciological data

%Weather_data - table of daily mean temperatures (*C) and total daily
%precipiation (m)

%AAD - Area Altitude Distribution: sq km of glacier surface area within
%regularly speced elevation bins (m)

%lapse_rate - provided lapse rate for modeling temperatures

%plot_ablation_model - switch for plotting ablation model (1) or not
%plotting (0)

%% Outputs: 

% ready - binary (0 or 1) flag stating that the finction is complete or
% failed 

%% general outline

%1) Use observed accumulation (swe) on the glacier at mass balance locations, modeled site
%temperature, and precipitation measured at a nearby weather station to
%calibrate precipitation ratios, aka
%site_precipitation/weather_station_precipitation

%2) Use measured ablation of snow and ice (seperately) to resolve inital
%Degree Day coefficents 

%% Begin function
dbstop if error
addpath functions

disp('%%%%%%%%%%%%%%%Inverting Melt and Precip Ratios%%%%%%%%%%%%%%%%%')
warning off MATLAB:rankDeficientMatrix
warning off MATLAB:ezplotfeval:NotVectorized

if isempty(years)
    years=unique(Glaciological_data.Year);
end
%% Convert dates of all glaciological observations from calender dates to datenumbers
dates=[];
for i=1:height(Glaciological_data)
    % spring date not allowed to be NaN but fall date is
    if strcmp(Glaciological_data.spring_date(i),'NaN')||strcmp(Glaciological_data.spring_date(i),'Nan')||strcmp(Glaciological_data.spring_date(i),'nan')
        dates=[dates;datenum(['3/1/',num2str(Glaciological_data.Year(i))])  NaN];
    elseif strcmp(Glaciological_data.fall_date(i),'NaN')||strcmp(Glaciological_data.fall_date(i),'Nan')||strcmp(Glaciological_data.fall_date(i),'nan')
        dates=[dates;datenum(Glaciological_data.spring_date(i))  NaN];
    else  
        dates=[dates;datenum(Glaciological_data.spring_date(i))  datenum(Glaciological_data.fall_date(i))];
    end
end
Glaciological_data.spring_date=dates(:,1);    
Glaciological_data.fall_date=dates(:,2); 
if exist(['data/',glacier,'/Input/Input_',glacier,'_SubSeasonal_Glaciological_Data.csv'])
    additional_data=readtable(['data/',glacier,'/Input/Input_',glacier,'_SubSeasonal_Glaciological_Data.csv']);
else
    additional_data=[];
end
    

all_site_names=unique(Glaciological_data.site_name);%get the name of every site ever measured on glacier selected

for i=1:length(all_site_names)%for each site create variables that will be populated as calibrations are derived
    Precipitation_data_for_all_sites{i}=[];
    site_snow_ablation_data{i}=[];
    site_ice_ablation_data{i}=[];
end


%% Set default degree day factors(DDFs) for glacier selected
% Degree_Day_Factors_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Degree_Day_Factors.csv'];
% if ~exist(Degree_Day_Factors_path) %if file doesn't exist then we need approximate DDFs to start
    ks=-0.005;              %approximate DDFs for snow and ice melt in m w.e./1^*C > 0*C
    ki=-0.005;


%% import precipitation ratios for glacier selected
% precipitation_ratios_path=['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Precipitation_Ratios.csv'];%path for precipitation ratios
all_site_names=unique(Glaciological_data.site_name);
precipitation_ratios=ones(length(all_site_names),1);
precipitation_rsq=ones(length(all_site_names),1);
precipitation_rmse=ones(length(all_site_names),1);
precipitation_ratios_table=table(all_site_names,precipitation_ratios,precipitation_rsq,precipitation_rmse,'VariableNames',{'site_name','precipitation_ratios','rsquared','RMSE'});


for i=1:length(all_site_names)
    site_elevation(i,1)=nanmean(Glaciological_data.elevation(strcmp(Glaciological_data.site_name,all_site_names(i))));
end
%% set time-system to stratigraphic to match observations
% all observations are in a floating date stratigraphic time-system
% Hence, we need to compare precipitation at wx stations and mass balance
% sites in during this time framce
time_system=1;

%% Calibrate site precipitation models
% this is done first since we need to make the least assumptions in order
% to get the data needed to calibrate

    for i=1:length(years)
        insitu_data=Glaciological_data(Glaciological_data.Year==years(i) & ~isnan(Glaciological_data.spring_date) & ~isnan(Glaciological_data.fall_date) & ~isnan(Glaciological_data.bw),:); %& ~isnan(Glaciological_data.ba),:);
        %We need to get the dates for eavery measured accumulation event. 
        %These accumulation events are between the mass minimum and the
        %measured winter balance at a location and any measured summer
        %accumulation
        if isempty(insitu_data)%there are no observations from this year
        else %there are observations from this year

            %construct dates for previous year to find mass minimum
            previous_fall_observation_dates=datenum(['9/30/',num2str(years(i)-1)],'mm/dd/yyyy').*ones(height(insitu_data),1);
            previous_spring_observation_dates=datenum(['5/1/',num2str(years(i)-1)],'mm/dd/yyyy').*ones(height(insitu_data),1);
            previous_year=(years(i)-1).*ones(height(insitu_data),1);
            insitu_data.fall_date=previous_fall_observation_dates;
            insitu_data.spring_date=previous_spring_observation_dates;
            insitu_data.Year=previous_year;
            integration_method=nan;
            integration_surface=nan;
            %get mass minimum date of previous year so we know what time-frame
            %to pull precipitation measured at the weather station to match the
            %stratigraphic winter balance measurements made in spring
            %since we are only getting the date of the mass minimum we onlly
            %need to change the dates
%             
            [~,~,~,previous_stratigraphic_mass_minimum_date] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);

            %now get mass minimum dates for each site of current year in case there was any
            %summer accumulation measured during fall trip
            insitu_data=Glaciological_data(Glaciological_data.Year==years(i) & ~isnan(Glaciological_data.spring_date) & ~isnan(Glaciological_data.fall_date) & ~isnan(Glaciological_data.bw),:);%& ~isnan(Glaciological_data.ba)
            [~,mass_maximum_date_numbers,~,stratigraphic_mass_minimum_date] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);

            %calculate precipitation measured at the gauge and compile it with
            %accumulation measured at the stake
            for site=1:height(insitu_data) %loop through each site
                if ~isnan(insitu_data.bw(site)) && (mass_maximum_date_numbers(site)>=insitu_data.spring_date(site)-30) %if there was a observation made in the spring
                    %using a precipitation ratio of one allows use to get the
                    %precipitation measured at the weather station and not
                    %the site... tricksy programmer

                    [Site_Weather]=Model_Site_Weather(Weather_data,insitu_data.elevation(site),previous_stratigraphic_mass_minimum_date(site),insitu_data.spring_date(site),glacier,lapse_rate,1);
                    if (sum(Site_Weather.T_flag)/height(Site_Weather))<.5 && (sum(Site_Weather.P_flag)/height(Site_Weather))<.5
                        Precipitation_data_for_all_sites{strcmp(insitu_data.site_name(site),all_site_names)}=[Precipitation_data_for_all_sites{strcmp(insitu_data.site_name(site),all_site_names)}; nansum(Site_Weather.Precipitation),insitu_data.bw(site) .5 insitu_data.Year(site)];
                    end
                    
                   
                end
                if insitu_data.summer_accumulation(site)>0 %if there was any summer accumulation measured during fall observations
                    if insitu_data.fall_date(site)<stratigraphic_mass_minimum_date(site)
                        disp(cell2mat(['Warning: Model predicting mass minimum data after fall measurements when new snow was recorded at site. Indicating the mass minimum date predicted is wrong in ',num2str(insitu_data.Year(site)),' at site ',insitu_data.site_name(site)]))
                    else
                    [Site_Weather]=Model_Site_Weather(Weather_data,insitu_data.elevation(site),stratigraphic_mass_minimum_date(site),insitu_data.fall_date(site),glacier,lapse_rate,1);
                        if (sum(Site_Weather.T_flag)/height(Site_Weather))<.5 && (sum(Site_Weather.P_flag)/height(Site_Weather))<.5
                            Precipitation_data_for_all_sites{strcmp(insitu_data.site_name(site),all_site_names)}=[Precipitation_data_for_all_sites{strcmp(insitu_data.site_name(site),all_site_names)}; nansum(Site_Weather.Precipitation),insitu_data.summer_accumulation(site) .5 insitu_data.Year(site)];
                        end
                    end
                end
            end        
        end
        if ~isempty(additional_data)
            bonus_insitu_data=additional_data(additional_data.Year==years(i) & additional_data.db_mwe>=0,:);
%             for i=1:height(bonus_insitu_data)
            
        end

    end
   
    %find maximum amount of precipitation for axis limit in calibration plot
    %below
    max_precipitation=0;
    for i=1:length(Precipitation_data_for_all_sites)
        if ~isempty(Precipitation_data_for_all_sites{i})
            if max(Precipitation_data_for_all_sites{i}(:,2)) > max_precipitation
                max_precipitation=max(Precipitation_data_for_all_sites{i}(:,2))+1;
            end
        end
    end

    % plot precipitation measured at the weather station to that measured on
    % the glacier
    all_precipitation_data=[];
    
    for i=1:length(Precipitation_data_for_all_sites)%for each site
        if isempty(Precipitation_data_for_all_sites{i}) %than no data during selected time-frame
            precipitation_ratios(i,1)=NaN; 
            precipitation_rsq(i,1)=NaN;
            precipitation_rmse(i,1)=NaN;
        else
            precipitation_ratios(i,1)=Precipitation_data_for_all_sites{i}(:,1)\Precipitation_data_for_all_sites{i}(:,2);
            [r,p] = corr(Precipitation_data_for_all_sites{i}(:,1),Precipitation_data_for_all_sites{i}(:,2),'type','Kendall');
            precipitation_rsq(i,1)=r;
            precipitation_rmse(i,1)=nanmean(abs((Precipitation_data_for_all_sites{i}(:,1).*precipitation_ratios(i,1))-Precipitation_data_for_all_sites{i}(:,2)));
        end


        if isempty(Precipitation_data_for_all_sites{i}) %then we had no data during the selected time period

        else %otherwise plot the data we do have
            Precipitation_data_for_all_sites{i}(:,5)=Precipitation_data_for_all_sites{i}(:,1)-Precipitation_data_for_all_sites{i}(:,2);
            all_precipitation_data=[all_precipitation_data;Precipitation_data_for_all_sites{i}(:,1) Precipitation_data_for_all_sites{i}(:,2)];
            if plot_calibration==1
                figure(1);hold on
                subplot(ceil(sqrt(length(all_site_names))),ceil(sqrt(length(all_site_names))),i)
                scatter(Precipitation_data_for_all_sites{i}(:,1),Precipitation_data_for_all_sites{i}(:,2),'filled');hold on
                y_est = Precipitation_data_for_all_sites{i}(:,1)*precipitation_ratios(i,1);
                plot(Precipitation_data_for_all_sites{i}(:,1),y_est,'-r')
                title(cell2mat(all_site_names(i)))
                set(gca, 'XColor', 'k', 'YColor', 'k','XAxisLocation','bottom','Xtick',0:300:1200,...
                'YAxisLocation','left','Ytick',-9:1:0,'fontname','arial ','fontsize',12,'TickLength',[0.025 0.025],'linewidth',2);
                xlim([0 max_precipitation])
                ylim([0 max_precipitation])
                xlabel('Weather Station (m)','fontweight','bold')
                ylabel('Site (m w.e.)','fontweight','bold')
                box on
                axis square
                set(gcf, 'PaperPositionMode', 'auto');
                print -depsc2 gates_epoch2.eps
            end
        end
    end
    for site=1:length(all_site_names)
        if isnan(precipitation_ratios(site))||isinf(precipitation_ratios(site))
            elevation_difference=abs(site_elevation(:,1)-site_elevation(site,1));
            [~,index]=min(elevation_difference(elevation_difference>0));
            precipitation_ratios(site,1)=precipitation_ratios(index,1);
        end
    end
    precipitation_ratios=[precipitation_ratios; all_precipitation_data(:,1)\all_precipitation_data(:,2)];

    [r,p] = corr(all_precipitation_data(:,1),all_precipitation_data(:,2),'type','Kendall');
    precipitation_rsq=[precipitation_rsq; r];
    precipitation_rmse=[precipitation_rmse; nanmean(abs((all_precipitation_data(:,1).*precipitation_ratios(end,1))-all_precipitation_data(:,2)))];
    
    %place all calibrated precipitation ratios and fit statistics into table
    %and save it
    precipitation_ratio_table=table([all_site_names;'All_sites'],[site_elevation; nan],precipitation_ratios,precipitation_rsq,precipitation_rmse,'VariableNames',{'site_name','elevation','precipitation_ratios','rsquared','RMSE'})
    precipitation_ratio_table=precipitation_ratio_table(~isnan(precipitation_ratio_table.precipitation_ratios),:);
%     writetable(precipitation_ratio_table,precipitation_ratios_path)
    
    %% Calibrate Degree Day Model For Snow 
    %since we have little ablation data solely measuring ice ablation, or
    % measurements which do not include the end of the accumulation season we
    % need to loop back through to use our calibrated precipitation ratios to
    % fill in unknown accumulation (difference from spring measurement to mass
    % maximum) and calibrate using the summer balance from the observed bw
    % value corrected to the mass maximum to either the observed ba or the mass
    % minimum date (which ever comes first)
difference_in_ks=[];
difference_in_ki=[];
    % precipratio=precipitation_ratio_table.precipitation_ratios; %start using calibrated precipitation ratios
while isempty(difference_in_ks)||sum(difference_in_ks>0.00001)||isempty(difference_in_ki)||sum(difference_in_ki>0.00001)
all_snow_ablation_data=[];
% clearvars site_snow_ablation_data all_snow_ablation_data
for i=1:length(all_site_names)%for each site create variables that will be populated as calibrations are derived
    site_snow_ablation_data{i}=[];
    site_ice_ablation_data{i}=[];
end
if plot_calibration==1 && ~isempty(difference_in_ks) && ~sum(difference_in_ks>0.00001)
    figure(1);clf
    figure(2);clf
    figure(3);clf
    figure(4);clf
    figure(5);clf
end

        for i=1:length(years)
        insitu_data=Glaciological_data(Glaciological_data.Year==years(i) & ~isnan(Glaciological_data.bw) & ~isnan(Glaciological_data.fall_date) & ~isnan(Glaciological_data.ba),:); %subset glaciological data for year

        %We need to get the dates for eavery measured ablation measurement. 
        %These will typically be between the stratigraphic mass maximum and the
        %measurement date in the fall
        if ~isempty(insitu_data)
        [stratigraphic_mass_maximum_correction,stratigraphic_mass_maximum_date,~,stratigraphic_mass_minimum_date] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);

        %loop through each site to calculate the Positive Degree Days for each
        %time-interval
            for site=1:height(insitu_data)
                if ~isnan(insitu_data.ba(site)) && ~isnan(insitu_data.bw(site)) && insitu_data.ba(site)>0 && (insitu_data.spring_date(site)>=stratigraphic_mass_maximum_date(site)) %make sure both bw and ba were measured and confirm melting only snow during measurement interval
                    if stratigraphic_mass_maximum_date(site) <= insitu_data.spring_date(site)       %if mass maximum date was before spring observation 
                        time_and_swe_1=[insitu_data.spring_date(site) insitu_data.bw(site)];        %use spring observation date as time one and apply no correction
                    elseif stratigraphic_mass_maximum_date(site) > insitu_data.spring_date(site)    %if mass maximum date was after spring observation 
                        time_and_swe_1=[stratigraphic_mass_maximum_date(site) (insitu_data.bw(site)+stratigraphic_mass_maximum_correction(site))]; %we need to correct the measurement to the mass maximum date                 
                    end
                    if stratigraphic_mass_minimum_date(site) >= insitu_data.fall_date(site)         %if mass minimum date was after fall observation 
                        time_and_swe_2=[insitu_data.fall_date(site) insitu_data.ba(site)];          %than use fall observation date
                    elseif stratigraphic_mass_minimum_date(site) < insitu_data.fall_date(site)      %if mass minimum date was before fall observation
                        time_and_swe_2=[stratigraphic_mass_minimum_date(site) insitu_data.ba(site)]; %use the mass minimum date                  
                    end
                    bs=time_and_swe_2(1,2)-time_and_swe_1(1,2);                                      %calculate the summer balance for site
                    precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,insitu_data.site_name(site)));
                    if isempty(precipitation_ratio)
                        precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,'All_sites'));
                    end
                    [Site_Weather]=Model_Site_Weather(Weather_data,insitu_data.elevation(site),time_and_swe_1(1,1),time_and_swe_2(1,1),glacier,lapse_rate,precipitation_ratio); %model site temperatures
                    if (sum(Site_Weather.T_flag)/height(Site_Weather))<.5 && (sum(Site_Weather.P_flag)/height(Site_Weather))<.5

                                accumulation=0;
                                Positive_Degree_Days=0;
                                
                                for day=1:height(Site_Weather)
                                   accumulation=accumulation+Site_Weather.Precipitation(day,1);
                                   if Site_Weather.Temperature(day,1)>0
                                        ablation=Site_Weather.Temperature(day,1).*ks;
                                   else
                                       ablation=0;
                                   end

                                   accumulation=accumulation+ablation;
                                   plot(day,accumulation)
                                       if accumulation<=0 && Site_Weather.Temperature(day,1)>0
                                            accumulation=0;
                                            Positive_Degree_Days=Positive_Degree_Days+Site_Weather.Temperature(day,1); 
                                       end
                                end
                     %get positive degree days (T > 0*C)
                    site_snow_ablation_data{strcmp(insitu_data.site_name(site),all_site_names)}=[site_snow_ablation_data{strcmp(insitu_data.site_name(site),all_site_names)}; Positive_Degree_Days bs time_and_swe_1(1,1) time_and_swe_2(1,1) insitu_data.Year(1)];
                
                    end
                end

                if ~isnan(insitu_data.winter_ablation(site)) && insitu_data.winter_ablation(site) < 0 %if any winter ablation was measured at site
                    previous_year_data_index=find(Glaciological_data.Year==years(i)-1);             %then we need to look at the previous years temps
                    if ~isempty(previous_year_data_index)                                           %if there is data/a measurements date from previous year
                        site_index=strcmp(Glaciological_data.site_name(previous_year_data_index),insitu_data.site_name(site)); %find data for site
                        site_index_for_previous_year=previous_year_data_index(site_index);
                        if ~isempty(site_index_for_previous_year)                                   %if we have data from previous year at the site we can get the PDDs for the measured winter ablation
                            if Glaciological_data.ba(site_index_for_previous_year)>0 && Glaciological_data.ba(site_index_for_previous_year) > abs(insitu_data.winter_ablation(site))
                                %get previous years mass minumum date
                                [~,~,~,previous_year_stratigraphic_mass_minimum_date] = Find_Mass_Maximum_and_Minimum_Adjustments(Glaciological_data(site_index_for_previous_year,:),Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);
                                %get the PDDs between observation and mass
                                %minimum date
                                precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,insitu_data.site_name(site))); 
                                [Site_Weather]=Model_Site_Weather(Weather_data,insitu_data.elevation(site),Glaciological_data.fall_date(site_index_for_previous_year),previous_year_stratigraphic_mass_minimum_date,glacier,lapse_rate,precipitation_ratio);
                                if (sum(Site_Weather.T_flag)/height(Site_Weather))<.5 && (sum(Site_Weather.P_flag)/height(Site_Weather))<.5
                                accumulation=0;
                                Positive_Degree_Days=0;
                                
                                for day=1:height(Site_Weather)
                                   accumulation=accumulation+Site_Weather.Precipitation(day,1);
                                   if Site_Weather.Temperature(day,1)>0
                                        ablation=Site_Weather.Temperature(day,1).*ks;
                                   else
                                       ablation=0;
                                   end
                                   accumulation=accumulation+ablation;
                                       if accumulation<=0 && Site_Weather.Temperature(day,1)>0
                                            accumulation=0;
                                            Positive_Degree_Days=Positive_Degree_Days+Site_Weather.Temperature(day,1); 
                                       end
                                end
                                site_snow_ablation_data{strcmp(insitu_data.site_name(site),all_site_names)}=[site_snow_ablation_data{strcmp(insitu_data.site_name(site),all_site_names)}; Positive_Degree_Days insitu_data.winter_ablation(site) Glaciological_data.fall_date(site_index_for_previous_year) previous_year_stratigraphic_mass_minimum_date insitu_data.Year(1)];
                                end
                            end
                        end
                    end
                end 
                
            end
        end
        if ~isempty(additional_data)
            bonus_insitu_data=additional_data(additional_data.Year==years(i) & additional_data.db_mwe<=0 & strcmp(additional_data.Surface1,'snow') & strcmp(additional_data.Surface2,'snow'),:);
             for i=1:height(bonus_insitu_data)
                 precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,bonus_insitu_data.Site_Name(i)));
                    if isempty(precipitation_ratio)
                        precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,'All_sites'));
                    end
                    time_1=datenum(bonus_insitu_data.Date1(i));
                    time_2=datenum(bonus_insitu_data.Date2(i));
                    [Site_Weather]=Model_Site_Weather(Weather_data,bonus_insitu_data.Elevation_m(i),time_1,time_2,glacier,lapse_rate,precipitation_ratio); %model site temperatures
                    if (sum(Site_Weather.T_flag)/height(Site_Weather))<.5 && (sum(Site_Weather.P_flag)/height(Site_Weather))<.5

                                accumulation=0;
                                Positive_Degree_Days=0;
                                
                                for day=1:height(Site_Weather)
                                   accumulation=accumulation+Site_Weather.Precipitation(day,1);
                                   if Site_Weather.Temperature(day,1)>0
                                        ablation=Site_Weather.Temperature(day,1).*ks;
                                   else
                                       ablation=0;
                                   end

                                   accumulation=accumulation+ablation;
                                   plot(day,accumulation)
                                       if accumulation<=0 && Site_Weather.Temperature(day,1)>0
    %                                         degree_day
                                            accumulation=0;
                                            Positive_Degree_Days=Positive_Degree_Days+Site_Weather.Temperature(day,1); 
                                       end
                                end
                        %get positive degree days (T > 0*C)
                        if sum(strcmp(precipitation_ratio_table.site_name,bonus_insitu_data.Site_Name(i)))==0
                            all_snow_ablation_data=[all_snow_ablation_data; Positive_Degree_Days bonus_insitu_data.db_mwe(i) time_1 time_2 bonus_insitu_data.Year(1)];
                        else
                            site_snow_ablation_data{strcmp(bonus_insitu_data.Site_Name(i),all_site_names)}=[site_snow_ablation_data{strcmp(bonus_insitu_data.Site_Name(i),all_site_names)};Positive_Degree_Days bonus_insitu_data.db_mwe(i) time_1 time_2 bonus_insitu_data.Year(1)];
                        end
                    end
                end
        end
        end

    %find max snow ablation for axis limit in plot below
    max_melt=0;
    for i=1:length(site_snow_ablation_data)
        if ~isempty(site_snow_ablation_data{i})
            if min(site_snow_ablation_data{i}(:,2)) < max_melt
                max_melt=min(site_snow_ablation_data{i}(:,2))-1;
            end
        end
    end

    
    
    for i=1:length(site_snow_ablation_data)
        if isempty(site_snow_ablation_data{i}) %We have no data for this site dring selected time-period

        elseif length(site_snow_ablation_data{i}(:,2))>=3%if we have more than two observations than we should calibrate using robust option to minimize the influence of outliners
%             site_snow_ablation_data{i}=[site_snow_ablation_data{i};0 0 NaN NaN 1];
        elseif length(site_snow_ablation_data{i}(:,2))==2 %if we only have two observations than we can only perform a standard least squares minimization
%             site_snow_ablation_data{i}=[site_snow_ablation_data{i};0 0 NaN NaN 1];
        else%if we only have one observation than we have to use the ratio of precipitation measured at the site verse that at the weather station
%             site_snow_ablation_data{i}=[site_snow_ablation_data{i};0 0 NaN NaN 1];
        end

        if isempty(site_snow_ablation_data{i}) %no data no DDFs and no stats
            Snow_Degree_Day_factor(i,1)=nan; 
            Snow_Degree_Day_factor_rsq(i,1)=nan;
            Snow_Degree_Day_factor_rmse(i,1)=nan;

        else %otherwise we have calibrated DDFs for snow for this site 
            all_snow_ablation_data=[all_snow_ablation_data;site_snow_ablation_data{i}];
            Snow_Degree_Day_factor(i,1)=site_snow_ablation_data{i}(:,1)\site_snow_ablation_data{i}(:,2);
            [r,p] = corr(site_snow_ablation_data{i}(:,1),site_snow_ablation_data{i}(:,2),'type','Kendall');
            Snow_Degree_Day_factor_rsq(i,1)=r;
            Snow_Degree_Day_factor_rmse(i,1)=nanmean(abs((site_snow_ablation_data{i}(:,1).*Snow_Degree_Day_factor(i,1))-site_snow_ablation_data{i}(:,2)));

        end

        if isempty(site_snow_ablation_data{i})
        else
            if plot_calibration==1 && ~isempty(difference_in_ks) && ~sum(difference_in_ks>0.00001)
            figure(2);hold on
            subplot(ceil(sqrt(length(all_site_names))),ceil(sqrt(length(all_site_names))),i)
            scatter(site_snow_ablation_data{i}(:,1),site_snow_ablation_data{i}(:,2),'filled');hold on
            y_est = site_snow_ablation_data{i}(:,1)*Snow_Degree_Day_factor(i,1);
             plot(site_snow_ablation_data{i}(:,1),y_est,'-r') 
            title(cell2mat(all_site_names(i)))
            set(gca, 'XColor', 'k', 'YColor', 'k','XAxisLocation','bottom','Xtick',0:300:1200,...
            'YAxisLocation','left','Ytick',-9:1:0,'fontname','arial ','fontsize',12,'TickLength',[0.025 0.025],'linewidth',2);
            xlim([0 (abs(max_melt).*300)])
            ylim([max_melt 0])
            xlabel('PDDs (^oC)','fontweight','bold')
            ylabel('A_{sfc} (m w.e.)','fontweight','bold')
            box on
            axis square
            set(gcf, 'PaperPositionMode', 'auto');
            print -depsc2 gates_epoch2.eps
            end

        end
    end
    %fit global DDF for snow using all sites data
    snow_ablation_model=all_snow_ablation_data(:,1)\all_snow_ablation_data(:,2);
    Snow_Degree_Day_factor(i+1,1)=(all_snow_ablation_data(:,1)\all_snow_ablation_data(:,2)); 
    % Snow_Degree_Day_factor_rsq=[Snow_Degree_Day_factor_rsq;snow_ablation_model.Rsquared.Ordinary];
    % Snow_Degree_Day_factor_rmse=[Snow_Degree_Day_factor_rmse;snow_ablation_model.RMSE];
    [r,p] = corr(all_snow_ablation_data(:,1),all_snow_ablation_data(:,2),'type','Kendall');
    Snow_Degree_Day_factor_rsq(i+1,1)=r;
    Snow_Degree_Day_factor_rmse(i+1,1)=nanmean(abs((all_snow_ablation_data(:,1).*Snow_Degree_Day_factor(end))-all_snow_ablation_data(:,2)));
    difference_in_ks=abs(ks-Snow_Degree_Day_factor(end));
    ks=Snow_Degree_Day_factor(end);
    if plot_calibration==1 && ~isempty(difference_in_ks) && ~sum(difference_in_ks>0.00001)    
        %plot global fit 
        figure(3);hold on
        scatter(all_snow_ablation_data(:,1),all_snow_ablation_data(:,2),'filled');hold on
        y_est = all_snow_ablation_data(:,1)*ks;
         plot(all_snow_ablation_data(:,1),y_est,'-r') 
         set(gca, 'XColor', 'k', 'YColor', 'k','XAxisLocation','bottom','Xtick',0:300:1200,...
        'YAxisLocation','left','Ytick',-9:1:0,'fontname','arial ','fontsize',12,'TickLength',[0.025 0.025],'linewidth',2);
        xlim([0 (abs(max_melt).*300)])
        ylim([max_melt 0])
        xlabel('PDDs (^oC)','fontweight','bold')
        ylabel('A_{sfc} (m w.e.)','fontweight','bold')
        box on
        axis square
        set(gcf, 'PaperPositionMode', 'auto');
        print -depsc2 gates_epoch2.eps
    end

    %% Calibrate Degree Day Model for Ice
    % since most of the measurements of ice ablation also inculde some snow
    % melt we need to use our calibrated precipitation and snow melt
    % coefficents to resolve this final parameter.
all_ice_ablation_data=[];
% site_ice_ablation_data=[];
    for i=1:length(years)
        insitu_data=Glaciological_data(Glaciological_data.Year==years(i) & ~isnan(Glaciological_data.spring_date) & ~isnan(Glaciological_data.fall_date) & ~isnan(Glaciological_data.bw) & ~isnan(Glaciological_data.ba),:);
        %We need to get the dates for eavery measured ablation measurement. 
        %These will typically be between the stratigraphic mass maximum and the
        %measurement date in the fall  
        if isempty(insitu_data)
        else
        [stratigraphic_mass_maximum_correction,stratigraphic_mass_maximum_date,stratigraphic_mass_minimum_correction,stratigraphic_mass_minimum_date] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);

            %loop through each site
            for site=1:height(insitu_data)
                if ~isnan(insitu_data.ba(site)) && ~isnan(insitu_data.bw(site)) && insitu_data.ba(site)<0 && (insitu_data.spring_date(site)>=stratigraphic_mass_maximum_date(site))%make sure both bw and ba were measured and confirm melting only snow during measurement interval

                    if stratigraphic_mass_maximum_date(site) <= insitu_data.spring_date(site)           %spring observation happend after mass maximum, so should ablating mass from now through fall
                        time_and_swe_1=[insitu_data.spring_date(site) insitu_data.bw(site)];

                    elseif stratigraphic_mass_maximum_date(site) > insitu_data.spring_date(site)    %spring bservation happend before the mass maximum so we need to account for any additional snow that would have fallen and then melted
                        time_and_swe_1=[stratigraphic_mass_maximum_date(site) (insitu_data.bw(site)+stratigraphic_mass_maximum_correction(site))];                  

                    end

                    if stratigraphic_mass_minimum_date(site) >= insitu_data.fall_date(site)%fall observation happend before mass minimum, so date and mass balance of observation are used
                        time_and_swe_2=[insitu_data.fall_date(site) insitu_data.ba(site)];

                    elseif stratigraphic_mass_minimum_date(site) < insitu_data.fall_date(site)%fall observation happend after mass minimum, so we need to use the end of the ablation year as predicted by the ablation model
                        time_and_swe_2=[stratigraphic_mass_minimum_date(site) insitu_data.ba(site)];                  

                    end

                    bs=time_and_swe_2(1,2);                         %Since we are just solving for ice ablation the majority of ba values shouldd be of ice ablation
                    precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,insitu_data.site_name(site)));
                    if isempty(precipitation_ratio)
                        precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,'All_sites'));
                    end
                    [Site_Weather]=Model_Site_Weather(Weather_data,insitu_data.elevation(site),time_and_swe_1(1,1),time_and_swe_2(1,1),glacier,lapse_rate,precipitation_ratio);
                    if (sum(Site_Weather.T_flag)/height(Site_Weather))<.5 && (sum(Site_Weather.P_flag)/height(Site_Weather))<.5
                    accumulation=0;
                                Positive_Degree_Days=0;
                                for day=1:height(Site_Weather)
                                   accumulation=accumulation+Site_Weather.Precipitation(day,1);
                                   ablation=Site_Weather.Temperature(day,1).*ks;
                                   accumulation=accumulation+ablation;
                                       if accumulation<=0 && Site_Weather.Temperature(day,1)>0
                                            accumulation=0;
                                            Positive_Degree_Days=Positive_Degree_Days+Site_Weather.Temperature(day,1); 
                                       end
                                 end
    %                 Positive_Degree_Days=sum(Site_Weather.Temperature(Site_Weather.Temperature>0,1));
                    Positive_Degree_Days = Positive_Degree_Days + (time_and_swe_1(1,2)/Snow_Degree_Day_factor(end,1));
                    site_ice_ablation_data{strcmp(insitu_data.site_name(site),all_site_names)}=[site_ice_ablation_data{strcmp(insitu_data.site_name(site),all_site_names)}; Positive_Degree_Days bs time_and_swe_1(1,1) time_and_swe_2(1,1) insitu_data.Year(1)];
                    end
                end

                if ~isnan(insitu_data.winter_ablation(site)) && insitu_data.winter_ablation(site) < 0
                    previous_year_data_index=find(Glaciological_data.Year==years(i)-1); 
                    site_index=strcmp(Glaciological_data.site_name(previous_year_data_index),insitu_data.site_name(site));
                    site_index_for_previous_year=previous_year_data_index(site_index);
                    if ~isempty(site_index_for_previous_year) && ~isnan(Glaciological_data.ba(site_index_for_previous_year)) && ~isnan(Glaciological_data.fall_date(site_index_for_previous_year))
                        if Glaciological_data.ba(site_index_for_previous_year)< 0 
                            [~,~,~,previous_year_stratigraphic_mass_minimum_date] = Find_Mass_Maximum_and_Minimum_Adjustments(Glaciological_data(site_index_for_previous_year,:),Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);
                            precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,insitu_data.site_name(site)));
                            if isempty(precipitation_ratio)
                                precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,'All_sites'));
                            end
                            [Site_Weather]=Model_Site_Weather(Weather_data,insitu_data.elevation(site),Glaciological_data.fall_date(site_index_for_previous_year),previous_year_stratigraphic_mass_minimum_date,glacier,lapse_rate,precipitation_ratio);
                            if (sum(Site_Weather.T_flag)/height(Site_Weather))<.5 && (sum(Site_Weather.P_flag)/height(Site_Weather))<.5
                            accumulation=0;
                                Positive_Degree_Days=0;
                                for day=1:height(Site_Weather)
                                   accumulation=accumulation+Site_Weather.Precipitation(day,1);
                                   ablation=Site_Weather.Temperature(day,1).*ks;
                                   accumulation=accumulation+ablation;
                                       if accumulation<=0 && Site_Weather.Temperature(day,1)>0
                                            accumulation=0;
                                            Positive_Degree_Days=Positive_Degree_Days+Site_Weather.Temperature(day,1); 
                                            
                                       end
                                end
                            site_ice_ablation_data{strcmp(Glaciological_data.site_name(site_index_for_previous_year),all_site_names)}=[site_ice_ablation_data{strcmp(Glaciological_data.site_name(site_index_for_previous_year),all_site_names)}; Positive_Degree_Days,Glaciological_data.winter_ablation(site_index_for_previous_year) Glaciological_data.fall_date(site_index_for_previous_year) previous_year_stratigraphic_mass_minimum_date insitu_data.Year(1)];                    
                            end
                        end
                    end
                end
            end
        end
        if ~isempty(additional_data)
            bonus_insitu_data=additional_data(additional_data.Year==years(i) & additional_data.db_mwe<=0 & strcmp(additional_data.Surface1,'ice') & strcmp(additional_data.Surface2,'ice'),:);
             for i=1:height(bonus_insitu_data)
                 precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,bonus_insitu_data.Site_Name(i)));
                    if isempty(precipitation_ratio)
                        precipitation_ratio=precipitation_ratio_table.precipitation_ratios(strcmp(precipitation_ratio_table.site_name,'All_sites'));
                    end
                    time_1=datenum(bonus_insitu_data.Date1(i));
                    time_2=datenum(bonus_insitu_data.Date2(i));
                    [Site_Weather]=Model_Site_Weather(Weather_data,bonus_insitu_data.Elevation_m(i),time_1,time_2,glacier,lapse_rate,precipitation_ratio); %model site temperatures
                    if (sum(Site_Weather.T_flag)/height(Site_Weather))<.5 && (sum(Site_Weather.P_flag)/height(Site_Weather))<.5

                                accumulation=0;
                                Positive_Degree_Days=0;
                                
                                for day=1:height(Site_Weather)
                                   accumulation=accumulation+Site_Weather.Precipitation(day,1);
                                   if Site_Weather.Temperature(day,1)>0
                                        ablation=Site_Weather.Temperature(day,1).*ks;
                                   else
                                       ablation=0;
                                   end

                                   accumulation=accumulation+ablation;
                                   plot(day,accumulation)
                                       if accumulation<=0 && Site_Weather.Temperature(day,1)>0
    %                                         degree_day
                                            accumulation=0;
                                            Positive_Degree_Days=Positive_Degree_Days+Site_Weather.Temperature(day,1); 
                                       end
                                end

    %                 Positive_Degree_Days=sum(Site_Weather.T(Site_Weather.T>0,1));  
                        %get positive degree days (T > 0*C)
                        if sum(strcmp(precipitation_ratio_table.site_name,bonus_insitu_data.Site_Name(i)))==0
                            all_ice_ablation_data=[all_ice_ablation_data; Positive_Degree_Days bonus_insitu_data.db_mwe(i) time_1 time_2 bonus_insitu_data.Year(1)];
                        else
                            site_ice_ablation_data{strcmp(bonus_insitu_data.Site_Name(i),all_site_names)}=[site_ice_ablation_data{strcmp(bonus_insitu_data.Site_Name(i),all_site_names)};Positive_Degree_Days bonus_insitu_data.db_mwe(i) time_1 time_2 bonus_insitu_data.Year(1)];
                        end
                    end
                end
        end
    end
    max_melt=0;
    
    for i=1:length(site_ice_ablation_data)
        if isempty(site_ice_ablation_data{i})
        elseif length(site_ice_ablation_data{i}(:,2))>=3
        elseif length(site_ice_ablation_data{i}(:,2))==2
        else
        end
        if isempty(site_ice_ablation_data{i})
            ice_Degree_Day_factor(i,1)=nan; 
            ice_Degree_Day_factor_rsq(i,1)=nan;
            ice_Degree_Day_factor_rmse(i,1)=nan;   
        else
            all_ice_ablation_data=[all_ice_ablation_data;site_ice_ablation_data{i}];
            ice_Degree_Day_factor(i,1)=site_ice_ablation_data{i}(~isnan(site_ice_ablation_data{i}(:,1))&~isnan(site_ice_ablation_data{i}(:,2)),1)\site_ice_ablation_data{i}(~isnan(site_ice_ablation_data{i}(:,1))&~isnan(site_ice_ablation_data{i}(:,2)),2);
            [r,p] = corr(site_ice_ablation_data{i}(~isnan(site_ice_ablation_data{i}(:,1))&~isnan(site_ice_ablation_data{i}(:,2)),1),site_ice_ablation_data{i}(~isnan(site_ice_ablation_data{i}(:,1))&~isnan(site_ice_ablation_data{i}(:,2)),2),'type','Kendall');
            ice_Degree_Day_factor_rsq(i,1)=r;
            ice_Degree_Day_factor_rmse(i,1)=nanmean(abs((site_ice_ablation_data{i}(:,1).*ice_Degree_Day_factor(i,1))-site_ice_ablation_data{i}(:,2)));
        end

        if max_melt==0
            max_melt=min(all_ice_ablation_data(:,2));
        end


        if isempty(site_ice_ablation_data{i})
        else
            if plot_calibration==1  && ~isempty(difference_in_ki) && ~sum(difference_in_ki>0.00001)
            figure(4);hold on
            subplot(ceil(sqrt(length(all_site_names))),ceil(sqrt(length(all_site_names))),i)
            scatter(site_ice_ablation_data{i}(:,1),site_ice_ablation_data{i}(:,2),'filled');hold on
            y_est = all_ice_ablation_data(:,1)*ki;
            plot(all_ice_ablation_data(:,1),y_est,'-r')
            title(cell2mat(all_site_names(i)))
            set(gca, 'XColor', 'k', 'YColor', 'k','XAxisLocation','bottom','Xtick',0:300:1200,...
            'YAxisLocation','left','Ytick',-9:1:0,'fontname','arial ','fontsize',12,'TickLength',[0.025 0.025],'linewidth',2);
            xlim([0 (abs(max_melt).*300)])
            ylim([max_melt 0])
            xlabel('PDDs (^oC)','fontweight','bold')
            ylabel('A_{sfc} (m w.e.)','fontweight','bold')
            box on
            axis square
            set(gcf, 'PaperPositionMode', 'auto');
            print -depsc2 gates_epoch2.eps
            end

        end
    end

    %fit global DDF for ice using all sites data
    ice_ablation_model=all_ice_ablation_data(~isnan(all_ice_ablation_data(:,1))&~isnan(all_ice_ablation_data(:,2)),1)\all_ice_ablation_data(~isnan(all_ice_ablation_data(:,1))&~isnan(all_ice_ablation_data(:,2)),2);
    ice_Degree_Day_factor(i+1,1)=(all_ice_ablation_data(~isnan(all_ice_ablation_data(:,1))&~isnan(all_ice_ablation_data(:,2)),1)\all_ice_ablation_data(~isnan(all_ice_ablation_data(:,1))&~isnan(all_ice_ablation_data(:,2)),2)); 
    [r,p] = corr(all_ice_ablation_data(~isnan(all_ice_ablation_data(:,1))&~isnan(all_ice_ablation_data(:,2)),1),all_ice_ablation_data(~isnan(all_ice_ablation_data(:,1))&~isnan(all_ice_ablation_data(:,2)),2),'type','Kendall');
    ice_Degree_Day_factor_rsq(i+1,1)=r;
    ice_Degree_Day_factor_rmse(i+1,1)=nanmean(abs((all_ice_ablation_data(:,1).*ice_Degree_Day_factor(end))-all_ice_ablation_data(:,2)));
    difference_in_ki=abs(ki-ice_Degree_Day_factor(end));
    ki=ice_Degree_Day_factor(end);

    if plot_calibration==1  && ~isempty(difference_in_ki) && ~sum(difference_in_ki>0.0001)
    figure(5);hold on
    scatter(all_ice_ablation_data(:,1),all_ice_ablation_data(:,2),'filled');hold on
    y_est = all_ice_ablation_data(:,1)*ki;
    plot(all_ice_ablation_data(:,1),y_est,'-r') 
     set(gca, 'XColor', 'k', 'YColor', 'k','XAxisLocation','bottom','Xtick',0:300:1200,...
    'YAxisLocation','left','Ytick',-9:1:0,'fontname','arial ','fontsize',12,'TickLength',[0.025 0.025],'linewidth',2);
    xlim([0 (abs(max_melt).*300)])
    ylim([max_melt 0])
    xlabel('PDDs (^oC)','fontweight','bold')
    ylabel('A_{sfc} (m w.e.)','fontweight','bold')
    box on
    axis square
    set(gcf, 'PaperPositionMode', 'auto');
    print -depsc2 gates_epoch2.eps
    end
end
Degree_Day_Factor_table=table([all_site_names; 'All_Sites'],Snow_Degree_Day_factor,ice_Degree_Day_factor,Snow_Degree_Day_factor_rsq,ice_Degree_Day_factor_rsq,Snow_Degree_Day_factor_rmse,ice_Degree_Day_factor_rmse,'VariableNames',{'site_name','ks','ki','Rsquared_ks','Rsquared_ki','RMSE_ks','RMSE_ki'})
Degree_Day_Factor_table=Degree_Day_Factor_table(~(isnan(Degree_Day_Factor_table.ks)&isnan(Degree_Day_Factor_table.ki)),:);
ready=1;
end
