function [mass_maximum_site_adjustments,mass_maximum_date_numbers,mass_minimum_site_adjustments,mass_minimum_date_numbers,Modeled_Mass_Balance] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model)
%% primary objective of this function is to perform short term accumulation
%% and ablation modeling to make adjustments to insitu measurements to within
%% three time-systems

    % 1) stratigraphic - site change in mass balance between field
    % observation and site specific mass maximum and mass minimum
    
    % 2) purely observed - mass balance model is not employed and any
    % additional mass change not recorded is considered negligible
    
    % 3) Fixed-Date - site change in mass balance between field
    % observation and May 1st (hydrologic year spring) and September 30th
    % (end of hydrologic year)
    
    % 4) Combined stratigraphic/Floating-Date - site change in mass balance between field
    % observation and glacier-wide mass maximum and mass minimum date
    
%% General Outline
       %1) Initalizes variables needed 
                % site/weather station elevation differences for lapse
                % rate, weather data for the year, glacier hypsometry, max
                % precip temperatures and minimum melt temperatures
                
        %2) Calculate site site weather and first orer mass blance time-series 
                % calculate daily precipitation and temperature at the
                % site, generate a approximate mass balance time-series for
                % the site using ks(snow ablation coefficent) to find
                % approximate mass max and minimum timing, and generate
                % approiate bounds for the time-series in order to refine
                % time-series on a site by site basis in step 3
        
        %3) Refine site mass balance time-series site-by-site
                % determine if a site specific mass balance time-series
                % needs to be refined. This occurs if the site may have
                % melted to ice. This occurs if ba < 0, if ba >0 but before
                % the mass minimum date, or ba is NaN.
        
        %4) Find the mass max and minimum date(s) and adjustements
                % the date of mass maximums and minimums and needed
                % water equivalent adjustements needed to put the site
                % observations into the users selected time-series
                
        %5) Plot site mass balance time-series
                % graphic display of modeled masss balance time-series with
                % the observations dates and mass max and minimums dates
                % (based on users selected time-system) highlighted
                
               
%% Input data
%       1) Table of insitu_data for a single year
%       2) Lapse Rate to be used to model mass balance
%       3) ks: calibrated snow melt coefficent for classic temperature
%       index model
%       4) ki: calibrated ice melt coefficent for classic temperature
%       index model

%% Output
%       1) mass maximum date numbers
%       2) mass maximum water equivalent adjustments
%       3) mass minimum date numbers
%       4) mass minimum water equivalent adjustments
%       5) modeled daily mass balance for each site
%% 1) Initalizes variables needed
dbstop if error
if length(unique(Weather_data.Elevation(~isnan(Weather_data.Elevation))))>1
      fprintf(1,'ERROR: You have date from two different weather stations in you daily weather observations. This code does not deal with that yet. Either change the input of fix the code!')
end

%Elevations
weather_station_altitude=Weather_data.Elevation(1); 
Site_elevation_differance = insitu_data.elevation - weather_station_altitude;

%Year
year=insitu_data.Year(1);
Years=str2num(datestr(Weather_data.Date,'YYYY'));

%Weather data for that year
wx_data_for_balance_year=Weather_data(Years==year,:);

% approximate mass balance extreme dates
spring_hydrologic_date_number=datenum(['5/1/',num2str(year)]);
fall_hydrologic_date_number=datenum(['9/30/',num2str(year)]);
% loop through each site and find mass max and minimum
water_equivalent_snow=zeros(height(wx_data_for_balance_year),height(insitu_data));%vector to fill if there is any accumulation
ablation_days=zeros(height(wx_data_for_balance_year),height(insitu_data));%vector to fill if there is any accumulation

%get glacier hypsometry for the year
bin_centers=AAD(1,2:end);
if year==AAD(2,1)-1
    glacier_hypsometry=AAD(AAD(:,1)==year+1,2:end);
else
    glacier_hypsometry=AAD(AAD(:,1)==year,2:end);
end

minimum_temperature=0; % minimum temp below which dry snow falls rather than rain
maximum_temperature=1.7; % maximum temp above which rain falls, in between wet snow falls
%% 2) Calculate site site weather and first order mass blance time-series
        %A) for each site model the temperature and precipitation using the lapse rate and precipitation ratios provided
        %B) Calculate first order mass balance time-series using ks for a
        %   degree day factor.
        %C) Find approximate stratigraphic mass maximum and minimum dates
        %   to create better bounds for the mass balance time-series. This
        %   later allows for a more detailed analysis of the mass balance
        %   during the mass maximum and minimum periods
        %D) Trim weather and mass balance vectors to the bounds determined
        %   in step C. Initalize for site-by-site mass balance time-series
        %`  refinement
        
        
%A) for each site model the temperature and precipitation using the lapse rate and precipitation ratios provided
dates=[];
for site=1:height(insitu_data)    %for each site
%     if isnan(precipitation_ratios_table.Precipitation_ratios(strcmp(precipitation_ratios_table.site_name,insitu_data.site_name(site))));
        
    Precipitation_ratio=precipitation_ratios_table.precipitation_ratios(strcmp(precipitation_ratios_table.site_name,insitu_data.site_name(site))); %get the site specific precipitation ratio
    if isempty(Precipitation_ratio)||isnan(Precipitation_ratio)%then no precipitation ratio could be calibrated for this site and we need to find one to use
        Precipitation_ratio=precipitation_ratios_table.precipitation_ratios(strcmp(precipitation_ratios_table.site_name,'All_sites'));%so we use a global fit
    end
    year_and_day=yearday(insitu_data.spring_date(site));    %get day number (1-365/366) for spring observation
    observation_day_index(site,1)=year_and_day(2);          %stash the day number in this matrix
    year_and_day=yearday(insitu_data.fall_date(site));      %get day number (1-365/366) for fall observation
    if year_and_day(1)~=year && insitu_data.summer_accumulation(site)>0 %some ba measurements werent made until mid winter, so we provide a approximate end of season date to constrain the model. The final adjustment from this site will be zero, but we 
        year_and_day=yearday(datenum(['9/30/',num2str(year)]));% so we provide a approximate end of season date to constrain the model.
    end
    observation_day_index(site,2)=year_and_day(2);          %stash the day number in this matrix  
    
    date_numbers=datenum(wx_data_for_balance_year.Date);    %date numbers to later compare to field observation dates
    site_temperature(:,site)=wx_data_for_balance_year.T + (lapse_rate*Site_elevation_differance(site)/1000); %model temperature at glaciological site
    site_precipitation(:,site)=wx_data_for_balance_year.P/1000*Precipitation_ratio; %model precipitation at site
    water_equivalent_snow(site_temperature(:,site)<maximum_temperature,site)=site_precipitation(site_temperature(:,site)<maximum_temperature,site); %pulling modeled daily site precipitation that fell at cold enough temperatures it should be snow
    ablation_days(site_temperature(:,site)>minimum_temperature,site)=site_temperature(site_temperature(:,site)>minimum_temperature,site); %pulling modeled daily temperatures that were warm enough to melt either snow or ice
    
    %B) Calculate first order mass balance time-series using ks for a
        %   degree day factor.
    site_modeled_mass_balance_time_series(:,site)=cumsum(water_equivalent_snow(:,site)+(ablation_days(:,site).*ks)); %assuming we only melted snow. This should be fine for spring adjustments, but we recalculate the time-series for fall if ice was melted
    
    %C) Find approximate stratigraphic mass maximum and minimum dates
    %   to create better bounds for the mass balance time-series. This
    %   later allows for a more detailed analysis of the mass balance
    %   during the mass maximum and minimum periods
    [xmax,imax,xmin,imin] = extrema(site_modeled_mass_balance_time_series(:,site)); %find approximate maximum dates to, trim the modeled time-series
    
    %find spring maximums closest to spring observation dates
    [minimum_difference,index]=min(abs(insitu_data.spring_date(site)-date_numbers(imax)));    
    approximate_mass_max_date(site)=imax(index); %approximate mass max. This is refined after we reduce the time-series to the appropriate subsect of time
    approximate_spring_dates=sort([spring_hydrologic_date_number date_numbers(approximate_mass_max_date(site)) insitu_data.spring_date(site)]); %sort the important dates for the spring to find the earliest date needed
    year_day=yearday(approximate_spring_dates(1)); %get day number(1-365/3666) of earliest spring dates of interest. This is our date bound, for when we need to be modeling the mass balance
    daily_balance_start_and_end(site,1)=year_day(2); %stash spring bounds for each site
    
    %find spring maximums closest to fall observation dates
    [minimum_difference,index]=min(abs(insitu_data.fall_date(site)-date_numbers(imin)));
    approximate_mass_min_date(site)=imin(index);%approximate mass min. This is refined after we reduce the time-series to the appropriate subsect of time
    approximate_fall_dates=sort([fall_hydrologic_date_number date_numbers(approximate_mass_min_date(site)) insitu_data.fall_date(site)]);%sort the important dates for the fall to find the latest date needed
    year_day=yearday(approximate_fall_dates(end)); %get day number(1-365/3666) of earliest fall dates of interest. This is our date bound, for when we need to be modeling the mass balance
    daily_balance_start_and_end(site,2)=year_day(2);%stash fall bounds for each site
end

bounds=[min(daily_balance_start_and_end(:,1))-30 max(daily_balance_start_and_end(:,2))+30];%find earliest and latest dates bounds and then tack on and extra month on either end to make sure we can find the mass max and min
%bounds=[30 max(daily_balance_start_and_end(:,2))+30];%find earliest and latest dates bounds and then tack on and extra month on either end to make sure we can find the mass max and min

if bounds(1)<1||bounds(1)>=30
    bounds(1)=1;
end
if bounds(2)>365
    bounds(2)=365;
end

%D) Trim weather and mass balance vectors to the bounds determined
%   in step C. Initalize for site-by-site mass balance time-series
%`  refinement
date_numbers=date_numbers(bounds(1):bounds(2),:);
site_temperature=site_temperature(bounds(1):bounds(2),:);
site_precipitation=site_precipitation(bounds(1):bounds(2),:);
site_modeled_mass_balance_time_series=site_modeled_mass_balance_time_series(bounds(1):bounds(2),:);
water_equivalent_snow=water_equivalent_snow(bounds(1):bounds(2),:);
ablation_days=ablation_days(bounds(1):bounds(2),:);
observation_day_index=(observation_day_index-bounds(1))+1;

bounds=bounds-bounds(1)+1;%reset bounds in reference to the site mass balance data water_equivalent_snow and ablation_days
bounds(2)=bounds(2);
ablation_season_mid_point=round(mean(bounds)-bounds(1));
[~,stratigraphic_mass_maximum_index]=max(site_modeled_mass_balance_time_series(1:ablation_season_mid_point,:));
stratigraphic_mass_maximum_index=stratigraphic_mass_maximum_index';
stratigraphic_mass_maximum_date_numbers=date_numbers(stratigraphic_mass_maximum_index);

[~,stratigraphic_mass_minimum_index]=min(site_modeled_mass_balance_time_series(ablation_season_mid_point:end,:));
stratigraphic_mass_minimum_index=stratigraphic_mass_minimum_index'+(ablation_season_mid_point-1);
stratigraphic_mass_minimum_date_numbers=date_numbers(stratigraphic_mass_minimum_index);

observation_day_index(observation_day_index(:,2)>stratigraphic_mass_minimum_index,2)=stratigraphic_mass_minimum_index(observation_day_index(:,2)>stratigraphic_mass_minimum_index,1);

Cumulative_ablation_season_PDDs=zeros(bounds(2),height(insitu_data));
daily_site_ablation=zeros(bounds(2),height(insitu_data));
mass_maximum_site_adjustments=nan*zeros(1,height(insitu_data));
mass_minimum_site_adjustments=nan*zeros(1,height(insitu_data));

%% 3) Refine site mass balance time-series site-by-site
        % A) Determine if a site might have melted ice at any point using
        %   the time-series of Positive Degree Days and provided Degree Day
        %   Coefficents
        %       i)  if ba < 0
        %       ii) if ba > 0 but stratigraphic mass minimum date was later
        %           than the observation date
        %       iii) if ba is a NaN then we use bw--if available
            
% loop through each site to refine the mass balance model if needed.
% this is done by determining if and when the site transistioned to ice
% find the Transient snowline date based on what we know of the
% relationship to snow and ice melt with reguards to temperature
 for site=1:height(insitu_data)
   TSL_index=[];
    % if we know ice melted at some point (ba<0) or there is a chance ice melted at some point (ba date < stratigraphic minimum date) 
    if insitu_data.ba(site)<0 || (insitu_data.ba(site)>0 && insitu_data.fall_date(site)<stratigraphic_mass_minimum_date_numbers(site)) || (isnan(insitu_data.ba(site)) && ~isnan(insitu_data.bw(site)))
        
        Cumulative_ablation_season_PDDs(:,site)=cumsum(ablation_days(:,site)); %get the cumulative positive degree days
        
        % if the site was ice during fall observations
        if insitu_data.ba(site)<0 && ~contains(insitu_data.site_name(site),'TSL') %if measured annual balance was less than zero, we know the site was ice
            
            %reverse the cumulative PPDs time-series to count back from ba measurement to find when snowline passed
            Cumulative_ablation_season_PDDs(:,site)=abs(Cumulative_ablation_season_PDDs(observation_day_index(site,2),site)-Cumulative_ablation_season_PDDs(:,site));
            PDD_Equivalent=abs(insitu_data.ba(site)./ki); %calculate number of PDDs it took to melt the measured ba
            TSL_index=find(Cumulative_ablation_season_PDDs(:,site)<=PDD_Equivalent,1,'first');% index for when to swithch DDFs from snow to ice

        % if the measured ba was positive, we know the site was still snow. But we need to make sure it didnt melt out to ice    
        elseif insitu_data.ba(site)>0 && insitu_data.fall_date(site)<stratigraphic_mass_minimum_date_numbers(site)
            %zero cumulative PDD time-series to observation date
%             Cumulative_ablation_season_PDDs(:,site)=abs(Cumulative_ablation_season_PDDs(observation_day_index(site,2),site)-Cumulative_ablation_season_PDDs(:,site));
            Cumulative_ablation_season_PDDs(:,site)=Cumulative_ablation_season_PDDs(:,site)-Cumulative_ablation_season_PDDs(observation_day_index(site,2),site);
            %dont reverse the cumulative PDDs time-series, since we are
            %counting forward 
            PDD_Equivalent=abs(insitu_data.ba(site)./ks); %calculate number of PDDs it took to melt the measured ba
            TSL_index=find(Cumulative_ablation_season_PDDs(:,site)>=PDD_Equivalent,1,'first');
        % if ba was not measured, we have to use bw instead, as long as it
        % was measured
        elseif isnan(insitu_data.ba(site)) && ~isnan(insitu_data.bw(site)) && ~contains(insitu_data.site_name(site),'TSL')
            if insitu_data.spring_date(site)>=stratigraphic_mass_maximum_date_numbers(site)
                Cumulative_ablation_season_PDDs(:,site)=Cumulative_ablation_season_PDDs(:,site)-Cumulative_ablation_season_PDDs(observation_day_index(site,1),site);
                PDD_Equivalent=abs(insitu_data.bw(site)./ks); %calculate number of PDDs it took to melt the measured bw
               
            elseif insitu_data.spring_date(site)<stratigraphic_mass_maximum_date_numbers(site)
                Cumulative_ablation_season_PDDs(:,site)=Cumulative_ablation_season_PDDs(:,site)-Cumulative_ablation_season_PDDs(stratigraphic_mass_maximum_index(site,1),site);
                mass_maximum_site_adjustments=site_modeled_mass_balance_time_series(stratigraphic_mass_maximum_index(site,1),site)-site_modeled_mass_balance_time_series(observation_day_index(site,1),site);
                PDD_Equivalent=abs((insitu_data.bw(site)+mass_maximum_site_adjustments)./ks); %calculate number of PDDs it took to melt the measured bw
            end
            %dont reverse the cumulative PDDs time-series, since we are
            %counting forward 
            TSL_index=find(Cumulative_ablation_season_PDDs(:,site)>=PDD_Equivalent,1,'first');
        elseif insitu_data.bw(site)==0 && contains(insitu_data.site_name(site),'TSL')
            TSL_index=find(date_numbers==insitu_data.spring_date(site),1,'first');
        end
        
        % so now that we have determined if and when the site melted to ice
        accumulation=0;
         for day=bounds(1):bounds(2)
             if  ~isempty(TSL_index) && day >= TSL_index% than the snowline has migrated past and we are melting ice
               %Keep checking site precipitation in case
               %there was some mid summer snow
               accumulation=accumulation+water_equivalent_snow(day,site); %build snowpack if precipitation occured
               if accumulation>0 %then we are still melting snow at site
                    ablation=ablation_days(day,site).*ks; %melt snowpack
                    accumulation=accumulation+ablation; %add accumulation and ablation until snow melts away again
                    daily_site_ablation(day,site)=ablation_days(day,site).*ks; 
                    
               else % otherwise we assume we are melting ice at this site since balance is < 0 m w.e.
                   ablation=0; %clear ablation incase it snows again
                   accumulation=0;
                   daily_site_ablation(day,site)=ablation_days(day,site).*ki;
           
               end
             else
                daily_site_ablation(day,site)=ablation_days(day,site).*ks; 
             end
         end 
    else
        daily_site_ablation(:,site)=ablation_days(:,site).*ks;
    end
        
        site_modeled_mass_balance_time_series(:,site)=cumsum(water_equivalent_snow(:,site)+daily_site_ablation(:,site));
end        
%% 4) Find the mass max and minimum date(s) and adjustements
    [~,stratigraphic_mass_minimum_index]=min(site_modeled_mass_balance_time_series(ablation_season_mid_point:end,:));%find minimum value in mass balance time-series
    stratigraphic_mass_minimum_index=stratigraphic_mass_minimum_index+(ablation_season_mid_point-1); %index of mass minimum
    stratigraphic_mass_minimum_index= stratigraphic_mass_minimum_index';
    stratigraphic_mass_minimum_date_numbers=date_numbers(stratigraphic_mass_minimum_index); %mass minimum date numbers
    
if ~isdatetime(time_system)
    if time_system==1 || time_system==2%stratigraphic time-system. Mass maximum and minimum determined on a site by site basis
        % spring max already computed above
        mass_maximum_index=stratigraphic_mass_maximum_index;
        mass_minimum_index=stratigraphic_mass_minimum_index;
        mass_maximum_date_numbers=stratigraphic_mass_maximum_date_numbers;
        mass_minimum_date_numbers=stratigraphic_mass_minimum_date_numbers;
    elseif time_system==3 %fixed-date hydrologic year time-system. Mass maximum and minimums are determined on May 1st and September 30th each year. This corresponds with the water year
        mass_maximum_index(1:height(insitu_data),1)=find(date_numbers==spring_hydrologic_date_number);
        mass_minimum_index(1:height(insitu_data),1)=find(date_numbers==fall_hydrologic_date_number);
        mass_maximum_date_numbers=ones(height(insitu_data),1).*spring_hydrologic_date_number;
        mass_minimum_date_numbers=ones(height(insitu_data),1).*fall_hydrologic_date_number;
    elseif time_system==4%combined floating-date/stratigraphic system. Determine the glacier-wide average mass maximum and minimum dates using the glaciers AAD.
        [~,mass_maximum_index]=max(site_modeled_mass_balance_time_series(1:ablation_season_mid_point,:));
        [~,mass_minimum_index]=min(site_modeled_mass_balance_time_series(ablation_season_mid_point:end,:));
        mass_minimum_index=mass_minimum_index'+(ablation_season_mid_point-1);
        mass_maximum_index=mass_maximum_index';
        mass_maximum_date_numbers=date_numbers(mass_maximum_index);
        mass_minimum_date_numbers=date_numbers(mass_minimum_index);
        mass_maximum_index(1:height(insitu_data),1)=find(date_numbers==mass_maximum_date_numbers(1));
        mass_minimum_index(1:height(insitu_data),1)=find(date_numbers==mass_minimum_date_numbers(1));
        daily_glacier_wide=[];

        for day=min(mass_maximum_date_numbers):max(mass_maximum_date_numbers)
            balances=[];
            balances=site_modeled_mass_balance_time_series(date_numbers==day,1:length(site_modeled_mass_balance_time_series(1,:)))';
            if ~isempty(balances)
                point_balances=table(insitu_data.Year,insitu_data.site_name,insitu_data.elevation,balances,nan*ones(length(balances),1),nan*ones(length(balances),1),day*ones(length(balances),1),mass_minimum_date_numbers,'VariableNames',{'Year' 'site_name' 'elevation' 'bw' 'bs' 'ba' 'mass_maximum_date' 'mass_minimum_date'});
                Glacier_Wide_solutions = integrate_point_balance(glacier,year,point_balances,AAD,time_system,integration_method,integration_surface,0,0);
                daily_glacier_wide=[daily_glacier_wide;date_numbers(date_numbers==day) Glacier_Wide_solutions.Bw_mwe];
            else
                daily_glacier_wide=[daily_glacier_wide;nan nan];
            end
        end
       [~,max_index]=max(daily_glacier_wide(:,2));
       mass_maximum_date_numbers(1:height(insitu_data),1)=date_numbers(date_numbers==daily_glacier_wide(max_index,1));
       mass_maximum_index(1:height(insitu_data),1)=find(date_numbers==mass_maximum_date_numbers(1));
       daily_glacier_wide=[];
       balances=[];
       for day=min(mass_minimum_date_numbers):max(mass_minimum_date_numbers)
            balances=site_modeled_mass_balance_time_series(date_numbers==day,1:length(site_modeled_mass_balance_time_series(1,:)))';
            if ~isempty(balances)
                point_balances=table(insitu_data.Year,insitu_data.site_name,insitu_data.elevation,nan*ones(length(balances),1),nan*ones(length(balances),1),balances,mass_maximum_date_numbers,day*ones(length(balances),1),'VariableNames',{'Year' 'site_name' 'elevation' 'bw' 'bs' 'ba' 'mass_maximum_date' 'mass_minimum_date'});
                Glacier_Wide_solutions = integrate_point_balance(glacier,year,point_balances,AAD,time_system,integration_method,integration_surface,0,0);
                daily_glacier_wide=[daily_glacier_wide;date_numbers(date_numbers==day) Glacier_Wide_solutions.Ba_mwe];
            else
                daily_glacier_wide=[daily_glacier_wide;nan nan];
            end
       end
       [~,min_index]=min(daily_glacier_wide(:,2));
       mass_minimum_date_numbers(1:height(insitu_data),1)=date_numbers(date_numbers==daily_glacier_wide(min_index,1));
       mass_minimum_index(1:height(insitu_data),1)=find(date_numbers==mass_minimum_date_numbers(1));
    end
elseif isdatetime(time_system)
    mass_maximum_index(1:height(insitu_data),1)=find(date_numbers==datenum(time_system));
    mass_minimum_index(1:height(insitu_data),1)=find(date_numbers==datenum(time_system));
    mass_maximum_date_numbers(1:height(insitu_data),1)=datenum(time_system);
    mass_minimum_date_numbers(1:height(insitu_data),1)=datenum(time_system);
end
for site=1:height(insitu_data)
    if ~isnan(insitu_data.winter_ablation(site))||insitu_data.summer_accumulation(site)>0
        insitu_data.ba(site)=insitu_data.ba(site)+insitu_data.winter_ablation(site);
        observation_day_index(site,2)=stratigraphic_mass_minimum_index(site);
    end
    mass_maximum_site_adjustments(site)=diff([site_modeled_mass_balance_time_series(observation_day_index(site,1),site) site_modeled_mass_balance_time_series(mass_maximum_index(site,1),site)]);
    if site_modeled_mass_balance_time_series(mass_minimum_index(site,1),site) < site_modeled_mass_balance_time_series(observation_day_index(site,2),site)
         if mass_minimum_index(site,1)>=observation_day_index(site,2)
            mass_minimum_site_adjustments(site)=diff([site_modeled_mass_balance_time_series(observation_day_index(site,2),site) site_modeled_mass_balance_time_series(mass_minimum_index(site,1),site)]);
            
         else
            mass_minimum_site_adjustments(site)=diff([site_modeled_mass_balance_time_series(mass_minimum_index(site,1),site) site_modeled_mass_balance_time_series(observation_day_index(site,2),site)]);
         end
    elseif site_modeled_mass_balance_time_series(mass_minimum_index(site,1),site) >= site_modeled_mass_balance_time_series(observation_day_index(site,2),site)
            mass_minimum_site_adjustments(site)=diff([site_modeled_mass_balance_time_series(observation_day_index(site,2),site) site_modeled_mass_balance_time_series(mass_minimum_index(site,1),site)]);
    end
end
 
Modeled_Mass_Balance=table(datestr(date_numbers),'VariableNames',{'Date'});
for site=1:height(insitu_data)
    Modeled_Mass_Balance=[Modeled_Mass_Balance table(site_modeled_mass_balance_time_series(:,site),'VariableNames',insitu_data.site_name(site))];
end

    
%% 5) Plot site mass balance time-series
if plot_ablation_model==1
    figure();hold on
    for site=1:height(insitu_data)
        if site==1
            mass_balance_model=plot(date_numbers,site_modeled_mass_balance_time_series(:,site),'k');hold on
        else
            plot(date_numbers,site_modeled_mass_balance_time_series(:,site),'k');hold on
        end
        
        if site==1
            observation=scatter([insitu_data.spring_date(site) insitu_data.fall_date(site)],[site_modeled_mass_balance_time_series(date_numbers==insitu_data.spring_date(site),site) site_modeled_mass_balance_time_series(date_numbers==insitu_data.fall_date(site),site)],100,'k','*');hold on
            gw_max=line([date_numbers(mass_maximum_index(site)) date_numbers(mass_maximum_index(site))], [min(site_modeled_mass_balance_time_series(:))-.5 max(site_modeled_mass_balance_time_series(:))+.5],'linewidth',2,'LineStyle','--','Color',[.5 .5 .5]);
            gw_min=line([date_numbers(mass_minimum_index(site)) date_numbers(mass_minimum_index(site))], [min(site_modeled_mass_balance_time_series(:))-.5 max(site_modeled_mass_balance_time_series(:))+.5],'linewidth',2,'LineStyle','--','Color',[.5 .5 .5]);
            strat_max=scatter(date_numbers(stratigraphic_mass_maximum_index(site)),site_modeled_mass_balance_time_series(stratigraphic_mass_maximum_index(site),site),100,[.5 0 .5],'d');hold on
            strat_min=scatter(date_numbers(stratigraphic_mass_minimum_index(site)),site_modeled_mass_balance_time_series(stratigraphic_mass_minimum_index(site),site),100,[1 .5 0],'d');hold on
        else
            scatter([observation_day_index(site,1) observation_day_index(site,2)],[site_modeled_mass_balance_time_series(observation_day_index(site,1),site) site_modeled_mass_balance_time_series(observation_day_index(site,2),site)],100,'k','*');hold on
            scatter(date_numbers(stratigraphic_mass_maximum_index(site)),site_modeled_mass_balance_time_series(stratigraphic_mass_maximum_index(site),site),100,[.5 0 .5],'d');hold on
            scatter(date_numbers(stratigraphic_mass_minimum_index(site)),site_modeled_mass_balance_time_series(stratigraphic_mass_minimum_index(site),site),100,[1 .5 0],'d');hold on
        end
        
        text(date_numbers(end)+1,site_modeled_mass_balance_time_series(end,site),[num2str(insitu_data.elevation(site)),'m'])
    end
    legend([mass_balance_model,observation,gw_max,strat_min],'Modeled Mass Balance','Observations','Glacier-wide Extrema','Site Extrema','Location','EO');hold on
    datetick('x','mm-dd')
    title(['Modeled Cumulative Site Balances: ',num2str(year)],'fontname','arial','fontweight','bold','fontsize',14)
    ylim([min(site_modeled_mass_balance_time_series(:))-.5 max(site_modeled_mass_balance_time_series(:))+.5])
    xlim([date_numbers(1) date_numbers(end)+30])
    xlabel('Date (mm-dd)','fontname','arial ','fontsize',14)
    ylabel('Mass Balance (m w.e.)','fontname','arial ','fontsize',14)
%     ylabel('{\Delta}m (m w.e.)','fontname','arial ','fontsize',14)
    set(gca,'fontname','arial ','fontsize',14,'TickLength',[0.025 0.025],'Ytick',-6:2:2,'YTickLabel',-6:2:2,'linewidth',2)
%     set(gca, 'YGrid','on','XGrid','on','Xtick',-4:2:2, 'XTickLabel',-4:2:2,'Ytick',-4:2:2,'YTickLabel',-4:2:2,...
%             'FontName','Arial','FontSize',14,'LineWidth',2,'FontWeight','Bold', 'Box', 'on');
    box on
    set(gcf, 'PaperPositionMode', 'auto');
    print -depsc2 gates_epoch2.eps
end 
end

