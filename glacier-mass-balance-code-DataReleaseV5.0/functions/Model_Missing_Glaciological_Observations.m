function [Filled_Glaciological_data] = Model_Missing_Glaciological_Observations(glacier,Glaciological_data,years,sites,Weather_data,AAD,ks,ki,lapse_rate,precipitation_ratios_table,integration_method,integration_surface,plot_ablation_model)
%% Model_Missing_Glaciological_Ovservations.m
%function uses a mass balance model to fill in missing point balances in
%the input file. The mass balance model is a combined temperature index and
%precipitation. Mass balance model coefficients are calibrated in Calibrate_Precipitation_and_Ablation Models.m 
%and used here to scale local weather station data to a given site on the glacier

%% General outline 
% 1) Date formatting and data sorting -- Convert all dates (missing or real) to datenumbers and contruct a
% table the conatins all sites the user defined to be modeled

% 2) Determine missing dates and point balances  -- Using the
% Find_Mass_Maximum_and_Minimum_Adjustments.m function predict the mass max
% and minimum dates for sites not observed. Then using the modeled mass
% balance timeseries from that function, predict the missing point balance
% values.

%% Input 

%glacier -- Name of glacier user is working on

%Glaciological_data -- Table of glaciological point balances. User can
%   insert nan values for any site they desire to be filled for a specific
%   year

%years -- List of years the user want data to be filled during

%sites -- List of site the user wants to be filled. If provided then this
%   function will model the listed sites for every year. If sites is empty
%   the function will only model a site balance if it contains a nan value

%Weather_data -- Table of daily average temperatures (*C) and total daily
%   precipitation (mm)

%ks -- Calibrated degree day coefficient for snow

%ki -- Calibrated degree day coefficient for ice

%lapse_rate -- Lapse rate

%precipitation_ratios_table -- Table of calibrated precipitation ratios
%   between weather station and glaciological sites

%integration_method -- Integration method (Index, Linear balance gradient,
%   Piecewise balance profile
%integration_surface -- Integration surface (Conventional, Reference
%   Surface, or specific year hypsometry)
%plot_ablation_model -- Flag to display mass balance model or not

%% Output
%Filled_Glaciological_data -- Table of glaciological data, with missing
%   point balances filled using mass balance model

% prime engines--all systems go captain
%% 1) Date formatting and data sorting 
dbstop if error
time_system=1; %set time system to stratigraphic for modeling missing point balances
dates=[];% empty date vector to be filled in the next ten lines

for i=1:height(Glaciological_data)% for each entry in the glaciological data table we need to deal with NaN values and convert all dates to datenumbers
    %if both dates are NaNs
    if strcmp(Glaciological_data.spring_date(i),'NaN')||strcmp(Glaciological_data.spring_date(i),'Nan')||strcmp(Glaciological_data.spring_date(i),'nan')&&strcmp(Glaciological_data.fall_date(i),'NaN')||strcmp(Glaciological_data.fall_date(i),'Nan')||strcmp(Glaciological_data.fall_date(i),'nan')
       dates=[dates;NaN NaN];
    %elseif fall dates are NaNs
    elseif strcmp(Glaciological_data.spring_date(i),'NaN')||strcmp(Glaciological_data.spring_date(i),'Nan')||strcmp(Glaciological_data.spring_date(i),'nan')
       dates=[dates;NaN datenum(Glaciological_data.fall_date(i))];
    %elseif spring dates are NaNs
    elseif strcmp(Glaciological_data.fall_date(i),'NaN')||strcmp(Glaciological_data.fall_date(i),'Nan')||strcmp(Glaciological_data.fall_date(i),'nan')
        dates=[dates;datenum(Glaciological_data.spring_date(i))  NaN];
    % else both spring and fall dates are entered
    else  
        dates=[dates;datenum(Glaciological_data.spring_date(i))  datenum(Glaciological_data.fall_date(i))];
    end
end
if isempty(years)
    years=unique(Glaciological_data.Year);
end
%add date numbers back into glaciological data table
Glaciological_data.spring_date=dates(:,1);    
Glaciological_data.fall_date=dates(:,2);

%create table to put filled glaciological data back into
Filled_Glaciological_data=table([],[],[],[],[],[],[],[],[],[],[],'VariableNames',Glaciological_data.Properties.VariableNames);

%liste of sites that are not transient snow line observations. TSLs can't
%be modeled because they are non-stationary
sites=sites(~contains(sites,'TSL'));

%for each year look at the glaciological table and determine and sites
%listed by user to be filled are in the table for a given year. If a site
%is not listed for a given year, add it to the table
for year=1:length(years)
    insitu_data=Glaciological_data(Glaciological_data.Year==years(year),:);%insitu data for the year
    if isempty(insitu_data)&&isempty(sites) %if there is no entry in the glaciological data for the year and not sites listed than nothing needs to be done.
        Filled_Glaciological_data=[Filled_Glaciological_data;insitu_data];
        
    else %otherwise we need to populate the table for the desired sites to be modeled
        if isempty(insitu_data)&&~isempty(sites) %no data in table, but sites are listed, so we need to add them to the table
            filled_year=ones(length(sites),1).*years(year); %years
            filled_sites=sites; %sites
            filler_nans=ones(length(sites),1).*nan; %vector of nans
            filler_dates=[ones(length(sites),1).*datenum(['5/1/',num2str(years(year))]),ones(length(sites),1).*datenum(['9/30/',num2str(years(year))])]; %approximate mass extreme dates
            filler_fill_tag = zeros(length(sites),1); %vector for flags stating filling method
            site_elevations=[];
            for site=1:length(sites)%we dont know the site elevations so we need to interpolate from the measured elevations
                site_elevations=[Glaciological_data.Year(strcmp(Glaciological_data.site_name,sites(site))) Glaciological_data.elevation(strcmp(Glaciological_data.site_name,sites(site)))];
                elevations(site,1)=interp1(site_elevations(:,1),site_elevations(:,2),years(year),'linear','extrap');
            end
            %created insitu data table to be modeled
            insitu_data=table(filled_year,filled_sites,filler_dates(:,1),filler_dates(:,2),elevations,filler_nans,filler_nans,filler_nans,filler_nans,filler_fill_tag,filler_fill_tag,'VariableNames',Glaciological_data.Properties.VariableNames);

        elseif ~isempty(sites)% else sites are listed and there is some data for the year
            for site=1:length(sites) %for each site listed
                if sum(strcmp(insitu_data.site_name,sites(site)))==0 %determine if the site is not listed in the glaciological data for the year
                    site_elevations=[Glaciological_data.Year(strcmp(Glaciological_data.site_name,sites(site))) Glaciological_data.elevation(strcmp(Glaciological_data.site_name,sites(site)))];
                    elevation=interp1(site_elevations(:,1),site_elevations(:,2),years(year),'linear','extrap'); %interpolate site elevation
                    insitu_data=[insitu_data;table(years(year),sites(site),datenum(['5/1/',num2str(years(year))]),datenum(['9/30/',num2str(years(year))]),elevation,nan,nan,nan,nan,0,0,'VariableNames',Glaciological_data.Properties.VariableNames)];

                end
            end
        end




    spring_dates=insitu_data.spring_date;
    fall_dates=insitu_data.fall_date;
    
    if sum(isnan(insitu_data.spring_date))>0
        spring_dates=insitu_data.spring_date;
        insitu_data.spring_date(isnan(insitu_data.spring_date))=datenum(['5/1/',num2str(years(year))]);
    end
    if sum(isnan(insitu_data.fall_date))>0
        fall_dates=insitu_data.fall_date;
        insitu_data.fall_date(isnan(insitu_data.fall_date))=datenum(['9/30/',num2str(years(year))]);
    end
    
    
%% 2) Determine missing dates and point balances    
    
    % Use the Find_Mass_Maximum_and_Minimum.m function to get mass extreme
    % dates and modeled mass balance time series for the year
    [~,mass_maximum_date_numbers,~,mass_minimum_date_numbers,Modeled_Mass_Balance] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);
    
    if sum(isnan(spring_dates))>0 %if there is a missing spring date
        insitu_data.spring_date(isnan(spring_dates))=mass_maximum_date_numbers(isnan(spring_dates)); %use spring mass extreme date
    end
    if sum(isnan(fall_dates))>0 %if there is a missing fall date
        insitu_data.fall_date(isnan(fall_dates))=mass_minimum_date_numbers(isnan(fall_dates)); %use fall mass extreme date
    end
    
    %convert dates in the models mass balance time-series for sites to date
    %number
    Modeled_Mass_Balance.Date=datenum(Modeled_Mass_Balance.Date); 
    
    %for each site in the insitu data determine if there is data that needs
    %to be filled
    for site=1:height(insitu_data)            
        
        %if there is a fall measurement but not a spring measurement
        if isnan(insitu_data.bw(site)) && ~isnan(insitu_data.ba(site))
            
            % if fall measurement happened before the mass minimum date
            if mass_minimum_date_numbers(site)>=insitu_data.fall_date(site)      
                %remove any bias between modeled site mass balance
                %time-series and measured fall balance 
                Modeled_Mass_Balance(:,strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site)))=array2table((insitu_data.ba(site)-table2array(Modeled_Mass_Balance(Modeled_Mass_Balance.Date==insitu_data.fall_date(site),....
                    strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site)))))+table2array(Modeled_Mass_Balance(:,strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site)))));
                
                %extract the winter balance from the modeled site mass
                %balance time-series
                insitu_data.bw(site)=table2array(Modeled_Mass_Balance(Modeled_Mass_Balance.Date==mass_maximum_date_numbers(site),....
                    strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))));
                
                %set spring date to spring mass extreme date
                insitu_data.spring_date(site)=mass_maximum_date_numbers(site);
                
                %flag this value has been filled using the mass balance
                %model
                insitu_data.bw_fill(site) = 2;
                
            % fall measurement date is after the mass minimum date
            else 
                %remove any bias between modeled site mass balance
                %time-series and measured fall balance 
                Modeled_Mass_Balance(:,strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site)))=array2table((insitu_data.ba(site)-table2array(Modeled_Mass_Balance(Modeled_Mass_Balance.Date==mass_minimum_date_numbers(site),....
                    strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site)))))+table2array(Modeled_Mass_Balance(:,strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site)))));
                
                %extract the winter balance from the modeled site mass
                %balance time-series
                insitu_data.bw(site)=table2array(Modeled_Mass_Balance(Modeled_Mass_Balance.Date==mass_maximum_date_numbers(site),....
                    strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))));
                
                %set spring date to spring mass extreme date
                insitu_data.spring_date(site)=mass_maximum_date_numbers(site);
                
                %flag this value has been filled using the mass balance
                %model
                insitu_data.bw_fill(site) = 2;
            end
            
        %elseif there is a spring measurement but not a fall measurement
        elseif ~isnan(insitu_data.bw(site)) && isnan(insitu_data.ba(site))                 
                %remove any bias between modeled site mass balance
                %time-series and measured spring balance
                Modeled_Mass_Balance(:,strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site)))=array2table((insitu_data.bw(site)-table2array(Modeled_Mass_Balance(Modeled_Mass_Balance.Date==insitu_data.spring_date(site),....
                    strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site)))))+table2array(Modeled_Mass_Balance(:,strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site)))));
                %extract the fall balance from the modeled site mass
                %balance time-series
                insitu_data.ba(site)=table2array(Modeled_Mass_Balance(Modeled_Mass_Balance.Date==mass_minimum_date_numbers(site),....
                    strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))));
                
                %set fall date to fall mass extreme date
                insitu_data.fall_date(site)=mass_minimum_date_numbers(site);
                
                %flag this value has been filled using the mass balance
                %model
                insitu_data.ba_fill(site) = 2;
        end
    end
    
    %if a site is missing both spring and fall measurements we need to
    %model the entire mass balance year from the previous years mass
    %minimum
    if sum(isnan(insitu_data.bw))>0 && sum(isnan(insitu_data.ba))>0
        
        %biuld synthetic data from previous year to get the mass extreme
        %dates
        previous_fall_observation_dates=datenum(['9/30/',num2str(years(year)-1)],'mm/dd/yyyy').*ones(height(insitu_data),1);
        previous_spring_observation_dates=datenum(['5/1/',num2str(years(year)-1)],'mm/dd/yyyy').*ones(height(insitu_data),1);
        previous_year=(years(year)-1).*ones(height(insitu_data),1);
        previous_insitu_data=insitu_data;
        previous_insitu_data.fall_date=previous_fall_observation_dates;
        previous_insitu_data.spring_date=previous_spring_observation_dates;
        previous_insitu_data.Year=previous_year;
        [~,~,~,previous_mass_minimum_date_numbers] = Find_Mass_Maximum_and_Minimum_Adjustments(previous_insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);
        
        %for each site in insitu data
        for site=1:height(insitu_data)
            %if we are missing both spring and fall measurements
            if isnan(insitu_data.bw(site)) && isnan(insitu_data.ba(site))
                precipitation_ratio=precipitation_ratios_table.precipitation_ratios(strcmp(precipitation_ratios_table.site_name,insitu_data.site_name(site))); %get precipitation ratio for site
                if isempty(precipitation_ratio)||isnan(precipitation_ratio)%then no precipitation ratio could be calibrated for this site and we need to find one to use
                    precipitation_ratio=precipitation_ratios_table.precipitation_ratios(strcmp(precipitation_ratios_table.site_name,'All_sites'));%so we use a global fit
                end
                
                %model site weather from the prvious years mass minimum and
                %the following years mass minimum
                [Site_Weather]=Model_Site_Weather(Weather_data,insitu_data.elevation(site),previous_mass_minimum_date_numbers(site),mass_minimum_date_numbers(site),glacier,lapse_rate,precipitation_ratio);
                
                site_balance=0;%since we are starting at the stratigraphic mass minimum, restart the clock
                
                %for each day of the mass balance year at the site
                for day=1:height(Site_Weather)
                    
                    %add any accumulation to the site balance
                    site_balance=site_balance+Site_Weather.Precipitation(day,1);
                    
                    %determine if any ablation occured, and if so, was is it snow or ice 
                    if Site_Weather.Temperature(day,1)>=0 && site_balance>0 %site balance greater than zero means we are melting snow
                        ablation=Site_Weather.Temperature(day,1).*ks;
                    elseif Site_Weather.Temperature(day,1)>=0 && site_balance<=0 %site balance less than zero, we are likely melting ice
                        ablation=Site_Weather.Temperature(day,1).*ki;
                    else
                        ablation=0;
                    end
                    %add ablation to site balance
                    site_balance=site_balance+ablation;
                    
                    %if we are at the spring mass max date, then we have
                    %our spring mass balance value
                   if Site_Weather.Date(day)==mass_maximum_date_numbers(site)
                       insitu_data.bw(site)=site_balance;
                       
                    %else if we are at the fall mass min date, then we have
                    %of fall mass balance value
                   elseif Site_Weather.Date(day)==mass_minimum_date_numbers(site)
                       insitu_data.ba(site)=site_balance;
                   end
                end
                
                %finalize spring and fall mass balance dates
                insitu_data.spring_date(site)=mass_maximum_date_numbers(site);
                insitu_data.fall_date(site)=mass_minimum_date_numbers(site);
                
                %flag that these values have been filled using the mass
                %balance model
                insitu_data.bw_fill(site) = 2;
                insitu_data.ba_fill(site) = 2;
            end
        end 
    end
    
    %if there are any sites with identical elevations, make them
    %unsignificantly different, so we don't run into errors later down the
    %road due to non-unique elevation values
    n=0;   
    for site=1:height(insitu_data)        
        index=find(insitu_data.elevation==insitu_data.elevation(site));
        if length(index)>1
            n=n+0.01;
            insitu_data.elevation(site)=insitu_data.elevation(site)+n;
        end
    end

    %append insitu data for the year to the final filled glaciological data table        
    Filled_Glaciological_data=[Filled_Glaciological_data;insitu_data];
    end
end
end

