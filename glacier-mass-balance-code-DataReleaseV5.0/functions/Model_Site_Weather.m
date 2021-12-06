function [Site_Weather]=Model_Site_Weather(Weather_data,Elevation,time_one,time_two,glacier,lapse_rate,Precipitation_ratio)
%% Model_Site_Weather.m 
%models site temperature and water equivalent amount of snow using a distal weather station
% and a provided lapse rate and precipiation ratio(site
%precipitation/weather station precipiation), during a given period of time

%% INPUTS:

%Weather_data - table of daily average temperatures(degrees C) and total
%precipitation(m) measured at a distal weather station. Must be longer than
%the period to model(>=time_one &&>=time_two)

%Elevation - Elevation of site to model

%time_one - start of period to model

%time_two - end of period to model

%glacier - name of glacier. currently unused variable

%lapse_rate - provided lapse rate (d*C km^-1)

%Precipiation_ratio - ratio of precipiation between weather station and
%site

%% Outputs:

% Site_Weather - table of modeled site weather(Date, mean daily
% temperature, total daily precipiation

%% begin function
dbstop if error
Site_elevation_differance = Elevation - Weather_data.Elevation(1); %change in elevation between weather station and site
wx_data=Weather_data(Weather_data.Date>=time_one & Weather_data.Date<=time_two,:);%weather data during period to model
Site_Weather=table(wx_data.Date,zeros(height(wx_data),1),zeros(height(wx_data),1),wx_data.T_flag,wx_data.P_flag,'VariableNames',{'Date','Temperature','Precipitation','T_flag','P_flag'});%create table to populate modeled data
minimum_temperature=0; % minimum temp below which dry snow falls rather than rain
maximum_temperature=1.7; % maximum temp above which rain falls, in between wet snow falls
Site_Weather.Temperature=wx_data.T + (Site_elevation_differance.*(lapse_rate/1000)); %model temperature at glaciological site
if sum(Site_Weather.Temperature<maximum_temperature)==0 %all precipiation fell as rain and no snow fell
else
    Site_Weather.Precipitation(Site_Weather.Temperature<maximum_temperature,1)=wx_data.P(Site_Weather.Temperature<maximum_temperature)/1000*Precipitation_ratio; %model water equivalent snow at site
end
end

