clear all
close all
addpath functions
glacier='LemonCreek';

glaciological_data=readtable(['data/',glacier,'/Input/Input_',glacier,'_Glaciological_Data.csv']);
weather_data=readtable(['data/',glacier,'/Intermediate/',glacier,'FilledWx.csv']);
weather_data.Date=datenum(weather_data.Date);
precipitation_ratios=readtable(['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Precipitation_Ratios.csv']);
Degree_day_factors=readtable(['data/',glacier,'/Intermediate/Calibrated_',glacier,'_Degree_Day_Factors.csv']);
ks=Degree_day_factors.ks(end);
current_year=str2num(datestr(date,'yyyy'));
insitu_data=glaciological_data(glaciological_data.Year==glaciological_data.Year(end),:);
current_date=datenum(date);

lapse_rate=-5.0;
snow_density=.5
modeled_site_balance=table([],[],[],[],'VariableNames',{'Site','Elevation_m','SWE_m','SnowDepth_m'});

%%
for site=1:height(insitu_data)
    precipitation_ratio=precipitation_ratios.precipitation_ratios(strcmp(precipitation_ratios.site_name,insitu_data.site_name(site)));
    [Site_Weather]=Model_Site_Weather(weather_data,insitu_data.elevation(site),datenum(insitu_data.fall_date(site)),datenum(date),glacier,lapse_rate,precipitation_ratio);
    ablation=zeros(height(Site_Weather),1);
    ablation(Site_Weather.Temperature>0,1)=Site_Weather.Temperature(Site_Weather.Temperature>0)*ks;
    site_cumulative_balance=cumsum(ablation+Site_Weather.Precipitation);
    [~,minimum_index]=min(site_cumulative_balance);
    site_cumulative_balance=site_cumulative_balance-site_cumulative_balance(minimum_index);
    modeled_site_balance=[modeled_site_balance;table(insitu_data.site_name(site),insitu_data.elevation(site),site_cumulative_balance(end),site_cumulative_balance(end)./snow_density,'VariableNames',{'Site','Elevation_m','SWE_m','SnowDepth_m'})];
    
    figure(site);hold on
    plot(Site_Weather.Date,site_cumulative_balance./snow_density,'linewidth',2);hold on
    plot(Site_Weather.Date,4)
    datetick('x','mm/dd')
    title(['Site:',insitu_data.site_name(site)])
    ylabel('Snow Depth (m)')
    xlabel('Date (mm/dd)')
    ylim([0 20])
    box on
    set(gcf, 'PaperPositionMode', 'auto');
    print -depsc2 gates_epoch2.eps
end
