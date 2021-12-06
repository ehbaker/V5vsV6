function [ready]=fillWxData(glacier)
% this function fixes the missing glacier temperature and precip data by
% using the closest city data and estimating regression functions between
% the city data and the glacier data.
%
dbstop if error
warning off MATLAB:table:RowsAddedExistingVars
warning off stats:statrobustfit:IterationLimit


%ins

%formatSecondaryWxData(glacier); %correct/reformat data from cities
%first format secondary source of weather data. This will be used to fill
%any data gaps in the primary weather station used in mass balance
%analysis.
    
Secondary = readtable(['data/',glacier,'/Input//Input_',glacier,'_Secondary_Weather_data.csv']);
Secondary.Maximum_Temperature(Secondary.Maximum_Temperature==9999)=nan; %replace no data value with nan
Secondary.Minimum_Temperature(Secondary.Minimum_Temperature==9999)=nan;
Secondary.Precipitation(Secondary.Precipitation==9999)=nan;
%We need the daily mean temperature at the secondary weather station,
%this is defined
%find any data in US units and convert to SCI units.
Temperature(strcmp(Secondary.Units,'US'),1) = ((Secondary.Maximum_Temperature(strcmp(Secondary.Units,'US')) - 32)*(5/9) + (Secondary.Minimum_Temperature(strcmp(Secondary.Units,'US')) - 32)*(5/9))/2;
Precipitation(strcmp(Secondary.Units,'US'),1) = Secondary.Precipitation(strcmp(Secondary.Units,'US')) * 25.4;
%find data in SCI units
Temperature(strcmp(Secondary.Units,'SI'),1) = (Secondary.Maximum_Temperature(strcmp(Secondary.Units,'SI')) + Secondary.Minimum_Temperature(strcmp(Secondary.Units,'SI')))/2;
Precipitation(strcmp(Secondary.Units,'SI'),1) = Secondary.Precipitation(strcmp(Secondary.Units,'SI'));

SecondaryWx = table(Secondary.Date,Temperature,Precipitation,Secondary.Source,'VariableNames',{'Date' 'Temperature' 'Precipitation' 'Source'});

PrimaryWx=readtable(['data/',glacier,'/Input/Input_',glacier,'_Daily_Weather.csv']); %Weather from the nearest Wx station
MissingPrecipitation=readtable(['data/',glacier,'/Input/Input_',glacier,'_Missing_Precipitation.csv']); % Manual Precip totals from data logger failure dates

PrimaryWx=table(PrimaryWx.Date,PrimaryWx.Temperature,PrimaryWx.Precipitation,PrimaryWx.Elevation,'VariableNames',{'Date','Temperature','Precipitation','Elevation'});
if datenum(PrimaryWx.Date(end))<datenum(['9/30/',num2str(str2num(datestr(date,'yyyy'))+2)])
    PrimaryWx.Date=datenum(PrimaryWx.Date);
    dates_needed(:,1)=PrimaryWx.Date(end)+1:datenum(['9/30/',num2str(str2num(datestr(date,'yyyy'))+2)]);
    extend_table=table(dates_needed,nan*ones(length(dates_needed),1),nan*ones(length(dates_needed),1),PrimaryWx.Elevation(1)*ones(length(dates_needed),1),'VariableNames',{'Date','Temperature','Precipitation','Elevation'});
    PrimaryWx=[PrimaryWx;extend_table];
end
SecondaryWx=SecondaryWx(datenum(SecondaryWx.Date)>=PrimaryWx.Date(1),:);
PrimaryWx.Date = datetime(datestr(PrimaryWx.Date));
SecondaryWx.Date = datetime(SecondaryWx.Date);
AllWx = outerjoin(PrimaryWx, SecondaryWx,'Keys','Date','MergeKeys',1); %make 1 big table organized by date
AllWx.Properties.VariableNames = {'Date' 'T_primary' 'P_primary' 'Elevation' 'T_secondary' 'P_secondary' 'Source'};


%%Fill short gaps based on adjacent temperatures
Tnancount(1) = sum(isnan(AllWx.T_primary));
filledT = interp1gap(AllWx.T_primary,3)';% linear interpolation of temperature through small (<3) gaps 
Tnancount(2) = sum(isnan(filledT));

%%Calculate monthly regressions and fill
MonthlyTemperatureRegressions = table;
figure (); hold on 
title(['Monthly temperature regressions ',glacier, 'Glacier'])
monthofyear = month(AllWx.Date,'monthofyear');
Secondary_Sources=unique(AllWx.Source(~strcmp(AllWx.Source,'')));

for source=1:length(Secondary_Sources)
    color=[0 0 0];
    color(1,source)=1;
    if source==4
        fprintf(1,'ERROR: We want to minimize the number of changes in secondary data to three maximum. Less is better!')
    end
    for m = 1:12
        %%iteratively reweighted least squares regression
        X = AllWx.T_secondary(monthofyear==m & strcmp(AllWx.Source,Secondary_Sources(source)));
        Y = AllWx.T_primary(monthofyear==m & strcmp(AllWx.Source,Secondary_Sources(source)));
        
        
        b = robustfit(X,Y);

        %%store the coefficients in a table
        MonthlyTemperatureRegressions.month(m,1) = m;
        MonthlyTemperatureRegressions.intercept(m,1) = b(1);
        MonthlyTemperatureRegressions.slope(m,1) = b(2);
        MonthlyTemperatureRegressions.rsquared(m,1) = round((corr(Y(~isnan(Y) & ~isnan(X)),b(1) + b(2)*X(~isnan(Y) & ~isnan(X)))^2)*100)/100;


        %%plot the regression
        subplot(3,4,m); hold on
        scatter(X,Y,1,color);hold on
        lsline

        text(-28,27, unique(month(AllWx.Date(monthofyear==m),'name')))
        text(-28,(27-source*10),[cell2mat(Secondary_Sources(source)),' r^2 = ',num2str(MonthlyTemperatureRegressions.rsquared(m))], 'color', color);hold on
        axis([-30 30 -30 30])

        %%fill gaps
        filledT(monthofyear==m & isnan(filledT) & strcmp(AllWx.Source,Secondary_Sources(source))) = b(1) + AllWx.T_secondary(monthofyear==m & isnan(filledT) & strcmp(AllWx.Source,Secondary_Sources(source))).*b(2);
    end
    writetable(MonthlyTemperatureRegressions,['data/',glacier,'/Intermediate/',glacier,'_',cell2mat(Secondary_Sources(source)),'_TemperatureRegressions.csv'])
end

Tnancount(3) = sum(isnan(filledT));


%%Fill the remaining NaNs with the average daily temperature
dayofyear = day(AllWx.Date,'dayofyear');
meanDailyT = nan(366,1);
for day_number = 1:366
    meanDailyT(day_number) = nanmean(AllWx.T_primary(dayofyear==day_number));
    filledT(dayofyear==day_number & isnan(filledT)) = meanDailyT(day_number);
end
Tnancount(4) = sum(isnan(filledT));
 
%%Precipitation 
filledP = AllWx.P_primary;
Pnancount(1) = sum(isnan(filledP));

%%Calculate monthly precipitation regressions and fill
MonthlyPrecipitationRegressions = table;
figure (); hold on 
title(['Monthly precipitation regressions ',glacier, 'Glacier'])
for source=1:length(Secondary_Sources)
    color=[0 0 0];
    color(1,source)=1;
    for m = 1:12
        %%iteratively reweighted least squares regression
        X = AllWx.P_secondary(monthofyear==m & AllWx.P_primary>=0.001 & strcmp(AllWx.Source,Secondary_Sources(source)));
        Y = AllWx.P_primary(monthofyear==m & AllWx.P_primary>=0.001 & strcmp(AllWx.Source,Secondary_Sources(source)));
        b = robustfit(X,Y,'bisquare',4.685,'off'); %removes the intercept term so that 0 precip secondary = 0 precip primary
        if b==0 b = NaN; end

       %%store the coefficients in a table
        MonthlyPrecipitationRegressions.month(m,1) = m;
        MonthlyPrecipitationRegressions.slope(m,1) = b;
        MonthlyPrecipitationRegressions.rsquared(m,1) = round((corr(Y(~isnan(Y) & ~isnan(X)),b*X(~isnan(Y) & ~isnan(X)))^2)*100)/100;

        %%plot the regression
        subplot(3,4,m); hold on
        scatter(X,Y,1,color);hold on
        lsline
        text(10,220,unique(month(AllWx.Date(monthofyear==m),'name')))
        text(10,(220-(40*source)),[cell2mat(Secondary_Sources(source)),' r^2 = ' num2str(MonthlyPrecipitationRegressions.rsquared(m),2)], 'color', color)
        axis([0 250 0 250])

        %%fill gaps

        filledP(monthofyear==m & isnan(filledP) & strcmp(AllWx.Source,Secondary_Sources(source))) = AllWx.P_secondary(monthofyear==m & isnan(filledP) & strcmp(AllWx.Source,Secondary_Sources(source))).*b;
    end
end

Pnancount(2) = sum(isnan(filledP));
writetable(MonthlyPrecipitationRegressions,['data/',glacier,'/Intermediate/',glacier,'PrecipitationRegressions.csv'])

%%Fill the remaining NaNs with the average daily precip
dayofyear = day(AllWx.Date,'dayofyear');
meanDailyP = nan(366,1);
for day_number = 1:366
    meanDailyP(day_number) = nanmean(AllWx.P_primary(dayofyear==day_number));
    filledP(dayofyear==day_number & isnan(filledP)) = meanDailyP(day_number);
end
Pnancount(3) = sum(isnan(filledP));

%use available data to correct some of the missing values
if isempty(MissingPrecipitation)
else   
MissingPrecipitation.start_date = datetime(MissingPrecipitation.start_date);
MissingPrecipitation.end_date = datetime(MissingPrecipitation.end_date);
for n = 1:height(MissingPrecipitation) 
    start = find(AllWx.Date==MissingPrecipitation.start_date(n));
    finish = find(AllWx.Date==MissingPrecipitation.end_date(n));
    estP = sum(filledP(start:finish));
    adj = MissingPrecipitation.P(n)/estP;
    filledP(start:finish) = filledP(start:finish).*adj;
end
end
FilledWx = table;
FilledWx.Date = AllWx.Date;
FilledWx.T = filledT;
FilledWx.P = filledP;
FilledWx.T_flag=isnan(AllWx.T_primary);
FilledWx.P_flag=isnan(AllWx.P_primary);
FilledWx.Elevation=AllWx.Elevation;
writetable(FilledWx,['data/',glacier,'/Intermediate/',glacier,'FilledWx.csv'])
% if Tnancount(4)==0 && Pnancount(3)==0
ready=1;
disp('      %%%%%%   Weather Data has been filled    %%%%%%%    ')
end
