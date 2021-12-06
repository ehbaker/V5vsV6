function [Filled_Glaciological_data] = Profile_Fill_Missing_Glaciological_Observations(glacier,Glaciological_data,years,normal_site,sites,AAD,plot_profiles)
%Profile_Fill_Missing_Glaciological_Observations 
%The first part of this function computes a normal balance profile
%balance values are normalized to a central site, which is designated in
%the first few lines of the function. The normal profile shape is then used
%to fill missing values.

dbstop if error

if isnan(normal_site.Year_1)
    year_1=Glaciological_data.Year(1);
else
    year_1=normal_site.Year_1;
end
if isnan(normal_site.Year_2)
    year_2=Glaciological_data.Year(end);
else
    year_2=normal_site.Year_2;
end
normal_site=normal_site.Normalization_Site;
disp(['%Filling missing data using the average mass balance profile noramalized to site ',....
                cell2mat(normal_site),' between ',....
                num2str(year_1),' and ', num2str(year_2)])

Profile_data = Glaciological_data(Glaciological_data.Year>=year_1 & Glaciological_data.Year<=year_2,:);
bin_centers=AAD(1,2:end);
    
    normalized_ba(:,1) = Profile_data.Year(strcmp(Profile_data.site_name,normal_site));
    normalized_ba(:,2) = Profile_data.ba(strcmp(Profile_data.site_name,normal_site));
    normalized_bw(:,1) = Profile_data.Year(strcmp(Profile_data.site_name,normal_site));
    normalized_bw(:,2) = Profile_data.bw(strcmp(Profile_data.site_name,normal_site));
    dates=[];
 
    for i = 1:height(Profile_data)
        if isbetween(Profile_data.spring_date(i),datetime(strcat('03/15/',num2str(Profile_data.Year(i)))),datetime(strcat('06/15/',num2str(Profile_data.Year(i)))))<1
            Profile_data.bw(i) = NaN;
        end
        if isbetween(Profile_data.fall_date(i),datetime(strcat('08/15/',num2str(Profile_data.Year(i)))),datetime(strcat('10/31/',num2str(Profile_data.Year(i)))))<1
            Profile_data.ba(i) = NaN;
        end
        
        Y = table2array(Profile_data(i,'Year'));
        if ismember(Y,normalized_bw(:,1))
            normalized(i,1) = normalized_bw(normalized_bw(:,1)==Y,2);
            normalized(i,2) = table2array(Profile_data(i, 'bw')) - normalized(i,1);
        else
            normalized(i,1) = NaN;
            normalized(i,2) = NaN;
        end
        if ismember(Y,normalized_ba(:,1))
            normalized(i,3) = normalized_ba(normalized_ba(:,1)==Y,2);
            normalized(i,4) = table2array(Profile_data(i, 'ba')) - normalized(i,3);
        else
            normalized(i,3) = NaN;
            normalized(i,4) = NaN;
        end
        if strcmp(Profile_data.spring_date(i),'NaN')||strcmp(Profile_data.spring_date(i),'Nan')||strcmp(Profile_data.spring_date(i),'nan')&&strcmp(Profile_data.fall_date(i),'NaN')||strcmp(Profile_data.fall_date(i),'Nan')||strcmp(Profile_data.fall_date(i),'nan')
            dates=[dates;NaN NaN];
        elseif strcmp(Profile_data.spring_date(i),'NaN')||strcmp(Profile_data.spring_date(i),'Nan')||strcmp(Profile_data.spring_date(i),'nan')
            dates=[dates;NaN datenum(Profile_data.fall_date(i))];
        elseif strcmp(Profile_data.fall_date(i),'NaN')||strcmp(Profile_data.fall_date(i),'Nan')||strcmp(Profile_data.fall_date(i),'nan')
            dates=[dates;datenum(Profile_data.spring_date(i))  NaN];
        else  
            dates=[dates;datenum(Profile_data.spring_date(i))  datenum(Profile_data.fall_date(i))];
        end
    end
    Profile_data.spring_date=dates(:,1);    
    Profile_data.fall_date=dates(:,2);
    clear dates

    ppw=splinefit(Profile_data.elevation,normalized(:,2),2,2);
    winter_balance_profile_values=ppval(ppw,bin_centers);
    ppa=splinefit(Profile_data.elevation,normalized(:,4),2,2);
    annual_balance_profile_values=ppval(ppa,bin_centers);
    Predicted = table;
    Predicted.site_name = categorical(Profile_data.site_name);
    Predicted.bw = ppval(ppw,Profile_data.elevation) + normalized(:,1);
    Predicted.offset_bw = Predicted.bw - Profile_data.bw;
    Predicted.ba = ppval(ppa,Profile_data.elevation) + normalized(:,3);
    Predicted.offset_ba = Predicted.ba - Profile_data.ba;
    
    Predicted.site_name = categorical(Predicted.site_name);
    fmean = @(x) mean(x,'omitnan');
    meanSiteOffset = varfun(fmean,Predicted,'GroupingVariables','site_name',...
                      'InputVariables',{'offset_bw','offset_ba'});
    
    if plot_profiles == 1
    
    figure();hold on
    subplot(1,2,1)
    scatter(Profile_data.elevation,normalized(:,2),'b','filled');hold on
    plot(bin_centers,winter_balance_profile_values,'LineWidth',1,'color','b');hold on
    title((glacier),'fontname','arial ','fontsize',14,'fontweight','bold')
    ylabel('Winter Balance (m w.e.)')
    xlabel('Elevation (m)')
    axis square
    box on
    set(gcf, 'PaperPosition', [0 0 3 6]);
    print -depsc2 gates_epoch2.eps
    
    subplot(1,2,2)
    scatter(Profile_data.elevation,normalized(:,4),'k','filled');hold on
    plot(bin_centers,annual_balance_profile_values,'LineWidth',1,'color','k')
    ylabel('Annual Balance (m w.e.)')
    xlabel('Elevation (m)')
    axis square
    box on
    set(gcf, 'PaperPosition', [0 0 3 6]);
    print -depsc2 gates_epoch2.eps
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Fill %%%%%%%%%%%%%%%%%
    dates = [];
    for i = 1:height(Glaciological_data)
        if strcmp(Glaciological_data.spring_date(i),'NaN')||strcmp(Glaciological_data.spring_date(i),'Nan')||strcmp(Glaciological_data.spring_date(i),'nan')&&strcmp(Glaciological_data.fall_date(i),'NaN')||strcmp(Glaciological_data.fall_date(i),'Nan')||strcmp(Glaciological_data.fall_date(i),'nan')
            dates=[dates;NaN NaN];
        elseif strcmp(Glaciological_data.spring_date(i),'NaN')||strcmp(Glaciological_data.spring_date(i),'Nan')||strcmp(Glaciological_data.spring_date(i),'nan')
            dates=[dates;NaN datenum(Glaciological_data.fall_date(i))];
        elseif strcmp(Glaciological_data.fall_date(i),'NaN')||strcmp(Glaciological_data.fall_date(i),'Nan')||strcmp(Glaciological_data.fall_date(i),'nan')
            dates=[dates;datenum(Glaciological_data.spring_date(i))  NaN];
        else  
            dates=[dates;datenum(Glaciological_data.spring_date(i))  datenum(Glaciological_data.fall_date(i))];
        end
    end
    Glaciological_data.spring_date=dates(:,1);    
    Glaciological_data.fall_date=dates(:,2);
               
Filled_Glaciological_data = table([],[],[],[],[],[],[],[],[],[],[],'VariableNames',Glaciological_data.Properties.VariableNames);    
   for year=1:length(years)
    insitu_data=Glaciological_data(Glaciological_data.Year==years(year),:);
    if isempty(insitu_data)&&isempty(sites)
        Filled_Glaciological_data=[Filled_Glaciological_data;insitu_data];
    else
    
    if ~isempty(insitu_data) && ~isempty(sites) %if we have some data and a list of sites
        for site=1:length(sites)
            if sum(strcmp(insitu_data.site_name,sites(site)))==0 
                site_elevations=[Glaciological_data.Year(strcmp(Glaciological_data.site_name,sites(site))) Glaciological_data.elevation(strcmp(Glaciological_data.site_name,sites(site)))];
                elevation=interp1(site_elevations(:,1),site_elevations(:,2),years(year),'linear','extrap');
                average_spring_date = nanmedian(insitu_data.spring_date); %since the profile will be fitted to the existing data
                average_fall_date = nanmedian(insitu_data.fall_date);
                insitu_data=[insitu_data;table(years(year),sites(site),average_spring_date,average_fall_date,elevation,nan,nan,nan,nan,0,0,'VariableNames',Glaciological_data.Properties.VariableNames)];   
            end
        end
    end
    
    if sum(isnan(insitu_data.bw))>0 && sum(~isnan(insitu_data.bw))>1%Now fill the balance values
        index = isnan(insitu_data.bw);
        subindex = and(abs(insitu_data.spring_date - nanmedian(insitu_data.spring_date(~index)))<10,~index);
        insitu_data.spring_date(index) = nanmedian(insitu_data.spring_date(subindex));
        if sum(subindex)>1 
            if max(insitu_data.spring_date(subindex)) - min(insitu_data.spring_date(subindex))<10  %if obs within 10 days, estimate with profile
                year_bw = insitu_data.bw(subindex);
                elevations = insitu_data.elevation(subindex);
                temp_bw = ppval(ppw,elevations);
                offset = mean(year_bw - temp_bw);
    
                insitu_data.bw(index) = ppval(ppw,insitu_data.elevation(index)) + offset + meanSiteOffset.Fun_offset_bw(meanSiteOffset.site_name==categorical(cellstr(normal_site)));
                insitu_data.bw_fill(index) = ones(sum(index),1);
            end
        end
    end
    clear index subindex
    
    if sum(isnan(insitu_data.ba))>0 && sum(~isnan(insitu_data.ba))>1
        index = isnan(insitu_data.ba);
        subindex = and(abs(insitu_data.fall_date - nanmedian(insitu_data.fall_date(~index)))<10,~index);
        insitu_data.fall_date(index) = nanmedian(insitu_data.fall_date(subindex));
        if sum(subindex)>1 
            if max(insitu_data.fall_date(subindex)) - min(insitu_data.fall_date(subindex))<10 
                year_ba = insitu_data.ba(subindex);
                elevations = insitu_data.elevation(subindex);
                temp_ba = ppval(ppa,elevations);
                offset = mean(year_ba - temp_ba);
    
                insitu_data.ba(index) = ppval(ppa,insitu_data.elevation(index)) + offset + meanSiteOffset.Fun_offset_ba(meanSiteOffset.site_name==categorical(cellstr(normal_site))); 
                insitu_data.ba_fill(index) = ones(sum(index),1);
            end
        end
    end
    
        Filled_Glaciological_data = [Filled_Glaciological_data;insitu_data];
    end 
   end
   
end

