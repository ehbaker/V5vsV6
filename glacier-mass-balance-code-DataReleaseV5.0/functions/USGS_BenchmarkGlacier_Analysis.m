function [final_point_balances,final_geodetic_data,Final_Glacier_Wide_solutions] = USGS_BenchmarkGlacier_Analysis(glacier,all_sites,years,Glaciological_data,Weather_data,AAD,ks,ki,precipitation_ratios_table,lapse_rate,time_system,integration_surface,integration_method,nan_incomplete_glaciological_data,Geodetic_Calibration_index,plot_integration,plot_ablation_model,Include_Previous_Glacierwide_Solutions)
%UNTITLED2 Summary of this function goes here,Include_Previous_Glacierwide_Solutions
%   Detailed explanation goes here
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
[final_point_balances,Glacier_Wide_solutions_table]=Calculate_GlacierWide_Balance(glacier,all_sites,years,Glaciological_data,Weather_data,AAD,ks,ki,precipitation_ratios_table,lapse_rate,time_system,integration_surface,integration_method,plot_integration,plot_ablation_model,nan_incomplete_glaciological_data)
Lapse_Rate_used=cell(height(final_point_balances),1);
Lapse_Rate_used(:,1)={lapse_rate};
time_systems_used=cell(height(final_point_balances),1);
time_systems_used(:,1)={time_system_string};
integration_method_used=cell(height(final_point_balances),1);
integration_method_used(:,1)={integration_method_string};
Integration_surface_used=cell(height(final_point_balances),1);
Integration_surface_used(:,1)={Integration_surface_string};

final_point_balances=[final_point_balances table(Lapse_Rate_used,time_systems_used,integration_method_used,Integration_surface_used,'VariableNames',{'Lapse_Rate','Time_System','Integration_Method','Integration_Surface'})];
writetable(final_point_balances,['data/',glacier,'/Output/','Output_',glacier,'_Adjusted_Point_Balances.csv']);

Lapse_Rate_used=cell(height(Glacier_Wide_solutions_table),1);
Lapse_Rate_used(:,1)={lapse_rate};
time_systems_used=cell(height(Glacier_Wide_solutions_table),1);
time_systems_used(:,1)={time_system_string};
integration_method_used=cell(height(Glacier_Wide_solutions_table),1);
integration_method_used(:,1)={integration_method_string};
Integration_surface_used=cell(height(Glacier_Wide_solutions_table),1);
Integration_surface_used(:,1)={Integration_surface_string};

if sum(isnan(Glacier_Wide_solutions_table.Ba_mwe))>=1
    for site=1:length(all_sites)
        number_of_observations(site,1)=sum(strcmp(final_point_balances.site_name,all_sites(site)));
    end
    [max_site_observations,max_site_observations_index]=max(number_of_observations(:,1));
    if max_site_observations==height(Glacier_Wide_solutions_table)
        bas=final_point_balances.ba(strcmp(final_point_balances.site_name,all_sites(max_site_observations_index)));
        bws=final_point_balances.bw(strcmp(final_point_balances.site_name,all_sites(max_site_observations_index)));
        ba_Ba_lm=fitlm(bas(~isnan(Glacier_Wide_solutions_table.Ba_mwe)),Glacier_Wide_solutions_table.Ba_mwe(~isnan(Glacier_Wide_solutions_table.Ba_mwe)))
        bw_Bw_lm=fitlm(bws(~isnan(Glacier_Wide_solutions_table.Bw_mwe)),Glacier_Wide_solutions_table.Bw_mwe(~isnan(Glacier_Wide_solutions_table.Bw_mwe)))
        ba_Ba_CI=coefCI(ba_Ba_lm);
        figure();hold on
        scatter(bas(~isnan(Glacier_Wide_solutions_table.Ba_mwe)),Glacier_Wide_solutions_table.Ba_mwe(~isnan(Glacier_Wide_solutions_table.Ba_mwe)),100,[1 .5 0],'filled');hold on
        lsline
        text(-3,1,['R^2=',num2str(round(ba_Ba_lm.Rsquared.Ordinary,2))],'fontname','arial ','fontsize',14,'fontweight','bold')
        set(gca, 'YGrid','on','XGrid','on','Xtick',-4:2:2, 'XTickLabel',-4:2:2,'Ytick',-4:2:2,'YTickLabel',-4:2:2,...
            'FontName','Arial','FontSize',14,'LineWidth',2,'FontWeight','Bold', 'Box', 'on');
        ylabel('B_a (m w.e.)','fontname','arial ','fontsize',14,'fontweight','bold')
        xlabel('b_a (m w.e.)','fontname','arial ','fontsize',14,'fontweight','bold');
        box on
        axis square
        set(gcf, 'PaperPositionMode', 'auto');
        print -depsc2 gates_epoch2.eps
    end
    for year=1:height(Glacier_Wide_solutions_table)         
         if isnan(Glacier_Wide_solutions_table.Ba_mwe(year))
            point_balance=final_point_balances.ba(final_point_balances.Year==Glacier_Wide_solutions_table.Year(year)&strcmp(final_point_balances.site_name,all_sites(max_site_observations_index)));
            Glacier_Wide_solutions_table.Ba_mwe(year)=predict(ba_Ba_lm,point_balance);
         end
         if isnan(Glacier_Wide_solutions_table.Bw_mwe(year))
            point_balance=final_point_balances.bw(final_point_balances.Year==Glacier_Wide_solutions_table.Year(year)&strcmp(final_point_balances.site_name,all_sites(max_site_observations_index)));
            Glacier_Wide_solutions_table.Bw_mwe(year)=predict(bw_Bw_lm,point_balance);   
         end
    end
end
Glacier_Wide_solutions_table=[Glacier_Wide_solutions_table table(Lapse_Rate_used,time_systems_used,integration_method_used,Integration_surface_used,'VariableNames',{'Lapse_Rate','Time_System','Integration_Method','Integration_Surface'})];
s=cell(height(Glacier_Wide_solutions_table),1);
s(1:height(Glacier_Wide_solutions_table),1)={'Reanalyzed'};
Glacier_Wide_solutions_table.Source=s;
Glacier_Wide_solutions_table.Bs_mwe=Glacier_Wide_solutions_table.Ba_mwe-Glacier_Wide_solutions_table.Bw_mwe
% writetable(Glacier_Wide_solutions_table,['data/',glacier,'/Output/','Output_',glacier,'_Glacier_Wide_solutions_uncalibrated.csv']);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   D) Calculated Glacier-wide Mass balance
%%      3) Geodetic Calibration of Glaciological Time Series
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%         a) Determine geodetic correction selected and import geodetic data
% Geodetic_Calibration_index=get(handles.Geodetic_Calibration_Selection,'Value');
if Geodetic_Calibration_index==4
    a=readtable(['data/',glacier,'/Input/Input_',glacier,'_Geodetic_BreakPoints.csv']);
    break_point=a.Breakpoints';
else
    break_point=NaN;
end
geodetic_data=readtable(['data/',glacier,'/Input/Input_',glacier,'_Geodetics.csv']);
density=0.85;
%%         b) Account for varying aquisition dates of geodetic measurements using mass balance model         
geodetic_data.dh_m=geodetic_data.dh_m.*density; %convert mean surface elevation change to mean mass change using assumed density by Huss et al., 2013
geodetic_data.Uncertainty=geodetic_data.Uncertainty.*density;
geodetic_data.Properties.VariableNames(strcmp(geodetic_data.Properties.VariableNames,'dh_m'))={'Mass_Change_mwe'};%change variable name to reflect change to m we
geodetic_data.Mass_Change_mwe=geodetic_data.Mass_Change_mwe(1)-geodetic_data.Mass_Change_mwe;
%%%%%%%%%%%
%find ablation correction for geodetic time-series
if Geodetic_Calibration_index ~= 1%a geodetic calibration has been selected
    Ablation_adjustments_for_geodetics=table([],[],[],[],'VariableNames',{'Aquistion_Date' 'Mass_Minimum_Date' 'Number_of_Days' 'Mass_Difference'});
    Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions=Glacier_Wide_solutions_table;
    all_geodetic_data=geodetic_data;
    geodetic_data=geodetic_data(strcmp(geodetic_data.Source,'DEM'),:);
    for measurement=1:height(geodetic_data)%for each geodetic measurement
        
            measurement_year=str2num(datestr(geodetic_data.Date(measurement),'yyyy'));%get the year of data
            aquisition_date=geodetic_data.Date(measurement);
            if datenum(geodetic_data.Date(measurement))==datenum('02/11/2000','mm/dd/yyyy')%Aquisition date for SRTM
                measurement_year=measurement_year-1;%Due to snow penetration of x-band SAR, we assume DEM surfce is representative of 1999 mass minimum surface
                aquisition_date=datetime('09/30/1999');
            end

            if ~isempty(Glaciological_data(Glaciological_data.Year==measurement_year,:))% && isempty(data)
                insitu_data=Glaciological_data(Glaciological_data.Year==measurement_year,:);
                balance_sites_for_year=insitu_data.site_name;
                data=[];
            elseif isempty(Glaciological_data(Glaciological_data.Year==measurement_year,:))% && isempty(data)%then no final data exists for this year
                %and we need to model glaciological data
                insitu_data = Model_Missing_Glaciological_Observations(glacier,Glaciological_data,measurement_year,all_sites(~strcmp(all_sites,'TSL')),Weather_data,AAD,ks,ki,lapse_rate,precipitation_ratios_table,integration_method,integration_surface,plot_ablation_model);
                data=[];
                balance_sites_for_year=all_sites(~contains(all_sites,'TSL'));
            end
%           1) Apply 'Find_Mass_Maximum_and_Minimum_Adjustments.m' to find
    %           the change in mass balance at glaciological site for year
    %           of aquisition        
                
                [mass_maximum_site_adjustments,mass_maximum_date_numbers,mass_minimum_site_adjustments,mass_minimum_date_numbers] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,aquisition_date,integration_method,integration_surface,0);
                ba=insitu_data.ba+mass_minimum_site_adjustments';
                data=[];
                for site=1:length(balance_sites_for_year)%need to great final point balances table for integrate point balances function
                    data=[data; table(insitu_data.Year(strcmp(insitu_data.site_name,balance_sites_for_year(site))),insitu_data.site_name(strcmp(insitu_data.site_name,balance_sites_for_year(site))),insitu_data.elevation(strcmp(insitu_data.site_name,balance_sites_for_year(site))),...
                        nan,nan,ba(strcmp(insitu_data.site_name,balance_sites_for_year(site))),nan,{datestr(mass_minimum_date_numbers(strcmp(insitu_data.site_name,balance_sites_for_year(site))))},'VariableNames',{'Year' 'site_name' 'elevation' 'bw' 'bs' 'ba' 'mass_maximum_date' 'mass_minimum_date'})];
                end
                data=data(data.elevation>=min(AAD(1,2:end)-(AAD(1,3)-AAD(1,2)))&data.elevation<max(AAD(1,2:end)+(AAD(1,3)-AAD(1,2))),:);
                data=sortrows(data,'elevation');
                Acquisition_Date_Glacier_Wide = integrate_point_balance(glacier,measurement_year,data,AAD,time_system,integration_method,integration_surface,0,0);
            
                [mass_maximum_site_adjustments,mass_maximum_date_numbers,mass_minimum_site_adjustments,mass_minimum_date_numbers] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);
                ba=insitu_data.ba+mass_minimum_site_adjustments';
                data=[];
                for site=1:length(balance_sites_for_year)%need to great final point balances table for integrate point balances function
                    data=[data; table(insitu_data.Year(strcmp(insitu_data.site_name,balance_sites_for_year(site))),insitu_data.site_name(strcmp(insitu_data.site_name,balance_sites_for_year(site))),insitu_data.elevation(strcmp(insitu_data.site_name,balance_sites_for_year(site))),...
                        nan,nan,ba(site),nan,{datestr(mass_minimum_date_numbers(strcmp(insitu_data.site_name,balance_sites_for_year(site))))},...
                        'VariableNames',{'Year' 'site_name' 'elevation' 'bw' 'bs' 'ba' 'mass_maximum_date' 'mass_minimum_date'})];
                end
                data=data(data.elevation>=min(AAD(1,2:end)-(AAD(1,3)-AAD(1,2)))&data.elevation<max(AAD(1,2:end)+(AAD(1,3)-AAD(1,2))),:);
                data=sortrows(data,'elevation');
                Mass_Minimum_Glacier_Wide = integrate_point_balance(glacier,measurement_year,data,AAD,time_system,integration_method,integration_surface,plot_integration,0);
            if datenum(geodetic_data.Date(measurement))>=datenum(['7/15/',num2str(measurement_year)]) && (geodetic_data.Date(measurement)~=datetime('02/11/2000'))%measurement dates to fall dates
                Ablation_adjustments_for_geodetics=[Ablation_adjustments_for_geodetics; table(geodetic_data.Date(measurement),{datestr(Mass_Minimum_Glacier_Wide.Ba_Date,'mm/dd/yyyy')},abs(Mass_Minimum_Glacier_Wide.Ba_Date-datenum(geodetic_data.Date(measurement))),Mass_Minimum_Glacier_Wide.Ba_mwe-Acquisition_Date_Glacier_Wide.Ba_mwe,'VariableNames',{'Aquistion_Date' 'Mass_Minimum_Date' 'Number_of_Days' 'Mass_Difference'})];             
            elseif geodetic_data.Date(measurement)==datetime('02/11/2000')%Aquisition date for SRTM
                Ablation_adjustments_for_geodetics=[Ablation_adjustments_for_geodetics; table(geodetic_data.Date(measurement),{datestr(Mass_Minimum_Glacier_Wide.Ba_Date,'mm/dd/yyyy')},0,0,'VariableNames',{'Aquistion_Date' 'Mass_Minimum_Date' 'Number_of_Days' 'Mass_Difference'})]; ; 
                geodetic_data.Date(measurement)=datestr(Mass_Minimum_Glacier_Wide.Ba_Date);     
            else %measurement dates to spring dates
                Ablation_adjustments_for_geodetics=[Ablation_adjustments_for_geodetics; table(geodetic_data.Date(measurement),{datestr(Mass_Minimum_Glacier_Wide.Ba_Date,'mm/dd/yyyy')},abs(Mass_Minimum_Glacier_Wide.Ba_Date-datenum(geodetic_data.Date(measurement))),Mass_Minimum_Glacier_Wide.Ba_mwe-Acquisition_Date_Glacier_Wide.Ba_mwe,'VariableNames',{'Aquistion_Date' 'Mass_Minimum_Date' 'Number_of_Days' 'Mass_Difference'})]; 
            end
            
                Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions.Ba_mwe(Glacier_Wide_solutions_table.Year==measurement_year)=Glacier_Wide_solutions_table.Ba_mwe(Glacier_Wide_solutions_table.Year==measurement_year)-Ablation_adjustments_for_geodetics.Mass_Difference(measurement);
                Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions.Ba_mwe(Glacier_Wide_solutions_table.Year==measurement_year+1)=Glacier_Wide_solutions_table.Ba_mwe(Glacier_Wide_solutions_table.Year==measurement_year+1)+Ablation_adjustments_for_geodetics.Mass_Difference(measurement);
    end
    
    [Calibrated_Glacier_Wide_solutions] = Geodetic_Calibration(glacier,geodetic_data,Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions,Geodetic_Calibration_index,break_point)
    Calibrated_Glacier_Wide_solutions.Ba_mwe(2:end)=Glacier_Wide_solutions_table.Ba_mwe-Calibrated_Glacier_Wide_solutions.Calibration(2:end);
    
    if Include_Previous_Glacierwide_Solutions==1 && exist(['data/',glacier,'/Input/Input_',glacier,'_Previous_Glacierwide_Solutions.csv']) 
        Previous_Glacier_wide_Solutions=readtable(['data/',glacier,'/Input/Input_',glacier,'_Previous_Glacierwide_Solutions.csv']);
            %create table to populate with previous solutions and append to
            %reanalyzed solutions
            s=cell(height(Glacier_Wide_solutions_table),1);
            s(1:height(Glacier_Wide_solutions_table),1)={'Reanalyzed'};
            Glacier_Wide_solutions_table.Source=s;
            table_to_append=table([],[],[],[],[],[],[],[],[],[],[],{},{},{},{},'VariableNames',Glacier_Wide_solutions_table.Properties.VariableNames);
            a=table_to_append;
            for year=Previous_Glacier_wide_Solutions.Year(1):Previous_Glacier_wide_Solutions.Year(end)
                blank_table=a;
                blank_table(1,1:4)=Previous_Glacier_wide_Solutions(Previous_Glacier_wide_Solutions.Year==year,1:4);
                blank_table(1,5)={nan};
                blank_table.Bw_Date='mm/dd/yyyy';
                blank_table.Ba_Date='mm/dd/yyyy';
                blank_table(1,8:11)={nan};
                blank_table.Source=Previous_Glacier_wide_Solutions.Source(Previous_Glacier_wide_Solutions.Year==year);
                blank_table.Lapse_Rate={lapse_rate};
                blank_table(1,12:14)={'nan'};
                table_to_append=[table_to_append;blank_table];
            end
            Previous_Glacier_wide_Solutions=table_to_append;
            Previous_Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions=Previous_Glacier_wide_Solutions;
        for measurement=1:height(geodetic_data)
            Previous_Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions.Ba_mwe(Previous_Glacier_wide_Solutions.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy')))=Previous_Glacier_wide_Solutions.Ba_mwe(Previous_Glacier_wide_Solutions.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy')))-Ablation_adjustments_for_geodetics.Mass_Difference(measurement);
            Previous_Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions.Ba_mwe(Previous_Glacier_wide_Solutions.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy'))+1)=Previous_Glacier_wide_Solutions.Ba_mwe(Previous_Glacier_wide_Solutions.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy'))+1)+Ablation_adjustments_for_geodetics.Mass_Difference(measurement);
        end
        
        [Calibrated_Previous_Glacier_Wide_solutions] = Geodetic_Calibration(glacier,geodetic_data,Previous_Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions,Geodetic_Calibration_index,break_point);
        Calibrated_Previous_Glacier_Wide_solutions.Ba_mwe(2:end)=Previous_Glacier_wide_Solutions.Ba_mwe-Calibrated_Previous_Glacier_Wide_solutions.Calibration(2:end);
        writetable(Calibrated_Previous_Glacier_Wide_solutions(2:end,:),['data/',glacier,'/Output/Output_',glacier,'_Previous_Glacier_Wide_solutions_calibrated.csv']);
    data=[];
    for year=2:height(Calibrated_Glacier_Wide_solutions)
        if ~isempty(Calibrated_Previous_Glacier_Wide_solutions.Ba_mwe(Calibrated_Previous_Glacier_Wide_solutions.Year==Calibrated_Glacier_Wide_solutions.Year(year)))
            data(year-1,1)=Calibrated_Glacier_Wide_solutions.Ba_mwe(year);
            data(year-1,2)=Calibrated_Previous_Glacier_Wide_solutions.Ba_mwe(Calibrated_Previous_Glacier_Wide_solutions.Year==Calibrated_Glacier_Wide_solutions.Year(year));
        end
    end
    [r,p] = corr(data(:,2),data(:,1),'type','Pearson');
    MAE_from_Previous_solutions=nanmean(abs(data(:,1)-data(:,2)));
    disp(['***** Correlation between Previous and Reanalyzed Glacier-wide Solutions: r=',num2str(round(r,2)),'; p-value=',num2str(round(p,2)),'; MAE=',num2str(round(MAE_from_Previous_solutions,2))])
    figure();hold on
    errorbarxy(data(:,2),data(:,1),0.38*ones(length(data),1),0.38*ones(length(data),1),{'ko',[1 .5 0],[1 .5 0]});hold on
    plot((-3:3),(-3:3),'LineWidth',2,'color','k')
    
    text(0,-1.75,['r= ',num2str(round(r,2))],'fontweight','bold');
    text(0,-2.1,['MAE= ',num2str(round(MAE_from_Previous_solutions,2)),' m w.e.'],'fontweight','bold');
%     text(0,-2.25,['std= ',num2str(JIRP_Reanalysis_Ba_model_std),' m w.e.'],'fontweight','bold');
    xlim([-3 3])
    ylim([-3 3])
    xlabel('Conventional (m w.e.)','fontweight','bold')
    ylabel('JIRP (m w.e.)','fontweight','bold')
%     title(glaciers(glacier_index),'fontweight','bold')
    set(gca,'fontname','arial ','fontsize',12,'TickLength',[0.025 0.025],'linewidth',2)
    axis square
    box on
    set(gcf, 'PaperPositionMode', 'auto');
    print -depsc2 gates_epoch2.eps
    
    Calibrated_Previous_Glacier_Wide_solutions.bw_gradient_mwe_km=num2cell(Calibrated_Previous_Glacier_Wide_solutions.bw_gradient_mwe_km)
    Calibrated_Previous_Glacier_Wide_solutions.bs_gradient_mwe_km=num2cell(Calibrated_Previous_Glacier_Wide_solutions.bs_gradient_mwe_km);
    Calibrated_Previous_Glacier_Wide_solutions.ba_gradient_mwe_km=num2cell(Calibrated_Previous_Glacier_Wide_solutions.ba_gradient_mwe_km);
    Merged_timeseries=[Calibrated_Previous_Glacier_Wide_solutions(Calibrated_Previous_Glacier_Wide_solutions.Year<Calibrated_Glacier_Wide_solutions.Year(2),:);Calibrated_Glacier_Wide_solutions(2:end,:)]
    Merged_timeseries=Merged_timeseries(2:end,:);
    Merged_Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions=Merged_timeseries(:,1:end-1);
    for measurement=1:height(geodetic_data)
        Merged_Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions.Ba_mwe(Merged_timeseries.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy')))=Merged_timeseries.Ba_mwe(Merged_timeseries.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy')))-Ablation_adjustments_for_geodetics.Mass_Difference(measurement);
        Merged_Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions.Ba_mwe(Merged_timeseries.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy'))+1)=Merged_timeseries.Ba_mwe(Merged_timeseries.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy'))+1)+Ablation_adjustments_for_geodetics.Mass_Difference(measurement);
    end
    
    [Calibrated_Merged_Glacier_Wide_solutions] = Geodetic_Calibration(glacier,geodetic_data,Merged_Glacier_wides_solutions_adjusted_to_Geodetic_aquisitions,Geodetic_Calibration_index,break_point);
    for measurement=1:height(geodetic_data)
        Calibrated_Merged_Glacier_Wide_solutions.Ba_mwe(Calibrated_Merged_Glacier_Wide_solutions.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy')))=Calibrated_Merged_Glacier_Wide_solutions.Ba_mwe(Calibrated_Merged_Glacier_Wide_solutions.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy')))+Ablation_adjustments_for_geodetics.Mass_Difference(measurement);
        Calibrated_Merged_Glacier_Wide_solutions.Ba_mwe(Calibrated_Merged_Glacier_Wide_solutions.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy'))+1)=Calibrated_Merged_Glacier_Wide_solutions.Ba_mwe(Calibrated_Merged_Glacier_Wide_solutions.Year==str2num(datestr(Ablation_adjustments_for_geodetics.Mass_Minimum_Date(measurement),'yyyy'))+1)-Ablation_adjustments_for_geodetics.Mass_Difference(measurement);
    end
    
    for year=1:height(Previous_Glacier_wide_Solutions)
       Calibrated_Merged_Glacier_Wide_solutions.Calibration(Calibrated_Merged_Glacier_Wide_solutions.Year==Previous_Glacier_wide_Solutions.Year(year))=Previous_Glacier_wide_Solutions.Ba_mwe(year)-Calibrated_Merged_Glacier_Wide_solutions.Ba_mwe(Calibrated_Merged_Glacier_Wide_solutions.Year==Previous_Glacier_wide_Solutions.Year(year));
    end
    for year=1:height(Glacier_Wide_solutions_table)
       Calibrated_Merged_Glacier_Wide_solutions.Calibration(Calibrated_Merged_Glacier_Wide_solutions.Year==Glacier_Wide_solutions_table.Year(year))=Glacier_Wide_solutions_table.Ba_mwe(year)-Calibrated_Merged_Glacier_Wide_solutions.Ba_mwe(Calibrated_Merged_Glacier_Wide_solutions.Year==Glacier_Wide_solutions_table.Year(year));
    end
    Calibrated_Glacier_Wide_solutions=Calibrated_Merged_Glacier_Wide_solutions
    end
    final_geodetic_data=[geodetic_data Ablation_adjustments_for_geodetics(:,2:4)];
    final_geodetic_data.Mass_Change_mwe=round(final_geodetic_data.Mass_Change_mwe,2);
    final_geodetic_data.Mass_Difference=round(final_geodetic_data.Mass_Difference,2);
    final_geodetic_data.Mass_Change_mwe=(final_geodetic_data.Mass_Change_mwe-final_geodetic_data.Mass_Change_mwe(end)).*-1;
else
    
    [Calibrated_Glacier_Wide_solutions] = Geodetic_Calibration(glacier,geodetic_data,Glacier_Wide_solutions_table,Geodetic_Calibration_index,break_point)
    Calibrated_Glacier_Wide_solutions.Ba_mwe(2:end)=Glacier_Wide_solutions_table.Ba_mwe-Calibrated_Glacier_Wide_solutions.Calibration(2:end);
    
    final_geodetic_data=geodetic_data;
    
end
Calibrated_Glacier_Wide_solutions.Ba_mwe(1)=0;
Calibrated_Glacier_Wide_solutions.Calibration=round(Calibrated_Glacier_Wide_solutions.Calibration,2).*-1; 
Calibrated_Glacier_Wide_solutions.Ba_mwe=round(Calibrated_Glacier_Wide_solutions .Ba_mwe,2);
Final_Glacier_Wide_solutions=Calibrated_Glacier_Wide_solutions 
Final_Glacier_Wide_solutions.Bs_mwe=Final_Glacier_Wide_solutions.Ba_mwe-Final_Glacier_Wide_solutions.Bw_mwe;

end

