function[final_point_balances,Glacier_Wide_solutions_table]=Calculate_GlacierWide_Balance(glacier,balance_sites,years,Glaciological_data,Weather_data,AAD,ks,ki,precipitation_ratios_table,lapse_rate,time_system,integration_surface,integration_method,plot_integration,plot_ablation_model,nan_incomplete_glaciological_data)
%% Main objective: Calculate glacier-wide Balance for each year selected
    %Glacier-wide balance calculated from glaciological data, using a mass
    %balance model based on temperature and precipitation to adjust
    %observations to the time-system selected by user. Finally, adjusted
    %glaciological data is integrated over the glacier hypsometry using
    %both the integration method (index or gradient) and surface
    %(Reference-Surface or Conventional) selected by user.
    %A) for each year run Find_Mass_Maximum_and_Minimum_Adjustments.m to
    %   get mass maximum and minimum point balance and dates.
    
    %   1) If user defines a stratigraphic balance,
    %       Find_Mass_Maximum_and_Minimum_Adjustments.m can be run to find mass
    %       extreme dates and winter ablation and summer accumulation if not
    %       measured. Otherwise, we use the observations instead ofthe model
    
    %   2) If user defines either fixed-date(hydrologic-year) or
    %       floating-date/stratigraphic Find_Mass_Maximum_and_Minimum_Adjustments.m is used to
    %       adjust every site to the mass extrema in the selected
    %       time-system
    
    %       a) Additionally, Find_Mass_Maximum_and_Minimum_Adjustments.m
    %       must be used to determine the difference from the mass minimum
    %       date in the selected time-system to the stratigraphic
    %       time-system, since this is the time-system of the input data
    
    %B) Integrate point balance using method and surface selected by user
dbstop if error
Glacier_Wide_solutions_table=table([],[],[],[],[],[],[],[],[],[],'VariableNames',{'Year' 'Bw_mwe','Bs_mwe','Ba_mwe','ELA_m','Bw_Date','Ba_Date','bw_gradient_mwe_km','bs_gradient_mwe_km','ba_gradient_mwe_km'});
final_point_balances=table([],[],[],[],[],[],[],[],'VariableNames',{'Year' 'site_name' 'elevation' 'bw' 'bs' 'ba' 'mass_maximum_date' 'mass_minimum_date'});    
for year=1:length(years)
        insitu_data=[];
        next_years_insitu_data=[];
        start_of_balance_year_adjustments=[];
        point_balances=[];
        ba=[];
        bw=[];
        bs=[];
        for site=1:length(balance_sites)
            %since additional ablation can occur after the fall
            %measurements, if it was measured the next spring, then we use
            %that to get out full on stratigraphic ba, otherwise we need to
            %use 'Find_Mass_Maximum_and_Minimum_Adjustments.m' to model any
            %additional ablation
            current_years_data=[];
             if ~strcmp(balance_sites(site),'TSL')
                current_years_data=Glaciological_data(Glaciological_data.Year==years(year) & strcmp(Glaciological_data.site_name,balance_sites(site)),:);
                if ~isempty(current_years_data)
                        winter_ablation=Glaciological_data.winter_ablation(Glaciological_data.Year==years(year)+1 & strcmp(Glaciological_data.site_name,balance_sites(site)));
                    if isempty(winter_ablation)
                        current_years_data.winter_ablation=NaN;
                    else
                        current_years_data.winter_ablation=winter_ablation;
                    end
                    insitu_data=[insitu_data;current_years_data];
                end
            elseif strcmp(balance_sites(site),'TSL')
                current_years_data=Glaciological_data(Glaciological_data.Year==years(year) & contains(Glaciological_data.site_name,'TSL'),:);
                insitu_data=[insitu_data;current_years_data];
            end
        end
        if ~isempty(insitu_data)
            insitu_data=sortrows(insitu_data,'elevation');
        end

        if time_system==1 %Stratigraphic time-system selected
            [mass_maximum_site_adjustments,mass_maximum_date_numbers,mass_minimum_site_adjustments,mass_minimum_date_numbers,Modeled_Mass_Balance] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);
            for site=1:height(insitu_data)
                if isnan(insitu_data.winter_ablation(site)) 
                     ba(site,1)=insitu_data.ba(site)+mass_minimum_site_adjustments(site);
                else
                    ba(site,1)=insitu_data.ba(site)+insitu_data.winter_ablation(site);
                end
            end
            bw=insitu_data.bw+mass_maximum_site_adjustments'; 
            bs=ba-bw;
            point_balances=table(insitu_data.Year,insitu_data.site_name,insitu_data.elevation,bw,bs,ba,mass_maximum_date_numbers,mass_minimum_date_numbers,'VariableNames',{'Year' 'site_name' 'elevation' 'bw' 'bs' 'ba' 'mass_maximum_date' 'mass_minimum_date'});   
        elseif time_system==2%stratigraphic measurements, no correction to mass minimum
            [~,mass_maximum_date_numbers,mass_minimum_site_adjustments,mass_minimum_date_numbers,Modeled_Mass_Balance] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,1,integration_method,integration_surface,plot_ablation_model);
            for site=1:height(insitu_data)
                if isnan(insitu_data.winter_ablation(site)) 
                     ba(site,1)=insitu_data.ba(site)+mass_minimum_site_adjustments(site);
                else
                    ba(site,1)=insitu_data.ba(site)+insitu_data.winter_ablation(site);
                end
            end
            bw=insitu_data.bw; 
            bs=ba-bw;
            point_balances=table(insitu_data.Year,insitu_data.site_name,insitu_data.elevation,bw,bs,ba,insitu_data.spring_date,mass_minimum_date_numbers,'VariableNames',{'Year' 'site_name' 'elevation' 'bw' 'bs' 'ba' 'mass_maximum_date' 'mass_minimum_date'});   
        
        elseif time_system==3%hydrologic year 
            previous_years_insitu_data=insitu_data;
            previous_years_insitu_data.Year=previous_years_insitu_data.Year-1;
            previous_years_insitu_data.spring_date=previous_years_insitu_data.spring_date-365;
            previous_years_insitu_data.fall_date=previous_years_insitu_data.fall_date-365;
            [~,~,~,previous_years_selected_time_system_mass_minimum_date_numbers,Modeled_Mass_Balance] = Find_Mass_Maximum_and_Minimum_Adjustments(previous_years_insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,0);
            [~,~,~,previous_years_stratigraphic_time_system_mass_minimum_date_numbers] = Find_Mass_Maximum_and_Minimum_Adjustments(previous_years_insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,1,integration_method,integration_surface,0);
             start_of_balance_year_adjustments=zeros(length(previous_years_selected_time_system_mass_minimum_date_numbers),1);
            for site=1:height(insitu_data)
                if previous_years_selected_time_system_mass_minimum_date_numbers(site) > previous_years_stratigraphic_time_system_mass_minimum_date_numbers(site)
                    start_of_balance_year_adjustments(site,1)=table2array(Modeled_Mass_Balance(datenum(Modeled_Mass_Balance.Date)==previous_years_selected_time_system_mass_minimum_date_numbers(site),strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))))-table2array(Modeled_Mass_Balance(datenum(Modeled_Mass_Balance.Date)==previous_years_stratigraphic_time_system_mass_minimum_date_numbers(site),strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))));
                else
                    start_of_balance_year_adjustments(site,1)=table2array(Modeled_Mass_Balance(datenum(Modeled_Mass_Balance.Date)==previous_years_stratigraphic_time_system_mass_minimum_date_numbers(site),strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))))-table2array(Modeled_Mass_Balance(datenum(Modeled_Mass_Balance.Date)==previous_years_selected_time_system_mass_minimum_date_numbers(site),strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))));
                end
            end
            [mass_maximum_site_adjustments,mass_maximum_date_numbers,mass_minimum_site_adjustments,mass_minimum_date_numbers,Modeled_Mass_Balance] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);
            for site=1:height(insitu_data)
                if isnan(insitu_data.winter_ablation(site)) 
                     ba(site,1)=insitu_data.ba(site)-abs(start_of_balance_year_adjustments(site,1))+mass_minimum_site_adjustments(site)';
                else
                    ba(site,1)=insitu_data.ba(site)+insitu_data.winter_ablation(site)-abs(start_of_balance_year_adjustments(site,1))+mass_minimum_site_adjustments(site)';
                end
            end
            bw=insitu_data.bw-abs(start_of_balance_year_adjustments)+mass_maximum_site_adjustments';
            bs=ba-bw;
            point_balances=table(insitu_data.Year,insitu_data.site_name,insitu_data.elevation,round(bw,2),round(bs,2),round(ba,2),mass_maximum_date_numbers,mass_minimum_date_numbers,'VariableNames',{'Year' 'site_name' 'elevation' 'bw' 'bs' 'ba' 'mass_maximum_date' 'mass_minimum_date'});
        elseif time_system==4%combined floating date/stratigraphic system
            %have to generate mock data for desired sites for the previous
            %year to find the mass minimum correction. We don't always 

            previous_years_insitu_data=insitu_data;
            previous_years_insitu_data.Year=previous_years_insitu_data.Year-1;
            previous_years_insitu_data.spring_date=previous_years_insitu_data.spring_date-365;
            previous_years_insitu_data.fall_date=previous_years_insitu_data.fall_date-365;
            if height(insitu_data)<=2
                %we can't attain glacier-wide min date because not enough
                %data exists. So we have to just calculate stratigraphic
                %balance adjustment
                
                [~,~,~,previous_years_selected_time_system_mass_minimum_date_numbers,Modeled_Mass_Balance] = Find_Mass_Maximum_and_Minimum_Adjustments(previous_years_insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,1,integration_method,integration_surface,0);
            else
                [~,~,~,previous_years_selected_time_system_mass_minimum_date_numbers,Modeled_Mass_Balance] = Find_Mass_Maximum_and_Minimum_Adjustments(previous_years_insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,0);
            end
            [~,~,~,previous_years_stratigraphic_time_system_mass_minimum_date_numbers] = Find_Mass_Maximum_and_Minimum_Adjustments(previous_years_insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,1,integration_method,integration_surface,0);
             start_of_balance_year_adjustments=zeros(length(previous_years_selected_time_system_mass_minimum_date_numbers),1);
            for site=1:height(insitu_data)
                %since all input point balances are in a stratigraphic
                %time-system, we need to make adjustment for mass minimum
                %dates from previous year
                if previous_years_selected_time_system_mass_minimum_date_numbers(site) > previous_years_stratigraphic_time_system_mass_minimum_date_numbers(site)
                    %desired time-systems mass minimum came after the
                    %stratigraphic mass minimum of the previous year
                    start_of_balance_year_adjustments(site,1)=table2array(Modeled_Mass_Balance(datenum(Modeled_Mass_Balance.Date)==previous_years_selected_time_system_mass_minimum_date_numbers(site),strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))))-table2array(Modeled_Mass_Balance(datenum(Modeled_Mass_Balance.Date)==previous_years_stratigraphic_time_system_mass_minimum_date_numbers(site),strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))));
                else
                    %desired time-systems mass minimum came before the
                    %stratigraphic mass minimum of the previous year
                    start_of_balance_year_adjustments(site,1)=table2array(Modeled_Mass_Balance(datenum(Modeled_Mass_Balance.Date)==previous_years_stratigraphic_time_system_mass_minimum_date_numbers(site),strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))))-table2array(Modeled_Mass_Balance(datenum(Modeled_Mass_Balance.Date)==previous_years_selected_time_system_mass_minimum_date_numbers(site),strcmp(Modeled_Mass_Balance.Properties.VariableNames,insitu_data.site_name(site))));
                end
            end
            if height(insitu_data)<=2
                %we can't attain glacier-wide min date because not enough
                %data exists.. So we have to just calculate stratigraphic
                %balance adjustment
                [mass_maximum_site_adjustments,mass_maximum_date_numbers,mass_minimum_site_adjustments,mass_minimum_date_numbers,Modeled_Mass_Balance] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,1,integration_method,integration_surface,plot_ablation_model);
            else
                %otherwise we are fine
                [mass_maximum_site_adjustments,mass_maximum_date_numbers,mass_minimum_site_adjustments,mass_minimum_date_numbers,Modeled_Mass_Balance] = Find_Mass_Maximum_and_Minimum_Adjustments(insitu_data,Weather_data,lapse_rate,ks,ki,precipitation_ratios_table,glacier,AAD,time_system,integration_method,integration_surface,plot_ablation_model);
            end
            
            for site=1:height(insitu_data)
                if isnan(insitu_data.winter_ablation(site)) 
                     ba(site,1)=insitu_data.ba(site)-abs(start_of_balance_year_adjustments(site,1))+mass_minimum_site_adjustments(site)';
                else
                    ba(site,1)=insitu_data.ba(site)+insitu_data.winter_ablation(site)-abs(start_of_balance_year_adjustments(site,1))+mass_minimum_site_adjustments(site)';
                end
            end
            bw=insitu_data.bw-abs(start_of_balance_year_adjustments)+mass_maximum_site_adjustments';
            bs=ba-bw;
            point_balances=table(insitu_data.Year,insitu_data.site_name,insitu_data.elevation,round(bw,2),round(bs,2),round(ba,2),mass_maximum_date_numbers,mass_minimum_date_numbers,'VariableNames',{'Year' 'site_name' 'elevation' 'bw' 'bs' 'ba' 'mass_maximum_date' 'mass_minimum_date'});
            %insert mass max and min for loop here
            
        end

        
        
        
    point_balances.mass_maximum_date=datestr(point_balances.mass_maximum_date);
    point_balances.mass_minimum_date=datestr(point_balances.mass_minimum_date);
    
    final_point_balances=[final_point_balances;point_balances];
    [Glacier_Wide_solutions] = integrate_point_balance(glacier,years(year),point_balances,AAD,time_system,integration_method,integration_surface,plot_integration,nan_incomplete_glaciological_data); 
    if isnan(Glacier_Wide_solutions.Bw_Date)
        Glacier_Wide_solutions.Bw_Date='mm/dd/yyyy';
    else
        Glacier_Wide_solutions.Bw_Date=datestr(Glacier_Wide_solutions.Bw_Date,'mm/dd/yyyy');
    end
    if isnan(Glacier_Wide_solutions.Ba_Date)
        Glacier_Wide_solutions.Ba_Date='mm/dd/yyyy';
    else
        Glacier_Wide_solutions.Ba_Date=datestr(Glacier_Wide_solutions.Ba_Date,'mm/dd/yyyy');
    end
    Glacier_Wide_solutions_table=[Glacier_Wide_solutions_table;Glacier_Wide_solutions];
    end


end

