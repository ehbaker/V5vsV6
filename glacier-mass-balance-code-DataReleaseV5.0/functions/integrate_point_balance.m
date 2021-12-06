function [Glacier_Wide_solutions,bals] = integrate_point_balance(glacier,year,point_balances,AAD,time_system,integration_method,integration_surface,plot_integration,nan_incomplete_glaciological_data)
%This function integrates point balances over the glacier hypsometry
%Integration can be done two ways
%The index method or the gradient method
dbstop if error

bin_centers=AAD(1,2:end);
if integration_surface(1,1)==1
    if year<AAD(2,1)
        disp('WARNING: The AAD provided does not cover the full range of glaciological data. First available hypsometry used. Updating AAD file recomended')
        glacier_hypsometry=AAD(2,2:end);
    elseif year>AAD(end,1)
        disp('WARNING: The AAD provided does not cover the full range of glaciological data. The most recent hypsometry is used instead. Updating AAD file recomended')
        glacier_hypsometry=AAD(end,2:end);
    else
        glacier_hypsometry=AAD(AAD(:,1)==year,2:end);
    end
elseif integration_surface(1,1)==2
    if year<AAD(2,1)
        disp('WARNING: The AAD provided does not cover the full range of glaciological data. Last available hypsometry used. Updating AAD file recomended')
        glacier_hypsometry=AAD(2,2:end);
    else
        glacier_hypsometry=AAD(2,2:end);
    end
elseif integration_surface(1,1)==3
    if integration_surface(1,2)<AAD(2,1)
        disp('WARNING: The AAD provided does not cover the full range of glaciological data. Last available hypsometry used. Updating AAD file recomended')
        glacier_hypsometry=AAD(2,2:end);
    else
        glacier_hypsometry=AAD(AAD(:,1)==integration_surface(1,2),2:end);
    end
end
if nan_incomplete_glaciological_data==1&&(height(point_balances)<=2||max(point_balances.elevation)<(sum(bin_centers.*(glacier_hypsometry./sum(glacier_hypsometry))))||min(point_balances.elevation)>(sum(bin_centers.*(glacier_hypsometry./sum(glacier_hypsometry)))))

        Glacier_Wide_solutions=table(year,nan,nan,nan,nan,NaN,NaN,nan,nan,nan,'VariableNames',{'Year' 'Bw_mwe','Bs_mwe','Ba_mwe','ELA_m','Bw_Date','Ba_Date','bw_gradient_mwe_km','bs_gradient_mwe_km','ba_gradient_mwe_km'}); 
%         disp(['Not enough data to calculate glacier wide balances in ',num2str(year)])
        if max(point_balances.elevation)<(sum(bin_centers.*(glacier_hypsometry./sum(glacier_hypsometry))))
            disp(['No data available for the upper portion of glacier in ',num2str(year),'. Unabale to calculate glacier-wide balances'])
        elseif min(point_balances.elevation)>(sum(bin_centers.*(glacier_hypsometry./sum(glacier_hypsometry))))
            disp(['No data available for the lower portion of glacier in ',num2str(year),'. Unabale to calculate glacier-wide balances'])
        else
            disp(['Only ',num2str(height(point_balances)),' are available for glacier in ',num2str(year),'. Unabale to calculate glacier-wide balances'])
        end
else
    if height(point_balances)<=2
        disp('Brace for impact!!')
        disp(['No data available for the lower portion of glacier in ',num2str(year),'. Unabale to calculate glacier-wide balances'])
        disp(['Suggest toggling on "NaN incomplete glaciological data" to model glacier-wide using stake to glacier-wide regression'])
    end
        
    if integration_method==1%Index Method
        %get area weights for stakes
        [zone_weights,split_elevation] = get_weights(bin_centers,glacier_hypsometry,point_balances.elevation);
        zone_weights(zone_weights==inf|zone_weights==-inf|isnan(zone_weights))=0;
        if plot_integration==1
            bin_area=glacier_hypsometry(glacier_hypsometry~=0);
            this_yr_bin_centers=bin_centers(glacier_hypsometry~=0);
            mean_zone_elevation=[];
            zone_width=[];
            for k=1:length(split_elevation)
                if k==1
                    mean_zone_elevation(k,1)=split_elevation(k,1)-(split_elevation(k,1)-(this_yr_bin_centers(1)-((this_yr_bin_centers(2)-this_yr_bin_centers(1))/2)))/2;
                    zone_width(k,1)=split_elevation(k,1)-(this_yr_bin_centers(1)-(this_yr_bin_centers(2)-this_yr_bin_centers(1))/2);
                elseif k==length(split_elevation)
                     zone_width(k,1)=(split_elevation(k,1)-split_elevation(k-1,1));
                     zone_width(k+1,1)=(((this_yr_bin_centers(end)+(this_yr_bin_centers(end)-this_yr_bin_centers(end-1))/2)-split_elevation(k,1)));
                     mean_zone_elevation(k,1)=split_elevation(k,1)-((split_elevation(k,1)-split_elevation(k-1,1))/2);
                     mean_zone_elevation(k+1,1)=split_elevation(k,1)+(((this_yr_bin_centers(end)+(this_yr_bin_centers(end)-this_yr_bin_centers(end-1))/2)-split_elevation(k,1))/2);
                else
                    zone_width(k,1)=(split_elevation(k,1)-split_elevation(k-1,1));
                    mean_zone_elevation(k,1)=mean([split_elevation(k,1) split_elevation(k-1,1)]);
                end
            end

                zone_area=sum(glacier_hypsometry).*zone_weights;
                
                figure(year);hold on 
                subplot(1,3,1)
                title('Winter Balance Integration')
                hold on
                yyaxis left
                areas=zone_area(~isnan(zone_area) & zone_area~=0,1);
                for j=1:length(mean_zone_elevation)
                    bar(mean_zone_elevation(j,1),areas(j,1),zone_width(j,1),'FaceColor',[.5 .5 .5])
                end
                ylabel('Area (km^2)')
                yyaxis right
                scatter(point_balances.elevation,point_balances.bw,200,'s','filled')
                ylabel('balance (m w.e.)')
                ylim([round(min(min(point_balances.bw)))-1 round(max(max(point_balances.bw)))+1])
                xlim([bin_centers(1)-50 bin_centers(end)+50])
                axis square
                box on

                subplot(1,3,2)
                title('Summer Balance Integration')
                hold on
                yyaxis left
                for j=1:length(mean_zone_elevation)
                    bar(mean_zone_elevation(j,1),areas(j,1),zone_width(j,1),'FaceColor',[.5 .5 .5])
                end
                ylabel('Area (km^2)')
                yyaxis right
                scatter(point_balances.elevation,point_balances.bs,200,'s','filled')
                ylabel('balance (m w.e.)')
                ylim([round(min(min(point_balances.bs)))-1 round(max(max(point_balances.bs)))+1])
                xlim([bin_centers(1)-50 bin_centers(end)+50])
                axis square
                box on

                subplot(1,3,3)
                title('Annual Balance Integration')
                hold on
                yyaxis left
                areas=zone_area(~isnan(zone_area) & zone_area>=0,1);
                for j=1:length(mean_zone_elevation)
                    bar(mean_zone_elevation(j,1),areas(j,1),zone_width(j,1),'FaceColor',[.5 .5 .5])
                end
                ylabel('Area (km^2)')
                yyaxis right
                scatter(point_balances.elevation,point_balances.ba,200,'s','filled')
                ylabel('balance (m w.e.)')
                ylim([round(min(min(point_balances.ba)))-1 round(max(max(point_balances.ba)))+1])
                xlim([bin_centers(1)-50 bin_centers(end)+50])
                axis square
                box on
                set(gcf, 'PaperPosition', [0 0 3 6]);
                print -depsc2 gates_epoch2.eps

        end
        if sum(isnan(point_balances.bw))>=sum(~isnan(point_balances.bw))
            Bw=NaN;
        else
            Bw=round(nansum(point_balances.bw.*zone_weights),2);
        end
        if sum(isnan(point_balances.ba))>=sum(~isnan(point_balances.bw))
            Bs=NaN;
        else
            Bs=round(nansum(point_balances.bs.*zone_weights),2);
        end

            Ba=nansum(point_balances.ba.*zone_weights);
        if time_system==1 %stratigraphic
            Bw_Date=round(nansum(datenum(point_balances.mass_maximum_date).*zone_weights));
            Ba_Date=round(nansum(datenum(point_balances.mass_minimum_date).*zone_weights));
        else
            Bw_Date=datenum(point_balances.mass_maximum_date(1,:));
            Ba_Date=datenum(point_balances.mass_minimum_date(1,:));
        end
        bw_gradient={NaN};
        bs_gradient={NaN};
        ba_gradient={NaN};
        if sum(~isnan(point_balances.ba))==height(point_balances)&&sum(point_balances.ba==0)~=height(point_balances)
            if length(unique(point_balances.ba))~=length(point_balances.ba)%Interp1 fails in the case there are identical point balances
                point_balances_edited=point_balances;%so we need to make identical values different, but not edit final values used to calculated glacier-wide
                n=0;
                for site=1:height(point_balances)
                    n=n+0.01;
                    index=find(point_balances_edited.ba==point_balances.ba(site));%find identical sites
                    if length(index)>1
                        point_balances_edited.ba(index(1))=point_balances.ba(index(1))+n;%add just enough to make point balances different without changing the final point balances used
                    end
                end
%                 if (height(point_balances)~=length(unique(point_balances.elevation)))||(height(point_balances)~=length(unique(point_balances.ba)))
%                     
%                 end
                ELA=round(interp1(point_balances_edited.ba,point_balances_edited.elevation,0,'linear','extrap'));
            elseif length(point_balances.ba)==1||sum(point_balances.ba==0)==height(point_balances)
                ELA=NaN;%there is not enough data to calculate ELA
            else %otherwise there is not issue
                ELA=round(interp1(point_balances.ba,point_balances.elevation,0,'linear','extrap'));
            end
        else
            ELA=nan;
        end
    Glacier_Wide_solutions=table(year,Bw,Bs,Ba,ELA,Bw_Date,Ba_Date,bw_gradient,bs_gradient,ba_gradient,'VariableNames',{'Year' 'Bw_mwe','Bs_mwe','Ba_mwe','ELA_m','Bw_Date','Ba_Date','bw_gradient_mwe_km','bs_gradient_mwe_km','ba_gradient_mwe_km'});       

    elseif integration_method==2 %Linear balance profile

        bin_weights(:,1)=(glacier_hypsometry(1:end)./sum(glacier_hypsometry(1:end)));
        if ~isnan(point_balances.bw)
            % get winter balance gradient and balance
            winter_lm=fitlm(point_balances.elevation,point_balances.bw);
            bw_gradient={round((double(cell2mat(table2cell(winter_lm.Coefficients(2,1))))*1000),1)};
            Bw=round(nansum(bin_weights.*predict(winter_lm,bin_centers')),2);
            Bw_Date=round(nansum(bin_weights.*predict(fitlm(point_balances.elevation,datenum(point_balances.mass_maximum_date)),bin_centers')));
        else
            bw_gradient={nan};
            Bw=nan;
            Bw_Date=nan;
        end
        if ~isnan(point_balances.bs)
            % get summer balance gradient and balance
            summer_lm=fitlm(point_balances.elevation,point_balances.bs);
            bs_gradient={round((double(cell2mat(table2cell(summer_lm.Coefficients(2,1))))*1000),1)};
            Bs=round(nansum(bin_weights.*predict(summer_lm,bin_centers')),2);    
        else
            bs_gradient={nan};
            Bs=nan;
        end
        if ~isnan(point_balances.ba)
            % get annual gradient and balance
            annual_lm=fitlm(point_balances.elevation,point_balances.ba);
            ba_gradient={round((table2array(annual_lm.Coefficients(2,1))*1000),1)};
            Ba=nansum(bin_weights.*predict(annual_lm,bin_centers'));
            Ba_Date=round(nansum(bin_weights.*predict(fitlm(point_balances.elevation,datenum(point_balances.mass_minimum_date)),bin_centers')));
            ELA=round(predict(fitlm(point_balances.ba,point_balances.elevation),0));
        else
            ba_gradient={nan};
            Ba=nan;
            Ba_Date=nan;
            ELA=nan;
        end

        Glacier_Wide_solutions=table(year,Bw,Bs,Ba,ELA,Bw_Date,Ba_Date,bw_gradient,bs_gradient,ba_gradient,'VariableNames',{'Year','Bw_mwe','Bs_mwe','Ba_mwe','ELA_m','Bw_Date','Ba_Date','bw_gradient_mwe_km','bs_gradient_mwe_km','ba_gradient_mwe_km'});       

        if plot_integration==1
            figure(year);hold on 
            subplot(1,3,1)
            title('Winter Balance Integration')
            hold on
            yyaxis left
            bar(bin_centers,glacier_hypsometry,1,'FaceColor',[.5 .5 .5])
            ylabel('Area (km^2)')
            yyaxis right
            scatter(point_balances.elevation,point_balances.bw,'b','fill');
            plot(bin_centers,predict(winter_lm,bin_centers'),'LineWidth',1,'color','b')
            ylabel('balance (m w.e.)')
            ylim([round(min(min(point_balances.bw)))-1 round(max(max(point_balances.bw)))+1])
            xlim([bin_centers(1)-50 bin_centers(end)+50])
            axis square
            box on

            subplot(1,3,2)
            title('Summer Balance Integration')
            hold on
            yyaxis left
            bar(bin_centers,glacier_hypsometry,1,'FaceColor',[.5 .5 .5])
            ylabel('Area (km^2)')
            yyaxis right
            scatter(point_balances.elevation,point_balances.bs,'r','fill');
            plot(bin_centers,predict(summer_lm,bin_centers'),'LineWidth',1,'color','b')
            ylabel('balance (m w.e.)')
            ylim([round(min(min(point_balances.bs)))-1 round(max(max(point_balances.bs)))+1])
            xlim([bin_centers(1)-50 bin_centers(end)+50])
            axis square
            box on

            subplot(1,3,3)
            title('Annual Balance Integration')
            hold on
            yyaxis left
            bar(bin_centers(1:end),glacier_hypsometry,1,'FaceColor',[.5 .5 .5])
            ylabel('Area (km^2)')
            yyaxis right
            scatter(point_balances.elevation,point_balances.ba,'o','fill');
            plot(bin_centers,predict(annual_lm,bin_centers'),'LineWidth',1,'color','b')
            ylabel('balance (m w.e.)')
            ylim([round(min(min(point_balances.ba)))-1 round(max(max(point_balances.ba)))+1])
            xlim([bin_centers(1)-50 bin_centers(end)+50])
            axis square
            box on
            set(gcf, 'PaperPosition', [0 0 3 6]);
            print -depsc2 gates_epoch2.eps
        end
    elseif integration_method==3% Piecewise Linear (free-knot 2 part linear spline)
    %     if exists
        bin_weights(:,1)=(glacier_hypsometry(1:end)./sum(glacier_hypsometry(1:end)));
        if ~isnan(point_balances.ba)
            % get annual gradient and balance

            pp=splinefit(point_balances.elevation,point_balances.ba,2,2);
            while length(pp.coefs(:,1))>1 && pp.coefs(2,1)<0
                breaks=[min(point_balances.elevation) point_balances.elevation(find(point_balances.elevation<pp.breaks(2),1,'last')) max(point_balances.elevation)];
                pp=splinefit(point_balances.elevation,point_balances.ba,breaks,2);
            end

            break_point_elevation=round(pp.breaks(1,2));
            break_point_balance=ppval(pp,break_point_elevation);
            annual_balance_profile_values=ppval(pp,bin_centers);
            integrated_ba=annual_balance_profile_values'.*bin_weights;
            glacier_wide=sum(integrated_ba);
            annual_lm=fitlm(point_balances.elevation,point_balances.ba);
            s1=polyfit(bin_centers(bin_centers<break_point_elevation),annual_balance_profile_values(bin_centers<break_point_elevation),1);
            s2=polyfit(bin_centers(bin_centers>=break_point_elevation),annual_balance_profile_values(bin_centers>=break_point_elevation),1);
            %ba_gradient={num2str([round(s1(1)*1000,2) round(s2(1)*1000,2)])};
            ba_gradient={num2str([s1(1) s1(2) break_point_elevation s2(1) s2(2)])};%to get lucus gradient info
            Ba=sum(integrated_ba);
            Ba_Date=round(nansum(bin_weights.*predict(fitlm(point_balances.elevation,datenum(point_balances.mass_minimum_date)),bin_centers')));
            ELA=interp1(annual_balance_profile_values,bin_centers,0, 'linear','extrap');
            bals=point_balances;
            point_balances=point_balances(~strcmp(point_balances.site_name,'bp'),:);
        else
            ba_gradient=nan;
            Ba=nan;
            Ba_Date=nan;
            ELA=nan;
        end
        if ~isnan(point_balances.bw)
            % get winter balance gradient and balance
%             point_balances.mass_maximum_date=datenum(point_balances.mass_maximum_date);
%             point_balances.mass_maximum_date=datestr(point_balances.mass_maximum_date);
            pp=splinefit(point_balances.elevation,point_balances.bw,2,2);
            break_point_elevation=round(pp.breaks(1,2));
            break_point_balance=ppval(pp,break_point_elevation);
            point_balances.mass_maximum_date=datenum(point_balances.mass_maximum_date);
            point_balances.mass_minimum_date=datenum(point_balances.mass_minimum_date);
            a=table(point_balances.Year(1),{'bp'},break_point_elevation,break_point_balance,nan,nan,round(interp1(point_balances.elevation,point_balances.mass_maximum_date,break_point_elevation,'spline')),...
                round(interp1(point_balances.elevation,point_balances.mass_minimum_date,break_point_elevation,'spline')),'VariableNames',point_balances.Properties.VariableNames);

            point_balances=[point_balances; a];
            point_balances.mass_maximum_date=datestr(point_balances.mass_maximum_date);
            point_balances.mass_minimum_date=datestr(point_balances.mass_minimum_date);
            point_balances=sortrows(point_balances,{'elevation'});
            pp=splinefit(point_balances.elevation,point_balances.bw,2,2);
            break_point_elevation=round(pp.breaks(1,2));
            break_point_balance=ppval(pp,break_point_elevation);
            winter_balance_profile_values=ppval(pp,bin_centers);
            integrated_bw=winter_balance_profile_values'.*bin_weights;
            
            s1=polyfit(bin_centers(bin_centers<break_point_elevation),winter_balance_profile_values(bin_centers<break_point_elevation),1);
            s2=polyfit(bin_centers(bin_centers>=break_point_elevation),winter_balance_profile_values(bin_centers>=break_point_elevation),1);
            bw_gradient={num2str([round(s1(1)*1000,2) round(s2(1)*1000,2)])};
            bw_gradient={num2str([s1(1) s1(2) break_point_elevation s2(1) s2(2)])};%to get lucus gradient info
            
            Bw=round(nansum(bin_weights.*winter_balance_profile_values'),2);
            Bw_Date=round(nansum(bin_weights.*predict(fitlm(point_balances.elevation,datenum(point_balances.mass_maximum_date)),bin_centers')));
            point_balances=point_balances(~strcmp(point_balances.site_name,'bp'),:);
        else
            bw_gradient=nan;
            Bw=nan;
            Bw_Date=nan;
        end
        if ~isnan(point_balances.bw)&~isnan(point_balances.ba)
            Bs=Ba-Bw;
            bs_gradient={nan};
        else 
            bs_gradient=nan;
            Bs=nan;
        end

        Glacier_Wide_solutions=table(year,Bw,Bs,Ba,ELA,Bw_Date,Ba_Date,bw_gradient,bs_gradient,ba_gradient,'VariableNames',{'Year','Bw_mwe','Bs_mwe','Ba_mwe','ELA_m','Bw_Date','Ba_Date','bw_gradient_mwe_km','bs_gradient_mwe_km','ba_gradient_mwe_km'});       

        if plot_integration==1
            figure(year);hold on 
            if ~isnan(point_balances.bw)
                subplot(1,2,1)
                title('Winter Balance Integration')
                hold on
                yyaxis left
                bar(bin_centers,glacier_hypsometry,1,'FaceColor',[.5 .5 .5])
                ylabel('Area (km^2)')
                yyaxis right
                scatter(point_balances.elevation,point_balances.bw,'b','fill');
                plot(bin_centers,winter_balance_profile_values,'LineWidth',1,'color','b')
                ylabel('balance (m w.e.)')
                ylim([round(min(min(point_balances.bw)))-1 round(max(max(point_balances.bw)))+1])
                xlim([bin_centers(1)-50 bin_centers(end)+50])
                axis square
                box on
            end
            if ~isnan(point_balances.ba)
                subplot(1,2,2)
                title('Annual Balance Integration')
                hold on
                yyaxis left
                bar(bin_centers(1:end),glacier_hypsometry,1,'FaceColor',[.5 .5 .5])
                ylabel('Area (km^2)')
                yyaxis right
                scatter(point_balances.elevation,point_balances.ba,'o','fill');
                plot(bin_centers,annual_balance_profile_values,'LineWidth',1,'color','b')
                ylabel('balance (m w.e.)')
                ylim([round(min(min(point_balances.ba)))-1 round(max(max(point_balances.ba)))+1])
                xlim([bin_centers(1)-50 bin_centers(end)+50])
                axis square
                box on
                set(gcf, 'PaperPosition', [0 0 3 6]);
                print -depsc2 gates_epoch2.eps
            end
        end

    end
end
end
