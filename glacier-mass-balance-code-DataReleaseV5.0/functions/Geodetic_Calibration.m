function [Calibrated_Glacier_Wide_solutions] = Geodetic_Calibration(glacier,geodetic_data,Glacier_Wide_solutions,Geodetic_Calibration_index,break_points,data_for_calibration)
%% Primary Objective: Geodetically Calibrate Glaciological Time-series to remove systematic bias in the glaciological time-series
%       A) Correct all geodetic measurements to the stratigraphic end of
%           the mass balance year from the date of aquisition
%           1) Apply 'Find_Mass_Maximum_and_Minimum_Adjustments.m' to find
%               the change in mass balance at glaciological site for year
%               of aquisition
%           2) Integrate over the glacier hypsometry.
%               This is done for consistency with glaciological
%               measurements.
%           3) Add change in mass balance (+/-) to geodetic measurement to
%               correct to the end of mass balance year in selected time_system
%       B) Calibrate Glacier-wide glaciological solutions to corrected
%           geodetic measurements to remove systematic bias from glaciological
%           data
%           1) if no calibration is select, no calibration is performed
%           2) if Best-fit is selected, a single best-fit line is applied 
%               to the residual difference of cumulative Glacier-wide Ba time-series 
%               and the Geodetic time-series as a function of time. The slope of the 
%               best-fit line is the linear rate of systematic bias per
%               year, which is then removed from the glaciological
%               time-series 
%           3) if Piece-wise is selected, an absolute geodetic calibration
%               is applied for each measurement interval (n:n+1) of the
%               geodetic data provided
%           4) if Break-point is selected user will have been prompted in
%               mb_GUI to enter in year(s) to break the glaciological
%               time-series at
%               
%% A) Calibrate Glacier-wide glaciological solutions to corrected geodetic measurements to remove systematic bias from glaciological data
%           1) if no calibration is select, no calibration is performed
%           2) if Best-fit is selected, a single best-fit line is applied 
%               to the residual difference of cumulative Glacier-wide Ba time-series 
%               and the Geodetic time-series as a function of time. The slope of the 
%               best-fit line is the linear rate of systematic bias per
%               year, which is then removed from the glaciological
%               time-series 
%           3) if Piece-wise is selected, an absolute geodetic calibration
%               is applied for each measurement interval (n:n+1) of the
%               geodetic data provided
%           4) if Break-point is selected user will have been prompted in
%               mb_GUI to enter in year(s) to break the glaciological
%               time-series at

Glacier_Wide_solutions=[Glacier_Wide_solutions(1,:); Glacier_Wide_solutions];
Glacier_Wide_solutions.Year(1)=Glacier_Wide_solutions.Year(1)-1;
Glacier_Wide_solutions.Bw_mwe(1)=NaN;
Glacier_Wide_solutions.Bs_mwe(1)=NaN;
Glacier_Wide_solutions.Ba_mwe(1)=0;
geodetic_data=geodetic_data(str2num(datestr(geodetic_data.Date,'yyyy'))>=Glacier_Wide_solutions.Year(1)&str2num(datestr(geodetic_data.Date,'yyyy'))<=Glacier_Wide_solutions.Year(end),:);
sources=unique(geodetic_data.Source);
if Geodetic_Calibration_index==1%no geodetic calibration applied
    Calibrated_Glacier_Wide_solutions=[Glacier_Wide_solutions table(zeros(height(Glacier_Wide_solutions),1),'VariableNames',{'Calibration'})];
    calibrated_Ba=Glacier_Wide_solutions.Ba_mwe;
    calibrated_Bw=Glacier_Wide_solutions.Bw_mwe;
    calibrated_Bs=Glacier_Wide_solutions.Bs_mwe;
    calibration=zeros(height(Glacier_Wide_solutions),1);
elseif Geodetic_Calibration_index==2 %best fit calibration based on all the geodetic data
%     calibration_table=nan*ones(length(sources),1);
    cumulative_Ba=cumsum(Glacier_Wide_solutions.Ba_mwe);
    for source=1:length(sources)
        if strcmp(sources(source),'DEM')
        residuals=[];
            data=[];
            data=geodetic_data(strcmp(geodetic_data.Source,sources(source)),:);
            if height(data)==2
                geodetic_uncertainty_weights=ones(height(data),1);
            else
                geodetic_uncertainty_weights = 1-(data.Uncertainty - min(data.Uncertainty)) / ( max(data.Uncertainty) - min(data.Uncertainty) );
                geodetic_uncertainty_weights(1)=1;
            end
            for measurement=1:height(data)
                residuals(measurement,1)=str2num(datestr(data.Date(measurement,:),'yyyy'));
                residuals(measurement,2)=data.Mass_Change_mwe(measurement);
                residuals(measurement,3)=cumulative_Ba(Glacier_Wide_solutions.Year==str2num(datestr(data.Date(measurement,:),'yyyy')),1);
                residuals(measurement,4)=residuals(measurement,2)-residuals(measurement,3);
            end
            residuals(:,1)=residuals(:,1)-residuals(1,1);
            residuals(:,4)=residuals(:,4)-residuals(1,4);
            residuals_model=fitlm(residuals(:,1),residuals(:,4),'linear','weights',geodetic_uncertainty_weights);
%             calibration_table(source,1)=residuals_model.Coefficients{2,1};  
    
    
            calibrated_Ba=Glacier_Wide_solutions.Ba_mwe+residuals_model.Coefficients{2,1};
            calibrated_Bw=Glacier_Wide_solutions.Bw_mwe+(residuals_model.Coefficients{2,1}/2);
            calibrated_Bs=Glacier_Wide_solutions.Bs_mwe+(residuals_model.Coefficients{2,1}/2);
        else
            calibrated_Ba=Glacier_Wide_solutions.Ba_mwe;
            calibrated_Bw=Glacier_Wide_solutions.Bw_mwe;
            calibrated_Bs=Glacier_Wide_solutions.Bs_mwe;
        end
    end
    calibration=(residuals_model.Coefficients{2,1}.*ones(height(Glacier_Wide_solutions),1)).*-1;
%     Calibrated_Glacier_Wide_solutions=[Glacier_Wide_solutions table(residuals_model.Coefficients{2,1}.*ones(height(Glacier_Wide_solutions),1),'VariableNames',{'Calibration'})]
elseif Geodetic_Calibration_index==3 %piece-wise calibration per measurement interval
    calibration=nan*ones(length(sources),1);
    for source=1:length(sources)
            residuals=[];
            data=[];
            data=geodetic_data(strcmp(geodetic_data.Source,sources(source)),:);
%             geodetic_uncertainty_weights=1./(data.Uncertainty+1e-10);
%             geodetic_uncertainty_weights = 1-(data.Uncertainty - min(data.Uncertainty)) / ( max(data.Uncertainty) - min(data.Uncertainty) );
%             geodetic_uncertainty_weights = 1-(data.Uncertainty - min(data.Uncertainty)) / (15 - min(data.Uncertainty) );
            geodetic_uncertainty_weights = min(data.Uncertainty)./1
            calibrated_Ba=[];
            calibrated_Bw=[];
            calibrated_Bs=[];
            if height(data)>=3 && strcmp(sources(source),'DEM')
                for measurement=2:height(data)
                    year_1=str2num(datestr(data.Date(measurement-1,:),'yyyy'));
                    year_2=str2num(datestr(data.Date(measurement,:),'yyyy'));
                    Geodetic_mass_change=data.Mass_Change_mwe(measurement)-data.Mass_Change_mwe(measurement-1);
                    glaciological_mass_change=sum(Glacier_Wide_solutions.Ba_mwe(Glacier_Wide_solutions.Year<=year_2&Glacier_Wide_solutions.Year>year_1));
                    period_calibration=(glaciological_mass_change-Geodetic_mass_change)/(year_2-year_1);
                    if measurement==2 && measurement~=height(data)
                        calibrated_Ba=Glacier_Wide_solutions.Ba_mwe(Glacier_Wide_solutions.Year<=year_2)-period_calibration;
                        calibrated_Bw=Glacier_Wide_solutions.Bw_mwe(Glacier_Wide_solutions.Year<=year_2)-(period_calibration/2);
                        calibrated_Bs=Glacier_Wide_solutions.Bs_mwe(Glacier_Wide_solutions.Year<=year_2)-(period_calibration/2);
                    elseif measurement>2 && measurement~=height(data)
                        calibrated_Ba=[calibrated_Ba;Glacier_Wide_solutions.Ba_mwe(Glacier_Wide_solutions.Year>year_1 & Glacier_Wide_solutions.Year<=year_2)-period_calibration];
                        calibrated_Bw=[calibrated_Bw;Glacier_Wide_solutions.Bw_mwe(Glacier_Wide_solutions.Year>year_1 & Glacier_Wide_solutions.Year<=year_2)-(period_calibration/2)];
                        calibrated_Bs=[calibrated_Bs;Glacier_Wide_solutions.Bs_mwe(Glacier_Wide_solutions.Year>year_1 & Glacier_Wide_solutions.Year<=year_2)-(period_calibration/2)];
                    elseif measurement>2 && measurement==height(data)
                        calibrated_Ba=[calibrated_Ba;Glacier_Wide_solutions.Ba_mwe(Glacier_Wide_solutions.Year>year_1)-period_calibration];
                        calibrated_Bw=[calibrated_Bw;Glacier_Wide_solutions.Bw_mwe(Glacier_Wide_solutions.Year>year_1)-(period_calibration/2)];
                        calibrated_Bs=[calibrated_Bs;Glacier_Wide_solutions.Bs_mwe(Glacier_Wide_solutions.Year>year_1)-(period_calibration/2)];
                    elseif measurement==2 && measurement==height(data)
                        calibrated_Ba=[calibrated_Ba;Glacier_Wide_solutions.Ba_mwe-period_calibration];
                        calibrated_Bw=[calibrated_Bw;Glacier_Wide_solutions.Bw_mwe-(period_calibration/2)];
                        calibrated_Bs=[calibrated_Bs;Glacier_Wide_solutions.Bs_mwe-(period_calibration/2)];
               
                    end
                end
            elseif height(data)<=2 && strcmp(sources(source),'DEM')
                print('WARNNING: Atleast two geodetic measurements are needed from a piece-wise calibration')
                cumulative_Ba=cumsum(Glacier_Wide_solutions.Ba_mwe);
                for source=1:length(sources)
                    if strcmp(sources(source),'DEM')
                    residuals=[];
                        data=[];
                        data=geodetic_data(strcmp(geodetic_data.Source,sources(source)),:);
                        if height(data)==2
                            geodetic_uncertainty_weights=ones(height(data),1);
                        else
                            geodetic_uncertainty_weights = 1-(data.Uncertainty - min(data.Uncertainty)) / ( max(data.Uncertainty) - min(data.Uncertainty) );
                            geodetic_uncertainty_weights(1)=1;
                        end
                        for measurement=1:height(data)
                            residuals(measurement,1)=str2num(datestr(data.Date(measurement,:),'yyyy'));
                            residuals(measurement,2)=data.Mass_Change_mwe(measurement);
                            residuals(measurement,3)=cumulative_Ba(Glacier_Wide_solutions.Year==str2num(datestr(data.Date(measurement,:),'yyyy')),1);
                            residuals(measurement,4)=residuals(measurement,2)-residuals(measurement,3);
                        end
                        residuals(:,1)=residuals(:,1)-residuals(1,1);
                        residuals(:,4)=residuals(:,4)-residuals(1,4);
                        residuals_model=fitlm(residuals(:,1),residuals(:,4),'linear','weights',geodetic_uncertainty_weights);
            %             calibration_table(source,1)=residuals_model.Coefficients{2,1};  


                        calibrated_Ba=Glacier_Wide_solutions.Ba_mwe+residuals_model.Coefficients{2,1};
                        calibrated_Bw=Glacier_Wide_solutions.Bw_mwe+(residuals_model.Coefficients{2,1}/2);
                        calibrated_Bs=Glacier_Wide_solutions.Bs_mwe+(residuals_model.Coefficients{2,1}/2);
                    else
                        calibrated_Ba=Glacier_Wide_solutions.Ba_mwe;
                        calibrated_Bw=Glacier_Wide_solutions.Bw_mwe;
                        calibrated_Bs=Glacier_Wide_solutions.Bs_mwe;
                    end
                end
            else
                
            end 
    end
    if length(sources)==1
        calibration=Glacier_Wide_solutions.Ba_mwe-calibrated_Ba;
    elseif length(sources)==2
        calibration=Glacier_Wide_solutions.Ba_mwe-calibrated_Ba;
    elseif isempty(sources)
        calibration=0;
    end
%     Calibrated_Glacier_Wide_solutions=[Glacier_Wide_solutions table(calibration,'VariableNames',{'Calibration'})];
elseif Geodetic_Calibration_index==4 %break fpoint calibration. User defined year(s) after which a new calibration should be applied
    calibration=nan*ones(length(sources),1);%vector of calibrations to be populated moving forward
    cumulative_Ba=cumsum(Glacier_Wide_solutions.Ba_mwe);
    for source=1:length(sources)
            residuals=[];
            data=[];
            data=geodetic_data(strcmp(geodetic_data.Source,sources(source)),:);
%             geodetic_uncertainty_weights=1./(data.Uncertainty+1e-10);
%             geodetic_uncertainty_weights = 1-(data.Uncertainty - min(data.Uncertainty)) / ( max(data.Uncertainty) - min(data.Uncertainty) );
%              geodetic_uncertainty_weights = 1-(data.Uncertainty - min(data.Uncertainty)) / (15 - min(data.Uncertainty) );
 %           geodetic_uncertainty_weights = 1./(data.Uncertainty).^2;
            geodetic_uncertainty_weights = min(data.Uncertainty(data.Uncertainty~=0))./data.Uncertainty
            geodetic_uncertainty_weights(geodetic_uncertainty_weights==Inf)=2;
%             geodetic_uncertainty_weights(1)=1;
            geodetic_years=str2num(datestr(data.Date,'yyyy'));
            calibrated_Ba=Glacier_Wide_solutions.Ba_mwe;
            calibrated_Bw=Glacier_Wide_solutions.Bw_mwe;
            calibrated_Bs=Glacier_Wide_solutions.Bs_mwe;
            if height(data)>=2 && strcmp(sources(source),'DEM')
%                 if length(break_points)==1
                    break_points=break_points(break_points>Glacier_Wide_solutions.Year(2));
                    break_points=[Glacier_Wide_solutions.Year(1);break_points';Glacier_Wide_solutions.Year(end)];
                    for index=2:length(break_points)
                        cumulative_Ba=cumsum(calibrated_Ba);%cumulative Ba time series to compare to geodetic solutions
                         residuals=[];
                         data.Mass_Change_mwe=data.Mass_Change_mwe-(data.Mass_Change_mwe(1)-cumulative_Ba(Glacier_Wide_solutions.Year==geodetic_years(1)));
                         residuals(:,1)=geodetic_years(geodetic_years>break_points(index-1)&geodetic_years<=break_points(index));
                         residuals(:,2)=data.Mass_Change_mwe(geodetic_years>break_points(index-1)&geodetic_years<=break_points(index));
                         for i=1:length(residuals(:,1))  
                             residuals(i,3)=cumulative_Ba(Glacier_Wide_solutions.Year==residuals(i));
                         end
                         residuals(:,4)=residuals(:,2)-residuals(:,3);                         
                         residuals(:,5)=geodetic_uncertainty_weights(geodetic_years>break_points(index-1)&geodetic_years<=break_points(index));
                        if index>=2
                             residuals=[break_points(index-1) cumulative_Ba(Glacier_Wide_solutions.Year==break_points(index-1)) 0 0 max(residuals(:,5));residuals];
                         end  
                         residuals(:,1)=residuals(:,1)-residuals(1,1);
                         if residuals(1,4)~=0
                            residuals(:,4)=residuals(:,4)-residuals(1,4);
                         end
                         residuals_model=fitlm(residuals(:,1),residuals(:,4),'linear','weights',residuals(:,5));
                         calibrated_Ba(Glacier_Wide_solutions.Year>break_points(index-1))=calibrated_Ba(Glacier_Wide_solutions.Year>break_points(index-1))+residuals_model.Coefficients{2,1};
                         calibrated_Bw(Glacier_Wide_solutions.Year>break_points(index-1))=calibrated_Bw(Glacier_Wide_solutions.Year>break_points(index-1))+(residuals_model.Coefficients{2,1}/2);
                         calibrated_Bs(Glacier_Wide_solutions.Year>break_points(index-1))=calibrated_Bs(Glacier_Wide_solutions.Year>break_points(index-1))+(residuals_model.Coefficients{2,1}/2);

                    end
            elseif height(data)<2 && strcmp(sources(source),'DEM')
                print('WARNNING!!! Atleast three geodetic measurements are needed from a break-point calibration')
            else
                
            end
    if length(sources)==1
        calibration=Glacier_Wide_solutions.Ba_mwe-calibrated_Ba;
    elseif length(sources)==2
        calibration=Glacier_Wide_solutions.Ba_mwe-calibrated_Ba;
    elseif isempty(sources)
        calibration=0;
    end
    end
    end
    if length(sources)==1
%         calibration=Glacier_Wide_solutions.Ba_mwe-calibrated_Ba;
        Glacier_Wide_solutions.Ba_mwe=round(calibrated_Ba,2);       
        Glacier_Wide_solutions.Bw_mwe=round(calibrated_Bw,2);
        Glacier_Wide_solutions.Bs_mwe=round(calibrated_Bs,2);
    elseif length(sources)==2
%         calibration=Glacier_Wide_solutions.Ba_mwe-calibrated_Ba;
        Glacier_Wide_solutions.Ba_mwe=round(calibrated_Ba,2);
        Glacier_Wide_solutions.Bw_mwe=round(calibrated_Bw,2);
        Glacier_Wide_solutions.Bs_mwe=round(calibrated_Bs,2);
    elseif isempty(sources)
%         calibration=0;
    end
    Glacier_Wide_solutions.Bw_mwe(1)=NaN;
    Glacier_Wide_solutions.Bs_mwe(1)=NaN;
    Glacier_Wide_solutions.Ba_mwe(1)=0;
    Calibrated_Glacier_Wide_solutions=[Glacier_Wide_solutions table(calibration,'VariableNames',{'Calibration'})];
% end
end


