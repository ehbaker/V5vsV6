clear 
close all

data_directory='data\';
glaciers={'Gulkana';'Wolverine';'LemonCreek';'SouthCascade';'Sperry'};                                    %list all processed (coregistered DEMs for directory/glacier

for glacier=1:length(glaciers)
    glaciological_solutions=readtable([data_directory,char(glaciers(glacier)),'\Output\Output_',char(glaciers(glacier)),'_Glacier_Wide_solutions_calibrated.csv']);
    geodetic_solutions=readtable([data_directory,char(glaciers(glacier)),'\Output\Output_',char(glaciers(glacier)),'_Geodetics_Adjusted_to_Mass_Minimum.csv']);
    geodetic_solutions=geodetic_solutions(strcmp(geodetic_solutions.Source,'DEM')&str2num(datestr(geodetic_solutions.Date,'yyyy'))>=glaciological_solutions.Year(1),:);
    geodetic_solutions.Mass_Change_mwe=geodetic_solutions.Mass_Change_mwe-geodetic_solutions.Mass_Change_mwe(1);
    geodetic_solutions.Mass_Change_mwe=geodetic_solutions.Mass_Change_mwe+sum(glaciological_solutions.Ba_mwe(glaciological_solutions.Year<=str2num(datestr(geodetic_solutions.Date(1,:),'yyyy'))));
figure(1);hold on
if strcmp(glaciers(glacier),'Gulkana')
subplot(5,1,1); hold on
axis([1950 2020 -40 10]); 
set(gca, 'YGrid','on', 'XTickLabel',[], 'Position',[.1 .81 .82 .18], 'Box', 'on');
text(1952,-12,'Gulkana Glacier', 'color', 'k'); hold on
plot(glaciological_solutions.Year,cumsum(glaciological_solutions.Ba_mwe), '-o', 'markersize', 4, 'linewidth', 1,'markerfacecolor', [1 .97 .8],'markeredgecolor',[.7 .7 .7]); hold on
errorbar(str2num(datestr(geodetic_solutions.Date,'yyyy')),geodetic_solutions.Mass_Change_mwe,geodetic_solutions.Uncertainty,'sk','linewidth', 1,'markersize', 8);hold on
legend('glaciological','geodetic');
elseif strcmp(glaciers(glacier),'Wolverine')
subplot(5,1,2);
plot(glaciological_solutions.Year,cumsum(glaciological_solutions.Ba_mwe), '-o', 'markersize', 4, 'linewidth', 1,'markerfacecolor', [1 .97 .8],'markeredgecolor',[.7 .7 .7]); hold on
errorbar(str2num(datestr(geodetic_solutions.Date,'yyyy')),geodetic_solutions.Mass_Change_mwe,geodetic_solutions.Uncertainty,'sk','linewidth', 1,'markersize', 8);
text(1954,-12,'Wolverine Glacier', 'color', 'k');
axis([1953 2020 -40 10]);  
set(gca, 'YGrid','on','XTickLabel',[], 'Position',[.1 .63 .82 .18]);
% 
elseif strcmp(glaciers(glacier),'LemonCreek')
subplot(5,1,3);
plot(glaciological_solutions.Year,cumsum(glaciological_solutions.Ba_mwe), '-o', 'markersize', 4, 'linewidth', 1,'markerfacecolor', [1 .97 .8],'markeredgecolor',[.7 .7 .7]); hold on
errorbar(str2num(datestr(geodetic_solutions.Date,'yyyy')),geodetic_solutions.Mass_Change_mwe,geodetic_solutions.Uncertainty,'sk','linewidth', 1,'markersize', 8);
text(1954,-12,'Lemon Creek Glacier', 'color', 'k');
axis([1953 2020 -40 10]);  
set(gca, 'YGrid','on','XTickLabel',[], 'Position',[.1 .45 .82 .18]);
ylabel('\Sigma B_a [m w.e.]');

elseif strcmp(glaciers(glacier),'SouthCascade')
subplot(5,1,4);
plot(glaciological_solutions.Year,cumsum(glaciological_solutions.Ba_mwe), '-o', 'markersize', 4, 'linewidth', 1,'markerfacecolor', [1 .97 .8],'markeredgecolor',[.7 .7 .7]); hold on
errorbar(str2num(datestr(geodetic_solutions.Date,'yyyy')),geodetic_solutions.Mass_Change_mwe,geodetic_solutions.Uncertainty,'sk','linewidth', 1,'markersize', 8);
text(1954,-12,'South Cascade Glacier', 'color', 'k');
axis([1953 2020 -40 10]);  
set(gca, 'YGrid','on','Position',[.1 .27 .82 .18]);

elseif strcmp(glaciers(glacier),'Sperry')
subplot(5,1,5);
plot(glaciological_solutions.Year,cumsum(glaciological_solutions.Ba_mwe), '-o', 'markersize', 4, 'linewidth', 1,'markerfacecolor', [1 .97 .8],'markeredgecolor',[.7 .7 .7]); hold on
errorbar(str2num(datestr(geodetic_solutions.Date,'yyyy')),geodetic_solutions.Mass_Change_mwe,geodetic_solutions.Uncertainty,'sk','linewidth', 1,'markersize', 8);
text(1954,-12,'Sperry Glacier', 'color', 'k');
axis([1953 2020 -40 10]);    
set(gca, 'YGrid','on','Position',[.1 .09 .82 .18]);

xlabel('Time [year]');
end
end