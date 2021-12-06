figure(5);hold on
scatter(all_ice_ablation_data(all_ice_ablation_data(:,5)~=2020,1),all_ice_ablation_data(all_ice_ablation_data(:,5)~=2020,2),'filled','r');hold on
scatter(all_ice_ablation_data(all_ice_ablation_data(:,5)==2020,1),all_ice_ablation_data(all_ice_ablation_data(:,5)==2020,2),'filled','b');hold on
y_est = all_ice_ablation_data(:,1)*ki;
plot(all_ice_ablation_data(:,1),y_est,'-r') 

 set(gca, 'XColor', 'k', 'YColor', 'k','XAxisLocation','bottom','Xtick',0:300:1200,...
'YAxisLocation','left','Ytick',-9:1:0,'fontname','arial ','fontsize',12,'TickLength',[0.025 0.025],'linewidth',2);
xlim([0 (abs(max_melt).*300)])
ylim([max_melt 0])
xlabel('PDDs (^oC)','fontweight','bold')
ylabel('A_{sfc} (m w.e.)','fontweight','bold')
box on
axis square
set(gcf, 'PaperPositionMode', 'auto');
print -depsc2 gates_epoch2.eps