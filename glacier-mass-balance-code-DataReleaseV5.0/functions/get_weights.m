function [zone_weights,meanElev] = get_weights(bin_centers,gl_area,zees)
% This function calculates the weight for each zone of a glacier with any
% number of mass balance stakes. Weights are then used to scale measured 
% balances, yielding the glacier wide average, or specific balance.
%
% inputs: Area Altitude Distribution information, in 100 m bins from
% terminus up to the head of the glacier
% (1) bin_centers - a vector of AAD bin median values in meters
% (2) gl_area - area of the glacier in each of the bins in (1), in km^2
% (3) zees - surface altitude of each index site for a given balance year
% if it is zero, means that shouldn't use that site in balance zones
%
% output: zone_weights - percentage of glacier in zones represented by each
% stake subsequently used to estimate specific balance.  
%%
dbstop if error
zone_weights=NaN*ones(length(zees),1); %initialize
iextra=find(zees==0|zees==nan); %stakes not using in balance 
iuse=find(zees>0); %stakes using
zees=zees(iuse);
%    
binWid0=100; %100 m for now, may change to 30 m
total_area = sum(gl_area);
gl_area = gl_area.';
bin_centers = bin_centers.';
binWid=binWid0*ones(length(bin_centers),1);
binWid(1)=(bin_centers(2)-binWid0/2-bin_centers(1))*2;
binWid(end)=(bin_centers(end)-binWid0/2-bin_centers(end-1))*2;
%number of zones need to deal with is one less than number of sites
%because last zone weight is computed from residual
numZones=length(zees)-1;
%mean elevation of each zone
meanElev= (zees(1:end-1)+zees(2:end))/2;
%parse boundary bins
numBins= length(bin_centers);
dist=NaN*ones(numBins,numZones);
for i=1:numZones
    %distance from midpoints to each bin center, one col per zone
    dist(:,i)= abs(meanElev(i) - bin_centers);
end
%Identify bins split between zones
% min is row vector containing the minimum element from each column
[minD,shareRow]= min(dist); %keep indices, shareRow vector numZonesX1
for i=1:numZones-1
    if shareRow(i)==shareRow(i+1) %should have smaller bins so divide problem bin
        s=shareRow(i);
        wid=binWid(s)/3;
        add=[bin_centers(s)-binWid(s)/2+wid;bin_centers(s)+binWid(s)/2-wid];
        bin_centers=[bin_centers(1:s-1);add;bin_centers(s+1:end)];
        gl_area = [gl_area(1:s-1);gl_area(s)/2;gl_area(s)/2;gl_area(s+1:end)];
        binWid=[binWid(1:s-1);2*wid;2*wid;binWid(s+1:end)];
        %parse boundary bins
        numBins= length(bin_centers);
        dist=NaN*ones(numBins,numZones);
        for j=1:numZones
            %distance from midpoints to each bin center, one col per zone
            dist(:,j)= abs(meanElev(j) - bin_centers);
        end
        %Identify bins split between zones
        % min is row vector containing the minimum element from each column
        [minD,shareRow]= min(dist); %keep indices, shareRow vector numZonesX1
    end
end 

%bounds for split bins, one row for each zone
split_bnd = [bin_centers(shareRow)-binWid(shareRow)/2  bin_centers(shareRow)+binWid(shareRow)/2];
%partitioning ratios of split bins
%stop(here)
give_above=[0;gl_area(shareRow).*((split_bnd(:,2) - meanElev)./binWid(shareRow))];
give_below=gl_area(shareRow) - give_above(2:end);
%
limits=[0,shareRow];
part=NaN*ones(1,numZones);
index =[(limits(1:end-1)+1).',(limits(2:end)-1).']; %index of bins in each zone
for i=1:numZones %area of zone
    part(i)=give_above(i)+sum(gl_area(index(i,1):index(i,2)))+give_below(i);
end
upp_part=total_area-sum(part);% residual gives upper zone
% calculate desired weights by scaling by total area
zone_weights0 = [part upp_part]./total_area; 
zone_weights(iuse)=zone_weights0;
zone_weights(iextra)=0; %zero extra stakes



