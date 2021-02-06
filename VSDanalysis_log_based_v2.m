clc
clear all


%Get the radius
%----------------------------
mix_info={'Arnesh_trial_22'} % list as many mixes as needed

for mix_id=1:numel(mix_info)
    mix=string(mix_info(mix_id)) % only change this to load the new samples

    area_sum=[];
    for sample=4:6.
        void_info=readtable(sprintf('%s-%d.csv',mix,sample));
        area_sum=[area_sum;void_info(1:end,2)];
    end
    area_sum=table2array(area_sum);
    radius=(area_sum/pi).^0.5;

    % Logarithmic binning scheme
    min_r=10; max_r=500; % limit of radius
    nb = 20;% number of bins
    wb = (log10(max_r)-log10(min_r))/nb;% wb is the width of bins in log scale
    for i=1:nb+1
        r(i) = 10^(log10(min_r)+(i-1)*wb);
    end
    ab = zeros(1,nb); % mean areas of the bins
    for i=1:nb
        ab(i)=(r(i)+r(i+1))/2;
    end

    % f2d is the number fractions in the histogram of observations
    F2d = histcounts(radius, r)
    figure(1)
    %f1=stairs(r,[F2d,0])
    f1=plot(ab,F2d,'LineWidth',1.5); % concrete scan
    hold on
end

%figure(1)
%f1=stairs(r,[F2d,0])
%f1=plot(ab,F2d,'LineWidth',1.5); % concrete scan
% f1=bar(ab,diag(F2d),1,'stack','edgecolor','black','FaceColor','flat'); % concrete scan
% for k = 1:size(ab,2)
%     f1(k).CData = k;
% end
% colormap(jet);
xlim([0 500]);
ylim([1 10000]);
xlabel(' 2D Apparent void radius, r [µm]');
ylabel('Void count');
set(gca, 'YGrid', 'on', 'XGrid', 'off')
set(gca, 'YScale', 'log')
set(gca,'FontSize',12)
set(gcf,'units','inches','position',[2,2,6,4])

% legend_info=[];
% for i=1:nb
%     legend_info{i}=[mat2str(round(ab(i)))];
% end
% legend(legend_info)
%legend(mix_info)
lgd=legend (mix_info,'Interpreter','none','Location','NorthEast','NumColumns',1);
title(lgd,'Mix name')
