clc
clear all


%Get the radius
%----------------------------
mix_info={'Arnesh_trial_21','Arnesh_trial_22'} % list as many mixes as needed
mix_info_2D=[]
mix_info_3D=[]

%% 2D VSD

for mix_id=1:numel(mix_info)
    mix=string(mix_info(mix_id)) % only change this to load the new samples

    area_sum=[];
    for sample=1:3
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
mix_info_2D=strcat(mix_info,'-2D VSD');

%% 3D VSD

for mix_id=1:numel(mix_info)
    mix=string(mix_info(mix_id)) % only change this to load the new samples

    area_sum=[];
    for sample=1:3
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
    
    p = zeros(1,nb); rm = zeros(1,nb); % ap = zeros(1,nb);
    NA = zeros(nb);NA_recover = zeros(nb); e = zeros(nb); rm(1:end) = r(1:end-1);
    % fA = ab.*f2d; fA = vf.*fA/sum(fA); % Convert to area-weighting
    for j = nb:-1:1
        p(1:j) = (sqrt(r(j+1).^2-rm(1:j).^2)-sqrt(r(j+1).^2-r(2:j+1).^2))/r(j+1);% Eq. 5
        %ap(1:j) = 0.5.*p(1:j).*(rm(1:j)+a(2:j+1)+(a(j+1).*p(1:j).^2)./3);% Eq. 8
        NA(j,1:j) = p(1:j)/sum(p(1:j));% Eq. 9
        F3d(j) = F2d(j)/NA(j,j);% Eq. 10

        % Error correction right after stripping each single bin
        if F3d(j)<0
            e(j,j+1:nb) = F2d(j)./NA(j+1:nb,j); % Eq. 12
            F3d(j) = 0.0;
        end

        F2d(j) = 0.0;%Eq. 11
        F2d = F2d-F3d(j)*NA(j,:);
        NA_recover(j,:) =F3d(j)*NA(j,:);
    end; clear p ap am j %NA  %fA
    F3d
    
    figure(1)
    %f1=stairs(r,[F2d,0])
    f2=plot(ab,F3d,'LineWidth',1.5); % concrete scan
    hold on
end
mix_info_3D=strcat(mix_info,'-3D VSD');

%% plotting

xlim([0 500]);
ylim([10 10000]);
xlabel('Void radius [um]');
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

mix_info_sum = [mix_info_2D, mix_info_3D];
lgd=legend (mix_info_sum,'Interpreter','none','Location','NorthEast','NumColumns',1);
title(lgd,'Mix name')