clear; clc

%% load the data
load data_Adams2020.mat

%% pull out data
% E = data.E_rate./1000;
% Eer = data.E_error./1000;
E = data5(:,1)./1000;
Eer = data5(:,2)./1000;

% precip
% P = data.mean_trmm;
% Per = data.std_trmm;
P = data5(:,3);
Per = data5(:,4);

% ksn
% ksn = data.mean_ksn;
% ksner = data.std_ksn;
ksn = data5(:,5);
ksner = data5(:,6);

% ksnq
% ksnq = data.mean_ksn_q;
% ksnqer = data.std_ksn_q;
ksnq = data5(:,7);
ksnqer = data5(:,8);

% slope
% slo = data.mean_gradient;
% sloer = data.std_gradient;
slo = data5(:,9);
sloer = data5(:,10);

%% bin data by precip and pack up
p_bins = [1.22,2.1,2.7];

% calculate mean precip for each bin
Pm1 = mean(P(P<p_bins(1)));
Pm2 = mean(P(P>=p_bins(1) & P<=p_bins(2)));
Pm3 = mean(P(P>=p_bins(2) & P<=p_bins(3)));
Pm4 = mean(P>p_bins(3));

% low precip
E1 = E(P<p_bins(1));
slo1 = slo(P<p_bins(1));
slo1er = sloer(P<p_bins(1));
ksn1 = ksn(P<p_bins(1));
ksn1er = ksner(P<p_bins(1));
ksnq1 = ksnq(P<p_bins(1));
ksnq1er = ksnqer(P<p_bins(1));
dat1 = [E1,slo1,slo1er,ksn1,ksn1er,ksnq1,ksnq1er];
col1 = [0.8 0.1 0.1];

% mid precip
E2 = E(P>=p_bins(1) & P<=p_bins(2));
slo2 = slo(P>=p_bins(1) & P<=p_bins(2));
slo2er = sloer(P>=p_bins(1) & P<=p_bins(2));
ksn2 = ksn(P>=p_bins(1) & P<=p_bins(2));
ksn2er = ksner(P>=p_bins(1) & P<=p_bins(2));
ksnq2 = ksnq(P>=p_bins(1) & P<=p_bins(2));
ksnq2er = ksnqer(P>=p_bins(1) & P<=p_bins(2));
dat2 = [E2,slo2,slo2er,ksn2,ksn2er,ksnq2,ksnq2er];
col2 = [0.95 0.75 0.0];

% mid2high precip
E3 = E(P>=p_bins(2) & P<=p_bins(3));
slo3 = slo(P>=p_bins(2) & P<=p_bins(3));
slo3er = sloer(P>=p_bins(2) & P<=p_bins(3));
ksn3 = ksn(P>=p_bins(2) & P<=p_bins(3));
ksn3er = ksner(P>=p_bins(2) & P<=p_bins(3));
ksnq3 = ksnq(P>=p_bins(2) & P<=p_bins(3));
ksnq3er = ksnqer(P>=p_bins(2) & P<=p_bins(3));
dat3 = [E3,slo3,slo3er,ksn3,ksn3er,ksnq3,ksnq3er];
col3 = [54/255,195/255,201/255];

% high precip
E4 = E(P>p_bins(3));
slo4 = slo(P>p_bins(3));
slo4er = sloer(P>p_bins(3));
ksn4 = ksn(P>p_bins(3));
ksn4er = ksner(P>p_bins(3));
ksnq4 = ksnq(P>p_bins(3));
ksnq4er = ksnqer(P>p_bins(3));
dat4 = [E4,slo4,slo4er,ksn4,ksn4er,ksnq4,ksnq4er];
col4 = [0.1 0.1 0.8];

%% define priors fr parameters interms of mean and std and pack up
K_p = [-25,-5];  %log10 K
%Kp_p = log10(10.^K_p./(mean(P)^(3*0.45))); % log10 Kp
n_p = [1,15];
Sc_p = tand([20,50]); % Sc in degrees
D_p = [-4,-0.1]; % log10 D
priors = [K_p;n_p;Sc_p;D_p];

%% define constants in the model
mn = 0.45; % concavity index used in adams et al. (2020)
pr = 2700; % rock density
ps = 2000; % soil density;
B = pr/ps;



%% run low precip inversion
% pick starting place for the parameters
%param_start = [mean(priors,2);mn;B];

% here near MAP based on ealier model runs
param_start = [-8.8395    2.224    0.602979   -1.659 mn B]';

burn_in = 5e3;
n_iter = 5e4;
[params1,like_data1,mMAP1,accrate1] = HM_MCMC_river_hillslope_ksn(dat1,param_start,priors,burn_in,n_iter);
save('lowP_params_ksn.mat','params1');
save('lowP_like_ksn.mat','like_data1');

param_start = [-8.61313    2.1229    0.5896   -1.6852 mn B]';
[params1q,like_data1q,mMAPq1q,accrate1q] = HM_MCMC_river_hillslope_ksnq(dat1,param_start,priors,burn_in,n_iter);
save('lowP_params_ksnq.mat','params1q');
save('lowP_like_ksnq.mat','like_data1q');

%% plot the data
param_names = {'K','n','Sc','D'};
% % plot the chains
% [r,c] = size(params);
% for i = 1:c
%     subplot(c,1,i)
%     if i == 1 || i ==4
%         plot(10.^params(:,i));
%     elseif i == 3
%         plot(atand(params(:,i)));
%     else
%         plot(params(:,i));
%     end
%     xlabel('iterations');
%     ylabel(param_names{i})
% end


f_num = 86;
plot_data(dat1,col1,f_num)
E = logspace(log10(0.001),log10(10000),500)./1e3;
paramss = [mMAP1,mn,B];
paramssq = [mMAPq1q,mn,B];
plot_model_ksn_ksnq(paramss,paramssq,E,col1,f_num)
disp('low P results:')
disp(param_names)
disp(mMAP)
disp(mMAPq)
disp(' ')

%% mid precip inversion

param_start = [-12.9741    4.1203    0.6088   -1.9064 mn  B]';

[params2,like_data2,mMAP2,accrate2] = HM_MCMC_river_hillslope_ksn(dat2,param_start,priors,burn_in,n_iter);
save('midP_params_ksn.mat','params2');
save('midP_like_ksn.mat','like_data2');

param_start = [-12.9692    4.0276    0.6152   -1.8509 mn  B]';
[params2q,like_data2q,mMAPq2q,accrate2q] = HM_MCMC_river_hillslope_ksnq(dat2,param_start,priors,burn_in,n_iter);
save('midP_params_ksnq.mat','params2q');
save('midP_like_ksnq.mat','like_data2q');

% plot results
f_num = 86;
plot_data(dat2,col2,f_num)
E = logspace(log10(0.001),log10(10000),500)./1e3;
paramss = [mMAP2,mn,B];
paramssq = [mMAPq2q,mn,B];
plot_model_ksn_ksnq(paramss,paramssq,E,col2,f_num)
disp('mid P results:')
disp(param_names)
disp(mMAP)
disp(mMAPq)
disp(' ')


%% mid2high precip inversion
param_start = [-15.5182    5.7144    0.5817   -1.8293 mn  B];
[params3,like_data3,mMAP3,accrate3] = HM_MCMC_river_hillslope_ksn(dat3,param_start,priors,burn_in,n_iter);
save('mid2highP_params_ksn.mat','params3');
save('mid2highP_like_ksn.mat','like_data3');

param_start = [-12.2276    3.7809    0.5776   -1.6665 mn  B];
[params3q,like_data3q,mMAPq3q,accrate3q] = HM_MCMC_river_hillslope_ksnq(dat3,param_start,priors,burn_in,n_iter);
save('mid2highP_params_ksnq.mat','params3q');
save('mid2highP_like_ksnq.mat','like_data3q');

% plot results
f_num = 86;
plot_data(dat3,col3,f_num)
E = logspace(log10(0.001),log10(10000),500)./1e3;
paramss = [mMAP3,mn,B];
paramssq = [mMAPq3q,mn,B];
plot_model_ksn_ksnq(paramss,paramssq,E,col3,f_num)
disp('high P results:')
disp(param_names)
disp(mMAP)
disp(mMAPq)
disp(' ')

%% high precip inversion
param_start = [-20.7927    8.1877    0.5700   -1.8242 mn  B];
[params4,like_data4,mMAP4,accrate4] = HM_MCMC_river_hillslope_ksn(dat4,param_start,priors,burn_in,n_iter);
save('highP_params_ksn.mat','params4');
save('highP_like_ksn.mat','like_data4');

param_start = [-15.7333    5.2526    0.5703   -1.5352 mn  B];
[params4q,like_data4q,mMAPq4q,accrate4q] = HM_MCMC_river_hillslope_ksnq(dat4,param_start,priors,burn_in,n_iter);
save('highP_params_ksnq.mat','params4q');
save('highP_like_ksnq.mat','like_data4q');

% plot results
f_num = 86;
plot_data(dat4,col4,f_num)
E = logspace(log10(0.001),log10(10000),500)./1e3;
paramss = [mMAP4,mn,B];
paramssq = [mMAPq4q,mn,B];
plot_model_ksn_ksnq(paramss,paramssq,E,col4,f_num)
disp('high P results:')
disp(param_names)
disp(mMAP)
disp(mMAPq)
disp(' ')


% set axes
subplot(3,2,1)
xlim([0 5]); ylim([15 40])

subplot(3,2,3)
xlim([0 5]); ylim([0 350])

subplot(3,2,4)
xlim([0 5]); ylim([0 600])

subplot(3,2,5)
xlim([0 500]); ylim([15 40])

subplot(3,2,6)
xlim([0 600]); ylim([15 40])

saveas(gcf,['D:\dream\E vers S\稿件\modified by Sean Gallen 2nd\RE_ RE_ RE_ invitation\results_20231115（用Erosion rate）\1.2,2.0\highPksnq_posterior_cross_plots20231116.pdf']);

%% take a look at the acceptance rates and chains
par_name = '_params_ksn';
like_name = '_like_ksn';

pre_names = {'lowP','midP','mid2highP','highP'};
cols =  [0.8 0.1 0.1; 0.95 0.75 0.0;54/255 195/255 201/255; 0.1 0.1 0.8];

for i = 1:length(pre_names)                              % !!!!!!!!
    
    figure(101)
    subplot(1,2,1)
    load([pre_names{i},like_name,'.mat'])
%     like_data = (['like_data',num2str(i)]);
    plot(like_data4(:,3),'-','color',cols(i,:)); hold on  % !!!!!!!!

    subplot(1,2,2)
    load([pre_names{i},like_name,'q.mat'])
    plot(like_data4q(:,3),'-','color',cols(i,:)); hold on  % !!!!!!!!

    figure(102)
    load([pre_names{i},par_name,'.mat']);
    n=0;
    for j = 1:4
        subplot(4,2,j+n)
        plot(params4(:,j),'-','color', cols(i,:)); hold on  % !!!!!!!!
        n = n +1;
    end
    
    load([pre_names{i},par_name,'q.mat']);
    n=1;
    for j = 1:4
        subplot(4,2,j+n)
        plot(params4(:,j),'-','color', cols(i,:)); hold on  % !!!!!!!!
        n = n+1;
    end
end

figure(101)
subplot(1,2,1)
xlabel('iteration')
ylabel('acceptance rate (%)')
subplot(1,2,2)
xlabel('iteration')

figure(102)
ylabs = {'log10(K)','n','Sc','log10(D)'};
    n=0;
    for j = 1:4
        subplot(4,2,j+n)
        ylabel(ylabs{j})
        n = n +1;
    end
    xlabel('iteration')
subplot(4,2,8)
xlabel('iteration')

saveas(gcf,['D:\dream\E vers S\稿件\modified by Sean Gallen 2nd\RE_ RE_ RE_ invitation\result_20231110(用Erosion rate)\1.1,1.6,2.8\num(S2)20231111.pdf']);

%% matrix plots and determine parameters and uncertainties

par_name = '_params_ksn';
like_name = '_like_ksn';

pre_names = {'lowP','midP','mid2highP','highP'};
cols =  [0.8 0.1 0.1; 0.95 0.75 0.0;54/255 195/255 201/255; 0.1 0.1 0.8];

ksn_tab = nan(3,8);
ksnq_tab = nan(3,8);

fign = 200;
for i = 1:length(pre_names)
    
    figure(fign)
    load([pre_names{i},par_name,'.mat']);
    [S,ax]=plotmatrix(params(burn_in+1:end,:));
    ax(1,1).YLabel.String='log10(K)';
    ax(2,1).YLabel.String='n';
    ax(3,1).YLabel.String='Sc';
    ax(4,1).YLabel.String='log10(D)';
    ax(4,1).XLabel.String='log10(K)';
    ax(4,2).XLabel.String='n';
    ax(4,3).XLabel.String='Sc';
    ax(4,4).XLabel.String='log10(D)';
    title([pre_names{i},' k_{sn}'])

    ksn_tab(i,1) = mean(params(burn_in+1:end,1));
    ksn_tab(i,2) = std(params(burn_in+1:end,1));
    ksn_tab(i,3) = mean(params(burn_in+1:end,2));
    ksn_tab(i,4) = std(params(burn_in+1:end,2));
    ksn_tab(i,5) = mean(params(burn_in+1:end,3));
    ksn_tab(i,6) = std(params(burn_in+1:end,3));
    ksn_tab(i,7) = mean(params(burn_in+1:end,4));
    ksn_tab(i,8) = std(params(burn_in+1:end,4));

    fign = fign+1;
    figure(fign)
    load([pre_names{i},par_name,'q.mat']);
    [S,ax]=plotmatrix(params(burn_in+1:end,:));
    ax(1,1).YLabel.String='log10(K)';
    ax(2,1).YLabel.String='n';
    ax(3,1).YLabel.String='Sc';
    ax(4,1).YLabel.String='log10(D)';
    ax(4,1).XLabel.String='log10(K)';
    ax(4,2).XLabel.String='n';
    ax(4,3).XLabel.String='Sc';
    ax(4,4).XLabel.String='log10(D)';
    title([pre_names{i},' k_{snQ}'])

    ksnq_tab(i,1) = mean(params(burn_in+1:end,1));
    ksnq_tab(i,2) = std(params(burn_in+1:end,1));
    ksnq_tab(i,3) = mean(params(burn_in+1:end,2));
    ksnq_tab(i,4) = std(params(burn_in+1:end,2));
    ksnq_tab(i,5) = mean(params(burn_in+1:end,3));
    ksnq_tab(i,6) = std(params(burn_in+1:end,3));
    ksnq_tab(i,7) = mean(params(burn_in+1:end,4));
    ksnq_tab(i,8) = std(params(burn_in+1:end,4));

    fign = fign+1;

end

% save tables with the estimated parameters

att = {'mean_log10(K)','std_log10(K)','mean_n','std_n',...
    'mean_Sc','std_Sc','mean_log10(D)','std_log10(D)'};
T = array2table(ksn_tab,'VariableNames',att);
writetable(T,'ksn_results.xlsx');

att = {'mean_log10(Kp)','std_log10(Kp)','mean_n','std_n',...
    'mean_Sc','std_Sc','mean_log10(D)','std_log10(D)'};
T = array2table(ksnq_tab,'VariableNames',att);
writetable(T,'ksnq_results.xlsx');