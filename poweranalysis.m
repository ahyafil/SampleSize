%% illustrate power studies (Figure 2 of Barbosa et al, ...)
clear; close all;

figure;
star_color = [1 0 0]; 

%% PANEL A: SIMPLE POWER FOR T-TEST (loosely based on Brysbaert 2019 JCognition figure 1)
d = .4; .6; % effect size (Cohen's d)
nMax = 100; 50; % max number of participants
alpha = 0.05; % significance criterion
z_alpha = -norminv(alpha/2);

PowerA = zeros(1,nMax);
for n=1:nMax
    s = d * sqrt(n); %z-score of measured effect distribution
    PowerA(n) = 100*normcdf( s - z_alpha); % probability of significant effect (in line with true effect)
end
n_star = find(PowerA>80, 1);

h = subplot(221); hold on;
h.Position(3) = .75*h.Position(3);
h.Position(1) = h.Position(1) + .1;
plot(1:nMax, PowerA, 'k');
plot([1 n_star], PowerA(n_star)*[1 1], 'k--');
xlabel('number of participants (n)');
ylabel('power (%)');
ylim([0 100]);
box off;


% inset with power at n
n_inset = [12 n_star];
x_inset = [8 68]; [8 28];
y_inset = [10 50];
dx_inset = 40; 20;
dy_inset = 20;
d_range = -2:.01:3;
for i=1:length(n_inset)
    n = n_inset(i);
    d_boundary = z_alpha/sqrt(n);
    
    n_drange = length(d_range);
    ppdf = normpdf(d_range, d, 1/sqrt(n));
    plot(x_inset(i) + [0 dx_inset], y_inset(i)*[1 1], 'color', .5*[1 1 1], 'linewidth',.5);
    plot((x_inset(i)+ dx_inset*find(d_range>0,1)/n_drange)*[1 1], y_inset(i)+[0 dy_inset], 'color', .5*[1 1 1], 'linewidth',.5); % vertical axis
    signif_points = d_range>d_boundary;
    x_signif = x_inset(i)+find(signif_points)/n_drange*dx_inset;
    
    % significant area
    fill( x_inset(i)+[find(signif_points) fliplr(find(signif_points))]/n_drange*dx_inset,...
        y_inset(i) + [ppdf(signif_points)/max(ppdf)*dy_inset zeros(1,sum(signif_points))], [.7 .7 1], 'linestyle','none');
    
    % distribution of d values
    plot(x_inset(i)+(1:n_drange)/n_drange*dx_inset, ppdf/max(ppdf)*dy_inset+y_inset(i), 'k');
    
    text(x_inset(i)+dx_inset/2, y_inset(i)-3, 'observed d', 'horizontalalignment','center');
    
    % arrow from point
    if n~=n_star
        plot(n, PowerA(n), 'k.','markersize',8);
    end
    arrow([n PowerA(n)], [x_inset(i)+dx_inset/3 y_inset(i)+dy_inset*.8], 6);
end

% significance point
plot(n_star, PowerA(n_star), 'p', 'markersize',10, 'color', star_color, 'markerfacecolor',star_color);
xlim([0 nMax]);
title('standard power analysis');

%%% add subplot with distribution of effect
subplot('Position', [.02 .76 .1 .1]); hold on;
d_range = -2:.01:3;
n_drange = length(d_range);
ppdf = normpdf(d_range, d, 1);
fill( [d_range fliplr(d_range)],[ppdf zeros(1,n_drange)], [.9 .9 .9], 'linestyle','none');
plot([min(d_range) max(d_range)], [0 0], 'color', .5*[1 1 1], 'linewidth',.5);
plot([0 0], [0 max(ppdf)],'color', .5*[1 1 1], 'linewidth',.5);
plot(d*[1 1], [0 max(ppdf)],'--','color', .5*[1 1 1], 'linewidth',.5);
arrow([d normpdf(1)],[d+1 normpdf(1)],4);
%arrow([d+.5 normpdf(.5)],[d-.5 normpdf(.5)],4);
text(d+.5, normpdf(1)-.05,  '$\sigma$', 'horizontalalignment','center','interpreter','latex');
plot(d_range, ppdf, 'k');
text(d,-max(ppdf)*.1,'$\mu$', 'horizontalalignment','center', 'interpreter','latex');
text(mean(xlim),-max(ppdf)*.4,'effect', 'horizontalalignment','center');
text(mean(xlim),-max(ppdf)*1,'effect size:', 'horizontalalignment','center');
text(mean(xlim),-max(ppdf)*1.3,'$d= \mu / \sigma$', 'horizontalalignment','center', 'interpreter','latex');
text(mean(xlim), max(ppdf)*1.4, {'participant-wise','variability'}, 'fontangle','italic', 'horizontalalignment','center');
axis off;

%% PANEL: B POWER ANALYSIS FOR SIMULATED MODEL (DDM)

T = readtable('DDMeffectsize_simulations.csv'); % load model fits performed with PyDDM (already provided, see power_analysis_ddm.py to run again)
x = T.Var2(2:end);
nTot = length(x);

d = .5; % effect size (Cohen's d)
nMax = 50; % max number of participants
alpha = 0.05; % significance criterion
nSub = 500; % number of subsample set per sample size
z_alpha = -norminv(alpha/2);

PowerB = nan(nSub,nMax);
tStat = nan(nSub,nMax);

for n=1:nMax
    for s=1:nSub
        [PowerB(s,n),~,~,stats] = ttest(x(randsample(nTot, n)),0,'Tail','right', 'Alpha',alpha/2); % ttest for subsample of size n of all simulated dataset (in line with true effect)
        tStat(s,n) = stats.tstat ;
    end
end
PowerB_mean = 100*mean(PowerB);
n_star = find(PowerB_mean>80, 1);

h = subplot(222); hold on;
h.Position(3) = .75*h.Position(3);
h.Position(1) = h.Position(1) + .15;
wu(100*PowerB, 'k', 'linewidth',1);
plot([1 n_star], PowerB_mean(n_star)*[1 1], 'k--');
xlabel('number of participants (n)');
ylabel('power (%)');
ylim([0 100]);
box off;


% inset with power at n
n_inset = [12 n_star];
x_inset = [13 29];
y_inset = [10 40];
dx_inset = 20;
dy_inset = 20;
d_range = -2:.01:3;
for i=1:length(n_inset)
    n = n_inset(i);
    d_boundary = z_alpha/sqrt(n);
    
    xx = (-2:7) + (min(tStat(PowerB(:,12)>0,12))-2); % small offset just to have nicer bars
    
    [hh] = hist(tStat(:,n),xx); % histogram of t-statistic for that sample size
    [hh_s] = hist(tStat(PowerB(:,n)>0,n),xx); % ... only for significant
    n_xx = length(xx);
    
    
    plot(x_inset(i) + [0 dx_inset], y_inset(i)*[1 1], 'color', .5*[1 1 1], 'linewidth',.5);
    plot((x_inset(i)+ dx_inset*find(d_range>0,1)/n_drange)*[1 1], y_inset(i)+[0 dy_inset], 'color', .5*[1 1 1], 'linewidth',.5); % vertical axis
    
    % significant area
    xx2 = repelem(1:n_xx+1, [1 2*ones(1,n_xx-1) 1]);
    yy2 = repelem(hh_s,2);
    
    fill( x_inset(i)+[xx2 n_xx 0]/n_xx*dx_inset,...
         y_inset(i) + [yy2/max(hh)*dy_inset 0 0], [.7 .7 1], 'linestyle','none');
    
    % distribution of d values
    stairs(x_inset(i)+(1:n_xx)/n_xx*dx_inset, hh/max(hh)*dy_inset+y_inset(i), 'k');
    
    text(x_inset(i)+dx_inset/2, y_inset(i)-3, 't-stat', 'horizontalalignment','center');
    
    % arrow from point
    if n~=n_star
        plot(n, PowerB_mean(n), 'k.','markersize',8);
    end
    arrow([n PowerB_mean(n)], [x_inset(i)+dx_inset/3 y_inset(i)+dy_inset*.8], 6);
end

% significance point
plot(n_star, PowerB_mean(n_star), 'p', 'markersize',10, 'color', star_color, 'markerfacecolor',star_color);
xlim([0 50]);
title('simulation-based power analysis');

%%% add subplot with distribution of effect
subplot('Position', [.52 .75 .1 .1]); hold on;
d_range = -2:.01:3;
true_d = 0.1;
[hc, xx] = hist(x,20);
fill( [xx fliplr(xx)],[hc zeros(1,length(xx))], [.9 .9 .9], 'linestyle','none');
plot([min(xx) max(xx)], [0 0], 'color', .5*[1 1 1], 'linewidth',.5);
plot([0 0], [0 max(hc)],'color', .5*[1 1 1], 'linewidth',.5);
plot(xx, hc, 'k');
fill(true_d + .05*[0 -1 1], max(hc)*[1.05 1.2 1.2], 'k');
text(mean(xlim),-max(hc)*.2,{'estimated effect'}, 'horizontalalignment','center');
text(mean(xlim),-max(hc)*.5,'$\mu_A - \mu_B$', 'interpreter','latex','horizontalalignment','center');
text(mean(xlim),-max(hc)*.8,'(Drift-Diffusion Model)', 'horizontalalignment','center');
text(mean(xlim), max(hc)*1.6, {'participant-wise','variability'}, 'fontangle','italic', 'horizontalalignment','center');
axis off;


%% PANEL C: BAYESIAN SEQUENTIAL ANALYSIS  (adapted from Keysers et al Figure 3c)
subplot(223); hold on;
title('sequential analysis');
mu = .1; % average log-likelihood in favor of H_A for each participant
sigma = .2; % standard deviation
nMax = 40; % max number of participants
bound =  6; 30;
yl =  10;
n_star = 0;

while n_star<20 % we want an exampe with at least 20 subejcts before reaching bound
    llh = normrnd(mu, sigma, 1, nMax);
    BF = cumsum(llh);
    n_star = find(abs(BF)>log(bound),1); % first sample size reaching bound
end
nMax = 10*ceil((n_star+2)/10);
plot([0 nMax],[1 1], '--', 'color',.5*[1 1 1]);
plot([0 nMax],bound*[1 1], 'k');
plot([0 nMax],1/bound*[1 1], 'k');
fill([0 nMax nMax 0], [bound bound yl yl], .85*[1 1 1], 'linestyle','none');
fill([0 nMax nMax 0], 1./[bound bound yl yl], .85*[1 1 1], 'linestyle','none');
text(1, bound, 'stopping criterion', 'verticalalignment','bottom');

plot(exp(BF(1:n_star)), 'ko', 'markerfacecolor',.5*[1 1 1], 'markersize',5);

% significance point
plot(n_star,1/yl, 'p', 'markersize',10, 'color', star_color, 'markerfacecolor',star_color);

set(gca, 'yscale','log'); ylim([1/yl yl]);
%yticks([1/30 1/10 1/3 1 3 10 30]); yticklabels({'1/30','1/10','1/3','1', '3','10','30'});
yticks([1/10 1/6 1/3 1 3 6 10]); yticklabels({'1/10','1/6','1/3','1', '3','6','10'});

ylabel('Bayes Factor');
xlabel('number of participants (n)');

text(3, exp(1), 'Evidence for H_A', 'verticalalignment','bottom');
arrow([2 exp(1.1)], [2 exp(1.6)],4);
text(3, exp(-1.6), 'Evidence for H_0', 'verticalalignment','bottom');
arrow([2 exp(-1.1)], [2 exp(-1.6)],4);

%% PANEL D: POWER COUNTOURS FOR SELECTING NUMBER OF participantS AND NUMBER OF TRIALS (Figure 1 h of Baker et al. Psy Methods 2020)
nMax = 60; 200; % max number of participants
kMax = 1500; 100; % max number of trials
M = 0.6; 1; % mean effect (equivalent to d before)
sigma_w = 20;30;  % within-participant dispersion
sigma_b = 0.6; % between-participant dispersion
alpha = 0.05; % significance criterion
cmap = [1 1 1]-linspace(0,1,64)'*[.5 .5 0];

PowerD = zeros(nMax, kMax);
for k=1:kMax
    sigma_s = sqrt(sigma_b^2 + sigma_w^2/k); % sample standard deviation
    for n=1:nMax
        s = M * sqrt(n)/sigma_s; %z-score of measured effect distribution
        PowerD(n,k) = 100*normcdf( s - z_alpha); % probability of significant effect (in line with true effect)
    end
end

h = subplot(2,2,4); hold on;
h.Position(3) = .8*h.Position(3);
h.Position(1) = h.Position(1) + .15;


% plot 2D map with contours
imagesc(PowerD');
colormap(cmap);

contour(PowerD', [20 40 60 90], 'color','k');
MM = contour(PowerD', 80*[1 1], 'color','k', 'linewidth',2);
xlabel('number of participants (n)');
ylabel('number of trials (k)');
axis xy;

% select best trade-off
[~,i_min] = min((MM(1,:)/nMax).^2+ (MM(2,:)/kMax).^2); % not sure this is what is used in original method
n_star = MM(1,i_min);
k_star = MM(2,i_min);
plot(n_star, k_star,'p', 'markersize',10, 'color', star_color, 'markerfacecolor',star_color);

hc = colorbar;
set(get(hc, 'Title'), 'String', 'power (%)');

title('Joint power analysis');

%%% add subplot with distribution of effect
subplot('Position', [.52 .26 .1 .1]); hold on;
d_range = -2*sigma_w:.1:sigma_w*2;
n_drange = length(d_range);
ppdf = normpdf(d_range, M, sigma_w);
fill( [d_range fliplr(d_range)],[ppdf zeros(1,n_drange)], [.9 .9 .9], 'linestyle','none');
plot([min(d_range) max(d_range)], [0 0], 'color', .5*[1 1 1], 'linewidth',.5);
plot([0 0], [0 max(ppdf)],'color', .5*[1 1 1], 'linewidth',.5);
plot(d_range, ppdf, 'k');
text(mean(xlim),-max(ppdf)*.2,'effect', 'horizontalalignment','center');
text(mean(xlim), max(ppdf)*1.4, {'trial-wise','variability'}, 'fontangle','italic', 'horizontalalignment','center');
axis off;

%% add panel letters and arrows

% create invisible panels
axes('position',[0 0 1,1],'visible','off');
text(.02, .95, 'a','fontweight','bold');
text(.52, .95, 'b','fontweight','bold');
text(.02, .45, 'c','fontweight','bold');
text(.52, .45, 'd','fontweight','bold');
xlim([0 1]);
ylim([0 1]);

% add arrows
% h_a = arrow([.13 .77], [.18 .75],8);
% h_a = arrow([.62 .77], [.67 .75],8);
% h_a = arrow([.62 .28], [.67 .25],8);


