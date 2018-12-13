%% find the best fit radius distribution given the neighbor distribution
clear all; close all; clc
global r_data n_data omega gene_index

load r2.csv
r = r2;
r(:,2:4) = r(:,2:4)/mean(diff(r(:,1)));
r_data = r;
load n.csv
n(:,2:4) = n(:,2:4)/mean(diff(n(:,1)));
n_data = n;
gene_index = 2;

gene={'','spp1','eye53-1','1020hh-1'};

omega = 3.68*pi;
options = optimset('Display','iter');


r_fitted = fminsearch(@radius,r(:,4),options);

r = r_fitted;
r(r<0)=0;
r = r/(sum(r)*mean(diff(r_data(:,1))));
r = r*3;


%% calculate the neighbor and contact distribution using the fitted radius distribution



clear all; close all; clc
load r.csv
r(:,2:4) = r(:,2:4)/mean(diff(r(:,1)));
r_original = r;

load fit_r_from_neighbor.csv
r = fit_r_from_neighbor;

load n.csv
n(:,2:4) = n(:,2:4)/mean(diff(n(:,1)));
load z.csv
z(:,2:4) = z(:,2:4)/mean(diff(z(:,1)));

gene={'','spp1','eye53-1','1020hh-1'};

omega = 3.68*pi;
for jj = 2:4
    subplot(3,3,(jj-2)*3+1); 
    plot(r(:,1), r(:,jj)); hold on; plot(r_original(:,1), r_original(:,jj),'.r'); legend('modified','original'); axis tight;
    subplot(3,3,(jj-2)*3+2);
    [neighbor_calc(jj-1,:) SSE] = neighbor(omega, n(:,[1,jj]),r(:,[1 jj]),1,gene{jj}); axis tight; hold on;
end

pc=[0 0.51 0.41 0.47];
for jj = 2:4
    subplot(3,3,(jj-2)*3+3);
    [contact_calc(jj-1,:) SSE] = contact(omega, z(:,[1 jj]),r(:,[1 jj]), pc(jj),1,gene{jj}); axis tight;hold on;
end

