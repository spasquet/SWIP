function [gauss_weight, gauss_center_new] = stack_weight_lam_new(f,v,maxwinsize,a)

%%% S. Pasquet - V17.08.28
% Calculate wavelength based weights to stack multi-window dispersion
%
% [gauss_weight, gauss_center] = stack_weight_lam_new(f,v,maxwinsize)

% a = 1.5; % increase a for larger gaussian (default 1.5)

[vv,ff] = meshgrid(v,f);
ll = vv./ff;

max_resamp_win = maxwinsize.*10.^(1./sqrt(0.5*maxwinsize));

% gauss_center = linspace(0,max(max_resamp_win),length(maxwinsize));
% diff_gauss_center = mean(diff(gauss_center));

% diff_gauss_center = mean(diff(max_resamp_win));
% gauss_center = max_resamp_win - [2*diff(max_resamp_win)/3 0];

% lam_min = min(maxwinsize);
lam_min = min(max_resamp_win);
lam_max = max(max_resamp_win);

diff_gauss_center = (lam_max - lam_min)/length(max_resamp_win);
sigma = 1/sqrt(2*pi);
% shift = 3*diff_gauss_center*sigma;

% Modif 4 sept 2020
shift=0;

% gauss_center = max_resamp_win;
% gauss_center = maxwinsize + 3*diff_gauss_center*sigma;

% Modif 4 sept 2020
gauss_center = maxwinsize*1.5 + shift;

%%
% close all
gauss_weight = zeros(size(ll,1),size(ll,2),length(gauss_center));
% cm = haxby(length(gauss_center));
for i = 1:length(gauss_center)
    ll_uni = unique(ll);
        
%     gauss_tst = exp(-0.5*((gauss_center(i)-ll_uni)./diff_gauss_center*a).^2);
    

    gauss_tst = 1/(sigma*sqrt(2*pi))*exp(-0.5*((ll_uni-gauss_center(i)+shift)/(sigma*diff_gauss_center*a)).^2);
    gauss_tmp = 1/(sigma*sqrt(2*pi))*exp(-0.5*((ll-gauss_center(i)+shift)/(sigma*diff_gauss_center*a)).^2);
    
    
%     ind_min = find(abs(gauss_tst(ll_uni>gauss_center(i)) - threshold) == min(abs(gauss_tst(ll_uni>gauss_center(i)) - threshold)),1);
%     ll_uni_up = ll_uni(ll_uni>gauss_center(i));
%     gauss_center(i) = gauss_center(i) - (ll_uni_up(ind_min) - max_resamp_win(i));
%     gauss_tmp = exp(-0.5*((gauss_center(i)-ll)./diff_gauss_center*a).^2);
    
    gauss_center_new(i) = gauss_center(i)-shift;
    if i == length(gauss_center)
        gauss_tmp(ll>gauss_center_new(i)) = 1;
        gauss_tst(ll_uni>gauss_center_new(i)) = 1;
    elseif i == 1
        gauss_tmp(ll<gauss_center_new(i)) = 1;
        gauss_tst(ll_uni<gauss_center_new(i)) = 1;
    end
    
    gauss_weight(:,:,i) = gauss_tmp;
    
%     plot_img([],vv,ll,gauss_weight(:,:,i),[],1,1,1,[],[],[],[],[],[0 max(max_resamp_win)]);
%     
%     figure(10); plot(ll_uni,gauss_tst,'-','linewidth',3);  hold on;
%     xlim([0 max(max_resamp_win)]);
%     line([max_resamp_win(i) max_resamp_win(i)],[0 1]);
end


% hold off;