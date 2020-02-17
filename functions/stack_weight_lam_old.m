function [gauss_weight, gauss_center] = stack_weight_lam_old(f,v,maxwinsize)

%%% S. Pasquet - V17.08.28
% Calculate wavelength based weights to stack multi-window dispersion
%
% gauss_weight = stack_weight_lam_new(f,v,nWvec,resampvec)

a = 2; % decrease a for larger gaussian

[vv,ff] = meshgrid(v,f);
ll = vv./ff;

max_resamp_win = maxwinsize.*10.^(1./sqrt(0.5*maxwinsize));
gauss_center = linspace(0,max(max_resamp_win),length(maxwinsize));
diff_gauss_center = mean(diff(gauss_center));

%%
% close all
gauss_weight = zeros(size(ll,1),size(ll,2),length(gauss_center));
for i = 1:length(gauss_center)
    gauss_tmp = exp(-0.5*((gauss_center(i)-ll)./diff_gauss_center*a).^2);
    if i == length(gauss_center)
        gauss_tmp(ll>max(max_resamp_win)) = 1;
    end
    gauss_weight(:,:,i) = gauss_tmp;
%     plot_img([],vv,ll,gauss_weight(:,:,i),[],1,1,1,[],[],[],[],[],[0 max(max_resamp_win)]);
%     
%     figure(10); plot(ll(:),gauss_tmp(:),'.');  hold on;
%     xlim([0 max(max_resamp_win)]);
end
% hold off;