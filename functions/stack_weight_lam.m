function [w_all,w_sum]=stack_weight_lam(f,v,nWvec,maxlam)

%%% S. Pasquet - V17.08.28
% Calculate wavelength based weights to stack multi-window dispersion
%
% [w_all,w_sum]=stack_weight_lam(f,v,nWvec,maxlam,maxlam_weight,stdgauss)
   
plot_weight=1;
if plot_weight==1
    close all;
end

maxlam_wght = 1;
std_wght = 1;
[vv,ff]=meshgrid(v,f);
ll=vv./ff;
maxlam = maxlam*maxlam_wght;
mid_i = (1+length(nWvec))/2;

if length(nWvec)>1
    w_all=zeros(length(f),length(v),length(nWvec));
    for i=1:length(nWvec)
        if i==1
            stdgauss = (maxlam(i+1)-maxlam(i))*std_wght;
            w_all_temp=(1/(stdgauss*sqrt(2*pi)))*exp(-0.5.*((ll-maxlam(i))/stdgauss).^2);
            w_all_temp=w_all_temp/max(max(w_all_temp));
            w_all_temp(ll<maxlam(i))=max(w_all_temp(:));
        elseif i==length(nWvec)
            stdgauss = (maxlam(i)-maxlam(i-1))*std_wght;
            w_all_temp=(1/(stdgauss*sqrt(2*pi)))*exp(-0.5.*((ll-maxlam(i))/stdgauss).^2);
            w_all_temp=w_all_temp/max(max(w_all_temp));
            w_all_temp(ll>maxlam(i))=max(w_all_temp(:));
        else
            stdgauss = (maxlam(i+1)-maxlam(i-1))*std_wght*0.5;
            w_all_temp=(1/(stdgauss*sqrt(2*pi)))*exp(-0.5.*((ll-maxlam(i))/stdgauss).^2);
            w_all_temp=w_all_temp/max(max(w_all_temp));
            %             stdgauss=stdgauss/(std_wght/abs(i-mid_i));
            %             stdgauss=stdgauss+std_wght*stdgauss*sqrt(abs(maxlam(i)-mean(maxlam))/mean(maxlam));
            %             stdgauss=stdgauss+std_wght*stdgauss/(1/(abs(maxlam(i)-mean(maxlam))/mean(maxlam)));
            if i<mid_i
                w_all_temp(ll<maxlam(i))=(1/(stdgauss*sqrt(2*pi)))*exp(-0.5.*((ll(ll<maxlam(i))-maxlam(i))/stdgauss).^2);
                w_all_temp(ll<maxlam(i))=w_all_temp(ll<maxlam(i))/max(max(w_all_temp(ll<maxlam(i))));
            elseif i>mid_i
                w_all_temp(ll>maxlam(i))=(1/(stdgauss*sqrt(2*pi)))*exp(-0.5.*((ll(ll>maxlam(i))-maxlam(i))/stdgauss).^2);
                w_all_temp(ll>maxlam(i))=w_all_temp(ll>maxlam(i))/max(max(w_all_temp(ll>maxlam(i))));
            end
        end
        w_all(:,:,i) = w_all_temp;
        if plot_weight==1 && length(nWvec)>1
            figure; surf(ff,vv,w_all_temp,'edgecolor','none'); view(0,90);
            figure(100); hold on; plot(vv(:,85)./ff(:,85),w_all_temp(:,85)); xlim([0 50]); hold off;
            figure(101); hold on; plot(ff(:,85),w_all_temp(:,85)); hold off;

        end
    end
    w_sum = sum(w_all,3);
else
    w_sum = ones(length(f),length(v));
    w_all = ones(length(f),length(v));
end

if plot_weight==1 && length(nWvec)>1
    figure; surf(ff,vv,w_sum,'edgecolor','none'); view(0,90); colorbar;
end

