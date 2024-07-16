function [xzvfinal, xi, zi, vi_filt] = vel2xzvtomo(xx,xzv,F2,min_vel,max_vel,showplot)

Xfinal=zeros(size(xzv));
Zfinal=Xfinal;
for ii = 1:length(xzv)
    Xfinal(ii) = xzv(ii,1).xz(1);
    Zfinal(ii) = xzv(ii,1).xz(2);
end

[mi,ni] = size(xx);
xi=reshape(Xfinal,mi-1,ni-1);
zi=reshape(Zfinal,mi-1,ni-1);
vi=F2(xi,zi);

if ~isempty(max_vel)
    vi(vi>max_vel) = max_vel;
end
if ~isempty(min_vel)
    vi(vi<min_vel) = min_vel;
end
vi_filt = median_filt_2D(vi',15,1);
vi_filt = vi_filt';

% for i=1:size(vi_filt,1)
%     ind_nan_first = find(~isnan(vi_filt(i,:)),1,'first') - 1;
%     ind_nan_last = find(~isnan(vi_filt(i,:)),1,'last') + 1;
%     if ind_nan_first >= 1
%         vi_filt(i,1:ind_nan_first) = vi_filt(i,ind_nan_first + 1);
%     end
%     if ind_nan_last <= size(vi_filt,2)
%         vi_filt(i,ind_nan_last:end) = vi_filt(i,ind_nan_last - 1);
%     end
% end

if ~isempty(max_vel) && ~isempty(min_vel)
    
    %     vel_mean = mean(vi_filt,2);
    %     vel_std = std(vi_filt,[],2);
    %     vel_min = min(vi_filt,[],2);
    %     vel_max= max(vi_filt,[],2);
    %
    %     a_max = 4;
    %     a_min = 1/(a_max/2);
    %     distrand=rand(1)-0.5;
    %     distrand=(distrand)/0.5;
    %     a_mean=0.5*(log10(a_min)+log10(a_max));
    %     a_diff=0.5*(log10(a_max)-log10(a_min));
    %     a_all=a_mean+distrand*a_diff;
    %     a_all = 10.^a_all;
    %
    %     b_max = std(vi_filt(:));
    %     b_min = min_vel;
    %     distrand=rand(1)-0.5;
    %     distrand=(distrand)/0.5;
    %     b_mean=0.5*(b_min+b_max);
    %     b_diff=0.5*(b_max-b_min);
    %     b_all=b_mean+distrand*b_diff;
    %     vi_filt_mean = b_all+a_all*repmat(vel_mean,1,size(vi_filt,2));

    vi_1D = mov_aver(nanmax(vi_filt,2)',5,1,size(vi_filt,1));
    vi_filt_mean = repmat(vi_1D,1,size(vi_filt,2));
else
    vi_filt_mean = vi_filt;
end

finalVel=reshape(vi_filt_mean,size(xi,1)*size(xi,2),1);
xzvfinal=xzv;
for ii = 1:length(xzvfinal)
    xzvfinal(ii,1).v = finalVel(ii);
    xzvfinal(ii,1).u = finalVel(ii)^-1;
end

if showplot == 1
    plot_img([],xi,zi,vi_filt,haxby(32),1,0,1,12,'X (m)',...
        'Altitude (m)','Velocity (m/s)',[],[],[min_vel max_vel],[],[],[],[],[],[],[1 20 24 12],[],1,1);
    
    plot_img([],xi,zi,vi_filt_mean,haxby(32),1,0,1,12,'X (m)',...
        'Altitude (m)','Velocity (m/s)',[],[],[min_vel max_vel],[],[],[],[],[],[],[24 20 24 12],[],1,1);
end