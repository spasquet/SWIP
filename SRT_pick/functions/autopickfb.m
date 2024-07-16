function [seismomat_test, tpick_auto, xpick_auto] = autopickfb(seismomat,xseis,tseis,sx,autopick,plot_fig)

%%
% close all
amp_lim = 10;
offmin = 5;

% noise_fac1 = 10;
% noise_fac2 = 250;
% pow1 = 0.2;
% pow2 = 1;

noise_fac1 = 5;
noise_fac2 = 100;

pow1 = 0.2;
pow2 = 1;

xpick_auto = xseis;

offset = abs(xpick_auto-sx);
noise_fac_law = logspace(log10(noise_fac1),log10(noise_fac2),length(offset))';
pow_law = logspace(log10(pow1),log10(pow2),length(offset))';

if ~isempty(find(tseis<0, 1))
    noise_std = std(abs(seismomat(:,tseis<0)),[],2);
    noise_std = median_filt(noise_std',7,1,length(noise_std));
else
    seismomat_tmp = seismomat;
%     noise_std = 0.002*offset.*std(abs(seismomat_tmp),[],2);
    noise_std = 0.05*mean(abs(seismomat_tmp),2);
    noise_std = median_filt(noise_std',7,1,length(noise_std));
end

% signal_std = std(abs(seismomat(:,tseis>0)),[],2);
% signal_std = median_filt(signal_std',7,1,length(signal_std));

seismomat_temp = seismomat;
seismomat_temp(:,tseis<=2*min(diff(tseis))) = 0;

noise_law = noise_fac_law.*offset.^-pow_law;
% noise_law(offset>0.5*max(offset)) = noise_fac2*offset(offset>0.5*max(offset)).^-pow2;
noise_law(isinf(noise_law)) = max(noise_fac_law);

cond = bsxfun(@gt,abs(seismomat_temp),noise_law.*noise_std);
% imagesc(cond')
%%

if autopick == 1
    source_pos = find(offset == min(offset),1);
    noise_threshold = min(log10(noise_std))+0.05*abs(median(log10(noise_std)));
    noise_threshold_ind_left = find(log10(noise_std)>=noise_threshold & xpick_auto<sx,1,'last');
    noise_threshold_ind_right = find(log10(noise_std)>=noise_threshold & xpick_auto>sx,1,'first');
    dist_left = abs(source_pos - noise_threshold_ind_left);
    dist_right = abs(source_pos - noise_threshold_ind_right);
    offset_min = max([dist_left dist_right]);
    if offset_min < offmin/unique(diff(xpick_auto))
        offset_min = offmin/unique(diff(xpick_auto));
    end
%     plot_fig = 1;

    cond_filt = despike_2D(cond,7,50,[1 length(xpick_auto)],[1 length(tseis)],source_pos,offset_min,plot_fig);
%     cond_filt = despike_2D(cond_filt,7,150,[1 length(xseis)],[1 length(tseis)],source_pos,offset_min,1);
else
    cond_filt = cond;
end

[ok,idx] = max(cond_filt,[],2);
% idx(idx <=4 & xseis~=sx) = 1;
norm_value = ok.*abs(seismomat_temp((1:numel(idx))'+numel(idx)*(idx-1)));
id_ok = find(norm_value~=0);
norm_value = interp1(xseis(id_ok),norm_value(id_ok),xseis);
seismomat_test = bsxfun(@rdivide,seismomat,norm_value);

% plot_wiggle(2,-1*seismomat_test',xseis,tseis)
% keyboard

seismomat_test(seismomat_test>amp_lim)=amp_lim;
seismomat_test(seismomat_test<-amp_lim)=-amp_lim;
seismomat_test = seismomat_test./amp_lim;
seismomat_test(isnan(mean(seismomat_test,2)),:) = 0;

if autopick == 1
    TSEIS = repmat(tseis,length(xpick_auto),1);
    TSEIS2 = TSEIS(:,idx);
    tpick_auto = TSEIS2(1,:);
    
%     if ismember(sx,xseis)
%         if sx~=min(xseis)
%             tpick_auto(xseis<sx) = median_filt(tpick_auto(xseis<sx),3,1,length(tpick_auto(xseis<sx)))';
%         end
%         if sx~=max(xseis)
%             tpick_auto(xseis>sx) = median_filt(tpick_auto(xseis>sx),3,1,length(tpick_auto(xseis>sx)))';
%         end
%     end

    xpick_auto(tpick_auto==0 & xseis'~=sx) = [];
    xpick_auto(tpick_auto>max(tseis)) = [];
    tpick_auto(tpick_auto==0 & xseis'~=sx) = [];
    tpick_auto(tpick_auto>max(tseis)) = [];
    
else
    xpick_auto = [];
    tpick_auto = [];
end