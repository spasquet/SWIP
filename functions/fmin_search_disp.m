function [flim,dspmat_all]=fmin_search_disp(dispstruct,dir_dat_xmid,f,specfmin)

%%% S. Pasquet - V17.05.23
%%% Look for minimum frequency where dispersion image reaches ...
% [flim,specfileOK]=fmin_search(specstruct,dir_dat_xmid,dt,specampmin,specfmin)

flim=NaN*zeros(length(dispstruct),1);
dspmat_all = [];

for ip=1:length(dispstruct)
    dispfile=dispstruct(ip).name;
    dspmat=dsp2dat(fullfile(dir_dat_xmid,[dispfile(1:end-3),'dsp']),0,0);  
    dspmat_all = [dspmat_all; dspmat];
end

dspamp_mean = nanmean(dspmat_all,2);
dspamp_std = nanstd(dspmat_all,2);
try
    ind = find(dspamp_mean==max(dspamp_mean) & f>=specfmin,1,'first');
catch
    ind = find(dspamp_mean==max(dspamp_mean) & f'>=specfmin,1,'first');
end

flim = f(ind);