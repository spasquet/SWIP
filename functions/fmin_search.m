function [flim,specfileOK,specmat_all]=fmin_search(av,specstruct,dir_dat_xmid,dt,specampmin,specfmin)

%%% S. Pasquet - V17.05.23
%%% Look for minimum frequency where spectrogram normalized amplitude reaches specampmin
% [flim,specfileOK]=fmin_search(specstruct,dir_dat_xmid,dt,specampmin,specfmin)
flim=NaN*zeros(length(specstruct),1);
specmat_all = [];

for ip=1:length(specstruct)
    specfile=specstruct(ip).name;
    specmat=dsp2dat(fullfile(dir_dat_xmid,[specfile(1:end-4),'spec']),0,0);
    fspecmax=1/dt;
    nfspec=size(specmat,2);
    dfspec=fspecmax/2/(nfspec-1);
    fspec=0:dfspec:fspecmax/2;
    nx=size(specmat,1);
    flimtmp=ones(nx,1)*max(fspec)*NaN;
    
    specmat_all = [specmat_all; specmat];
    
    for jj=1:nx
        ind=find(specmat(jj,:)>=specampmin & fspec>specfmin,1,'first');
        if isempty(ind)==0
            flimtmp(jj)=fspec(ind);
        else
            flimtmp(jj)=0;
        end
    end
    flim(ip)=mean(flimtmp(isnan(flimtmp)==0))-std(flimtmp(isnan(flimtmp)==0));
end

specamp_mean = mean(specmat_all);
specamp_std = std(specmat_all);

if av == 1
    if isempty(specstruct)==0
        specfileOK=specstruct(flim==min(flim)).name;
        flim=mean(flim);
    else
        flim=NaN;
        specfileOK=[];
    end
elseif av == 2 || av == 0
%     ind = find(specamp_mean>=specampmin & fspec>specfmin,1,'first');
    ind = find(specamp_mean>=specamp_mean(2)+specamp_std(2) & fspec>specfmin,1,'first');
    
    if isempty(specstruct)==0
        specfileOK=specstruct(flim==min(flim)).name;
        flim = fspec(ind-1);
    else
        flim=NaN;
        specfileOK=[];
    end
end

%%
% plot_curv(1,fspec,specmat_tmp,[],'.-','r',2,0,0,0,16,'Frequency (Hz)','Mean norm. amplitude',[],...
%     [0 100],[],[],[],[],[],specampmin,flim,[],[],[],[],[]);
