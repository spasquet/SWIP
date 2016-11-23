function [flim,specfileOK]=fmin_search(specstruct,dir_dat_xmid,dt,specampmin,specfmin)

%%% S. Pasquet - V16.11.18
%%% Look for minimum frequency where spectrogram normalized amplitude reaches specampmin
% [flim,specfileOK]=fmin_search(specstruct,dir_dat_xmid,dt,specampmin,specfmin)

flim=NaN*zeros(length(specstruct),1); 
for ip=1:length(specstruct)
    specfile=specstruct(ip).name;
    specmat=dsp2dat(fullfile(dir_dat_xmid,[specfile(1:end-4),'spec']),0,0);
    fspecmax=1/dt;
    nfspec=size(specmat,2);
    dfspec=fspecmax/2/(nfspec-1);
    fspec=0:dfspec:fspecmax/2;
    nx=size(specmat,1);
    flimtmp=ones(nx,1)*max(fspec)*NaN;
    
    for jj=1:nx
        ind=find(specmat(jj,:)>=specampmin & fspec>specfmin,1,'first');
        if isempty(ind)==0
            flimtmp(jj)=fspec(ind);
        else
            flimtmp(jj)=0;
        end
    end
    flim(ip)=mean(flimtmp(isnan(flimtmp)==0));
end
if isempty(specstruct)==0
    specfileOK=specstruct(flim==min(flim)).name;
    flim=mean(flim);
else
    flim=NaN;
    specfileOK=[];
end
end
