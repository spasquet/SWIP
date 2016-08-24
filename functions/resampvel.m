function [freqresamp,vresamp,deltaresamp]=resampvel(freq,vel,deltac,resampvec,sampling,withnan)

%%% S. Pasquet - V16.5.17
% Resample phase velocity in lambda or frequency

if nargin<6
    withnan=0;
end
% Resample in lambda or frequency
if sampling==1
    lambda=round((vel./freq)*1e10)/1e10;
    %     dlam=mean(unique(diff(resampvec)));
    %     dl=unique(diff(lambda));
    %     if mean(abs(diff(dl)))<1e-10 || isnan(mean(abs(diff(dl))))==1
    %         lambda=round(lambda/dlam)*dlam;
    %     end
    if length(unique(lambda))<length(lambda)
        [lambda,I]=unique(lambda,'first');
        vel=vel(I);        
    end
    vresamp=interp1(lambda,vel,resampvec,'linear');
    freqresamp=vresamp./resampvec;
    deltaresamp=interp1(lambda,deltac,resampvec,'linear');
    if exist('flim','var')==1
        vresamp(freqresamp<flim)=NaN;
    end
    if withnan==0
        freqresamp=fliplr(freqresamp(isnan(vresamp)==0));
        deltaresamp=fliplr(deltaresamp(isnan(vresamp)==0));
        vresamp=fliplr(vresamp(isnan(vresamp)==0));
        [freqresamp,idx]=sort(freqresamp);
        deltaresamp=deltaresamp(idx);
        vresamp=vresamp(idx);
    end
else
    if exist('flim','var')==1
        resampvec=resampvec(resampvec>flim);
    end
    vresamp=interp1(freq,vel,resampvec,'linear');
    freqresamp=resampvec;
    deltaresamp=interp1(freq,deltac,resampvec,'linear');
    if withnan==0
        freqresamp=freqresamp(isnan(vresamp)==0);
        deltaresamp=deltaresamp(isnan(vresamp)==0);
        vresamp=vresamp(isnan(vresamp)==0);
        [freqresamp,idx]=sort(freqresamp);
        deltaresamp=deltaresamp(idx);
        vresamp=vresamp(idx);
    end

end