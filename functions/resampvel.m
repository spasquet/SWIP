function [freqresamp,vresamp,deltaresamp]=resampvel(freq,vel,deltac,resampvec,sampling,withnan)

%%% S. Pasquet - V16.11.18
% Resample phase velocity in lambda or frequency
% [freqresamp,vresamp,deltaresamp]=resampvel(freq,vel,deltac,resampvec,sampling,withnan)

if nargin<6
    withnan=0;
end
% Resample in lambda or frequency
if sampling==1
    lambda=round((vel./freq)*1e10)/1e10;
    resampvec = round(resampvec*1e10)/1e10;
%     lambda = vel./freq;
    if length(unique(lambda))<length(lambda)
        [lambda,I]=unique(lambda,'first');
        vel=vel(I);    
        deltac=deltac(I);
    end
    vresamp=interp1(lambda,vel,resampvec,'linear');
    freqresamp=vresamp./resampvec;
    deltaresamp=interp1(lambda,deltac,resampvec,'linear');
    if exist('flim','var')==1
        vresamp(freqresamp<flim)=NaN;
    end
    
%     if max(vresamp./freqresamp) > 123 && max(vresamp./freqresamp) < 130
%        keyboard 
%     end
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