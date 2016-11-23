function [vmaxamp,fmaxamp]=findpeak(dspmat,f,v,fpick,vpick,wl)

% S. Pasquet - V16.11.18
% findpeak.m looks at dispersion image maximum in specified window at the picked frequency
% [vmaxamp,fmaxamp]=findpeak(dspmat,f,v,fpick,vpick,wl)

nbpick=length(fpick);
fmaxamp=[]; vmaxamp=[];
for i=1:nbpick
    Fvec=repmat(fpick(i),size(f));
    ind=find(abs(Fvec-f)==min(abs(Fvec-f)),1,'first');
    fmaxamp(i)=f(ind);
    if isnan(vpick(i))~=1
        if wl(i)>1
            condV=(v>=vpick(i)-wl(i) & v<=vpick(i)+wl(i));
            ind2=find(dspmat(ind,condV)==max(dspmat(ind,condV)));
            nbind2=length(ind2);
            if nbind2>1
                if mod(nbind2,2)==1
                    ind2=ind2(1+(nbind2-1)/2);
                else
                    ind2=ind2(nbind2/2);
                end
            end
            vwin=v(condV);
            vmaxamp(i)=vwin(ind2);
        else
            vmaxamp(i)=vpick(i);
        end
    else
        vmaxamp(i)=NaN;
    end
end

end