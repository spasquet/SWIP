function [vmaxamp,fmaxamp]=findpeak(dspmat,f,v,fpick,vpick,wl)

% S. Pasquet - V18.05.23
% findpeak.m looks at dispersion image maximum in specified window at the picked frequency
% [vmaxamp,fmaxamp]=findpeak(dspmat,f,v,fpick,vpick,wl)

nbpick=length(fpick);
fmaxamp=[]; vmaxamp=[];
for i=1:nbpick
    if fpick(i)>max(f) || fpick(i)<min(f)
        continue
    end
    Fvec=repmat(fpick(i),size(f));
    ind=find(abs(Fvec-f)==min(abs(Fvec-f)),1,'first');
    fmaxamp(i)=f(ind);
    if vpick(i)>max(v)
        vpick(i)=max(v);
    end
    if isnan(vpick(i))~=1
        if wl(i)>1
            vmin=vpick(i)-wl(i);
            vmax=vpick(i)+wl(i);
            condV=(v>=vmin & v<=vmax);
            if isempty(find(condV==1, 1))
                ind_v=find(abs(vpick(i)-v)==min(abs(vpick(i)-v)),1,'first');
                condV(ind_v)=1;
            end
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
            if ~isempty(ind2)
                vmaxamp(i)=vwin(ind2);
            else
                vmaxamp(i)=NaN;
            end
        else
            vmaxamp(i)=vpick(i);
        end
    else
        vmaxamp(i)=NaN;
    end
end

end