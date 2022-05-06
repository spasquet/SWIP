function [vpout,vpoutmin,vpoutmax,vsout,vpoutstd]=velresamp(zin,vpin,zin2,vsin2,poisMIN,verbose,plot)

%%% S. Pasquet - V17.05.10
% Resample velocity model according to another one
% [vpout,vpoutmin,vpoutmax,vsout]=velresamp(zin,vpin,zin2,vsin2,poisMIN,verbose,plot)

if nargin>3
    modifpois=1;
    if exist('vsin2','var')==0
        modifpois=0;
        vsout=[];
    end
else
    modifpois=0;
    verbose=0;
    vsout=[];
    plot=0;
end
vpoutstd = [];
nlay=length(zin2)-1;

if nlay<length(zin)
    vpout=zeros(1,nlay);
    vpoutmin=vpout;
    vpoutmax=vpout;
    vpoutstd=vpout;
    for ll=1:nlay
        if isempty(vpin(zin>=zin2(ll) & zin<=zin2(ll+1)))==1 && ll>1
            vpout(ll)=vpout(ll-1);
            vpoutmin(ll)=vpoutmin(ll-1);
            vpoutmax(ll)=vpoutmax(ll-1);
            vpoutstd(ll)=vpoutstd(ll-1);
        else
            vpout(ll)=mean(vpin(zin>=zin2(ll) & zin<=zin2(ll+1)));
            vpoutmin(ll)=min(vpin(zin>=zin2(ll) & zin<=zin2(ll+1)));
            vpoutmax(ll)=max(vpin(zin>=zin2(ll) & zin<=zin2(ll+1)));
            vpoutstd(ll)=std(vpin(zin>=zin2(ll) & zin<=zin2(ll+1)));
        end
        if isnan(vpout(ll))==1
            vpout(ll)=vpout(ll-1);
            vpoutmin(ll)=vpoutmin(ll-1);
            vpoutmax(ll)=vpoutmax(ll-1);
            if verbose==1
                fprintf(['\n  NaN value for layer ',num2str(ll),'\n']);
            end
        end
        if modifpois==1
            modifvs(ll)=0;
            while poisson(vpout(ll),vsin2(ll))<=poisMIN
                vpout(ll)=vpout(ll)+1;
                vsin2(ll)=vsin2(ll)-1;
                modifvs(ll)=modifvs(ll)+1;
            end
            if modifvs(ll)>0 && verbose==1
                fprintf(['\n  Vs corrected of ',num2str(modifvs(ll)),...
                    'm/s in layer ',num2str(ll),'\n']);
            end
        end
    end
    if plot==1
        figure(2);
        stairs(zin,vpin(2:end),'g');
        view(90,90)
        hold on
        stairs(zin2,[vpout';vpout(end)],'r');
        drawnow
        hold off
    end
    
else
    ll=1; vpoutmin=[]; vpoutmax=[];
    for ll2=2:length(zin)
        while zin(ll2)>zin2(ll) && ll<length(zin2)
            vpout(ll)=vpin(ll2-1);
            if isnan(vpout(ll))==1
                vpout(ll)=vpout(ll-1);
                if verbose==1
                    fprintf(['\n  NaN value for layer ',num2str(ll),'\n']);
                end
            end
            if modifpois==1
                modifvs(ll)=0;
                while poisson(vpout(ll),vsin2(ll))<=poisMIN
                    vpout(ll)=vpout(ll)+1;
                    vsin2(ll)=vsin2(ll)-1;
                    modifvs(ll)=modifvs(ll)+1;
                end
                if modifvs(ll)>0 && verbose==1
                    fprintf(['\n  Vs corrected of ',num2str(modifvs(ll)),...
                        'm/s in layer ',num2str(ll),'\n']);
                end
            end
            ll=ll+1;
        end
    end
    if plot==1
        figure(2);
        stairs(zin,vpin,'b');
        view(90,90)
        hold on
        stairs(zin2,[vpout';vpout(end)],'r');
        drawnow
        hold off
    end
end
vpout=vpout';
if modifpois==1
    vsout=vsin2;
end
end