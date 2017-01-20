function [fsave,Vsave,deltacsave,modenext,closefig,xmidprev]=matpickamp(dspmat,f,v,filepick,pickstyle,modeinit,...
    err,smoothpick,nW,dx,fac,maxerr,minvelerr,sigma)

%%% S. Pasquet - V17.01.16
% Pick amplitudes of dispersion image
% [fsave,Vsave,deltacsave,modenext,closefig,xmidprev]=matpickamp(dspmat,f,v,filepick,...
%    pickstyle,modeinit,err,smoothpick,nW,dx,fac,maxerr,minvelerr,sigma)

fprintf(['\n  Picking mode ',num2str(modeinit),'\n']);
hold on;
h=get(gca,'children');
button=1;
i=0; closefig=0;
modenext=[]; xmidprev=0;
fsave=[];Vsave=[];deltacsave=[];

if exist(filepick,'file')==2
    Vprev=load(filepick);
    Vi=interp1qr(Vprev(:,1),Vprev(:,2),f')';
    deltac=interp1qr(Vprev(:,1),Vprev(:,3),f')';
%     vmaxamp=Vi(isnan(Vi)==0);
%     fmaxamp=f(isnan(Vi)==0);
    fmaxamp=Vprev(:,1);
    vmaxamp=Vprev(:,2);
    deltacamp=Vprev(:,3);
    i=length(Vi);
    hold on;
    h0=plot(fmaxamp,vmaxamp,'w.');
else
    Vi=[]; deltac=[];  fmaxamp=[]; vmaxamp=[]; deltacamp=[];
end

fprintf('\n  ENTER or close window : Save current picks and go to next Xmid');
fprintf('\n  BACKSPACE : Save current picks and go to previous Xmid');
fprintf('\n  H : Save current picks and go to next (higher) mode');
fprintf('\n  L : Save current picks and go to previous (lower) mode');
fprintf('\n  W : Save current picks and stop script');
fprintf('\n  A : Discard current picks and go to next Xmid (keep previous picks)');
fprintf('\n  Z : Discard current picks and go to previous Xmid (keep previous picks)');
fprintf('\n  N : Discard current picks and go to next (higher) mode (keep previous picks)');
fprintf('\n  P : Discard current picks and go to previous (lower) mode (keep previous picks)');
fprintf('\n  X : Discard current picks and stop script (keep previous picks)');
fprintf('\n  M : Switch between manual and semi-automatic picking');
fprintf('\n  S : Switch between smooth and regular picking');
fprintf('\n  D : Delete one or several points');
fprintf('\n  R : Reset all picks');
fprintf('\n  E : Change error style (no error, percentage or lorentz)');
fprintf('\n  C : Open colormap editor\n');

while button==1
    i=i+1;
    try
    [Fpick,Vpick,button]=ginput(1);
    catch
        closefig=1;
        if isempty(Vi)==1
            Vsave=0;
            fprintf('\n  No picks - Go to next Xmid\n');
        end
        if pickstyle==1
            fsave=f; Vsave=Vi; deltacsave=deltac;
        else
            fsave=fmaxamp'; Vsave=vmaxamp'; deltacsave=deltacamp';
        end
        fprintf(['\n  Mode ',num2str(modeinit),' saved - Go to next Xmid\n']);
        return
    end
    if isempty(button)==1 % Press ENTER to end picking, save pick file and go to next Xmid
        if isempty(Vi)==1
            Vsave=0;
            fprintf('\n  No picks - Go to next Xmid\n');
        end
        if pickstyle==1
            fsave=f; Vsave=Vi; deltacsave=deltac;
        else
            fsave=fmaxamp'; Vsave=vmaxamp'; deltacsave=deltacamp';
        end
        fprintf(['\n  Mode ',num2str(modeinit),' saved - Go to next Xmid\n']);
        break
    elseif button==8 % Press BACKSPACE to end picking, save pick file and go to previous Xmid
        if isempty(Vi)==1
            Vsave=0;
            fprintf('\n  No picks - Go to previous Xmid\n');
        end
        if pickstyle==1
            fsave=f; Vsave=Vi; deltacsave=deltac;
        else
            fsave=fmaxamp'; Vsave=vmaxamp'; deltacsave=deltacamp';
        end
        fprintf(['\n  Mode ',num2str(modeinit),' saved - Go to previous Xmid\n']);
        xmidprev=1;
        break
    elseif button==119 % Press W to end picking, save pick file and stop script
        if isempty(Vi)==1
            Vsave=0;
        end
        if pickstyle==1
            fsave=f; Vsave=Vi; deltacsave=deltac;
        else
            fsave=fmaxamp'; Vsave=vmaxamp'; deltacsave=deltacamp';
        end
        fprintf(['\n  Mode ',num2str(modeinit),' saved - Stopping script\n']);
        xmidprev=-1;
        break
    elseif button==104 % Press H to end picking, save pick file and go to higher mode
        if isempty(Vi)==1
            Vsave=0;
            fprintf('\n  No picks - Go to higher mode\n');
        end
        if pickstyle==1
            fsave=f; Vsave=Vi; deltacsave=deltac;
        else
            fsave=fmaxamp'; Vsave=vmaxamp'; deltacsave=deltacamp';
        end
        modenext=modeinit+1;
        fprintf(['\n  Mode ',num2str(modeinit),' saved - Go to higher mode\n']);
        break
    elseif button==108 % Press L to end picking, save pick file and go to lower mode
        if isempty(Vi)==1
            Vsave=0;
            fprintf('\n  No picks - Go to lower mode\n');
        end
        if pickstyle==1
            fsave=f; Vsave=Vi; deltacsave=deltac;
        else
            fsave=fmaxamp'; Vsave=vmaxamp'; deltacsave=deltacamp';
        end
        if modeinit>0
            modenext=modeinit-1;
        else
            modenext=0;
        end
        fprintf(['\n  Mode ',num2str(modeinit),' saved - Go to lower mode\n']);
        break
    elseif button==97 % Press A to abort, keep existing pick file and go to next Xmid
        fprintf('\n  Abort and keep existing picks - Go to next Xmid\n');
        Vsave=0;
        break
    elseif button==120 % Press X to abort, keep existing pick file and stop script
        fprintf('\n  Abort and keep existing picks - Stopping script\n');
        Vsave=0; xmidprev=-1;
        break
    elseif button==122 % Press Z to abort, keep existing pick file and go to previous Xmid
        fprintf('\n  Abort and keep existing picks - Go to previous Xmid\n');
        Vsave=0; xmidprev=1;
        break
    elseif button==110 % Press N to abort and go to higher (next) mode
        modenext=modeinit+1;
        fprintf('\n  Abort and keep existing picks - Go to next mode\n');
        Vsave=0;
        break
    elseif button==112 % Press P to abort and go to lower (previous) mode
        if modeinit>0
            modenext=modeinit-1;
        else
            modenext=0;
        end
        fprintf('\n  Abort and keep existing picks - Go to previous mode\n');
        Vsave=0;
        break
    elseif button==99 % Press C to display colorbar editor
        fprintf('\n  Change colormap\n');
        fprintf('\n  Go in Tools/Standard Colormaps or edit current one then press OK\n');
        colormapeditor;
        i=i-1; button=1;
        continue
    elseif button==109 % Press M to select manual or semi-auto picking 
        if pickstyle==1
            pickstyle=0;
            style='manual';
            if exist('h0','var')==1
                delete(h0); hold on;
            end
            h0=plot(fmaxamp,vmaxamp,'w.');
        else
            pickstyle=1;
            style='semi-automatic';
            if exist('h0','var')==1
                delete(h0); hold on;
            end
            if length(vmaxamp)>1
                h0=plot(f,Vi,'w.');
            else
                h0=plot(fmaxamp,vmaxamp,'w.');
            end
        end
        fprintf(['\n  Switch picking style -> ',style,' picking\n']);
        i=i-1; button=1;
        continue
    elseif button==115 % Press S to select smooth or regular picking
        if smoothpick==1
            smoothpick=0;
            style2='regular';
        else
            smoothpick=1;
            style2='smooth';
        end
        fprintf(['\n  Switch picking style -> ',style2,' picking\n']);
        i=i-1; button=1;
        continue
    elseif button==101 % Press E to change error style
        if err==1
            err=2;
            errstyle='percentage';
        elseif err==2
            err=0;
            errstyle='no';
        elseif err==0
            err=1;
            errstyle='Lorentz';
        end
        fprintf(['\n  Change error style -> ',errstyle,' error\n']);
        i=i-1; button=1;
        continue
    elseif button==100 % Press D to delete points (clic outside of data to end)
        if i==1
            i=i-1; button=1;
            continue
        end
        indselec=1;
        fprintf('\n  Delete points');
        fprintf('\n  Keep the mouse button pushed while moving the pointer to delete points');
        fprintf('\n  Press ENTER to resume picking or close the window to go to the next Xmid\n');
        while isempty(indselec)==0
            try
                indselec=selectdata('action','delete',...
                'selectionmode','brush','ignore',h);
            catch
                fprintf(['\n  Mode ',num2str(modeinit),' saved - Go to next Xmid\n']);
                closefig=1;
                return
            end
            if isempty(indselec)==0
                if pickstyle==1
                    Vi(indselec)=NaN;
                else
                    vmaxamp(indselec)=NaN;
                    fmaxamp(indselec)=NaN;
                    deltacamp(indselec)=NaN;
                    vmaxamp=vmaxamp(isnan(vmaxamp)~=1);
                    fmaxamp=fmaxamp(isnan(fmaxamp)~=1);
                    deltacamp=deltacamp(isnan(deltacamp)~=1);
                end
            end
        end
        button=1; fprintf('\n  Back to picking\n');
%         if isempty(f(isnan(Vi)~=1))==0 && length(f(isnan(Vi)~=1))>1
%             vmaxamp=interp1(f(isnan(Vi)~=1),Vi(isnan(Vi)~=1),fmaxamp,'linear');
%         else
%             vmaxamp=fmaxamp*NaN;
%         end
%         vmaxamp=vmaxamp(isnan(vmaxamp)~=1);
%         fmaxamp=fmaxamp(isnan(vmaxamp)~=1);
        continue
    elseif button==114 % Press R to reset picks
        fprintf('\n  Reset all picks\n');
        Vi=[]; deltac=[]; fmaxamp=[]; vmaxamp=[];
        button=1; i=1;
        if exist('h0','var')==1
            delete(h0); clear h0;
        end
        continue
    elseif button~=1;
        button=1; i=i-1;
        continue
    end
    if pickstyle==1
        if err==1
            wl=lorentzerr(Vpick,Vpick/Fpick,nW,dx,fac,maxerr,minvelerr);
        elseif err==2
            wl=0.01*sigma*Vpick;
        end
    else
        wl=1;
    end
    if i>1
        [vmaxamp(end+1),fmaxamp(end+1)]=findpeak(dspmat,f,v,Fpick,Vpick,wl);
        if ismember(fmaxamp(end),fmaxamp(1:end-1))==1
            vmaxamp(fmaxamp==fmaxamp(end))=vmaxamp(end);
            vmaxamp=vmaxamp(1:end-1);
            fmaxamp=fmaxamp(1:end-1);
        end
        if length(vmaxamp)>1
            Vi=interp1(fmaxamp,vmaxamp,f,'linear');
            if pickstyle==1
                if err==1
                    wl=lorentzerr(Vi,Vi./f,nW,dx,fac,maxerr,minvelerr);
                    wl(wl==0)=NaN;
                elseif err==2
                    wl=0.01*sigma*Vi;
                    wl(wl==0)=NaN;
                end
            else
                wl=ones(size(Vi));
            end
            if err==1
                deltac=lorentzerr(Vi,Vi./f,nW,dx,fac,maxerr,minvelerr);
                deltac(deltac==0)=NaN;
                if size(fmaxamp,1)==1 && size(fmaxamp,2)==2
                    fmaxamp=fmaxamp'; vmaxamp=vmaxamp';
                end
                deltacamp=lorentzerr(vmaxamp',vmaxamp'./fmaxamp',nW,dx,fac,maxerr,minvelerr);
                deltacamp(deltac==0)=NaN;
            elseif err==2
                deltac=0.01*sigma*Vi;
                deltac(deltac==0)=NaN;
                deltacamp=0.01*sigma*vmaxamp;
                deltacamp(deltac==0)=NaN;
            else
                deltac=zeros(size(Vi));
                deltacamp=zeros(size(vmaxamp));
            end
            if size(deltac,1)~=size(Vi,1)
                deltac=deltac';
            end
            if size(deltacamp,1)~=size(vmaxamp,1)
                deltacamp=deltacamp';
            end
            Vi=findpeak(dspmat,f,v,f,Vi,wl);
            if smoothpick==1
                Vi=median_filt(Vi,5,1,length(Vi));
                Vi=mov_aver(Vi',3,1,length(Vi));
                Vi=Vi';
            end
        else
            [vmaxamp,fmaxamp]=findpeak(dspmat,f,v,Fpick,Vpick,wl);
        end
    else
        [vmaxamp,fmaxamp]=findpeak(dspmat,f,v,Fpick,Vpick,wl);
    end
    [fmaxamp,I]=sort(fmaxamp);
    vmaxamp=vmaxamp(I);
    hold on
    if exist('h0','var')==1
        delete(h0);
    end
    if length(vmaxamp)>1 && pickstyle==1
        h0=plot(f,Vi,'w.');
    else
        h0=plot(fmaxamp,vmaxamp,'w.');
    end
end
end
