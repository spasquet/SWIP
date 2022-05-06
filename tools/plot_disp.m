function plot_disp(dspfile,targetfile,Dlogscale,fMIN,fMAX,imgform,imgres,flip,eb)

%%% S. Pasquet - V16.11.22
% Quick plot of dispersion image and curves
% plot_disp(dspfile,targetfile,Dlogscale,fMIN,fMAX,imgform,imgres,flip,eb)

if exist('dspfile','var')==0 || isempty(dspfile)==1
    [dspfile,dsppath]=uigetfile('*.dsp','Select dispersion file');
    if dspfile==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a dispersion file');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
else
    dsppath=[];
end

if exist('targetfile','var')==0 || isempty(targetfile)==1
        [targetfile,targetpath]=uigetfile('*.target','Select dispersion curve or cancel');
else
    targetpath=[];
end

run('SWIP_defaultsettings')

[dspmat,f,v]=dsp2dat(fullfile(dsppath,dspfile),flip,0);

if Dlogscale==0
    fig=plot_img([],f,v,-dspmat',map0,0,axerev,0,fs,...
        'Frequency (Hz)','Phase velocity (m/s)',[],[fMIN fMAX],[],[],...
        [],[],[],[],[],[],[1 1 24 18],[]);
else
    dspmat=1./(1-dspmat);
    dspmat(isinf(dspmat))=max(max(dspmat(isinf(dspmat)==0)));
    fig=plot_img_log([],f,v,dspmat',flipud(map0),0,axerev,0,fs,...
        'Frequency (Hz)','Phase velocity (m/s)',[],[fMIN fMAX],[],[1 length(map0)],...
        [],[],[],[],[],[],[1 1 24 18],[]);
end

if targetfile~=0
    % Read target file to get picked dispersion curves
    [freqresamp,vresamp,deltaresamp,modes]=targ2pvc(fullfile(targetpath,targetfile));
    npvc=length(modes);
    for ip=1:npvc
        % Resample in lambda or frequency
        if length(freqresamp{modes(ip)+1})>1
            [freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1}]=...
                resampvel(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                deltaresamp{modes(ip)+1},resampvec,sampling,1);
        end
    end
    
    for ip=1:npvc
        hold on;
        col='w';
        if eb==1
            errorbar(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},'.','Color',col,...
                'linewidth',1.5,'markersize',8);
        elseif eb==0
            plot(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},'.','Color',col,...
                'linewidth',1.5,'markersize',10);
        end
    end
end
hold off

file1=[dspfile(1:end-3),'disp.',imgform];
save_fig(fig,file1,imgform,imgres,1);

end