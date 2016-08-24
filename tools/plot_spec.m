function plot_spec(specfile,fMIN,fMAX,imgform,imgres)
%%% S. Pasquet - V16.5.3
% Quick plot of spectrogram

if exist('specfile','var')==0 || isempty(specfile)==1
    [specfile,specpath]=uigetfile('*.spec','Select spectrogram file');
    if specfile==0
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   Please select a spectrogram file');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        return
    end
else
    specpath=[];
end

run('SWIP_defaultsettings')

[specmat,fspec,xspec]=spec2dat(fullfile(specpath,specfile),0);

fig=plot_img([],fspec,xspec,-specmat,map0,0,axerev,1,fs,...
    'Frequency (Hz)','Gx (m)',[],[fMIN fMAX],[],[-10e5 0],...
    [],[],[],[],[],[1 1 24 18],[]);

sizeax=get(gca,'Position');
set(gca,'ActivePositionProperty','Position');
set(gca,'position',[sizeax(1),sizeax(2),sizeax(3),sizeax(4)/3]);
% set(gcf,'PaperPosition',[0 0 24 18],'PaperSize',[24 9.5]);

file1=[specfile,'.',imgform];
save_fig(fig,file1,imgform,imgres,1);

end