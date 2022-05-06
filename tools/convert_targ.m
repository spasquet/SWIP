function convert_targ(resampvec)

%%% S. Pasquet - V16.11.22
% Convert target files in pvc files
% convert_targ(resampvec)

if exist('resampvec','var')==0 || isempty(resampvec)==1
   resampvec = 0:0.5:200; 
end

dir_targ=uigetdir('./','Select folder containing target files to convert');
targstruct=dir(fullfile(dir_targ,'*.target'));
dir_pvc=fullfile(dir_targ,'file.pvc');
if exist(dir_pvc,'dir')~=7
    mkdir(dir_pvc);
end

for i=1:length(targstruct)
    nametarg=targstruct(i).name;
    Xmid=str2double(nametarg(1:end-7));
%     pvcfile=[num2str(Xmid),targfile(end-6:end)];
    [freqresamp,vresamp,deltaresamp,modes]=targ2pvc(fullfile(dir_targ,nametarg));
    % Read target file to get picked dispersion curves
    npvc=length(modes);
    for ip=1:npvc
        % Resample in lambda or frequency
        if length(freqresamp{modes(ip)+1})>1
            plot_curv(modes(ip),freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                deltaresamp{modes(ip)+1},'.-',[0 0 1],[],0,0,...
                0,16,'Frequency (Hz)','Phase velocity (m/s)',...
                [],[],[],[],[],[],[],[],[],[0 0 24 18],[]);
            hold on

            [freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1}]=...
                resampvel(freqresamp{modes(ip)+1},vresamp{modes(ip)+1},deltaresamp{modes(ip)+1},resampvec,0,0);
            
            namepvc=fullfile(dir_pvc,[num2str(Xmid),'.M',num2str(modes(ip)),'.txt']);
            dlmwrite(namepvc,[freqresamp{modes(ip)+1}',vresamp{modes(ip)+1}',deltaresamp{modes(ip)+1}'],'\t');

            plot_curv(modes(ip),freqresamp{modes(ip)+1},vresamp{modes(ip)+1},...
                deltaresamp{modes(ip)+1},'.-',[1 0 0],[],0,0,...
                0,16,'Frequency (Hz)','Phase velocity (m/s)',...
                [],[],[],[],[],[],[],[],[],[0 0 24 18],[]);
        end
    end
end
