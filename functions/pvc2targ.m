function lmaxpick=pvc2targ(pvcstruct,dir_pick,nametarg,wave,sampling,resampvec,flim,maxerr)

%%% S. Pasquet - V16.6.28
% Convert .pvc ASCII file in .target dinver file

if nargin<7
    flim=0;
end
if nargin<8
    maxerr=0.5;
end

if strcmp(wave,'R')==1
    wave='Rayleigh';
elseif strcmp(wave,'L')==1
    wave='Love';
end

% if lininv==1
%     named=[num2str(XmidT(ix),['%10.',num2str(bcscale),'f']),...
%         '.contents.d'];
%     fid2=fopen(named,'w');
% end

nmode=length(pvcstruct);
lmaxpick=zeros(1,nmode);

fid0=fopen('contents.xml','w');
fprintf(fid0,'%s \n','<Dinver>');
fprintf(fid0,'%s \n','<pluginTag>DispersionCurve</pluginTag>');
fprintf(fid0,'%s \n','<pluginTitle>Surface Wave Inversion</pluginTitle>');
fprintf(fid0,'%s \n','<TargetList>');
fprintf(fid0,'%s \n',' <ModalCurveTarget type="dispersion">');
fprintf(fid0,'%s \n','<selected>true</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_Normalized</misfitType>');
for ip=1:nmode
    pvcfile=pvcstruct(ip).name;
    m=str2double(pvcfile(end-4)); % Mode number
    Vprev=load(fullfile(dir_pick,pvcfile));
    % Resample in lambda or frequency
    [freqresamp,vresamp,deltaresamp]=resampvel(Vprev(:,1),Vprev(:,2),...
        Vprev(:,3),resampvec,sampling);
    vresamp=vresamp(freqresamp>flim);
    deltaresamp=deltaresamp(freqresamp>flim);
    freqresamp=freqresamp(freqresamp>flim);
    if isempty(freqresamp)==0 && length(freqresamp)>1
        lmaxpick(ip)=max(vresamp./freqresamp);
    else
        lmaxpick(ip)=NaN;
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
        fprintf('\n   No sample in the range of resampvec - Target file empty');
        fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
        fclose(fid0);
        unix(['rm -rf ',nametarg]);
        return
    end
    nsp=length(vresamp);
    fprintf(fid0,'%s \n','<ModalCurve>');
    fprintf(fid0,'%s \n',['<name>',pvcfile,'</name>']);
    fprintf(fid0,'%s \n',['<log>',num2str(nsp),...
        ' samples loaded from file ',pvcfile,'</log>']);
    fprintf(fid0,'%s \n','<Mode>');
    fprintf(fid0,'%s \n','<slowness>Phase</slowness>');
    fprintf(fid0,'%s \n',['<polarisation>',...
        wave,'</polarisation>']);
    fprintf(fid0,'%s \n','<ringIndex>0</ringIndex>');
    fprintf(fid0,'%s %i %s \n','<index>',m,'</index>');
    fprintf(fid0,'%s \n','</Mode>');
    for ic=1:nsp
        %         stdev=1/(vresamp(ic)-deltaresamp(ic))-...
        %             1/(vresamp(ic)+deltaresamp(ic));
        
        stdev=deltaresamp(ic)/((vresamp(ic))^2-deltaresamp(ic)^2);
        if stdev>maxerr*(1/vresamp(ic))
            stdev=maxerr*(1/vresamp(ic));
        end
        
        fprintf(fid0,'%s \n','<StatPoint>');
        fprintf(fid0,'%s %4.18f %s\n','<x>',...
            freqresamp(ic),'</x>');
        fprintf(fid0,'%s %4.18f %s\n',' <mean>',...
            1/vresamp(ic),'</mean>');
        fprintf(fid0,'%s %4.18f %s\n','<stddev>',...
            stdev,'</stddev>');
        fprintf(fid0,'%s \n','<weight>1</weight>');
        fprintf(fid0,'%s \n','<valid>true</valid>');
        fprintf(fid0,'%s \n',' </StatPoint>');
        
        %         if lininv==1
        %             fprintf(fid2,'%s %d %4.6f %4.6f %4.6f\n',...
        %                 'SURF96 R C X',m,1/freqresamp(ic),...
        %                 vresamp(ic)/1000,deltaresamp(ic)/1000);
        %         end
    end
    fprintf(fid0,'%s \n','</ModalCurve>');
end
fprintf(fid0,'%s \n','</ModalCurveTarget>');
fprintf(fid0,'%s \n','<AutocorrTarget>');
fprintf(fid0,'%s \n','<selected>false</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_NormalizedBySigmaOnly</misfitType>');
fprintf(fid0,'%s \n','<AutocorrCurves>');
fprintf(fid0,'%s \n','</AutocorrCurves>');
fprintf(fid0,'%s \n','</AutocorrTarget>');
fprintf(fid0,'%s \n','<ModalCurveTarget type="ellipticity">');
fprintf(fid0,'%s \n','<selected>false</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_LogNormalized</misfitType>');
fprintf(fid0,'%s \n','</ModalCurveTarget>');
fprintf(fid0,'%s \n','<ValueTarget type="ellipticity peak">');
fprintf(fid0,'%s \n','<selected>false</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_Normalized</misfitType>');
fprintf(fid0,'%s \n','<StatValue>');
fprintf(fid0,'%s \n','<mean>0</mean>');
fprintf(fid0,'%s \n','<stddev>0</stddev>');
fprintf(fid0,'%s \n','<weight>1</weight>');
fprintf(fid0,'%s \n','<valid>false</valid>');
fprintf(fid0,'%s \n','</StatValue>');
fprintf(fid0,'%s \n','</ValueTarget>');
fprintf(fid0,'%s \n','<RefractionTarget type="Vp">');
fprintf(fid0,'%s \n','<selected>false</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_Normalized</misfitType>');
fprintf(fid0,'%s \n','</RefractionTarget>');
fprintf(fid0,'%s \n','<RefractionTarget type="Vs">');
fprintf(fid0,'%s \n','<selected>false</selected>');
fprintf(fid0,'%s \n','<misfitWeight>1</misfitWeight>');
fprintf(fid0,'%s \n','<minimumMisfit>0</minimumMisfit>');
fprintf(fid0,'%s \n','<misfitType>L2_Normalized</misfitType>');
fprintf(fid0,'%s \n','</RefractionTarget>');
fprintf(fid0,'%s \n','</TargetList>');
fprintf(fid0,'%s \n','</Dinver>');

fclose(fid0);
% if lininv==1
%     fclose(fid2);
% end
% Conversion in target
unix(['tar czf ',nametarg,' contents.xml']);
delete('contents.xml');

nametarg=strrep(nametarg,'\','\\');
fprintf(['\n  Target file saved as ',nametarg,'\n']);

end
