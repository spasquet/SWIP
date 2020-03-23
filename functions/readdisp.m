function [D,dispselec]=readdisp(file,nmaxmod,runnb,nbest,marg,freq,vphmin,vphmax)

%%% S. Pasquet - V16.11.18
% Read dispersion computed by gpdc (single forward model) or by gpdcreport
% (dinver inversion results). Look for nbest models, or for the models
% included within the error bars
% [D,dispselec]=readdisp(file,nmaxmod,runnb,nbest,marg,freq,vphmin,vphmax)

version = geopsy_version();

if exist('runnb','var')==0 && exist('nbest','var')==0
    runnb=1;
    fwd=1;
else
    fwd=0;
    if exist('freq','var')==1 && isempty(freq)==0
        freq=round(10000*freq)/10000;
    end
end

fileID=fopen(file);
% Read velocity model
formatSpec = '%f %f';
ligne=fgets(fileID); % Read line
if ligne==-1
%     fprintf('\n  Empty dispersion file - Re-run forward calculation\n');
    D=[]; dispselec=[];
else
%     if version > 2
%         keyboard
%         ind_slash = strfind(ligne,'/');
%         nmaxmod = str2double(ligne(3:ind_slash-1));
%     end
    
    D=cell(nmaxmod,2);
    modnum=zeros(nmaxmod,1)*NaN;
    modok=modnum;
    misround=modnum;
    nfreqsample=modnum;
    
    if fwd==1
        fgets(fileID); % Read line
    end
    nn=0;
    for n=1+(runnb-1)*nmaxmod:(runnb-1)*nmaxmod+nmaxmod
        nn=nn+1;
        ligne=fgets(fileID); % Read line
        if length(ligne)>1
            if fwd~=1
                if version == 2
                    modnum(nn)=sscanf(ligne,'%*s %*s %*s %*s %*s %*s %*s %*s %d')+(runnb-1)*nmaxmod;
                    misround(nn)=sscanf(ligne,'%*s %*s %*s %*s %*s %*s %*s %*s %*s %*6s %f');
                else
                    ind_parenthese = strfind(ligne,')');
                    ind_2points = strfind(ligne,':');
                    modnum(nn) = str2double(ligne(ind_parenthese+1:ind_2points-1))+(runnb-1)*nmaxmod;
                    
                    ind_value = strfind(ligne,'value=');
                    misround(nn) = str2double(ligne(ind_value+6:end));
                end
            end
            D(nn,1:2)=textscan(fileID,formatSpec);
            freqcal=round(10000*cat(2,D{nn,1}))/10000;
            veltmp=1./cat(2,D{nn,2}); % store inversion results in temp variable
            nfreqsample(nn)=length(freqcal);
            if isempty(freqcal)==1
                D(nn,1:2)={NaN};
            end
            % Skip models if calculated dispersion curves have less points
            % than picked dispersion curves
            if fwd~=1 && (isempty(freqcal)==0 && (min(freqcal)>min(freq) || max(freqcal)<max(freq)))
%                 plot(freq,vphmin,'k-');
%                 hold on;
%                 plot(freq,vphmax,'k-');
%                 plot(freqcalnew,veltmpnew,'rx');
%                 drawnow;
                veltmp=[];
            end
            if fwd~=1 && isempty(veltmp)==0
                if nbest==0 % Select models within the error bars
                    if nn==1
                        freqcalprev=freqcal*0;
                    end
                    if isequal(freqcal,freqcalprev)==0
                        vphminnew=interp1qr(freq',vphmin',freqcal);
                        vphmaxnew=interp1qr(freq',vphmax',freqcal);
                    end
                    veltmpnew=veltmp(isnan(vphminnew)==0);
                    freqcalnew=freqcal(isnan(vphminnew)==0);
                    vphmaxnew2=vphmaxnew(isnan(vphminnew)==0);
                    vphminnew2=vphminnew(isnan(vphminnew)==0);
                    ref=zeros(size(freqcalnew));
                    for kk=1:length(freqcalnew)
                        ref(kk)=veltmpnew(kk)>=vphminnew2(kk) && veltmpnew(kk)<=vphmaxnew2(kk);
                    end
%                     plot(freqcalnew,veltmpnew,'rx');
%                     hold on;
%                     plot(freqcalnew,vphminnew2,'k-');
%                     plot(freqcalnew,vphmaxnew2,'k-');
%                     hold off;
%                     drawnow
                    if sum(ref(ref==1))>=length(ref)-marg
                        modok(nn)=1;
                    else
                        modok(nn)=0;
                    end
                else % Select nbest models
                    if nn<=nbest
                        modok(nn)=1;
                    else
                        modok(nn)=0;
                    end
                end
            else
                modok(nn)=0;
            end
        else
            if fwd==0
                dispselec=struct('modnum',modnum,'modok',modok,...
                    'misround',misround,'nfreqsample',nfreqsample);
            end
            fclose(fileID);
            return
        end
        freqcalprev=freqcal;
    end
    fclose(fileID);
    if fwd==0
        dispselec=struct('modnum',modnum,'modok',modok,...
            'misround',misround,'nfreqsample',nfreqsample);
    end
end
end