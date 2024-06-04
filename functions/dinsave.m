function dinsave(filename,thk,vp,vs,rho)

%%% S. Pasquet - V22.05.04
% Save velocity model in dinver format for forward calculation
% dinsave(filename,thk,vp,vs,rho)

wsl = ispc_wsl;

filename_unix = unix_wsl_path(filename,wsl);

unix_cmd(['rm -rf -- ',filename_unix],wsl);
nl=length(thk);
if length(vp)==1
    vp=repmat(vp,size(thk));
end
if length(vs)==1
    vs=repmat(vs,size(thk));
end
if length(rho)==1
    rho=repmat(rho,size(thk));
end

fid0=fopen(filename,'w'); % Open file for writing
fprintf(fid0,'%d \n',nl); % Write nb of layers
for i=1:1:nl-1 % Loop over nb of layers - 1
    % Write thickness and velocity in dinver format
    fprintf(fid0,'%f %f %f %f \n',[thk(i) vp(i) vs(i) rho(i)]);
end
% Write last line
fprintf(fid0,'%f %f %f %f \n',[0 vp(i+1) vs(i+1) rho(i+1)]);
fclose(fid0); % Close file
