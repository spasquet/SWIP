function status=matdinver(targetfile,paramfile,nrun,itmax,ns0,ns,nr,dir_out,verbose)

% S. Pasquet - V20.03.23
% matdinver execute dinver inversion through matlab

% status=matdinver(targetfile,paramfile,nrun,itmax,ns0,ns,nr,dir_out,verbose)
%
% targetfile = name of the file containing dispersion curves
% paramfile = name of the file containing the parameter space
% nrun = Nb of run
% itmax = Nb of iteration per run (not for dinver > 2.0)
% ns0 = Nb of starting models
% ns = Nb of models created at each iterations
% nr = Nb of previous models used to build new sub-parameter space for the
% next ns models
% Handles dinver > 2

version = geopsy_version();

for j=1:nrun
    fprintf(['\n      Run ',num2str(j),'\n']);
    dlmwrite(fullfile(dir_out,'input.txt'),1);
        
    if version == 2
        com1=['dinver -i DispersionCurve -optimization -target ',targetfile,...
            ' -param ',paramfile,' -itmax ',num2str(itmax),' -ns0 ',...
            num2str(ns0),' -ns ',num2str(ns),' -nr ',num2str(nr),' -f -nobugreport -o ',...
            fullfile(dir_out,['run_0',num2str(j),'.report']),' < ',fullfile(dir_out,'input.txt')];
    else
        com1=['dinver -i DispersionCurve -optimization -target ',targetfile,...
            ' -param ',paramfile,' -ns0 ',num2str(ns0),' -ns ',num2str(ns*itmax+ns0),' -nr ',num2str(nr),' -f -nobugreport -o ',...
            fullfile(dir_out,['run_0',num2str(j),'.report']),' < ',fullfile(dir_out,'input.txt')];
    end
    
    if verbose==0
        [status,~]=unix(com1);
    else
        [status]=unix(com1);
    end
    if status~=0
        [status,din_err]=unix(com1);
        if isempty(strfind(din_err,'Cannot open file for writing'))==0
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n      Invalid working directory name');
            fprintf('\n   Remove special characters and spaces');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n\n');
    elseif isempty(strfind(din_err,'parameters'))==0
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Invalid .param file');
            fprintf('\n     Go to next Xmid');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!\n\n');
        elseif isempty(strfind(din_err,'target'))==0
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!');
            fprintf('\n   Invalid .target file');
            fprintf('\n     Go to next Xmid');
            fprintf('\n  !!!!!!!!!!!!!!!!!!!!!!\n\n');   
        end
        break
    end
    delete(fullfile(dir_out,'input.txt'));
end