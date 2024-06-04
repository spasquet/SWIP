function matgpdcreport(dir_in,exportopt,nrun,nmode,nmaxmod,wave,nf,fmin,fmax)

%%% S. Pasquet - V22.05.04
% Read dinver report files
% matgpdcreport(dir_in,exportopt,nrun,nmode,nmaxmod,wave,nf,fmin,fmax)

if nargin>6
    resamp=1;
else
    resamp=0;
end

for j=1:nrun

    for k=1:nmode
        if exportopt==0 || exportopt==2
            com1=['gpdcreport ',fullfile(dir_in,['run_0',num2str(j),'.report']),...
                ' -p',wave,' ',num2str(k-1),' -best ',num2str(nmaxmod),' > ',...
                fullfile(dir_in,['best',num2str(j),'.M',num2str(k-1),'.txt'])];
            unix(com1);
            if resamp==1
                com1=['mv ',fullfile(dir_in,['best',num2str(j),'.M',num2str(k-1),'.txt']),...
                    ' ',fullfile(dir_in,'best_old.txt')];
                unix(com1);
                com1=['gpcurve ',fullfile(dir_in,'best_old.txt'),' -resample ',...
                    num2str(nf),' -min ',num2str(fmin),' -max ',num2str(fmax),...
                    ' > ',fullfile(dir_in,['best',num2str(j),'.M',num2str(k-1),'.txt'])];
                unix(com1);
                delete(fullfile(dir_in,'best_old.txt'));
            end
        end
    end
    if exportopt==1 || exportopt==2
        com1=['gpdcreport ',fullfile(dir_in,['run_0',num2str(j),'.report']),...
            ' -vs -best ',num2str(nmaxmod),' > ',...
            fullfile(dir_in,['vs',num2str(j),'.txt'])];
        unix(com1);
        com1=['gpdcreport ',fullfile(dir_in,['run_0',num2str(j),'.report']),...
            ' -vp -best ',num2str(nmaxmod),' > ',...
            fullfile(dir_in,['vp',num2str(j),'.txt'])];
        unix(com1);
        com1=['gpdcreport ',fullfile(dir_in,['run_0',num2str(j),'.report']),...
            ' -rho -best ',num2str(nmaxmod),' > ',...
            fullfile(dir_in,['rho',num2str(j),'.txt'])];
        unix(com1);
    end
end