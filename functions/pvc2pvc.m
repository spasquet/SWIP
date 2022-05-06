function pvc2pvc(filein,nW,dx,incerr,maxerr,minvelerr,fileout)

%%% S. Pasquet - V16.11.18
% Convert old pvc files (2 columns) in new pvc files (3 columns)
% pvc2pvc(filein,nW,dx,incerr,maxerr,minvelerr,fileout)

olddata=dlmread(filein,'',1,0);
freq=olddata(:,1);
vs=olddata(:,2);
deltac=lorentzerr(vs',vs'./freq',nW,dx,incerr,maxerr,minvelerr);
dlmwrite(fileout,[freq,vs,deltac],'delimiter','\t');

end
