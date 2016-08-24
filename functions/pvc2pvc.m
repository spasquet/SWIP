function pvc2pvc(filein,nW,dx,incerr,maxerr,minvelerr,fileout)

%%% S. Pasquet - V16.2.17
% Convert old pvc files (2 columns) in new pvc files (3 columns)

olddata=dlmread(filein,'',1,0);
freq=olddata(:,1);
vs=olddata(:,2);
deltac=lorentzerr(vs',vs'./freq',nW,dx,incerr,maxerr,minvelerr);
dlmwrite(fileout,[freq,vs,deltac],'delimiter','\t');

end
