function deltac=lorentzerr(vel,lam,nW,dx,fac,maxerr,minvelerr,a)

%%% S. Pasquet - V17.05.23
%%% Compute Lorentzian error according to O'Neill's workflow
% deltac=lorentzerr(vel,lam,nW,dx,fac,maxerr,minvelerr,a)

if exist('fac','var')==0 || isempty(fac)==1    
    fac=10.^(1./sqrt((nW-1)*dx)); % factor to adapt error depending on window size     
end

if exist('maxerr','var')==0 || isempty(maxerr)==1    
    maxerr=0.5; % Maximum error ratio (0.5 means that the error won't be higher than 50% of the velocity)
end

if exist('minvelerr','var')==0 || isempty(minvelerr)==1    
    minvelerr=5; % Minimum (rror (in m/s)
end

if exist('a','var')==0 || isempty(a)==1
    a=0.75; % Default a parameter (0.5 recommended by A. O'Neill, 0.75 seems to fit better)
end

A=zeros(length(vel),1);
B=A; DELTAc=A; deltac=A;
% Error calculation
A(vel>0)=(1./(vel(vel>0)))-1./((2*vel(vel>0)./lam(vel>0))*((nW-1)*fac)*dx);
B(vel>0)=(1./(vel(vel>0)))+1./((2*vel(vel>0)./lam(vel>0))*((nW-1)*fac)*dx);
DELTAc(vel>0)=abs((1./A(vel>0))-(1./B(vel>0)));
deltac(vel>0)=(10.^(-a)).*DELTAc(vel>0);
% Error higher limit
deltac(vel'>0 & deltac>(maxerr*vel)')=maxerr*vel(vel'>0 & deltac>(maxerr*vel)')';
% Error lower limit
deltac(vel'>0 & deltac<minvelerr)=minvelerr;
