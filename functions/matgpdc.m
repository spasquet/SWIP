function status=matgpdc(filevel,nmodemax,wave,nf,fmin,fmax,sampling,filedisp)

% S. Pasquet - V22.05.04
% Forward dispersion curve calculation using gpdc
% status=matgpdc(filevel,nmodemax,wave,nf,fmin,fmax,sampling,filedisp)

% Sampling type
% period = regular sampling in period
% frequency = regular sampling in frequency
% log = regular sampling in log (frequency)
if sampling==0
    sampling='period';
elseif sampling==1
    sampling='frequency';
else
    sampling='log';
end

com1=['gpdc ',filevel,' -' wave ' ' num2str(nmodemax) ...
    ' -n ' num2str(nf) ' -min ' num2str(fmin) ' -max ' ...
    num2str(fmax) ' -s ' sampling ' > ' filedisp];
[status,~]=unix(com1);
end