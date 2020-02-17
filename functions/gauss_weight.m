function weight = gauss_weight(misfit)

% keyboard
nbin = 10;

t = linspace(min(misfit),max(misfit),nbin);
gauss=(1/(4*std(misfit)*sqrt(2*pi)))*exp(-0.5.*((t-min(misfit))/(4*std(misfit))).^2);
coef=gauss/sum(gauss);

weight=zeros(size(misfit));
for gg=1:length(t)-1
    flag1=find(misfit>=t(gg) & misfit<t(gg+1));
    weight(flag1)=coef(gg)/length(flag1);
end
weight=weight/sum(weight);
