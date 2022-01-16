clear;
clc
close all

load ('vv.mat') % then we obtain veloicty model c

global M dt h v;
M=7;

AA=zeros(M,M);
b=zeros(M,1);
dt=0.0025;



h=20;


coeffJune28=zeros(floor(max(max(c))- min(min(c)) )+1, M+1);
tic;
iii=1;

for v= floor(min(min(c))):floor(max(max(c)))
    x0=0.001*ones(1,M+1);
    options = optimset('TolFun',10^-20,'TolX',10^-20,'MaxFunEvals',8000,'MaxIter',200);
    
    lb=-5*ones(M+1,1);
    ub=5*ones(M+1,1);
    [x,fval,out,iteration]= fmincon(@myfun,x0,[],[],[],[],lb,ub,[],options);
    coeffJune28(iii,1:length(x)-1)=x(1:end-1);
    coeffJune28(iii,8)=x(end);
    iii=iii+1;
end
toc
save('figure5bCoeff.mat','coeffJune28','dt')