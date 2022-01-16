%E-SGFD
% ��ʱ 563.031032 �롣
clear
clc %%%%%%%
close all
nt=1609;    % number of time steps
eps=.6;     % stability
isnap=20;    % snapshot sampling
load('vv')

c1=flipud(c);

v=c1;
nx=800;
nx=nx+45*2;
nz=475;
nz=nz+45*2;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90
        vv(ii+45,jj+45)=v(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45
        vv(ii+45,jj)=v(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45:nx
        vv(ii+45,jj)=v(ii,800);
    end
end


for ii=1:45  %%top
    for jj=1:nx
        vv(ii,jj)=vv(46,jj);
    end
end

for ii=nz-44:nz  %%bottom
    for jj=1:nx
        vv(ii,jj)=vv(nz-45,jj);
    end
end


clear v
v=vv;



% vv(25:nz-25,25:nx-25)=v;


dx=20;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.0025; % calculate time step from stability criterion
tau=dt;


f0=40;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^2*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=(diff(src))/dx^2;				% time derivative to obtain gaussian

taper=ones(nz,nx);
for i=1:50
    for j=1:nx
        taper(i,j)=0.5-0.5*cos(pi*(i-1)/(50-1));
        taper(nz-i+1,j)=taper(i,j);
    end
end
for i=1:nz
    for j=1:50
        taper(i,j)=taper(i,j)*(0.5-0.5*cos(pi*(j-1)/(50-1)));
        taper(i,nx-j+1)=taper(i,j);
    end
end

zs=46;
xs=600-150+25;

seis_record=zeros(nt,nx);
p=zeros([nz nx]); pnew=p; pold=p;
d2px=zeros(nz,nx,7);
d1px=zeros(nz,nx,7);
d1pz=zeros(nz,nx,7);
VxzAddedpoint=p;
d1pxz=p;
AddedpointPx=p;
AddedpointPz=p;
r=v*dt/h;

load ('figure5aCoeff.mat')
coeff=zeros(nz,nx,8);
for ii=1:nz
    for jj=1:nx
        coeff(ii,jj,:)=coeffJune(floor((v(ii,jj)-1486)/10)+1,:);
    end
end
tic

Vx=zeros(nz,nx);
Vz=zeros(nz,nx);
M=7;
for it=1:nt-2,
%     it
    
    % it is very time consuming to assign  different coefficients for different
    % velocities
    
    
    d2px(:,:,1)=((circshift(Vx,[0 0])-circshift(Vx,[0 1]))) +(circshift(Vz,[0])-circshift(Vz,[1]));
    d2px(:,:,2)=(circshift(Vx,[0 -1])-circshift(Vx,[0 2]))+(circshift(Vz,[-1])-circshift(Vz,[2])) ;
    d2px(:,:,3)=(circshift(Vx,[0 -2])-circshift(Vx,[0 3])) +(circshift(Vz,[-2])-circshift(Vz,[3]));
    d2px(:,:,4)=(circshift(Vx,[0 -3])-circshift(Vx,[0 4])) +(circshift(Vz,[-3])-circshift(Vz,[4]));
    d2px(:,:,5)=(circshift(Vx,[0 -4])-circshift(Vx,[0 5])) +(circshift(Vz,[-4])-circshift(Vz,[5]));
    d2px(:,:,6)=(circshift(Vx,[0 -5])-circshift(Vx,[0 6])) +(circshift(Vz,[-5])-circshift(Vz,[6]));
    d2px(:,:,7)=(circshift(Vx,[0 -6])-circshift(Vx,[0 7])) +(circshift(Vz,[-6])-circshift(Vz,[7]));
    
    
    VxzAddedpoint=(circshift(Vx,[-1 0])-circshift(Vx,[-1 1]) +circshift(Vx,[1 0])-circshift(Vx,[1 1])   ) ...
        +(circshift(Vz,[0 -1])-circshift(Vz,[1 -1]) +circshift(Vz,[0 1])-circshift(Vz,[1 1])   );
    
    
    %%whether important reduce in compuratation time
    d1pxz=coeff(:,:,1).*(d2px(:,:,1));
    
    for m=2:M
        d1pxz=d1pxz+coeff(:,:,m).*(d2px(:,:,m));
    end
    
    d1pxz=d1pxz+coeff(:,:,M+1).*(VxzAddedpoint);
    
    
    
    p=p-dt*v.^2.*d1pxz/h;
    p(zs,xs)= p(zs,xs)+src(it);
    
    d1px(:,:,1)=(circshift(p,[0 -1])-circshift(p,[0 0]));
    d1px(:,:,2)=(circshift(p,[0 -2])-circshift(p,[0 1]));
    d1px(:,:,3)=(circshift(p,[0 -3])-circshift(p,[0 2]));
    d1px(:,:,4)=(circshift(p,[0 -4])-circshift(p,[0 3]));
    d1px(:,:,5)=(circshift(p,[0 -5])-circshift(p,[0 4]));
    d1px(:,:,6)=(circshift(p,[0 -6])-circshift(p,[0 5]));
    d1px(:,:,7)=(circshift(p,[0 -7])-circshift(p,[0 6]));
    
    
    
    d1pz(:,:,1)=(circshift(p,[-1])-circshift(p,[0]));
    d1pz(:,:,2)=(circshift(p,[-2])-circshift(p,[1]));
    d1pz(:,:,3)=(circshift(p,[-3])-circshift(p,[2]));
    d1pz(:,:,4)=(circshift(p,[-4])-circshift(p,[3]));
    d1pz(:,:,5)=(circshift(p,[-5])-circshift(p,[4]));
    d1pz(:,:,6)=(circshift(p,[-6])-circshift(p,[5]));
    d1pz(:,:,7)=(circshift(p,[-7])-circshift(p,[6]));
    
    
    AddedpointPx=(circshift(p,[-1 -1])-circshift(p,[-1 0]) +circshift(p,[1 -1])-circshift(p,[1 0])   ) ;
    AddedpointPz= (circshift(p,[-1 -1])-circshift(p,[0 -1]) +circshift(p,[-1 1])-circshift(p,[0 1])   );
    
        
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    d1pxSum=coeff(:,:,1).*(d1px(:,:,1));
    
    for m=2:M
        d1pxSum=d1pxSum+coeff(:,:,m).*(d1px(:,:,m));
    end
    
    d1pxSum=d1pxSum+coeff(:,:,M+1).*(AddedpointPx);
    
    d1pzSum=coeff(:,:,1).*(d1pz(:,:,1));
    
    for m=2:M
        d1pzSum=d1pzSum+coeff(:,:,m).*(d1pz(:,:,m));
    end
    
    d1pzSum=d1pzSum+coeff(:,:,M+1).*(AddedpointPz);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    Vx=Vx-dt*d1pxSum/h;
    Vz=Vz-dt*d1pzSum/h;
    
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
%     Vx=taper.*Vx;
%     Vz=Vz.*taper; % for reviewing
    if rem(it,isnap)== 0,
        imagesc(x,z,p), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(p))))
        drawnow
    end
    
    
    seis_record(it,:)=p(zs,:);
      if it==500
        pp1=p;
    elseif it==1000
        pp2=p;
    elseif it==1500
        pp3=p;
    elseif it==2000
        pp4=p;
    elseif it==2500
        pp5=p;
    end
    
    if it==500
        Vx1=Vx;
    elseif it==1000
        Vx2=Vx;
    elseif it==1500
        Vx3=Vx;
    elseif it==2000
        Vx4=Vx;
    elseif it==2500
        Vx5=Vx;
    end
end
toc
save('figure5a40Hz.mat')
figure;imagesc(v(45:end-45,45:end-45))
xlabel('x/dx')
ylabel('z/dz')
h = colorbar;
set(get(h,'title'),'string','m/s');
hold on ;plot(xs,zs-45,'*r')