function [xx_cs,Qy]=least_Square(time,mass,QQ)
% input data:
% time: the study months in year
% mass: size: 3717*months (3717: the numbers of mass SHCs (degree 2-60)
% QQ: monthly full covarice matrix of SHCs (3717*3717*months)

% Output data:
% xx_cs: Estimated parameter
% Qx: Covariance matrix of parameter

%Initialization

% Converts geoid coefficients (gc) to mass coefficients (mc)
T=diag(v_gc2mc(60));
IM=eye(3717);
N0=zeros(7*3717,7*3717);
b0=zeros(7*3717,1);
t0=time(1);
tt=time;

for i=1:157

x0(:,i)=mass(:,i);

%Observation matrix
A0=[1*IM (tt(i)-t0)*IM ((tt(i)-t0).^2)*IM cos(2*pi*(tt(i)-t0))*IM sin(2*pi*(tt(i)-t0))*IM cos(4*pi*(tt(i)-t0))*IM sin(4*pi*(tt(i)-t0))*IM];

NNP=inv(T'*100)*inv(QQ(:,:,i))*inv(T*100); 

N0=N0+A0'*NNP*A0;
b0=b0+A0'*NNP*x0(:,i);
end

xx_cs(:,1)=((N0)\eye(size(N0)))*b0;
Qy=(N0)\eye(size(N0));
