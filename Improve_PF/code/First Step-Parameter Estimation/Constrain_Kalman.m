function [P11,X11,fomg,PP]=Constrain_Kalman(time,mass,t,sig,F,xx_cs,Qx,QQ)
% input data:
% time: the study months in year
% mass: size: 3717*months (3717: the numbers of mass SHCs (degree 2-60)
% t: the correlation time of FOGM
% sig: the variance of FOGM
% F: constrain matrix of Fisher-test  
% F(:,1): trend test; F(:,2): annual signals test; F(:,3): Semi-annual signals test
% xx_cs: initial parameter (usually obtain from leaat-square)
% Qx: initial covariance matrix of parameter (usually obtain from leaat-square)
% QQ: monthly full covarice matrix of SHCs (3717*3717*months)

% Output data:
% X11: The parameters of determined component in the final months 7m*1
% P11: The covariance matrice of X11 7m*7m
% fogm: The parameters of random component in all months 3717*months
% PP: The covariance matrices of fogm 3717*3717*months

% Converts geoid coefficients (gc) to mass coefficients (mc)
T=diag(v_gc2mc(60));

R=T*QQ(:,:,1)*T'; 
m=3717;
IM=eye(m);
Q=diag(ones(1,8*m),0)*0.000;   %process noise
I=diag(ones(1,8*m),0);
t0=time(1);
%Observation matrix
H=[1*IM (time(1)-t0)*IM ((time(1)-t0)^2)*IM cos(2*pi*(time(1)-t0))*IM sin(2*pi*(time(1)-t0))*IM cos(4*pi*(time(1)-t0))*IM sin(4*pi*(time(1)-t0))*IM 1*IM];

% Constrain Matrix using Fisher-test
HI=blkdiag(IM,diag(F(:,1)),diag(F(:,1)),diag(F(:,2)),diag(F(:,2)),diag(F(:,3)),diag(F(:,3)),IM);
C=eye(8*3717)-HI;
C(all(C==0,2),:)=[];

%Initialized month

%prediction
X01=[xx_cs;zeros(m,1)];            %%
P01=blkdiag(Qy,eye(m)*10.^(14)); 

%Gain matrix
K=P01*H'*inv(H*P01*H'+R(:,:,1));
%Update
v=pos(:,1)-H*X01;
X11=X01+K*v;
P11=(I-K*H)*P01;

%Fisher-test Constrain
X11_=X11-P11*C'*inv(C*P11*C')*C*X11;
P11_=P11-P11*C'*inv(C*P11*C')*C*P11;
XX(:,1)=X11_;
PP(:,1)=diag(P11_(7*m+1:8*m,7*m+1:8*m));
fomg(:,1)=X11_(7*m+1:8*m,1);

% The rest months
for i=2:length(time)
    i
   X00=X11;
   P00=P11;
   H=[1*IM (time(i)-t0)*IM ((time(i)-t0)^2)*IM (cos(2*pi*(time(i)-t0)))*IM (sin(2*pi*(time(i)-t0)))*IM (cos(4*pi*(time(i)-t0)))*IM (sin(4*pi*(time(i)-t0)))*IM 1*IM];
   dt=(time(i)-time(i-1))*365.25;
   B=blkdiag(eye(7*m),eye(m)*diag(exp(-dt/t)));
   Q(7*m+1:8*m,7*m+1:8*m)=eye(m)*sig*(1-exp(-2*dt/t)); 
   %prediction
   X01=B*X00;
   P01=B*P00*B'+Q;
   
   %Gain matrix
    R=T*QQ(:,:,i)*T';
    K=P01*H'*inv(H*P01*H'+R);
   
   %Update
   v=pos(:,i)-H*X01;
   X11=X01+K*v;%8m*1
   P11=(I-K*H)*P01;%8m*8m

   %Fisher-test Constrain
   X11_=X11-P11*C'*inv(C*P11*C')*C*X11;
   P11_=P11-P11*C'*inv(C*P11*C')*C*P11;


   PP(:,:,i)=P11_(7*m+1:8*m,7*m+1:8*m);
   XX(:,i)=X11_;
   fomg(:,i)=X11_(7*m+1:8*m,1);
end
