function zz=filter_C(xx,Qx,Ymn1,Sy01,lamda)
% input data:
% xx：Unfiltered Parameters 7m*1
% Qx: The covariance matrix of parameters
% Ymn1: The matrix of transformation from grid to SHC. Grid is 180*360 [89.5:-89.5 0.5:359.5]
% lamda: balance factor for the covariance matrices of signals and noise

% Output data:
% zz: Filtered parameter

%Initialization
m=3717;
gnum=64800; 
gm0=zeros(gnum,1);
gm=zeros(gnum,7);
num=1;
ee(num)=1;
invN=Qy;

%Scale Matrix
D=get_scaleM(60);

while ee(num)>=0.002
num=num+1;
mmk=num-1 
gm0=gm(:,2);

Syn=Sy01+lamda*invN; 
invsyn=(Syn)\eye(size(Syn));
z(:,mmk)=Sy01*invsyn*xx_cs;

for j=1:7
    csm=vector2cs(z((j-1)*m+1:j*m,mmk),60);
    grid=gmt_cs2grid(csm,0,1);
    gg0=reshape(grid',180*360,1);
    gg(:,j)=gg0.^2;%m
end
gm(:,:)=gg; 

for jj=1:7
Sg=gm(:,jj);
Sg1=Sg.^(1/2);
C=D*Ymn1;
S1=C.*Sg1';
Sy0(:,:,jj)=S1*S1';
end

Sy01=blkdiag(Sy0(:,:,1),Sy0(:,:,2),Sy0(:,:,3),Sy0(:,:,4),Sy0(:,:,5),Sy0(:,:,6),Sy0(:,:,7));
ee(num)=(abs(mean(gm(:,2))-mean(gm0)))/(mean(gm(:,2)))
end
zz=z;
