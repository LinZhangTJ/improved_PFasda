function zz=filter_res(fomg,PP,Ymn1,Sy01,D,lamda)
% input data:
% fomg：Unfiltered Parameters 3717*months
% PP: The covariance matrix of parameters 3717*months
% Ymn1: The matrix of transformation from grid to SHC. Grid is 180*360 [89.5:-89.5 0.5:359.5]
% lamda: balance factor for the covariance matrices of signals and noise
% Sy01: Initial signal covariance matrix

% Output data:
% zz: Filtered parameter

%Initialization
m=3717;
gnum=64800; 
gm0=zeros(gnum,1);
gm=zeros(gnum,1);
num=1;
ee(1)=1; 

%Scale Matrix
D=get_scaleM(60);

while ee(num)>=0.002
num=num+1;
mmk=num-1
gm0=gm(:,1); 
for ii=1:157
ii

invN=diag(PP(:,ii));

Syn=Sy01+lamda*invN; 
invsyn=(Syn)\eye(size(Syn));
z(:,ii,mmk)=Sy01*invsyn*fomg(:,ii);
csm=vector2cs(z(:,ii,mmk),60);
grid(:,:,ii)=gmt_cs2grid(csm,0,1);

end
 
%signal variance
%rms^2
for i=1:180
    for j=1:360
        ss(:,:)=grid(i,j,:);
        rmss(i,j)=sqrt(sum(ss.^2)/157);
    end
end
    gg0=reshape(rmss',180*360,1);
    gm=gg0.^2;
    Sg=gm;
Sg1=Sg.^(1/2);
C=D*Ymn1;
S1=C.*Sg1';
Sy01=S1*S1';
ee(num)=(abs(mean(gm(:,1))-mean(gm0)))/(mean(gm(:,1)))
end
zz=z;
end

