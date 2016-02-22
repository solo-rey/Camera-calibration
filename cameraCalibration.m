function [calibrationMat, rotateMat, translateVect, alfa_u, beta_v, u0, v0, reprojectMat, avgErroru, avgErrorv, totaltime] = cameraCalibration(File3D, File2D)
tic;
[matrix3D]=textread('3d.txt');
matrix3D(1,:)=[]; 
[nR3D, nC3D]=size(matrix3D);

if(nC3D~=3)
error('Less than 3 columns in 3D point matrix!!');
end

X = matrix3D(:,1);
Y = matrix3D(:,2);
Z = matrix3D(:,3);

figure;
plot3(matrix3D(:,1),matrix3D(:,2),matrix3D(:,3),'.');
[matrix2D]=textread('2d.txt');
matrix2D(1,:)=[]; 
[nR2D, nC2D]=size(matrix2D);

if(nC2D~=2)
error('Less than 2 columns in 2D point matrix!!');
end

u_vals=matrix2D(:,1);
v_vals=matrix2D(:,2);

figure;
plot(matrix2D(:,1),matrix2D(:,2),'.');

if(nR3D~=nR2D)
error('2D and 3D points are not same!!');
end

o = ones(size(u_vals));
z = zeros(size(u_vals));
AoddRows  = [ X Y Z o z z z z -u_vals.*X -u_vals.* Y -u_vals.*Z -u_vals ];
AevenRows = [ z z z z X Y Z o -v_vals.*X -v_vals.* Y -v_vals.*Z -v_vals ];
A=[AoddRows; AevenRows];
[U, S, V] = svd(A,0);
m = V(:,end);
display(m);
M = reshape(m,4,3)';
display(M);
calibrationMat=M;
abs_lambda=sqrt(M(3,1)^2 + M(3,2)^2 + M(3,3)^2);
M = M / abs_lambda;
fronT=1;
if fronT
s = sign(M(3,4));
else
s = -sign(M(3,4));
end

T(3) = s*M(3,4);
R = zeros(3,3);
R(3,:)=s*M(3,1:3);
m1 = M(1,1:3)';
m2 = M(2,1:3)';
m3 = M(3,1:3)';
m4 = M(1:3,4);
u0 = m1'*m3;
v0 = m2'*m3;
alfa_u=sqrt( m1'*m1 - u0^2 );
beta_v=sqrt( m2'*m2 - v0^2 );
R(1,:) = s*(u0*M(3,1:3) - M(1,1:3) ) / alfa_u;
R(2,:) = s*(v0*M(3,1:3) - M(2,1:3) ) / beta_v;
T(1) = s*(u0*M(3,4) - M(1,4) ) / alfa_u;
T(2) = s*(v0*M(3,4) - M(2,4) ) / beta_v;
T = T';

translateVect=T;
display(R);

[U,S,V] = svd(R)
R = U*S*V';
             
rotateMat=R;
newMatrix2D=zeros(nR3D,2);

for i=1:1:nR3D
    NumberU=M(1,1)* matrix3D(i,1) + M(1,2)* matrix3D(i,2) + M(1,3)* matrix3D(i,3) + M(1,4);
    NumberV=M(2,1)* matrix3D(i,1) + M(2,2)* matrix3D(i,2) + M(2,3)* matrix3D(i,3) + M(2,4);
    DenSity=M(3,1)* matrix3D(i,1) + M(3,2)* matrix3D(i,2) + M(3,3)* matrix3D(i,3) + M(3,4);
    newMatrix2D(i,1)=NumberU/DenSity;
    newMatrix2D(i,2)=NumberV/DenSity;
end

reprojectMat=newMatrix2D;
errorDiff=reprojectMat-matrix2D;
avgErroru=mean(errorDiff(:,1));
avgErrorv=mean(errorDiff(:,2));
totaltime=toc;