clear all
%In this code, we use the combination of radial basis functions (RBFs) and polynomials
%to derive compact finite difference (FD) formulas. This method is applied to obtain 
%numerical solution of Poisson equations on regular or irregular domain.The
%method is used as meshless method
global Ns
global mp
global ep 

h=0.025; % 
Ns=25; %number of supports points Ns<=25
mp=Ns-1; 

%Exact solution. The boundary condition are taken from xact solutions.
u=@(x,y) sin(2.*(x - 0.1).^2).* cos((x - 0.3).^2) + (sin((2.* (y - 0.5).^2).^1))./(1 +2 .*x.^2 + y.^2);

uxx=@(x,y)  4.*cos((x - 3./10).^2).*cos(2.*(x - 1./10).^2) - (4.*sin(2.*(y - 1./2).^2))./(2.*x.^2 + y.^2 + 1).^2 - 2.*sin(2.*(x - 1./10).^2).*sin((x - 3./10).^2) +...
(32.*x.^2.*sin(2.*(y - 1./2).^2))./(2.*x.^2 + y.^2 + 1).^3 - sin(2.*(x - 1./10).^2).*cos((x - 3./10).^2).*(2.*x - 3./5).^2 - sin(2.*(x - 1./10).^2).*cos((x - 3./10).^2).*(4.*x - 2./5).^2 -... 
2.*sin((x - 3./10).^2).*cos(2.*(x - 1./10).^2).*(2.*x - 3./5).*(4.*x - 2./5);
uyy=@(x,y) (4.*cos(2.*(y - 1./2).^2))./(2.*x.^2 + y.^2 + 1) - (2.*sin(2.*(y - 1./2).^2))./(2.*x.^2 + y.^2 + 1).^2 + ...
(8.*y.^2.*sin(2.*(y - 1./2).^2))./(2.*x.^2 + y.^2 + 1).^3 - (sin(2.*(y - 1./2).^2).*(4.*y - 2).^2)./(2.*x.^2 + y.^2 + 1) - ...
(4.*y.*cos(2.*(y - 1./2).^2).*(4.*y - 2))./(2.*x.^2 + y.^2 + 1).^2;
 

uxxyy=@(x,y) -((16.*(-0.5 + y).* y.* cos(2 .*(-0.5 + y).^2))./(1 + 2.* x.^2 + y.^2).^2) -...
16.* (-0.3 + x).*(-0.1 + x) .*cos( 2.* (-0.1 + x).^2).* sin((-0.3 +  x).^2) + (-4 .*(-0.3 + x).^2 .*cos((-0.3 + x).^2) - ...
    2.* sin((-0.3 + x).^2)) .*sin(2 .*(-0.1 + x).^2) +  cos((-0.3 + x).^2).* (4.* cos(2 .*(-0.1 + x).^2) - ...
    16.* (-0.1 + x).^2 .*sin(2 .*(-0.1 + x).^2)) + (( 32.* x.^2)./(1 + 2 .*x.^2 + y.^2).^3 - ...
    4./(1 + 2.* x.^2 + y.^2).^2).* sin( 2.* (-0.5 + y).^2) + ((8 .*y.^2)./(1 + 2.* x.^2 + y.^2).^3 - ...
    2./(1 + 2 .*x.^2 + y.^2).^2).*sin(2.*(-0.5 + y).^2) +...
(4.* cos(2.* (-0.5 + y).^2) - 16 .*(-0.5 + y).^2 .*sin(2.* (-0.5 + y).^2))./( 1 + 2.* x.^2 + y.^2);

% Star domain 
x1=-2:h:2;
[X Y]=meshgrid(x1);
ctr0=[X(:) Y(:)];
[m n]=size(ctr0);
teta=0:h:2*pi-h;
mm=1;
for i=1:m
     [teta1,r1] = cart2pol(ctr0(i,1),ctr0(i,2));
     rho=1+0.1*(sin(7*teta1)+sin(teta1).^2);
    if r1<=rho-h/2
        ctr(mm,1)=ctr0(i,1);
        ctr(mm,2)=ctr0(i,2);
            mm=mm+1;

    end
end
teta=0:h:2*pi-h/2;
r0=1+0.1*(sin(7*teta)+sin(teta).^2);
bdctrx1=r0.*cos(teta);
bdctry1=r0.*sin(teta);
bdctrs=[bdctrx1' bdctry1'];
% To plot the domain uncomment the following two lines

% plot(ctr(:,1),ctr(:,2),'.')
% hold on
% plot(bdctrs(:,1),bdctrs(:,2),'r*')

[NI N]=size(ctr);
[Nb N]=size(bdctrs);
AFu=zeros(NI,NI);
BFu=zeros(NI,NI);
LU=zeros(NI,1);
bF=zeros(NI,1);
ctrs=[ctr;bdctrs];
for jc=1:NI
ep=1; %shape parameter
Au=zeros(Ns+mp,Ns+mp);
Bu=zeros(Ns+mp,Ns+mp);
Lb=zeros(Ns+mp,1);
%For each interior point x_jc, find the its (Ns-1) neighboring points.
dis=zeros(NI+Nb,2);
dis(:,1)=sqrt((ctrs(:,1)-ctrs(jc,1)).^2+(ctrs(:,2)-ctrs(jc,2)).^2);
[m2 n2]=size(dis);
 for i=1:m2
     dis(i,2)=i;
 end
sortctrs= sortrows(dis,1);
for i=1:Ns
Supportctrs(i,1)=ctrs(sortctrs(i,2),1);
Supportctrs(i,2)=ctrs(sortctrs(i,2),2);
end

 [cxxyy,shapeP]=Cxxyy(Supportctrs);
 Shape(jc)=shapeP;
 
 %	Put the obtained weights in matrix AFu (A_h). 
     for j=1:Ns
        if sortctrs(j,2)<= NI
            AFu(jc,sortctrs(j,2))=cxxyy(j);
        else
          bF(jc)= bF(jc)-cxxyy(j)*u(ctrs(sortctrs(j,2),1),ctrs(sortctrs(j,2),2));
        end
     end
    
     
         for j=Ns+1:Ns+mp
              bF(jc)= bF(jc)-cxxyy(j)*uxxyy(ctrs(sortctrs(j-Ns+1,2),1),ctrs(sortctrs(j-Ns+1,2),2));
         end
        bF(jc)= bF(jc)+uxxyy(ctrs(jc,1),ctrs(jc,2));
       
 end
Uaproximation=AFu\bF; 
norm2=norm(inv(AFu))
cond2=cond(AFu)

Uexact=u(ctr(:,1),ctr(:,2));
ERROR=abs(Uaproximation-Uexact);

%plot the absolute error
 plot3(ctr(:,1),ctr(:,2),ERROR,'.') 
%plot shape parameter
% plot3(ctr(:,1),ctr(:,2),Shape,'.') ( uncomment the following line)
AbsoluteERROR=max(abs(ERROR))

%The sparsity pattern of matrix AFu (A_h)( uncomment the following line)

%spy(AFu,'r.')
