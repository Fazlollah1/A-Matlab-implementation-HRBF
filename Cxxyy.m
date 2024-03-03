function [cxxyy,shapeP] = Cxxyy(Supportctrs)
global Ns
global mp
global ep 
Au=zeros(Ns+mp,Ns+mp);
Lb=zeros(Ns+mp,1);
%constracting distance matrix to compute the HRFD wights
for i=1:Ns
    for j=1:Ns
         r=sqrt((Supportctrs(i,1)-Supportctrs(j,1)).^2+(Supportctrs(i,2)-Supportctrs(j,2)).^2);
         Au(i,j)=sqrt(ep^2+r^2);
        if (i>1 & mp>0 & i+Ns-1<Ns+mp+1)

           r=sqrt((Supportctrs(i,1)-Supportctrs(j,1)).^2+(Supportctrs(i,2)-Supportctrs(j,2)).^2);
           Au(i+Ns-1,j)=(2/sqrt(ep^2+r^2))-(r^2/(ep^2+r^2)^(3/2));

        end
  end
     r=sqrt((Supportctrs(1,1)-Supportctrs(i,1)).^2+(Supportctrs(1,2)-Supportctrs(i,2)).^2);
     Lb(i)=(2/sqrt(ep^2+r^2))-(r^2/(ep^2+r^2)^(3/2));
end

for i=1:Ns
    
x=Supportctrs(i,1);
y=Supportctrs(i,2);
Rj(1)=1;
Rj(2)=x;
Rj(3)=y;
Rj(4)=x^2;
Rj(5)=x*y;
Rj(6)=y^2;
Rj(7)=x^3;
Rj(8)=(x^2)*y;
Rj(9)=x*y^2;
Rj(10)=y^3;
Rj(11)=x^4;
Rj(12)=(x^3)*y;
Rj(13)=(x^2)*y^2;
Rj(14)=x*y^3;
Rj(15)=y^4;
Rj(16)=x^5;              
Rj(17)=(x^4)*(y^1);      
Rj(18)=(x^3)*(y^2);       
Rj(19)=(x^2)*(y^3);     
Rj(20)=(x^1)*(y^4);     
Rj(21)=y^5;              
Rj(22)=x^6;              
Rj(23)=(x^5)*(y^1);      
Rj(24)=(x^4)*(y^2);     


LRj(1)=0;
LRj(2)=0;
LRj(3)=0;
LRj(4)=2;
LRj(5)=0;
LRj(6)=2;
LRj(7)=6*x;
LRj(8)=2*y;
LRj(9)=x*2;
LRj(10)=6*y;
LRj(11)=12*x^2;
LRj(12)=6*x*y;
LRj(13)=2*(y^2)+2*(x^2);
LRj(14)=6*x*y;
LRj(15)=12*y^2;
LRj(16)=20*x^3;
LRj(17)=12*(x^2)*(y^1);
LRj(18)=6*(x^1)*(y^2)+2*x^3;
LRj(19)=2*y^3+6*(x^2)*(y^1);
LRj(20)=12*(y^2)*(x^1);
LRj(21)=20*y^3;
LRj(22)=30*x^4;
LRj(23)=20*(x^3)*(y^1);
LRj(24)=12*(x^2)*(y^2)+2*x^4;
 
    for j=Ns+1:Ns+mp
       Au(i,j)=Rj(j-Ns);
      
       if (i>1  & mp>0 & i+Ns-1<Ns+mp+1 )
        Au(i+Ns-1,j)=LRj(j-Ns); 
        end
    end
end
if mp>0
 Lb(Ns+1)=0;
 Lb(Ns+2)=0;
end

 for i=Ns+3:Ns+mp
x=Supportctrs(1,1);
y=Supportctrs(1,2);

Rj(1)=1;
Rj(2)=x;
Rj(3)=y;
Rj(4)=x^2;
Rj(5)=x*y;
Rj(6)=y^2;
Rj(7)=x^3;
Rj(8)=(x^2)*y;
Rj(9)=x*y^2;
Rj(10)=y^3;
Rj(11)=x^4;
Rj(12)=(x^3)*y;
Rj(13)=(x^2)*y^2;
Rj(14)=x*y^3;
Rj(15)=y^4;
Rj(16)=x^5;              
Rj(17)=(x^4)*(y^1);      
Rj(18)=(x^3)*(y^2);       
Rj(19)=(x^2)*(y^3);     
Rj(20)=(x^1)*(y^4);     
Rj(21)=y^5;              
Rj(22)=x^6;              
Rj(23)=(x^5)*(y^1);      
Rj(24)=(x^4)*(y^2);     


LRj(1)=0;
LRj(2)=0;
LRj(3)=0;
LRj(4)=2;
LRj(5)=0;
LRj(6)=2;
LRj(7)=6*x;
LRj(8)=2*y;
LRj(9)=x*2;
LRj(10)=6*y;
LRj(11)=12*x^2;
LRj(12)=6*x*y;
LRj(13)=2*(y^2)+2*(x^2);
LRj(14)=6*x*y;
LRj(15)=12*y^2;
LRj(16)=20*x^3;
LRj(17)=12*(x^2)*(y^1);
LRj(18)=6*(x^1)*(y^2)+2*x^3;
LRj(19)=2*y^3+6*(x^2)*(y^1);
LRj(20)=12*(y^2)*(x^1);
LRj(21)=20*y^3;
LRj(22)=30*x^4;
LRj(23)=20*(x^3)*(y^1);
LRj(24)=12*(x^2)*(y^2)+2*x^4;
 


Lb(i)=LRj(i-Ns);
 end 

 Bu=Au';
 %*******************************
 while (cond(Bu)> (5e+14)) && (ep >0.1)
     ep=0.9*ep;
     %constracting distance matrix to compute the HRFD wights
for i=1:Ns
    for j=1:Ns
         r=sqrt((Supportctrs(i,1)-Supportctrs(j,1)).^2+(Supportctrs(i,2)-Supportctrs(j,2)).^2);
         Au(i,j)=sqrt(ep^2+r^2);
        if (i>1 & mp>0 & i+Ns-1<Ns+mp+1)

           r=sqrt((Supportctrs(i,1)-Supportctrs(j,1)).^2+(Supportctrs(i,2)-Supportctrs(j,2)).^2);
           Au(i+Ns-1,j)=(2/sqrt(ep^2+r^2))-(r^2/(ep^2+r^2)^(3/2));

        end
  end
     r=sqrt((Supportctrs(1,1)-Supportctrs(i,1)).^2+(Supportctrs(1,2)-Supportctrs(i,2)).^2);
     Lb(i)=(2/sqrt(ep^2+r^2))-(r^2/(ep^2+r^2)^(3/2));
end

for i=1:Ns
    
x=Supportctrs(i,1);
y=Supportctrs(i,2);
Rj(1)=1;
Rj(2)=x;
Rj(3)=y;
Rj(4)=x^2;
Rj(5)=x*y;
Rj(6)=y^2;
Rj(7)=x^3;
Rj(8)=(x^2)*y;
Rj(9)=x*y^2;
Rj(10)=y^3;
Rj(11)=x^4;
Rj(12)=(x^3)*y;
Rj(13)=(x^2)*y^2;
Rj(14)=x*y^3;
Rj(15)=y^4;
Rj(16)=x^5;              
Rj(17)=(x^4)*(y^1);      
Rj(18)=(x^3)*(y^2);       
Rj(19)=(x^2)*(y^3);     
Rj(20)=(x^1)*(y^4);     
Rj(21)=y^5;              
Rj(22)=x^6;              
Rj(23)=(x^5)*(y^1);      
Rj(24)=(x^4)*(y^2);     


LRj(1)=0;
LRj(2)=0;
LRj(3)=0;
LRj(4)=2;
LRj(5)=0;
LRj(6)=2;
LRj(7)=6*x;
LRj(8)=2*y;
LRj(9)=x*2;
LRj(10)=6*y;
LRj(11)=12*x^2;
LRj(12)=6*x*y;
LRj(13)=2*(y^2)+2*(x^2);
LRj(14)=6*x*y;
LRj(15)=12*y^2;
LRj(16)=20*x^3;
LRj(17)=12*(x^2)*(y^1);
LRj(18)=6*(x^1)*(y^2)+2*x^3;
LRj(19)=2*y^3+6*(x^2)*(y^1);
LRj(20)=12*(y^2)*(x^1);
LRj(21)=20*y^3;
LRj(22)=30*x^4;
LRj(23)=20*(x^3)*(y^1);
LRj(24)=12*(x^2)*(y^2)+2*x^4;
 
    for j=Ns+1:Ns+mp
       Au(i,j)=Rj(j-Ns);
      
       if (i>1  & mp>0 & i+Ns-1<Ns+mp+1 )
        Au(i+Ns-1,j)=LRj(j-Ns); 
        end
    end
end
if mp>0
 Lb(Ns+1)=0;
 Lb(Ns+2)=0;
end

 for i=Ns+3:Ns+mp
x=Supportctrs(1,1);
y=Supportctrs(1,2);

Rj(1)=1;
Rj(2)=x;
Rj(3)=y;
Rj(4)=x^2;
Rj(5)=x*y;
Rj(6)=y^2;
Rj(7)=x^3;
Rj(8)=(x^2)*y;
Rj(9)=x*y^2;
Rj(10)=y^3;
Rj(11)=x^4;
Rj(12)=(x^3)*y;
Rj(13)=(x^2)*y^2;
Rj(14)=x*y^3;
Rj(15)=y^4;
Rj(16)=x^5;              
Rj(17)=(x^4)*(y^1);      
Rj(18)=(x^3)*(y^2);       
Rj(19)=(x^2)*(y^3);     
Rj(20)=(x^1)*(y^4);     
Rj(21)=y^5;              
Rj(22)=x^6;              
Rj(23)=(x^5)*(y^1);      
Rj(24)=(x^4)*(y^2);     


LRj(1)=0;
LRj(2)=0;
LRj(3)=0;
LRj(4)=2;
LRj(5)=0;
LRj(6)=2;
LRj(7)=6*x;
LRj(8)=2*y;
LRj(9)=x*2;
LRj(10)=6*y;
LRj(11)=12*x^2;
LRj(12)=6*x*y;
LRj(13)=2*(y^2)+2*(x^2);
LRj(14)=6*x*y;
LRj(15)=12*y^2;
LRj(16)=20*x^3;
LRj(17)=12*(x^2)*(y^1);
LRj(18)=6*(x^1)*(y^2)+2*x^3;
LRj(19)=2*y^3+6*(x^2)*(y^1);
LRj(20)=12*(y^2)*(x^1);
LRj(21)=20*y^3;
LRj(22)=30*x^4;
LRj(23)=20*(x^3)*(y^1);
LRj(24)=12*(x^2)*(y^2)+2*x^4;
 


Lb(i)=LRj(i-Ns);
 end
 Bu=Au';
 end
  %*******************************

 cxxyy=Bu\Lb;
shapeP=ep;
end