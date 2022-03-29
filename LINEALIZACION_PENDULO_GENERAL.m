clc, clear all, close all;

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas
Se pueden obtener las fi genéricas del sistema, derivarlas genéricamente y
obtener expresiones que al ser valuadas en el punto de operación entregan A
y B

Duda/conclusión, el sistema solo se comporta parecido al no lineal si se
linealiza en el equilibrio estable o inestable
-------------------------------------------------------------------------
%}

m=0.1; F=0.1; l=0.6; g=9.8; M=0.5; 
X0=[0 0 pi 0]; u0=[0];
x1=X0(1);
x2=X0(2);
x3=X0(3);
x4=X0(4);

h1=M+m-m*cos(x3)^2;
b=-F/h1;
B=u0-F*x2-m*l*x4^2*sin(x3)-m*g*sin(x3)*cos(x3);
h2=m*l*x4^2*cos(x3)-m*g*(-sin(x3)^2+cos(x3)^2);
h3=2*m*cos(x3)*sin(x3);
h=(h1*h2-B*h3)/(h1^2);
a=2*m*l*x4*sin(x3)/h1;
c=b*cos(x3)/l;
w=g*cos(x3)/l-(1/l)*(h*cos(x3)-B*sin(x3)/h1);
z=-a*cos(x3)/l;
p=1/h1;
q=-p*cos(x3)/l;

A=[0 1 0 0; 0 b h a; 0 0 0 1; 0 c w z];
B=[0; p; 0; q];   %Para cualquier valor X0,u0, estas matrices quedan con el modelo linealizado

%Simulacion
deltat=10^-4;
ts=10;
u=zeros(1,ts/deltat);
t=0:deltat:(ts-deltat);

valoreslin=zeros(4,ts/deltat);
X=[0 0 pi-0.1 0];

for i=1:1:(ts/deltat)
    valoreslin(1,i)=X(1,1);
    valoreslin(2,i)=X(1,2);
    valoreslin(3,i)=X(1,3);
    valoreslin(4,i)=X(1,4);
    xp=A*(X-X0)'+B*u(1,i);
    X=X+deltat*xp';
end  

plot(t,valoreslin(3,:));