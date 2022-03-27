clear all, close all;

%{
-------------------------------------------------------------------------
                    Comentarios/conclusiones/dudas

-------------------------------------------------------------------------
%}

m=0.1; F=0.1; l=0.6; g=9.8; M=0.5; 

X01=[0 0 -0.01 0];
deltat=10^-4;
ts=10;
u=zeros(1,ts/deltat);

t=0:deltat:(ts-deltat);
valores=zeros(4,ts/deltat);
valores(1,1)=X01(1);
valores(2,1)=X01(2);
valores(3,1)=X01(3);
valores(4,1)=X01(4);
delta_dd=0;
phi_dd=0;

for i=1:1:(ts/deltat-1)
       delta_dd=(u(1,i)-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*l*phi_dd*cos(valores(3,i)))/(M+m);
       phi_dd=(g*sin(valores(3,i))-delta_dd*cos(valores(3,i)))/l;
       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
       valores(2,i+1)=valores(2,i)+deltat*delta_dd;
       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
       valores(4,i+1)=valores(4,i)+deltat*phi_dd;
%       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
%       valores(2,i+1)=valores(2,i)+deltat*(u-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*g*cos(valores(3,i))*sin(valores(3,i)))/(M+m-m*cos(valores(3,i))*cos(valores(3,i)));
%       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
%       valores(4,i+1)=valores(4,i)+deltat*(g*sin(valores(3,i))/l-valores(2,i+1)*cos(valores(3,i))/l);
%esta version de Euler no funciona, no se por qué, del sistema no lineal
end

figure
plot(t,valores(3,:))
xlabel('Tiempo');
ylabel('Ángulo respecto a la normal');
title('Ángulo inicial = -0.01, m=0.1');
hold on;



m=0.1; F=0.1; l=0.6; g=9.8; M=0.5; 

X01=[0 0 3.01 0];
deltat=10^-4;
ts=10;
u=zeros(1,ts/deltat);

t=0:deltat:(ts-deltat);
valores=zeros(4,ts/deltat);
valores(1,1)=X01(1);
valores(2,1)=X01(2);
valores(3,1)=X01(3);
valores(4,1)=X01(4);
delta_dd=0;
phi_dd=0;

for i=1:1:(ts/deltat-1)
       delta_dd=(u(1,i)-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*l*phi_dd*cos(valores(3,i)))/(M+m);
       phi_dd=(g*sin(valores(3,i))-delta_dd*cos(valores(3,i)))/l;
       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
       valores(2,i+1)=valores(2,i)+deltat*delta_dd;
       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
       valores(4,i+1)=valores(4,i)+deltat*phi_dd;
%       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
%       valores(2,i+1)=valores(2,i)+deltat*(u-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*g*cos(valores(3,i))*sin(valores(3,i)))/(M+m-m*cos(valores(3,i))*cos(valores(3,i)));
%       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
%       valores(4,i+1)=valores(4,i)+deltat*(g*sin(valores(3,i))/l-valores(2,i+1)*cos(valores(3,i))/l);
%esta version de Euler no funciona, no se por qué, del sistema no lineal
end

figure
plot(t,valores(3,:))
xlabel('Tiempo');
ylabel('Ángulo respecto a la normal');
title('Ángulo inicial = 3.01, m=0.1');
hold on;










%Repeticion con el doble de masa
m=0.2; F=0.1; l=0.6; g=9.8; M=0.5; 

X01=[0 0 -0.01 0];
deltat=10^-4;
ts=10;
u=zeros(1,ts/deltat);

t=0:deltat:(ts-deltat);
valores=zeros(4,ts/deltat);
valores(1,1)=X01(1);
valores(2,1)=X01(2);
valores(3,1)=X01(3);
valores(4,1)=X01(4);
delta_dd=0;
phi_dd=0;

for i=1:1:(ts/deltat-1)
       delta_dd=(u(1,i)-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*l*phi_dd*cos(valores(3,i)))/(M+m);
       phi_dd=(g*sin(valores(3,i))-delta_dd*cos(valores(3,i)))/l;
       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
       valores(2,i+1)=valores(2,i)+deltat*delta_dd;
       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
       valores(4,i+1)=valores(4,i)+deltat*phi_dd;
%       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
%       valores(2,i+1)=valores(2,i)+deltat*(u-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*g*cos(valores(3,i))*sin(valores(3,i)))/(M+m-m*cos(valores(3,i))*cos(valores(3,i)));
%       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
%       valores(4,i+1)=valores(4,i)+deltat*(g*sin(valores(3,i))/l-valores(2,i+1)*cos(valores(3,i))/l);
%esta version de Euler no funciona, no se por qué, del sistema no lineal
end

figure
plot(t,valores(3,:))
xlabel('Tiempo');
ylabel('Ángulo respecto a la normal');
title('Ángulo inicial = -0.01, m=0.2');
hold on;



m=0.2; F=0.1; l=0.6; g=9.8; M=0.5; 

X01=[0 0 3.01 0];
deltat=10^-4;
ts=10;
u=zeros(1,ts/deltat);

t=0:deltat:(ts-deltat);
valores=zeros(4,ts/deltat);
valores(1,1)=X01(1);
valores(2,1)=X01(2);
valores(3,1)=X01(3);
valores(4,1)=X01(4);
delta_dd=0;
phi_dd=0;

for i=1:1:(ts/deltat-1)
       delta_dd=(u(1,i)-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*l*phi_dd*cos(valores(3,i)))/(M+m);
       phi_dd=(g*sin(valores(3,i))-delta_dd*cos(valores(3,i)))/l;
       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
       valores(2,i+1)=valores(2,i)+deltat*delta_dd;
       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
       valores(4,i+1)=valores(4,i)+deltat*phi_dd;
%       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
%       valores(2,i+1)=valores(2,i)+deltat*(u-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*g*cos(valores(3,i))*sin(valores(3,i)))/(M+m-m*cos(valores(3,i))*cos(valores(3,i)));
%       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
%       valores(4,i+1)=valores(4,i)+deltat*(g*sin(valores(3,i))/l-valores(2,i+1)*cos(valores(3,i))/l);
%esta version de Euler no funciona, no se por qué, del sistema no lineal
end

figure
plot(t,valores(3,:))
xlabel('Tiempo');
ylabel('Ángulo respecto a la normal');
title('Ángulo inicial = 3.01, m=0.2');
hold on;




%Evaluación lineal vs no lineal equilibrio inestable
m=0.01; F=0.1; l=1.2; g=9.8; M=0.5; 
X0=[0 0 0 0]; u0=[0];
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
ts=2.5;
u=zeros(1,ts/deltat);
t=0:deltat:(ts-deltat);

valoreslin=zeros(4,ts/deltat);
X=[0 0 -0.01 0];

for i=1:1:(ts/deltat)
    valoreslin(1,i)=X(1,1);
    valoreslin(2,i)=X(1,2);
    valoreslin(3,i)=X(1,3);
    valoreslin(4,i)=X(1,4);
    xp=A*(X-X0)'+B*u(1,i);
    X=X+deltat*xp';
end  

valores=zeros(4,ts/deltat);
X01=[0 0 -0.01 0];
valores(1,1)=X01(1);
valores(2,1)=X01(2);
valores(3,1)=X01(3);
valores(4,1)=X01(4);
delta_dd=0;
phi_dd=0;

for i=1:1:(ts/deltat-1)
       delta_dd=(u(1,i)-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*l*phi_dd*cos(valores(3,i)))/(M+m);
       phi_dd=(g*sin(valores(3,i))-delta_dd*cos(valores(3,i)))/l;
       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
       valores(2,i+1)=valores(2,i)+deltat*delta_dd;
       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
       valores(4,i+1)=valores(4,i)+deltat*phi_dd;
%       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
%       valores(2,i+1)=valores(2,i)+deltat*(u-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*g*cos(valores(3,i))*sin(valores(3,i)))/(M+m-m*cos(valores(3,i))*cos(valores(3,i)));
%       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
%       valores(4,i+1)=valores(4,i)+deltat*(g*sin(valores(3,i))/l-valores(2,i+1)*cos(valores(3,i))/l);
%esta version de Euler no funciona, no se por qué, del sistema no lineal
end

figure
plot(t,valoreslin(3,:));
hold on;
plot(t,valores(3,:));
legend({'Linealizado','No lineal'},'Location','southwest')
xlabel('Tiempo');
ylabel('Ángulo respecto a la normal');
title('Ángulo inicial = -0.01, equilibrio inestable, comparación');
%se observa la divergencia del sistema linealizado





%Evaluación lineal vs no lineal equilibrio estable
m=0.5; F=0.1; l=12; g=9.8; M=0.5; 
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
X=[0 0 pi-0.8 0];

for i=1:1:(ts/deltat)
    valoreslin(1,i)=X(1,1);
    valoreslin(2,i)=X(1,2);
    valoreslin(3,i)=X(1,3);
    valoreslin(4,i)=X(1,4);
    xp=A*(X-X0)'+B*u(1,i);
    X=X+deltat*xp';
end  

valores=zeros(4,ts/deltat);
X01=[0 0 pi-0.8 0];
valores(1,1)=X01(1);
valores(2,1)=X01(2);
valores(3,1)=X01(3);
valores(4,1)=X01(4);
delta_dd=0;
phi_dd=0;

for i=1:1:(ts/deltat-1)
       delta_dd=(u(1,i)-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*l*phi_dd*cos(valores(3,i)))/(M+m);
       phi_dd=(g*sin(valores(3,i))-delta_dd*cos(valores(3,i)))/l;
       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
       valores(2,i+1)=valores(2,i)+deltat*delta_dd;
       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
       valores(4,i+1)=valores(4,i)+deltat*phi_dd;
%       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
%       valores(2,i+1)=valores(2,i)+deltat*(u-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*g*cos(valores(3,i))*sin(valores(3,i)))/(M+m-m*cos(valores(3,i))*cos(valores(3,i)));
%       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
%       valores(4,i+1)=valores(4,i)+deltat*(g*sin(valores(3,i))/l-valores(2,i+1)*cos(valores(3,i))/l);
%esta version de Euler no funciona, no se por qué, del sistema no lineal
end

figure
plot(t,valoreslin(3,:));
hold on;
plot(t,valores(3,:));
legend({'Linealizado','No lineal'},'Location','southwest')
xlabel('Tiempo');
ylabel('Ángulo respecto a la normal');
title('Ángulo inicial = pi-0.8, equilibrio estable, comparación');
%se observa la divergencia del sistema linealizado




%Si m y l se mantienen bajos:
m=0.01; F=0.1; l=0.6; g=9.8; M=0.5; 
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
X=[0 0 pi-0.8 0];

for i=1:1:(ts/deltat)
    valoreslin(1,i)=X(1,1);
    valoreslin(2,i)=X(1,2);
    valoreslin(3,i)=X(1,3);
    valoreslin(4,i)=X(1,4);
    xp=A*(X-X0)'+B*u(1,i);
    X=X+deltat*xp';
end  

valores=zeros(4,ts/deltat);
X01=[0 0 pi-0.8 0];
valores(1,1)=X01(1);
valores(2,1)=X01(2);
valores(3,1)=X01(3);
valores(4,1)=X01(4);
delta_dd=0;
phi_dd=0;

for i=1:1:(ts/deltat-1)
       delta_dd=(u(1,i)-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*l*phi_dd*cos(valores(3,i)))/(M+m);
       phi_dd=(g*sin(valores(3,i))-delta_dd*cos(valores(3,i)))/l;
       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
       valores(2,i+1)=valores(2,i)+deltat*delta_dd;
       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
       valores(4,i+1)=valores(4,i)+deltat*phi_dd;
%       valores(1,i+1)=valores(1,i)+deltat*valores(2,i);
%       valores(2,i+1)=valores(2,i)+deltat*(u-F*valores(2,i)+m*l*valores(4,i)*valores(4,i)*sin(valores(3,i))-m*g*cos(valores(3,i))*sin(valores(3,i)))/(M+m-m*cos(valores(3,i))*cos(valores(3,i)));
%       valores(3,i+1)=valores(3,i)+deltat*valores(4,i);
%       valores(4,i+1)=valores(4,i)+deltat*(g*sin(valores(3,i))/l-valores(2,i+1)*cos(valores(3,i))/l);
%esta version de Euler no funciona, no se por qué, del sistema no lineal
end

figure
plot(t,valoreslin(3,:));
hold on;
plot(t,valores(3,:));
legend({'Linealizado','No lineal'},'Location','southwest')
xlabel('Tiempo');
ylabel('Ángulo respecto a la normal');
title('Ángulo inicial = pi-0.8, equilibrio estable, comparación');
%se observa la divergencia del sistema linealizado, aunque es mas lenta que
%si l y m fueran mas grandes







%{
%Sistema lineal
A=[0 1 0 0; 0 -F/M -m*g/M 0; 0 0 0 1; 0 F/(l*M) g*(1+m/M)/l 0];
B=[0; 1/M; 0; -1/(M*l)];
C=[0 0 1 0];
D=0;

valoreslin=zeros(4,ts/deltat);
% valoreslin(1,1)=X01(1);
% valoreslin(2,1)=X01(2);
% valoreslin(3,1)=X01(3);
% valoreslin(4,1)=X01(4);
X=[0 0 -0.01 0];
aux=[0 0 0 0];  %-> equilibrio inestable

for i=1:1:(ts/deltat)
    valoreslin(1,i)=X(1,1);
    valoreslin(2,i)=X(1,2);
    valoreslin(3,i)=X(1,3);
    valoreslin(4,i)=X(1,4);
    xp=A*(X-aux)'+B*u(1,i);
    X=X+deltat*xp';
end    
% valoreslin(1,1)=X01(1);
% valoreslin(2,1)=X01(2);
% valoreslin(3,1)=X01(3);
% valoreslin(4,1)=X01(4);
% 
% for i=2:1:(ts/deltat)
%     valoreslin(1,i)=valoreslin(1,i-1)+deltat*valoreslin(2,i-1);
%     valoreslin(2,i)=valoreslin(2,i-1)+deltat*(-F*valoreslin(2,i-1)/M-m*g*valoreslin(3,i-1)/M);
%     valoreslin(3,i)=valoreslin(3,i-1)+deltat*valoreslin(4,i-1);
%     valoreslin(4,i)=valoreslin(4,i-1)+deltat*(F*valoreslin(2,i-1)/(l*M)+(1+m/M)*g*valoreslin(3,i-1)/l);
% %esta forma de Euler no anda, no se por qué, del sistema linealizado
% end

% figure
%plot(t,valoreslin(3,:));
%}
