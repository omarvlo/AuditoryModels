clear; clc; close all;

Nx=251;
L=3.5;		%--- cochlear length
Nf = 220;
fmax = 25000;
x=linspace(0,L,Nx); 
f=logspace(1.3,log10(fmax),Nf)';

%% Membragrama Dr. Mario

f_mario = (1/2*pi).*sqrt(((1e9*exp(-2*x))/0.15)-((200^2)/(2*(0.15^2))));
f_mario = f_mario/10;

x_mario = (log(((((2*pi.*f).^2)*0.15)+((40000)/(2*0.15)))/10^9))/-2;
x_mario = (log(((4*0.15*pi^2)*(f.^2)+((200^2)/(2*0.15))^2)/10^9))/(-2)


%% Membragrama Omar

Nx=251;
%L=3.5;		%--- cochlear length
Nf = 220;
fmax = 25000;
% x=linspace(0,L,Nx); 
% f=logspace(1.3,log10(fmax),Nf)';
% 
% f_omar = 3.695*1e4*exp(-1.485*x);
% 
% plot(x,f_omar,'b','LineWidth',1)
% hold on
L=0.035;
x=linspace(0,L,Nx); 
%f_omar = 3.455*1e4*exp(-1.394*x);
f_omar = 2.003*1e4*exp(-141.2*x)
%x = x*100;
% plot(x,f_omar,'m','LineWidth',1)
% 
% legend('Antigua','Propuesta')
% lgd = legend;
% lgd.FontSize = 20;
% 
% axis([0 3.5 0 40000])

% a = 2.021e4
% b = -142.2
x_omar = (log(f/(3.695*1e4)))/-1.485;

%% Membragrama Greenwood

Nx=251;
L=35;		
x=linspace(0,L,Nx);


A = 165.4;
a = 0.06;

f_greenwood = A*((10.^(a*(L-x)))-1);
x_greenwood = L-(16.7*log10((0.006046.*(f))+1));
x_greenwood = x_greenwood/10;
%% Ploteo de frecuenca vs distancia
x = x/10;
figure(1)
plot(x,f_mario,'k','LineWidth',1)
hold on

plot(x,f_omar,'r','LineWidth',1)
hold on

plot(x,f_greenwood,'b','LineWidth',1);
xlabel('Posición en la coclea [cm]'), ylabel('Frecuencia [Hz]')

ax = gca;
ax.FontSize = 20;

legend('Jimenez M.','Propuesta','Greenwood')
lgd = legend;
lgd.FontSize = 20;
axis([0 3.5 0 fmax])

%% Ploteo de distancia vs frecuencia

figure(2)

plot(f,x_mario,'k--','LineWidth',0.8) 
hold on

% plot(f,x_omar,'r','LineWidth',1)
% hold on

plot(f,x_greenwood,'b--','LineWidth',0.8);
%axis([0 fmax 0 5.5])
xlabel('Frequency [Hz]'), ylabel('Position along the cochlea [cm]')