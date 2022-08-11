% Neely and Kim's (1986) cochlear model 
% Computes Pressures and Velocities
% Ref: S. Neely, D. Kim, "A Model for Active Elements in Cochlear Biomechanics", 
%J. Acoust. Soc. Am. 79(5), pp. 1472-1480, May 1986.


clear
clc
close all

disp(' Neely and Kim�s model (1986)');
disp(' Impedance with Nf=220 freq. points between 100 Hz and 50kHz and Nx=500 sections.');

Nf=220;
Nx=500;
L=0.035;		%--- cochlear length

f=logspace(2,log10(25000),Nf)';
x=linspace(0,L,Nx);
s=j*2*pi*f;

%[Zp,Zpmod,fQ,Mef]=NKim86Z(f,x);	%--- default gama=1
%fZa=fQ(1,:);fHz=fQ(2,:);fHp=fQ(3,:);
%QZa=fQ(4,:);QHz=fQ(5,:);QHp=fQ(6,:);
%plot(x,log(fZa),x,log(fHz),x,log(fHp))

%or simply: Zp=NKim86Z(f,x);
%--------------------------

%--- discrete computations ------
%================================

gama=10;
g=1;
b=0.4;
k1=1.1e10*exp(-4*x);
c1=200+15000*exp(-2*x);
m1=3e-2;
k2=7e7*exp(-4.4*x);
c2=100*exp(-2.2*x);
m2=0.5e-2;
k3=1e8*exp(-4*x);
c3=20*exp(-0.8*x);
k4=6.15e9*exp(-4*x);
c4=10400*exp(-2*x);

figure(1)
subplot(5,1,1)
plot(x,k1)
title('Par�metro k1'), xlabel('x'), ylabel('k1')
subplot(5,1,2)
plot(x,c1)
title('Par�metro c1'), xlabel('x'), ylabel('c1')
subplot(5,1,3)
plot(x,m1)
title('Par�metro m1'), xlabel('x'), ylabel('m1')
subplot(5,1,4)
plot(x,k2)
title('Par�metro k2'), xlabel('x'), ylabel('k2')
subplot(5,1,5)
plot(x,c2)
title('Par�metro c2'), xlabel('x'), ylabel('c2')

figure(2)
subplot(5,1,1)
plot(x,m2)
title('Par�metro m2'), xlabel('x'), ylabel('m2')
subplot(5,1,2)
plot(x,k3)
title('Par�metro k3'), xlabel('x'), ylabel('k3')
subplot(5,1,3)
plot(x,c3)
title('Par�metro c3'), xlabel('x'), ylabel('c3')
subplot(5,1,4)
plot(x,k4)
title('Par�metro k4'), xlabel('x'), ylabel('k4')
subplot(5,1,5)
plot(x,c4)
title('Par�metro c4'), xlabel('x'), ylabel('c4')



invs=1 ./s;
umf=ones(Nf,1);
umx=ones(1,Nx);

Z1=invs*k1 + c1(umf,:) + m1*s(:,umx);   %--- Z1=k1/s+c1+s*m1;
Z2=invs*k2 + c2(umf,:) + m2*s(:,umx);	%--- Z2=k2/s+c2+s*m2;
Z3=invs*k3 + c3(umf,:);			%--- Z3=k3/s+c3;
Z4=invs*k4 + c4(umf,:);			%--- Z4=k4/s+c4;
Hc=Z2./(Z2+Z3);

Zp=(g/b)*( Z1 + Hc.*(Z3-gama*Z4) );

%---------ZBM=(g/b)*Z1;
%---------ZOC=(g/b)*Hc.*(Z3-gama*Z4);



%======================================================================

disp(' Middle-ear connection')

%-------- middle ear -------------
km=2.1e6;
cm=4000;
mm=45e-2;
As=0.000001;
Am=0.000035;
gm=0.5;
Zm = s*mm + cm + km./s;
Pi = Am/gm/As;
%---------------------------------


disp(' Transmission Line Model')

Dx=L/Nx;
h=0.001;
rho=10;
zs=2*s*rho*Dx;
zp=h*Zp/Dx;
zp_prueba = zp;
zp(:,Nx)=sqrt(2*s*rho*h.*Zp(:,Nx))-zs;

[P,Ps,Zin]=TLineME(Pi,Zm,zs,zp);
VBM = -P./Zp;
VC=Hc.*VBM;

freq_re = db(P);

disp(' Ploting Responses every 100 sections...')
% 
% figure(1)
% P = P';
% 
% r = 1117;
% frecuencia = f(r)
% semilogx(x,db(P(:,r)))
% hold on
% 
% r = 2232;
% frecuencia = f(r)
% semilogx(x,db(P(:,r)))
% hold on
% 
% r = 3347;
% frecuencia = f(r)
% semilogx(x,db(P(:,r)))
% hold on
% 
% r = 4443;
% frecuencia = f(r)
% semilogx(x,db(P(:,r)))
% hold on

figure(4)
semilogx(f,db(P(:,1:100:Nx)))
title('|P(f)| (pressure)'), xlabel('f [Hz]'), ylabel('dB'), %pause
% semilogx(f,fase(P(:,1:100:Nx)))
% title('arg{P(f)}'), xlabel('f [Hz]'), ylabel('x pi'), pause
% 
% 
% 
% semilogx(f,db(VBM(:,1:100:Nx)))
% title('|VBM(f)| (BM velocity)'), xlabel('f [Hz]'), ylabel('dB'), pause
% semilogx(f,fase(VBM(:,1:100:Nx)))
% title('arg{VBM(f)}'), xlabel('f [Hz]'), ylabel('x pi'), pause
% 
% semilogx(f,db(VC(:,1:100:Nx)))
% title('|VC(f)| (cilia velocity)'), xlabel('f [Hz]'), ylabel('dB'), pause
% semilogx(f,fase(VC(:,1:100:Nx)))
% title('arg{VC(f)}'), xlabel('f [Hz]'), ylabel('x pi'), pause
% 
% 
% %---------------POWER------------------------
% clc
% disp(' Ratio of power absorbed by partition and input power for 2, 4 and 8 kHz')
% 
% 
% Potin=real(Ps.*conj(Ps./Zin))/2;
% PotBM =real(P.*conj(P./zp))/2;
% 
% for k=[155,131,107];
%  plot(x/L,PotBM(k,:)/Potin(k)),hold on
% end
% hold off,  xlabel('x/L'),ylabel('dB'),title('Pot/Potin')

%=======================================================

f_round = f;

P_abs = db(P);
for i=1:1:Nx
vector_p = P_abs(:,i);
[M,maximo] = max(vector_p);
distancia = i;
frecuencia(i) = f(maximo);
a = 1;
end

f_axis = zeros(1,50e3);

figure(5)
plot(x,frecuencia)
title('Relaci�n frecuencia-distancia'), xlabel('mm a partir de la ventana oval'), ylabel('Frecuencia (Hz)'),

figure(6)
plot(frecuencia,x)
title('Relaci�n distancia-frecuencia'), xlabel('Frecuencia (Hz)'), ylabel('mm a partir de la ventana oval'),