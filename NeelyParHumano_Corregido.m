% Neely and Kim's (1986) cochlear model 
% Parámetros para el ser humano

clear
clc
close all

disp(' Neely and Kim´s model (1986)');
disp(' Impedance with Nf=220 freq. points between 100 Hz and 50kHz and Nx=500 sections.');

Nf=220;
Nx=251;
L=0.035; 		%--- cochlear length %en metros
fmax = 20000;
fmin = 20;
f=logspace(log10(fmin),log10(fmax),Nf)';
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

gama=1;
g=1;
b=0.4;
k1=4.95e9*exp(-320*(x+0.00375));
c1=1+19700*exp(-179*(x+0.00375));
m1=1.35e-2;
k2=3.15e7*exp(-352*(x+0.00375));
c2=113*exp(-176*(x+0.00375));
m2=2.3e-3;
k3=4.5e7*exp(-320*(x+0.00375));
c3=22.5*exp(-64*(x+0.00375));
k4=2.82e9*exp(-320*(x+0.00375));
c4=9650*exp(-164*(x+0.00375));

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
km=2.63e8;
cm=2.8e4;
mm=2.96e-2;
As=3.2e-6;
Am=5.236e-5; %area efectiva de la membrana tímpanica (2/3) en m2
gm=0.5;
Zm = s*mm + cm + km./s;
Pi = Am/gm/As;
%---------------------------------

disp(' Transmission Line Model')

Dx=L/Nx;
h=0.001;
rho=1000;
zs=2*s*rho*Dx;
zp=h*Zp/Dx;
zp_prueba = zp;
zp(:,Nx)=sqrt(2*s*rho*h.*Zp(:,Nx))-zs;

[P,Ps,Zin]=TLineME(Pi,Zm,zs,zp);
VBM = -P./Zp;
VC=Hc.*VBM;

freq_re = db(P);

disp(' Ploting Responses every 100 sections...')

semilogx(f,db(P(:,1:50:Nx)),'LineWidth',2)
title('Relación frecuencia-presión sobre la membrana basilar'), xlabel('Frecuencia [Hz]'), ylabel('Presión [dB]'), %pause
axis([20 20000 -200 50])
ax = gca;
ax.FontSize = 18;
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

%frecuencia(1) = frecuencia(2);

%f_axis = zeros(1,50e3);

figure(2)
frecuencia_sm = smooth(frecuencia)
x = x*1000;
plot(x,frecuencia_sm,'k','LineWidth',2) 
title('Relación frecuencia-distancia'), xlabel('mm a partir de la ventana oval'), ylabel('Frecuencia (Hz)'),
ax = gca;
ax.FontSize = 18;

figure(3)
%x = x*10;
plot(frecuencia_sm,x,'k','LineWidth',2)
title('Relación distancia-frecuencia'), xlabel('Frecuencia (Hz)'), ylabel('mm a partir de la ventana oval'),
ax = gca;
ax.FontSize = 18;
% Probando la función f a partir del helicotrema

for i=1:1:(length(frecuencia)-1)
    diferencias(i) = ((frecuencia(i+1))./(frecuencia(i)));
end

cftool