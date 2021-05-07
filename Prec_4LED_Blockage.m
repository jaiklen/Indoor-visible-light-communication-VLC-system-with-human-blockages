clc;
clear all;
close all;
radius = 5;
area = pi*radius^2;
R_B = 0.2;
h_B = 0.85;
b = 10;
lambdaM = b./25; %intesity
%%MATLAB Codes to Calculate the Optical Power Distribution of LOS Link at Receiving Plane for a Typical Room
theta = 90;
FOV = 90;
%FOV of a receiver
% Noise-bandwidth factor %
I2 = 0.562;
% Data rate (Bit per second)
Rb = 115200;
% Ambient light power (Ampere) %
Iamb = 7E-8;
% Photodiode responsivity (A/W )%
R = 0.55;
% Electron charge (C)
q = 1.60E-19;
% Amplifier bandwidth (Hz)%
Ba = 4.5E6;
% Amplifier noise density (Ampere/Hz^0.5)%
Iamf = 5e-12 ;
% semi-angle at half power
ml=-log10(2)/log10(cosd(theta));
%Lambertian order of emission
P_LED=15;
%transmitted optical power by individual LED
nLED=60;
% number of LED array nLED*nLED
P_total=nLED*nLED*P_LED;
%Total transmitted power
Adet=10^-4;
%detector physical area of a PD
Ts=1;
%gain of an optical filter; ignore if no filter is used
index=1.5;
%refractive index of a lens at a PD; ignore if no lens is used

G_Con=(index^2)/(sind(FOV).^2);
%gain of an optical concentrator; ignore if no lens is used
%%
lx=5; ly=5; lz=3;
% room dimension in meter
h_T=2.15;
%the distance between source and receiver plane
%led position
x0lc=0;
y0lc=0;
rlc=2;
teta=-pi:0.01:pi;
xlc=rlc*cos(teta)+x0lc;
ylc=rlc*sin(teta)+y0lc;
%hold on
nlc=8;
tet=linspace(-pi,pi,nlc+1);
xilc=rlc*cos(tet)+x0lc;
yilc=rlc*sin(tet)+y0lc;
% %the distance between source and receiver plane

[XT,YT]=meshgrid([-lx/4 lx/4],[-ly/4 ly/4]);
% position of LED; it is assumed all LEDs are located at samepoint for
% faster simulation
% for one LED simulation located at the central of the room, use XT=0 and YT=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx=lx*5; Ny=ly*5;
% number of grid in the receiver plane
x=linspace(-lx/2,lx/2,Nx);
y=linspace(-ly/2,ly/2,Ny);
[XR,YR]=meshgrid(x,y);
D1=sqrt((XR-xilc(1,1)).^2+(YR-yilc(1,1)).^2+h_T^2);
D2=sqrt((XR-xilc(1,2)).^2+(YR-yilc(1,2)).^2+h_T^2);
D3=sqrt((XR-xilc(1,3)).^2+(YR-yilc(1,3)).^2+h_T^2);
D4=sqrt((XR-xilc(1,4)).^2+(YR-yilc(1,4)).^2+h_T^2);
D5=sqrt((XR-xilc(1,5)).^2+(YR-yilc(1,5)).^2+h_T^2);
D6=sqrt((XR-xilc(1,6)).^2+(YR-yilc(1,6)).^2+h_T^2);
D7=sqrt((XR-xilc(1,7)).^2+(YR-yilc(1,7)).^2+h_T^2);
D8=sqrt((XR-xilc(1,8)).^2+(YR-yilc(1,8)).^2+h_T^2);
% distance vector from source 1
cosphi_A1=h_T./D1;
cosphi_A2=h_T./D2;
cosphi_A3=h_T./D3;
cosphi_A4=h_T./D4;
cosphi_A5=h_T./D1;
cosphi_A6=h_T./D2;
cosphi_A7=h_T./D3;
cosphi_A8=h_T./D4;
% angle vector
receiver_angle_1=acosd(cosphi_A1);
receiver_angle_2=acosd(cosphi_A2);
receiver_angle_3=acosd(cosphi_A3);
receiver_angle_4=acosd(cosphi_A4);
receiver_angle_5=acosd(cosphi_A5);
receiver_angle_6=acosd(cosphi_A6);
receiver_angle_7=acosd(cosphi_A7);
receiver_angle_8=acosd(cosphi_A8);
% alternative methods to calculate angle, more accurate if theangle are
% negatives
% nr=[0 0 1];
% RT=[1.25 1.25]; % transmitter location
% for r=1:length(x)
% for c=1:length(y)
%
% angleA12=atan(sqrt((x(r)?1.25).^2+(y(c)?1.25).^2)./h);
% costheta(r,c)=cos(angleA12);
% end
% end
%
%%
for n = 1:10000
%% Blockage modelling using MHCP
x0 = 10; y0 = 10; %center of disc
lambdaB = (1-exp(-lambdaM*pi*R_B^2*radius))./(pi*radius*R_B^2);
pois_pointM = poissrnd(lambdaB*area); %generating poisson random number with mean lambda*area
thetaM=2*pi*(rand(1,pois_pointM)); %theta coordinates
rhoM=radius*sqrt(rand(1,pois_pointM)); %rho coordinates
[xM, yM] = pol2cart(thetaM,rhoM);
xM=xM+x0;yM=yM+y0; %Shift centre of disk to (x0,y0)
xM1= ceil(xM)./0.2; 
yM1= ceil(yM)./0.2;
xM2 = xM1+1;
yM2 = yM1+1;
%%
% Blockage distance calculation from source
d1_1 = sqrt(((xM-x0)-xilc(1,1)).^2 + ((yM-y0)-yilc(1,1)).^2);
d1_2 = sqrt(((xM-x0)-xilc(1,2)).^2 + ((yM-y0)-xilc(1,2)).^2);
d1_3 = sqrt(((xM-x0)-xilc(1,3)).^2 + ((yM-y0)-xilc(1,3)).^2);
d1_4 = sqrt(((xM-x0)-xilc(1,4)).^2 + ((yM-y0)-xilc(1,4)).^2);
d1_5 = sqrt(((xM-x0)-xilc(1,5)).^2 + ((yM-y0)-yilc(1,5)).^2);
d1_6 = sqrt(((xM-x0)-xilc(1,6)).^2 + ((yM-y0)-xilc(1,6)).^2);
d1_7 = sqrt(((xM-x0)-xilc(1,7)).^2 + ((yM-y0)-xilc(1,7)).^2);
d1_8 = sqrt(((xM-x0)-xilc(1,8)).^2 + ((yM-y0)-xilc(1,8)).^2);

r_B_1 = h_B.*d1_1./(h_T-h_B);
r_B_2 = h_B.*d1_2./(h_T-h_B);
r_B_3 = h_B.*d1_3./(h_T-h_B);
r_B_4 = h_B.*d1_4./(h_T-h_B);
r_B_5 = h_B.*d1_5./(h_T-h_B);
r_B_6 = h_B.*d1_6./(h_T-h_B);
r_B_7 = h_B.*d1_7./(h_T-h_B);
r_B_8 = h_B.*d1_8./(h_T-h_B);

d_B1 = ceil(r_B_1./0.1414);
d_B2 = ceil(r_B_2./0.1414);
d_B3 = ceil(r_B_3./0.1414);
d_B4 = ceil(r_B_4./0.1414);
d_B5 = ceil(r_B_5./0.1414);
d_B6 = ceil(r_B_6./0.1414);
d_B7 = ceil(r_B_7./0.1414);
d_B8 = ceil(r_B_8./0.1414);

d_B = [d_B1 d_B2 d_B3 d_B4 d_B5 d_B6 d_B7 d_B8];
%%
H_A1=(ml+1)*Adet.*cosphi_A1.^(ml+1)./(2*pi.*D1.^2);
H_A2=(ml+1)*Adet.*cosphi_A2.^(ml+1)./(2*pi.*D2.^2);
H_A3=(ml+1)*Adet.*cosphi_A3.^(ml+1)./(2*pi.*D3.^2);
H_A4=(ml+1)*Adet.*cosphi_A4.^(ml+1)./(2*pi.*D4.^2);
H_A5=(ml+1)*Adet.*cosphi_A5.^(ml+1)./(2*pi.*D5.^2);
H_A6=(ml+1)*Adet.*cosphi_A6.^(ml+1)./(2*pi.*D6.^2);
H_A7=(ml+1)*Adet.*cosphi_A7.^(ml+1)./(2*pi.*D7.^2);
H_A8=(ml+1)*Adet.*cosphi_A8.^(ml+1)./(2*pi.*D8.^2);

% channel DC gain for source 1,2,3 and 4
for i=1:size(xM1,2)
   
   if xM1(i) == pi*100 || xM1(i) == 1 && yM1(i) == pi*100 || yM1(i) == 1
       H_A1(xM1(i),yM1(i))=0;
       H_A2(xM1(i),yM1(i))=0;
       H_A3(xM1(i),yM1(i))=0;
       H_A4(xM1(i),yM1(i))=0;
       H_A5(xM1(i),yM1(i))=0;
       H_A6(xM1(i),yM1(i))=0;
       H_A7(xM1(i),yM1(i))=0;
       H_A8(xM1(i),yM1(i))=0;
   else
   H_A1(xM1(i),yM1(i))=0;
   H_A1(xM2(i),yM1(i))=0;
   H_A1(xM1(i),yM2(i))=0;
   H_A1(xM2(i),yM2(i))=0;
   H_A2(xM1(i),yM1(i))=0;
   H_A2(xM2(i),yM1(i))=0;
   H_A2(xM1(i),yM2(i))=0;
   H_A2(xM2(i),yM2(i))=0;
   H_A3(xM1(i),yM1(i))=0;
   H_A3(xM2(i),yM1(i))=0;
   H_A3(xM1(i),yM2(i))=0;
   H_A3(xM2(i),yM2(i))=0;
   H_A4(xM1(i),yM1(i))=0;
   H_A4(xM2(i),yM1(i))=0;
   H_A4(xM1(i),yM2(i))=0;
   H_A4(xM2(i),yM2(i))=0;
   
   H_A5(xM1(i),yM1(i))=0;
   H_A5(xM2(i),yM1(i))=0;
   H_A5(xM1(i),yM2(i))=0;
   H_A5(xM2(i),yM2(i))=0;
   H_A6(xM1(i),yM1(i))=0;
   H_A6(xM2(i),yM1(i))=0;
   H_A6(xM1(i),yM2(i))=0;
   H_A6(xM2(i),yM2(i))=0;
   H_A7(xM1(i),yM1(i))=0;
   H_A7(xM2(i),yM1(i))=0;
   H_A7(xM1(i),yM2(i))=0;
   H_A7(xM2(i),yM2(i))=0;
   H_A8(xM1(i),yM1(i))=0;
   H_A8(xM2(i),yM1(i))=0;
   H_A8(xM1(i),yM2(i))=0;
   H_A8(xM2(i),yM2(i))=0;
   end
  
end

for l = 1:length(d_B )
    for j = 1:length(d_B (1))
        for i = 1:size(xM1,2)
        if xM1(i) >= pi*25 || yM1(i)>=pi*25
            pk = 1;
        else
            xM3(j) = xM1(i)+j;
            yM3(j) = yM1(i)+j;
             
        end
       
        end
       xM4(l) = xM3;
       yM4(l) = yM3;
    end
   
end
 for g = 1:size(xM4,2)
   H_A1(xM4(g),yM4(g))=0; 
   H_A2(xM4(g),yM4(g))=0;
   H_A3(xM4(g),yM4(g))=0;
   H_A4(xM4(g),yM4(g))=0;
   H_A5(xM4(g),yM4(g))=0; 
   H_A6(xM4(g),yM4(g))=0;
   H_A7(xM4(g),yM4(g))=0;
   H_A8(xM4(g),yM4(g))=0;
 end
 
P_rec_A1=P_total.*H_A1.*Ts.*G_Con;
P_rec_A2=P_total.*H_A2.*Ts.*G_Con;
P_rec_A3=P_total.*H_A3.*Ts.*G_Con;
P_rec_A4=P_total.*H_A4.*Ts.*G_Con;
P_rec_A5=P_total.*H_A5.*Ts.*G_Con;
P_rec_A6=P_total.*H_A6.*Ts.*G_Con;
P_rec_A7=P_total.*H_A7.*Ts.*G_Con;
P_rec_A8=P_total.*H_A8.*Ts.*G_Con;
% received power from source 1;
P_rec_A1(find(abs(receiver_angle_1)>FOV))=0;
P_rec_A2(find(abs(receiver_angle_2)>FOV))=0;
P_rec_A3(find(abs(receiver_angle_3)>FOV))=0;
P_rec_A4(find(abs(receiver_angle_4)>FOV))=0;
P_rec_A5(find(abs(receiver_angle_5)>FOV))=0;
P_rec_A6(find(abs(receiver_angle_6)>FOV))=0;
P_rec_A7(find(abs(receiver_angle_7)>FOV))=0;
P_rec_A8(find(abs(receiver_angle_8)>FOV))=0;
% if the anlge of arrival is greater than FOV, no current is generated at % the photodiode.
P_rec_total=P_rec_A1+P_rec_A2+P_rec_A3+P_rec_A4+P_rec_A5+P_rec_A6+P_rec_A7+P_rec_A8;
P_rec_dBm=10*log10(P_rec_total);
end
%%
% % 5. For SNR %
% 
% Bn = I2 * Rb; % Noise-bandwidth (Sec^-1)%
% Pamb = Iamb / R; % Ambient light power (W) %
% % Shot-noise variance ( Ampere^2 )%
% omega_shot = 2 * q * R * (P_rec_total + Pamb) * Bn; 
% % Amplifier noise variance ( Ampere^2 )%
% omega_amplifier = Iamf^2 * Ba; 
% %Thermal noise variance
% omega_thermal = (8*pi*295*112E-8*1E-3*.562*1E6*1.38E-23)+((16*pi^2*1.38E-23*295*1.5*(112E-8)^2*1E-8*.56281E12)/.03);
% % Total noise variance ( Ampere^2 )%
% omega_total = omega_amplifier + omega_shot+omega_thermal; 
% % SNR %
% SNR = (( R * P_rec_total )^2)./ omega_total;
% SNRdB = 10*log(SNR);

surfc(x,y,P_rec_dBm);
xlabel('width (m)');
ylabel('Length (m)');
zlabel('Received Power (dBm)'); 
grid on