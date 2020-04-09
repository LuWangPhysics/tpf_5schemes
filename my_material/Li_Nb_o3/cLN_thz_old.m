function  [n_thz,alpha_thz]= cLN_thz_old(f_in,T)
%**********************************************
%here is For 80k only!!!
%***********************************************

f_t=f_in;
f = (f_in)/3e10; %Converting to /cm units

%For 6% CLN
a_c = 2.706*exp(-0.0008344*T) + 1.993*exp(0.001173*T);
b_c = -1.613e-13*T^3 + 9.679e-11*T^2 + 2.254e-9*T + 3.197e-5;
c_c = 3e-10;

%Terahertz Refractive Index
n_thz1 = 1.035*(a_c+b_c*f.^2 + c_c*f.^4);
n_thz1(f_t>100*1e12)=n_thz1(sum(f_t<100*1e12));
%(Valid between 10K and 300K, and 100/cm)

f_thz = (0.4:0.1:1.9)*1e12/3e10;
f0 = (0:.01:200)*1e12/3e10;


alpha_thz_300           = [0.2*6 0.3*7.5 0.6*8 0.9*10 1.1*13 1.2*15 1.4*16 1.5*18 1.6*20 1.7*23 1.8*25 1.9*29 1.9*33 1.9*36 2*40 2*45]*1e2;
alpha_thz_300_s         = smooth(alpha_thz_300,5);
alpha_thz_300_interp1   = interp1(f_thz,alpha_thz_300_s,f0,'linear','extrap');
alpha_thz_77            = [0.5*1.2 0.5*1.4 0.6*1.5 0.8*1.6 0.8*1.8 0.9*2 0.9*2.1 1.1*2.25 1.2*2.6 1.3*2.9 1.4*3.25 1.5*3.6 1.5*4.1 1.6*4.5 1.6*5 1.7*5.75]*1e2;
alpha_thz_77_s          = smooth(alpha_thz_77,5);
alpha_thz_77_interp1    = interp1(f_thz,alpha_thz_77_s,f0,'linear','extrap');

f_10 = [0.5:0.1:1,1.4,1.6,1.7,1.8,1.9]*1e12/3e10;
alpha_thz_10            = [0.25:.02:0.35,0.8*0.5, 1.2*0.6, 1.2*0.7,1.2*0.8,1.1*1.5]*1e2;
alpha_thz_10_s          = smooth(alpha_thz_10,5);
alpha_thz_10_interp1    = interp1(f_10,alpha_thz_10_s,f0,'linear','extrap');

T1 =[10;100;300];
%in unit /m
alpha_thz1 =10*interp2(f0,T1,[alpha_thz_10_interp1;alpha_thz_77_interp1;alpha_thz_300_interp1],f,T,'linear')./3;
alpha_thz1(f_t>100*1e12)=alpha_thz1(sum(f_t<100*1e12));

alpha_thz=alpha_thz1;
n_thz=n_thz1;
 end