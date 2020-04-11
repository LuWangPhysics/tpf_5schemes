function [h_r]= raman(omega_b,omega_0)

c=3e8;
f_tt=((omega_b-omega_0)./(2*pi))/c/100;                                                  
eps_infinity=0;%3.3;
h_r=ones(length(f_tt),1).*eps_infinity;
f_1=[130,248,274,307,628,692];                                             %wave number   
gamma_1=[75,21,14,25,34,49];
s_1=[4,16,1,0.16,2.55,0.13];
for k=1:length(f_1)
h_r=h_r+s_1(k)*f_1(k)^2./(f_1(k)^2-f_tt.^2+1i*gamma_1(k).*f_tt);
end




df=(omega_b(2)-omega_b(1))/2/pi;

h_t=fftshift(ifft(ifftshift(h_r))).*(length(omega_b)*df);
dt=1/(length(omega_b)*df);
%h_r(center)=1 automatically satisfied!
h_r = h_r./sum(h_t*dt);%-1;    


end