function [n_thz,alpha_thz]= cLN_thz(f_in,~)
%only valid for room temperature!! 300k!!
%--------------------------------------------------------
%calculate phonon resonance to be the same as stimulated raman
%--------------------------------------------------------
c=3e8;
f_tt=(f_in)/c/100;                                                  
eps_infinity=3.3;
chi_1=ones(length(f_tt),1).*eps_infinity;
f_1=[130,248,274,307,628,692];                                             %wave number   
gamma_1=[75,21,14,25,34,49];
s_1=[4,16,1,0.16,2.55,0.13];
for k=1:length(f_1)
chi_1=chi_1+s_1(k)*f_1(k)^2./(f_1(k)^2-f_tt.^2+1i*gamma_1(k).*f_tt);
end
alpha_resonance=-2*pi*f_in.*imag(sqrt(chi_1))./c;


        [n_ravi,alpha_thz1]= cLN_thz_old(f_in,80);
        %alpha 0.3THz around 7/cm reference
        %https://aip.scitation.org/doi/pdf/10.1063/1.1929859?class=pdf
        %https://link.springer.com/content/pdf/10.1007%2Fs10762-015-0165-5.pdf
        f_choice=3.24e12;%3.24e12;
        alpha_thz=5.*alpha_thz1;
        alpha_thz(f_in>f_choice)=alpha_resonance(f_in>f_choice); 
        n_thz=real(sqrt(chi_1));

   
end