function [n_ir_0,n_ir_b,n_thz_c,n_thz,alpha_thz,n_g,chi_2_0,n2_0]=LiNb03(omega_0,omega_b,omega_t,T,f_c_t)
       
        n_ir_0    =  cLN_ir(omega_0/(2*pi),T);                    % Refractive index at center wavelength 
        n_ir_b    = cLN_ir(omega_b/(2*pi),T); 
        
        [n_thz,alpha] =  cLN_thz(omega_t/(2*pi),T);
        alpha_thz=-alpha';
        [n_thz_c,~]   =  cLN_thz(f_c_t,T);
        
        [n_g,~]=Gro_v(omega_0,omega_b,n_ir_b);   
        
        
         chi_2_0=  2*168e-12;                  % 2nd order nonlinearity (d_eff= 168 pm/V, chi_2 = 2*d_eff)                                     %2nd order nonlinearity 
         n2_0  =  1.25e-19;                    %m^2/W
end