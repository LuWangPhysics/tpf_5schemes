function [n_gaussian,LN,e_f]=input_beam(my_input,deltaw,material_name)

        omega_0   =2*pi*my_c.c/my_input.lambda_0;
                                                                
        [LN]=f_t_discrit2D(omega_0,my_input,deltaw,material_name);          
                                                                          
        
        n_gaussian=1;                                                      %spatial profile
        %----------------------------------------------------------------
        %input peak field strength in time domain,gaussian input
        %----------------------------------------------------------------
        tau = my_input.tau_fwhm/(sqrt(2*log(2)));  
        e0 = sqrt(2^(1/n_gaussian)*sqrt(2)*my_input.energy/(sqrt(pi)*pi*my_input.sigma_in*my_input.sigma_homo*my_c.c*my_c.eps*LN.n_ir_0*tau*gamma((n_gaussian+1)/n_gaussian)));
        phi_c_b=exp(-1i.*(0.5.*my_input.GDD.*LN.delta_omega.^2+(1/6).*my_input.TOD.*LN.delta_omega.^3));
        e_f  = tau*e0*sqrt(pi).*exp(-LN.delta_omega.^2*tau^2/4).*phi_c_b;


   
end