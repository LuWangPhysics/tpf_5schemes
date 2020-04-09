function [k_vector,e_h,x0,LN,my_input]=tpf_spatial_temporal(my_input,LN,x0,e_f,alpha_tpf_aim)
        %-------------------------------------------------------------------------
        %define the spatial and temporal chirp
        %-------------------------------------------------------------------------

        tau = my_input.tau_fwhm/(sqrt(2*log(2)));  
        tan_tpf_air=(tan(alpha_tpf_aim)*LN.n_g);
        my_input.sigma_crystal=my_input.sigma_in;
        phi_2=1.5e-25;
        v=tan_tpf_air/my_c.c/phi_2;
        %total bandwidth FWHM
        bw=my_input.sigma_in.*sqrt(2*log(2))*sqrt(v^2*tau^2/2+1/my_input.sigma_in.^2)/tau/pi;
        
        
        %-------------------------------------------------------------------------
        %change the coordinate to the rotated coordinate.
        %-------------------------------------------------------------------------

        alpha_rotate=LN.alpha_tpf_aim;
        e0 = sqrt(2*sqrt(2)*my_input.energy/(sqrt(pi)*pi*my_input.sigma_in*my_input.sigma_homo*my_c.c*my_c.eps*LN.n_ir_0*tau*gamma(2)));
  
        %x in the ir frame
        x=x0.*cos(alpha_tpf_aim);   
        scale_x=(1/my_input.sigma_in^2+v^2*tau^2/2)^(1/4)*sqrt(my_input.sigma_in);
        scale_w=sqrt(1-tau^4*v^2/2/(1/my_input.sigma_in^2+v^2*tau^2/2));
        e_coe=tau*e0*sqrt(pi)*scale_x*scale_w;
        e_h = e_coe.*exp(-x.^2./(2*my_input.sigma_in.^2)).*exp(-tau^2.*(LN.delta_omega-v.*x).^2./4).*exp(-0.5i.*phi_2.*LN.delta_omega.^2);

        %----------------------------------------------------------------------
        %k vector in rotated coordinates with respect to propagation direction
        %----------------------------------------------------------------------
        k_vector=struct;
        
        k_ir_z=(LN.omega_b.*LN.n_ir_b./my_c.c);    
        k_vector.kx_ir=-k_ir_z.*sin(alpha_rotate);
        k_vector.kz_ir=cos(alpha_rotate).*k_ir_z;
        

        F_pump_ir=0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,LN.n_ir_b),1).*LN.df;
        %make sure when beam size change the fluence stays the same
        e_h=sqrt(my_input.energy/(pi*my_input.sigma_homo*my_input.sigma_crystal)).*e_h./sqrt(max(F_pump_ir));
        LN.F_peak_irframe=max(0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,LN.n_ir_b),1).*LN.df);
        F_pump_x=0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,k_vector.kz_ir./(LN.omega_b./my_c.c)),1).*LN.df;
        %cos theta roate is the projection of the energy flux to propagation direction.
        LN.F_pump=sum(F_pump_x)*(x0(2)-x0(1));
        

        %----------------------------------------------------------------
        %remove the part where kx>kz 
        %---------------------------------------------------------------
     
        LN.ir_index=logical(zeros(size(e_h)));
        %k vector for thz
        theta_thz=(alpha_rotate-alpha_tpf_aim);
        LN.theta_thz=theta_thz;
        LN.alpha_rotate=alpha_rotate;

        k_vector.kz_thz=cos(theta_thz).*LN.omega_t.*LN.n_thz/my_c.c; 
        k_vector.kx_thz=-sin(theta_thz).*LN.omega_t.*LN.n_thz/my_c.c; 

        dx=x0(2)-x0(1);
        k_vector.kx = 2*pi*(-(length(x0)-1)/(2*length(x0).*dx):1/(length(x0)*dx):(length(x0)-1)/(2*length(x0)*dx));
        k_vector.dkx=k_vector.kx(2)-k_vector.kx(1);
        %---------------------------------------------------
        %remove thz resonance peak numerical error
        %--------------------------------------------
        index_thz_prop=abs(abs(k_vector.kx))>(LN.omega_t.*LN.n_thz/my_c.c);
        index_thz_resonant=repmat((LN.n_thz<1),1,length(x0));
        LN.index_thz=index_thz_prop|index_thz_resonant;
        LN.os_magnify=1;
        LN.image_d=0;
        LN.theta_i=0;
        LN.bs_crystal_thz=my_input.sigma_in*abs(LN.os_magnify)/cos(LN.alpha_tpf_aim);
 
end
