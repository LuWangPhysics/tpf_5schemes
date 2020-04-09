function [k_vector,e_h,x0,LN,my_input]=tpf_grating(my_input,LN,x0,e_f,alpha_tpf_aim)

        %-------------------------------------------------------------------------
        %define the grating and imaging system
        %-------------------------------------------------------------------------

      
        %p=my_input.lambda_0*sqrt(1+2*LN.n_ir_0/LN.n_g/tan(alpha_tpf_aim)^2);
        p           = my_input.p;                                                                       

        g_order=1;                                                         %grating order
        coef_a=(LN.n_ir_0.^2*LN.n_g*p/my_input.lambda_0/2)*sqrt(my_input.lambda_0^2/(LN.n_g^2*p^2*tan(alpha_tpf_aim)^4)+4/(LN.n_ir_0^2))-LN.n_ir_0^2/2/tan(alpha_tpf_aim)^2;
        theta_i=asin(my_input.lambda_0*(1-coef_a/LN.n_ir_0/LN.n_g)/p);
        LN.theta_i=theta_i;        
        theta_o = asin(g_order*my_input.lambda_0/p - sin(LN.theta_i));  
        my_input.f2=(my_input.f1*my_input.lambda_0/(p*cos(theta_o)*tan(alpha_tpf_aim)*LN.n_g));     
        %-------------------------------------------------------------------------
 
 
 
        %grating bandwidth reduce the spread in kx of broadened ir spectrum
        grating_band_width=exp(-(LN.delta_omega.^2.*my_input.tau_fwhm.^2/(4*pi^2)).^4);
        F0_o = -2*pi*my_c.c*LN.delta_omega./(p*(LN.omega_0).^2.*cos(theta_o));
        F1_o = F0_o.*(- LN.delta_omega./LN.omega_0+0.5.*tan(theta_o).*F0_o);
        F_o =(F0_o+grating_band_width.*F1_o);
            
 
        LN.os_magnify=cos(theta_o)*my_input.f2/(cos(LN.theta_i)*my_input.f1);
        my_input.sigma_crystal=my_input.sigma_in*abs(LN.os_magnify);
        [lens_M,M_optical,c_23,LN.image_d]=optical_system(my_input,LN);  
        m23=c_23.*F_o;
        LN.bs_crystal_thz=my_input.sigma_in*abs(LN.os_magnify)/cos(LN.alpha_tpf_aim);

        %-------------------------------------------------------------------------
        %change the coordinate to the rotated coordinate.
        %-------------------------------------------------------------------------

        alpha_rotate=alpha_tpf_aim;


        e_h=zeros(length(LN.omega_b),length(x0));

        for k = 1:length(LN.omega_b)   

                    % matrix of grating
                    m_g_o = [cos(theta_o)./cos(LN.theta_i) 0 0;...
                            0 cos(LN.theta_i)./cos(theta_o) F_o(k);...
                            0 0 1];
                    m_out_o =M_optical*m_g_o;

                    e_h(k,:)  = field_homo(alpha_rotate,LN.n_g,m_out_o,my_input.sigma_in,x0,LN.omega_b(k),e_f(k),LN.omega_0);


        end


        %----------------------------------------------------------------------
        %k vector in rotated coordinates with respect to propagation direction
        %----------------------------------------------------------------------
        %k vector for ir. kz_ir is k vector in propagation direction, kx is
        %included in tpf angle term
        k_vector=struct;
        
        kx_tpf=LN.omega_b.*m23./my_c.c;

        k_ir_z=sqrt((LN.omega_b.*LN.n_ir_b./my_c.c).^2-kx_tpf.^2);    
        k_vector.kx_ir=cos(alpha_rotate).*kx_tpf-k_ir_z.*sin(alpha_rotate);
        k_vector.kz_ir=sin(alpha_rotate).*kx_tpf+cos(alpha_rotate).*k_ir_z;
        
       
        LN.F_peak_irframe=max(0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,LN.n_ir_b),1).*LN.df);

 
        %cos theta roate is the projection of the energy flux to propagation direction.
        F_pump_x=0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,k_vector.kz_ir./(LN.omega_b./my_c.c)),1).*LN.df;

        LN.F_pump=sum(F_pump_x)*(x0(2)-x0(1));

        %----------------------------------------------------------------
        %remove the part where kx>kz 
        %---------------------------------------------------------------
      
       ir_index=k_vector.kz_ir<(cos(alpha_rotate)*LN.omega_0.*LN.n_g./my_c.c/20);
        LN.ir_index=repmat(ir_index,length(x0),1);
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
        index_thz_resonant=repmat((LN.n_thz<1),1,length(x0))|repmat((LN.f_t>10e12),1,length(x0));
        LN.index_thz=index_thz_prop|index_thz_resonant;
       
end
