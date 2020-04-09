function [k_vector,e_h,x0,LN,my_input]=tpf_NLES(my_input,LN,x0,e_f,alpha_tpf_aim)

%-------------------------------------------------------------------------
%define the grating and imaging system
%-------------------------------------------------------------------------
                                        
        p           = my_input.p;                                                                                                                                                     
        f1          = my_input.f1;
  
        coe_b=LN.n_ir_0.*my_input.lambda_0/p/tan(alpha_tpf_aim)^2/(LN.n_ir_0-1);    
        LN.theta_o=asin((-coe_b+sqrt(coe_b^2+4))/2); 
        f2=f1*my_input.lambda_0/(tan(alpha_tpf_aim)*p*cos(LN.theta_o));    
        LN.theta_i=asin(my_input.lambda_0/p - sin(LN.theta_o));            %grating incidence angle  
        
        echelon_w =150e-6;                                                 %beamlet size
        LN.echelon_w=echelon_w;
        echelon_h  = echelon_w*(tan(alpha_tpf_aim));       
%--------------------------------------------------------------------------

        grating_band_width=exp(-(LN.delta_omega.^2.*my_input.tau_fwhm.^2/(4*pi^2)/4).^4);
        F0_o = -2*pi*my_c.c*LN.delta_omega./(p*(LN.omega_0).^2.*cos(LN.theta_o));
        F1_o = F0_o.*(- LN.delta_omega./LN.omega_0+0.5.*tan(LN.theta_o).*F0_o);
        F_o = (F0_o+grating_band_width.*F1_o);
       
        LN.os_magnify=cos(LN.theta_o)/cos(LN.theta_i)*(f2/f1);
        my_input.sigma_crystal=my_input.sigma_in*abs(LN.os_magnify);
%--------------------------------------------------------------------------
       %optical system
%--------------------------------------------------------------------------
        %for nles there is no entering the prism process. it directly start
        %to generate thz
        % s should be zero
        s=0;

         s1          = f1;                                                 %Grating to lens distance
         s2          = f1+f2;                                              %distance between two lens
         s3          = f2-s/LN.n_ir_0;                                     %Lens to crystal distance

        m_s1 = [1 s1 0;0 1 0;0 0 1];                                       %free propagation for s1 in vancumm
        m_l1  = [1 0 0;-1./f1 1 0;0 0 1];                                  %focusing
        m_s2 = [1 s2 0;0 1 0;0 0 1];                                       %free propagation for s2 in vancumm
        m_l2  = [1 0 0;-1./f2 1 0;0 0 1];                                  %focusing	
        m_s3 = [1 s3 0;0 1 0;0 0 1];                                       %free propagation for s2 in vancumm
        
        M_optical=m_s3*m_l2*m_s2*m_l1*m_s1;             
        LN.image_d=0;
        LN.bs_crystal_thz=0;



%-------------------------------------------------------------------------
%change the coordinate to the rotated coordinate.
%-------------------------------------------------------------------------
        %same as theta_tpf_in_air/n
        alpha_rotate=alpha_tpf_aim;
        LN.alpha_tpf_aim=alpha_tpf_aim;

        e_h=zeros(length(LN.omega_b),length(x0));

        for k = 1:length(LN.omega_b)   

                    % matrix of grating
                    m_g_o = [cos(LN.theta_o)./cos(LN.theta_i) 0 0;...
                            0 cos(LN.theta_i)./cos(LN.theta_o) F_o(k);...
                            0 0 1];
                    m_out_o =M_optical*m_g_o;

                    e_h(k,:)  = field_homo(alpha_rotate,LN.n_g,m_out_o,my_input.sigma_in,x0,LN.omega_b(k),e_f(k),LN.omega_0);


        end
         m23=(1-s2/f2+s1*(-(1/f2) - (1 - s2/f2)/f1)).*F_o;
        %--------------------------------------------------------------------------
        %setup zigzag LN x0 is tpf frame, x is IR frame
        %-------------------------------------------------------------------------
        x=x0.*cos(alpha_rotate);                                                                   
        center_x=echelon_w.*floor((x+echelon_w/2)./echelon_w);
        echelon_distance=echelon_h.*center_x./echelon_w;
        %create a spatial intensity modulation 
        x_module=exp(-((x-center_x).^2./(echelon_w/2.3).^2).^8);
        phase_echelon=(LN.delta_omega./my_c.c).*echelon_distance; 
        %opposite from the echelon distance phase
        %cancel out the phase brought in by the coordinate transform LN part
        e_h=x_module.*e_h.*exp(1i.*phase_echelon*(1-LN.n_g));
        
        %----------------------------------------------------------------------
        %in the thz frame
        %----------------------------------------------------------------------
        k_vector=struct;
        kx_tpf=LN.omega_b.*m23./my_c.c;

        k_ir_z=sqrt((LN.omega_b.*LN.n_ir_b./my_c.c).^2-kx_tpf.^2);    
        k_vector.kx_ir=cos(alpha_rotate).*kx_tpf-k_ir_z.*sin(alpha_rotate);
        k_vector.kz_ir=sin(alpha_rotate).*kx_tpf+cos(alpha_rotate).*k_ir_z;
               
        LN.F_peak_irframe=max(0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,LN.n_ir_b),1).*LN.df);
        F_pump_x=0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,k_vector.kz_ir./(LN.omega_b./my_c.c)),1).*LN.df;
        LN.F_pump=sum(F_pump_x)*(x0(2)-x0(1));
       
        %----------------------------------------------------------------
        %remove the part where kx>kz 
        %---------------------------------------------------------------
      
        ir_index=LN.f_b<262e12;
        LN.ir_index=repmat(ir_index,length(x0),1);
        %k vector for thz
        theta_thz=-(alpha_rotate-alpha_tpf_aim);
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
       

end
