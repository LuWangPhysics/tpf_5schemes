function [k_vector,e_h,x0,LN,my_input]=tpf_silicon(my_input,LN,x0,e_f,alpha_tpf_aim)

        alpha_rotate=LN.alpha_tpf_aim;
        %--------------------------------------------------------------------------
        %same as echelon, just separated beamlets without diffraction
        %-------------------------------------------------------------------------
        
        echelon_w =150e-6;
        n_silica=silica_n(LN.f_b);
        n_silica_0=1.45;
        echelon_h  = echelon_w*(tan(alpha_tpf_aim).*LN.n_g)/(n_silica_0-1);       
                                                                      
        %------------------------------------------------------------------------
        %construct diffraction for each section
        %------------------------------------------------------------------------
       
        x=x0.*cos(LN.alpha_tpf_aim);
        dxx=x(2)-x(1);
        kxx = 2*pi*(-(length(x)-1)/(2*length(x).*dxx):1/(length(x)*dxx):(length(x)-1)/(2*length(x)*dxx));
        dkxx=kxx(2)-kxx(1);
        
                                                   
        center_x=echelon_w.*floor((x+echelon_w/2)./echelon_w);
        my_input.sigma_crystal=my_input.sigma_in;
        
        [n_center,p_center]=hist(floor((x+echelon_w/2)./echelon_w),unique(floor((x+echelon_w/2)./echelon_w)));

        
       e_h=0;
        edge_side=fix((x(end)-4*my_input.sigma_in)/echelon_w);
        for nn=edge_side:(length(p_center)-edge_side)
            L_air=my_input.sigma_in*tan(alpha_tpf_aim)*LN.n_g/(n_silica_0-1)-p_center(nn)*echelon_h;
            L_silica=my_input.sigma_in*tan(alpha_tpf_aim)*LN.n_g/(n_silica_0-1)+p_center(nn)*echelon_h;
            
            e_temp=e_f.*exp(-((x-p_center(nn)*echelon_w).^2./(echelon_w/2.4).^2).^8).*exp(-x.^2./(2*my_input.sigma_in.^2));
            e_kx=fftshift(ifft(ifftshift(e_temp,2),[],2),2)./(length(kxx)*dkxx/(2*pi));
            
            e_kx=e_kx.*exp(1i.*L_air.*kxx.^2./(2.*LN.omega_b./my_c.c)).*exp(1i.*L_silica.*kxx.^2./(2.*n_silica.*LN.omega_b./my_c.c));
            e_temp=fftshift(fft(ifftshift(e_kx,2),[],2),2).*(length(kxx)*dkxx/(2*pi));
             
            
            e_h=e_h+e_temp;
        end
        
        
      

        phase_echelon=((LN.omega_b-LN.omega_0)./my_c.c).*center_x.*tan(alpha_tpf_aim).*LN.n_g;   
        e_h=e_h.*exp(-1i.*phase_echelon);


        %----------------------------------------------------------------------
        %k vector in rotated coordinates with respect to propagation direction
        %----------------------------------------------------------------------
        %k vector for ir. kz_ir is k vector in propagation direction, kx is
        %included in tpf angle term
        k_vector=struct;
        dx=x0(2)-x0(1);
        k_vector.kx = 2*pi*(-(length(x0)-1)/(2*length(x0).*dx):1/(length(x0)*dx):(length(x0)-1)/(2*length(x0)*dx));
        k_vector.dkx=k_vector.kx(2)-k_vector.kx(1);

        
        %adjust the initial kx distribution
        kx_tpf=LN.delta_omega.*tan(alpha_tpf_aim)*LN.n_g./my_c.c;
        e_h=e_h.*exp(1i.*kx_tpf.*x);
        
        k_ir_z=sqrt((LN.omega_b.*LN.n_ir_b./my_c.c).^2-kx_tpf.^2);    
        k_vector.kx_ir=cos(alpha_rotate).*kx_tpf-k_ir_z.*sin(alpha_rotate);
        k_vector.kz_ir=sin(alpha_rotate).*kx_tpf+cos(alpha_rotate).*k_ir_z;
        
        
        
        %IN ir frame
        LN.F_peak_irframe=max(0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,LN.n_ir_b),1).*LN.df);
        %in THz frame
        F_pump_x=0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,k_vector.kz_ir./(LN.omega_b./my_c.c)),1).*LN.df;

        LN.F_pump=sum(F_pump_x)*(x0(2)-x0(1));
        

        %----------------------------------------------------------------
        %remove the part where kx>kz 
        %---------------------------------------------------------------
        %kz_1(abs(kx_1)>kz_ir)=-1i.* kz_1(abs(kx_1)>kz_ir);
        LN.ir_index=ones(size(LN.omega_b));
       % LN.ir_index(abs(kx_1)>kz_ir)=0;
        %k vector for thz
        theta_thz=-(alpha_rotate-alpha_tpf_aim);
        LN.theta_thz=theta_thz;
        LN.alpha_rotate=alpha_rotate;

        k_vector.kz_thz=cos(theta_thz).*LN.omega_t.*LN.n_thz/my_c.c; 
        k_vector.kx_thz=-sin(theta_thz).*LN.omega_t.*LN.n_thz/my_c.c; 
        

        LN.theta_shg=LN.alpha_rotate-LN.shg_pm;
        k_vector.kx_shg=-LN.k_shg.*sin( LN.theta_shg);
        k_vector.kz_shg=LN.k_shg.*cos(LN.theta_shg);
        %---------------------------------------------------
        %remove thz resonance peak numerical error
        %--------------------------------------------
        index_thz_prop=abs(abs(k_vector.kx))>(LN.omega_t.*LN.n_thz/my_c.c);
        index_thz_resonant=repmat((LN.n_thz<1),1,length(x0));
        LN.index_thz=index_thz_prop|index_thz_resonant;
        %just for saving stirng
        LN.os_magnify=1;
        LN.alpha_tpf_aim=alpha_tpf_aim;
        LN.bs_crystal_thz=my_input.sigma_in*abs(LN.os_magnify)/cos(LN.alpha_tpf_aim);
        LN.image_d=0;
        LN.theta_i=0;
      
end