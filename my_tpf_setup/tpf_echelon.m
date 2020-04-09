function [k_vector,e_h,x0,LN,my_input]=tpf_echelon(my_input,LN,x0,e_f,alpha_tpf_aim)

        %--------------------------------------------------------------------------
        %setup optical imaging system
        %-------------------------------------------------------------------------

        f1=my_input.f1;
        f2=f1*sqrt(2/LN.n_ir_0/LN.n_g);
        LN.os_magnify=(-f2/my_input.f1);
        LN.image_d=my_input.image_distance*my_input.sigma_in*tan(alpha_tpf_aim)*abs(LN.os_magnify);
        LN.bs_crystal_thz=my_input.sigma_in*abs(LN.os_magnify)/cos(alpha_tpf_aim);
        my_input.sigma_crystal=my_input.sigma_in*abs(LN.os_magnify);
        
        echelon_w =150e-6;   %150e-6                                       %beamlet size width
        echelon_h = echelon_w*(tan(alpha_tpf_aim)*LN.n_g)*abs(LN.os_magnify)/2;  %depth                                                                     %smooth level of the step, large the a, less smooth it is
             
        %------------------------------------------------------------------------
        %construct echelon 
        %------------------------------------------------------------------------

        x=x0.*cos(alpha_tpf_aim);
        dxx=x(2)-x(1);
        kxx = 2*pi*(-(length(x)-1)/(2*length(x).*dxx):1/(length(x)*dxx):(length(x)-1)/(2*length(x)*dxx));
        dkxx=kxx(2)-kxx(1);
 
        alpha_rotate=alpha_tpf_aim;
        %------------------------------------------------------------------------
        %construct diffraction for each section
        %------------------------------------------------------------------------
      
        k_w0=LN.omega_b/my_c.c;
                                                                            
        echelon_w_final=(f2/f1)*echelon_w;
        center_x=echelon_w_final*floor((x+echelon_w_final/2)./echelon_w_final);
        e_h=sqrt(f1/f2).*e_f.*exp(-((x-center_x).^2.*f1^2./(f2*echelon_w/2.4).^2).^8).*exp(-x.^2.*f1^2./(2*my_input.sigma_in.^2.*f2^2));
        echelon_distance=-2*echelon_h*center_x./echelon_w_final;  
        phase_echelon=((LN.omega_b-LN.omega_0)./my_c.c).*echelon_distance;  
        e_h=e_h.*exp(1i.*phase_echelon);
              
        e_kx=fftshift(ifft(ifftshift(e_h,2),[],2),2)./(length(kxx)*dkxx/(2*pi));
        d_start=-LN.image_d;
        e_kx=e_kx.*exp(1i*kxx.^2.*d_start./(2.*k_w0));
        e_h=fftshift(fft(ifftshift(e_kx,2),[],2),2).*(length(kxx)*dkxx/(2*pi));

        %----------------------------------------------------------------------
        %k vector in rotated coordinates with respect to propagation direction
        %----------------------------------------------------------------------
        %k vector for ir. kz_ir is k vector in propagation direction, kx is
        %included in tpf angle term
        k_vector=struct;
        dx=x0(2)-x0(1);
        k_vector.kx = 2*pi*(-(length(x0)-1)/(2*length(x0).*dx):1/(length(x0)*dx):(length(x0)-1)/(2*length(x0)*dx));
        k_vector.dkx=k_vector.kx(2)-k_vector.kx(1);

        
        %adjust the initial kx distribution, the kx_tpf draw the tpf back
        %to the thz coordinate and this kx_tpf successfully predict the new
        %generated kx such that the kx domain don't need to be large
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

        LN.ir_index=ones(size(LN.omega_b));

        %k vector for thz
        theta_thz=-(alpha_rotate-alpha_tpf_aim);
        LN.theta_thz=theta_thz;
        LN.alpha_rotate=alpha_rotate;

        k_vector.kz_thz=cos(theta_thz).*LN.omega_t.*LN.n_thz/my_c.c; 
        k_vector.kx_thz=-sin(theta_thz).*LN.omega_t.*LN.n_thz/my_c.c; 
        

        %---------------------------------------------------
        %remove thz resonance peak numerical error
        %--------------------------------------------
        index_thz_prop=abs(abs(k_vector.kx))>(LN.omega_t.*LN.n_thz/my_c.c);
        index_thz_resonant=repmat((LN.n_thz<1),1,length(x0));
        LN.index_thz=index_thz_prop|index_thz_resonant;


end
