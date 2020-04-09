function [k_vector,e_h,x0,LN]=plain_input(my_input,LN,x0,e_f,alpha_tpf_aim,n_scan)
frame={'IR_frame','THz_frame'};
% 0 in IR frame 1 in thz frame
case_n=0; 

LN.frame=frame{case_n+1};
alpha_rotate=case_n*alpha_tpf_aim;

x=x0.*cos(alpha_rotate);


 
e_h=exp(-x.^2./(2*my_input.sigma_in.^2)).*e_f;
e_t=fftshift(ifft(ifftshift(e_h,1),[],1),1).*length(LN.omega_b).*LN.df;
slope=angle_test(x,e_t,LN.t0,my_c.c./LN.n_g)
%---------------------------------
%change the coordinate to the rotated coordinate.
%-------------------------------------------------------------------------






        %----------------------------------------------------------------------
        %k vector in rotated coordinates with respect to propagation direction
        %----------------------------------------------------------------------
        %k vector for ir. kz_ir is k vector in propagation direction, kx is
        %included in tpf angle term
        k_vector=struct;
        dx=x0(2)-x0(1);
        k_vector.kx = 2*pi*(-(length(x0)-1)/(2*length(x0).*dx):1/(length(x0)*dx):(length(x0)-1)/(2*length(x0)*dx));
        k_vector.dkx=k_vector.kx(2)-k_vector.kx(1);
        
        k_ir_z=(LN.omega_b.*LN.n_ir_b./my_c.c);    
        k_vector.kx_ir=-k_ir_z.*sin(alpha_rotate);
        k_vector.kz_ir=cos(alpha_rotate).*k_ir_z;
       

        %IN ir FRAME
        LN.F_peak_irframe=max(0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,LN.n_ir_b),1).*LN.df);
        %e_h=sqrt(707).*e_h./sqrt(LN.F_peak_irframe);
        %energy=(x0(2)-x0(1))*cos(alpha_rotate)*sum(0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,LN.n_ir_b),1).*LN.df)
        %cos theta roate is the projection of the energy flux to
        %  %in thz frame
        F_pump_x=0.5*my_c.c*my_c.eps*sum(bsxfun(@times,abs(e_h).^2,k_vector.kz_ir./(LN.omega_b./my_c.c)),1).*LN.df;

        LN.F_pump=sum(F_pump_x)*(x0(2)-x0(1));
        
        
        
        e_t=fftshift(ifft(ifftshift(e_h.*exp(-1i.*k_vector.kx_ir.*x0),1),[],1),1).*length(LN.omega_b).*LN.df;
        slope=angle_test(x0,e_t,LN.t0,my_c.c./LN.n_g)

        %----------------------------------------------------------------
        %remove the part where kx>kz 
        %---------------------------------------------------------------
        %kz_1(abs(kx_1)>kz_ir)=-1i.* kz_1(abs(kx_1)>kz_ir);
        LN.ir_index=ones(size(LN.omega_b));
       % LN.ir_index(abs(kx_1)>kz_ir)=0;
        %k vector for thz
        theta_thz=0*(alpha_tpf_aim-alpha_rotate);
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
        LN.image_d=my_input.sigma_in*tan(LN.alpha_crystal)*LN.os_magnify;
        LN.theta_i=0;
        %num_tpf_LN=atan(n_air*cos(theta_1).*(tan(theta_1)+tan(slope*pi/180))./(LN.n_ir_0.*cos(theta_2))-tan(theta_2));
end