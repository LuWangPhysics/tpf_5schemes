classdef initialize_OP

    methods
        
        function OP=initialize_OP
       
        end
        
        function phase_cor=phase_corr_prism(OP,x0,z0,dz,LN)
        %---------------------------------------------------------------
        %compensate for the phase of LN prism 
        %---------------------------------------------------------------
        x_edge=-LN.bs_crystal_thz+z0/tan(LN.alpha_tpf_aim);
        x_array=zeros(size(x0));
        x_array(x0>x_edge)=dz;
        phase_l=LN.prism_comp.*sin(LN.alpha_tpf_aim).*x_array.*(1+tanh(5e4*(x0+x_edge)))/2;
        phase_cor=exp(1i.*phase_l.*LN.omega_b.*(LN.n_ir_b-LN.n_g)./my_c.c);
         
        end
        
        
        function phase=propagation_phase(OP,kx0,x_array,kz0,z0)
            phase=exp(-1i.*kz0.*z0).*exp(-1i.*kx0.*x_array);
        end
        
        function phase=diffraction_phase(OP,kx0,kx,kz0,z0)
        phase=exp(1i.*(kx.^2-2*kx0.*kx).*z0./(2*kz0));
        end


        function e_x=fft2_kxtox(OP,phi_ir,k_vector)
        e_x=fftshift(ifft(ifftshift(phi_ir,2),[],2),2).*(length(k_vector.kx)*k_vector.dkx/(2*pi));
        end

        function e_kx=fft2_xtokx(OP,phi_ir,k_vector)
        e_kx=fftshift(fft(ifftshift(phi_ir,2),[],2),2)./(length(k_vector.kx)*k_vector.dkx/(2*pi));
        end

        function e_f=fft2_ttof(OP,e_t,LN)
        e_f=fftshift(fft(ifftshift(e_t,1),[],1),1)./(length(e_t(:,1)).*LN.df); 
        end

        function e_t=fft2_ftot(OP,e_f,LN)
        e_t=fftshift(ifft(ifftshift(e_f,1),[],1),1).*(length(e_f(:,1)).*LN.df); 
        end


        function P_spm_wx_ir=spm_wx_ir(OP,LN,e_ir_tx,kx_ir,kz_ir,x0,z0)
                   P_spm_wx=OP.fft2_ttof(exp(-1i.*LN.phase_modify).*abs(e_ir_tx).^2.*e_ir_tx,LN); 
                   %omega dimension
                   P_spm_wx_ir=bsxfun(@times, P_spm_wx,-LN.omega_b.^2.*LN.n_ir_b.^2.*my_c.eps./my_c.c).*conj(OP.propagation_phase(kx_ir,x0,kz_ir,z0));
        end

        function P_spm_wx_shg=spm_wx_shg(OP,LN,e_shg_tx,kx_shg,kz_shg,x0,z0)
                   P_spm_wx=OP.fft2_ttof(abs(e_shg_tx).^2.*e_shg_tx,LN); 
                   %omega dimension
                   P_spm_wx_shg=bsxfun(@times, P_spm_wx,-LN.omega_shg.^2.*LN.n_shg.^2.*my_c.eps./my_c.c).*conj(OP.propagation_phase(kx_shg,x0,kz_shg,z0));
        end

        function P_srs_wx=srs_wx(OP,LN,e_ir_tx,kx_ir,kz_ir,x0,z0)
                   h_r = raman(LN.omega_b,LN.omega_0);
                   I_ir_wx2=OP.fft2_ttof(exp(-1i.*LN.phase_modify).*abs(e_ir_tx).^2,LN);             
                   P_srs_tx=OP.fft2_ftot(bsxfun(@times, I_ir_wx2,h_r),LN); 
                   P_srs_wx=OP.fft2_ttof(P_srs_tx.*e_ir_tx,LN); 
                   P_srs_wx=bsxfun(@times, P_srs_wx,-LN.omega_b.^2.*LN.n_ir_b.^2.*my_c.eps./my_c.c).*conj(OP.propagation_phase(kx_ir,x0,kz_ir,z0));

        end
        

    end
end
