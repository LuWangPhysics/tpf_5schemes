function [RK] =pol_slow_vary(crystal_tilt,phi,LN,k_vector,x0,z0,dz,RK,OP) 
                
%--------------------------------------------------------------------
%RK updating
%-------------------------------------------------------------------
     
            phi_ir=phi.ir+RK.array(RK.n-1).*RK.ir_s{RK.n-1}.*dz;
            phi_thz=phi.thz+RK.array(RK.n-1).*RK.thz_s{RK.n-1}.*dz;
         

            e_ir_wx=OP.propagation_phase(k_vector.kx_ir,x0,k_vector.kz_ir,z0).*OP.fft2_kxtox(phi_ir.*OP.diffraction_phase(k_vector.kx_ir,k_vector.kx,k_vector.kz_ir,z0),k_vector);
            % phase compensation 0 for NLES and 1 for all other prism shapes
            e_ir_wx=e_ir_wx.*OP.phase_corr_prism(x0,z0,dz,LN);
            
            e_thz_wx=OP.propagation_phase(k_vector.kx_thz,x0,k_vector.kz_thz,z0).*OP.fft2_kxtox(phi_thz.*OP.diffraction_phase(k_vector.kx_thz,k_vector.kx,k_vector.kz_thz,z0),k_vector); 
            e_thz_wx_full=[conj(e_thz_wx(end:-1:1,:));zeros(1,length(x0));e_thz_wx];
            %reconstruct electric field  e(t,x)  
            e_ir_tx=exp(1i.*LN.phase_modify).*OP.fft2_ftot(e_ir_wx,LN);
     
            %get e_thz with frequency dependence
            e_thz_wx_full_chi=([LN.n_thz(end:-1:1);LN.n_thz_c;LN.n_thz]./LN.n_thz_c).^2.*e_thz_wx_full;
            e_thz_tx=OP.fft2_ftot(e_thz_wx_full_chi,LN);


%-------------------------------------------------------------------------
%calculate nonlinear terms for ir P_dfg(omega,x)   ,
%-----------------------------------------------------------------------
           %x dimensino of chi_2
           P_dfg_ir_wx=crystal_tilt.chi_2(z0).*OP.fft2_ttof(exp(-1i.*LN.phase_modify).*e_ir_tx.*e_thz_tx,LN); 
           P_dfg_ir_wx=-(LN.omega_b.^2 ./my_c.c^2).*P_dfg_ir_wx.*conj(OP.propagation_phase(k_vector.kx_ir,x0,k_vector.kz_ir,z0));
           %x dimension
           %stimulated raman         
           P_srs_wx=crystal_tilt.n_2(z0).*OP.srs_wx(LN,e_ir_tx,k_vector.kx_ir,k_vector.kz_ir,x0,z0);
           %P_spm_wx=crystal_tilt.n_2(z0).*OP.spm_wx_ir(LN,e_ir_tx,k_vector.kx_ir,k_vector.kz_ir,x0,z0);

           P_ir_w_kx=OP.fft2_xtokx((P_dfg_ir_wx+P_srs_wx),k_vector).*(0.5i./(k_vector.kz_ir));
           P_ir_w_kx(LN.ir_index)=0;
           RK.ir_s{RK.n}=P_ir_w_kx.*conj(OP.diffraction_phase(k_vector.kx_ir,k_vector.kx,k_vector.kz_ir,z0));
                 
%-------------------------------------------------------------------          
%Terahertz polarization
%----------------------------------------------------------------------  
            I_ir_wx=OP.fft2_ttof(abs(e_ir_tx).^2,LN);
            P_dfg_thz_wx=(LN.n_thz.^2./LN.n_thz_c.^2).*bsxfun(@times, I_ir_wx((length(LN.omega_b)+3)/2:end,:).*crystal_tilt.chi_2(z0),-0.5i.*LN.omega_t.^2 ./(my_c.c^2.*k_vector.kz_thz)).*conj(OP.propagation_phase(k_vector.kx_thz,x0,k_vector.kz_thz,z0));
            P_dfg_thz_w_kx=OP.fft2_xtokx(P_dfg_thz_wx,k_vector);
            RK.thz_s{RK.n}=(P_dfg_thz_w_kx).*conj(OP.diffraction_phase(k_vector.kx_thz,k_vector.kx,k_vector.kz_thz,z0))+bsxfun(@times,phi_thz,0.5.*LN.alpha_thz);
            RK.thz_s{RK.n}(LN.index_thz)=0;


            RK.n=RK.n+1;


end
