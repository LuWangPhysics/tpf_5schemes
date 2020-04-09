function [phi]=initialization_input_field(e_h,k_vector,LN,OP)
phi=struct;
%E(x,w)=e_h(x,w)exp(-ikx*x)exp(-ikz*z)
        e_ir_k=fftshift(fft(ifftshift(e_h,2),[],2),2)./(length(k_vector.kx)*k_vector.dkx/(2*pi));
        phi.ir=e_ir_k;  
    
        phi.thz=zeros(length(LN.omega_t),length(k_vector.kx));
        phi.dir=ones(size(phi.ir)).*0;

        phi.dthz=zeros(size(phi.thz));


e_ir_wx=OP.fft2_kxtox(e_ir_k,k_vector);
%construct initial condition for phi.dthz
e_ir_tx=exp(1i.*LN.phase_modify).*OP.fft2_ftot(e_ir_wx,LN);
 I_ir_wx=OP.fft2_ttof(abs(e_ir_tx).^2,LN);      
 %-------------------------------------------------------------------          
%Terahertz polarization
%----------------------------------------------------------------------          
           P_dfg_thz_wx=bsxfun(@times, I_ir_wx((length(LN.omega_b)+3)/2:end,:).*LN.chi_2_0',-LN.omega_t.^2 ./(my_c.c^2));
	       P_dfg_thz_w_kx=OP.fft2_xtokx(P_dfg_thz_wx,k_vector);
	       phi.dthz=0*P_dfg_thz_w_kx./(-2i*k_vector.kz_thz);   
       






end