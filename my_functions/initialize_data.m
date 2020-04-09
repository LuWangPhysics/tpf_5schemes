function data_out=initialize_data(z_mesh)
        data_out=struct;
        data_out.eff=zeros(z_mesh.N,1);
        data_out.energy_cons=zeros(z_mesh.N,1);
        data_out.ir_energy=zeros(z_mesh.N,1);
        data_out.shg_energy=zeros(z_mesh.N,1);
        data_out.spectrum=@spectrum;
        data_out.spatial_ir=@spatial_ir;
        data_out.spatial_thz=@spatial_thz;
        
        
function spec=spectrum(phi_ir,k_z,k0,k_vector)
spec=(k_z./k0).*sum(0.5*my_c.c*my_c.eps*abs(phi_ir).^2,2)*k_vector.dkx/(2*pi);
end
        
function spat=spatial_ir(OP,phi,LN,k_vector,z0)
 k0=LN.omega_b./my_c.c;
u_wx=fftshift(ifft(ifftshift(phi.ir.*OP.diffraction_phase(k_vector.kx_ir,k_vector.kx,k_vector.kz_ir,z0),2),[],2),2).*(length(k_vector.kx)*k_vector.dkx/(2*pi));
spat=sum(0.5*my_c.c*my_c.eps*(k_vector.kz_ir./k0).*abs(u_wx).^2,1)*LN.df;
end
        
function spat=spatial_thz(OP,phi,LN,k_vector,z0)
 k0=LN.omega_t./my_c.c;
u_wx=fftshift(ifft(ifftshift(phi.thz.*OP.diffraction_phase(k_vector.kx_thz,k_vector.kx,k_vector.kz_thz,z0),2),[],2),2).*(length(k_vector.kx)*k_vector.dkx/(2*pi));
spat=sum(0.5*my_c.c*my_c.eps*(k_vector.kz_thz./k0).*abs(u_wx).^2,1)*LN.df;
end
end