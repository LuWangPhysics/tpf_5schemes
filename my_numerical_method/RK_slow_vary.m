function [data_out,phi,z_mesh]=RK_slow_vary(crystal_tilt,data_out,phi,LN,k_vector,x0,z_mesh,RK,m,OP)
        RK.n=2;
    
       for rk_loop=1:4
        [RK]=pol_slow_vary(crystal_tilt,phi,LN,k_vector,x0,z_mesh.z0,z_mesh.dz,RK,OP);  
       end

       dz=z_mesh.dz;
        %Solve for phi term of ir and thz in (omega, kx) dimension
	 RK_sum=@(x) x{2}+2*x{3}+2*x{4} +x{5};
      

        if(z_mesh.thz_vanish==m)
            phi.thz=zeros(size(phi.thz));
        end
        
        phi.thz = phi.thz+ dz*(RK_sum(RK.thz_s))/6;

        phi.ir= phi.ir+ dz.*(RK_sum(RK.ir_s))/6;	



        z_mesh.z0=z_mesh.z0+z_mesh.dz;
        data_out=output_calculation(data_out,phi,k_vector,LN,m);
end
