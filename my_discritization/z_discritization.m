function [z_mesh]=z_discritization(my_input,name_material)
        z_mesh=struct;
       switch name_material
           case 'LN'
        z_mesh.dz=8e-7;%1e-7 SHG
           case 'KTP'
             z_mesh.dz=3e-7;
       end
             z_mesh.z0 =z_mesh.dz;%
       % dz_min=4e-7;
       
    
       
        z_mesh.L_thz=(my_input.d_beam_peak); 
        z_mesh.L_total=z_mesh.L_thz;

        z_mesh.N=fix(z_mesh.L_total/z_mesh.dz);
        z_mesh.thz_vanish=(fix(z_mesh.L_thz/z_mesh.dz)+2);



        
end
