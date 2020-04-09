%Lu Wang's awosome code 2019.12.30
        clear all
        tic
        restoredefaultpath; 
        initialize_my_path;       
        my_input=input_Lu();

%--------------------------------------------------------------------------
%Pump parameters
%--------------------------------------------------------------------------
        
        my_input.lambda_0=1018e-9;                                         %center wavelength
        my_input.f_c_t=0.3e12;                                             %target PM frequency
        my_input.tau_fwhm= 500e-15;                                        %Input pulse durationg FWHM
        my_input.sigma_in=0.5e-3;                                          %beamsize in TPF plane before all the optical elements
        my_input.sigma_homo=(3.5)*1e-3;                                    %beamsize in vertical plane
               
        name_material='LN';                                                %Choose nonlinear material LN
%------------------------------------------------------------------------
%choose "my_setup_name" between 5 different schemes:
%'grating', 'tpf_silicon', 'echelon_reflect' 'tpf_NLES'  %'spatial_temporal' 
%------------------------------------------------------------------------
        my_setup_name='grating';
        my_input.f1=300e-3;                                                %focal length of first lense
        my_input.p=1e-3/1500;                                              %grating density
        %damage energy with a gaussian spatial profile
        my_input.energy=2*pi*my_input.sigma_in*my_input.sigma_homo*sqrt(my_input.tau_fwhm/1e-18)/2;              
       
        my_input=my_input.input_Lu_construct(my_setup_name,name_material);  
        [n_gaussian,LN,e_f]=input_beam(my_input,my_input.deltaw,name_material);     %initialize input beam properties 
           
        name_addition=['_notes_'];                                           %extra notation in saving
               
%--------------------------------------------------------------------------
%crystal parameterscrystal_tilt
%--------------------------------------------------------------------------               
        LN.alpha_tpf_aim=acos(LN.n_g/LN.n_thz_c);                             %tpf pm angle
        LN.alpha_crystal=LN.alpha_tpf_aim;

%------------------------------------------------------------------------
%crystal coordinate setup
%------------------------------------------------------------------------
        x0=x_discritization(my_input,my_setup_name);
%------------------------------------------------------------------------  
%configure the imaging and input electric field according to the scheme
%------------------------------------------------------------------------  
        [k_vector,e_h,x0,LN,my_input]=configure_my_setup(my_setup_name,my_input,LN,x0,e_f,LN.alpha_tpf_aim);   
        z_mesh=z_discritization(my_input,name_material);

          
%------------------------------------------------------------------------  
%manage the saving form+data_out.shg_energy(m);

%------------------------------------------------------------------------                
        file_name=[my_setup_name name_addition  name_material];
        file_name=[file_name '_Temperature' num2str(my_c.T) 'trans_bs' num2str(my_input.sigma_homo*1e3) 'mm' 'bstpf' num2str(my_input.sigma_in*1e3) 'mm' '_d' num2str(LN.image_d)  'F_peak' num2str(LN.F_peak_irframe,'%1.0f') ];
        mkdir(['my_output/' file_name]);
        save_string=['my_output/' file_name '/' my_input.name '_' num2str(my_input.lambda_0.*1e6) 'um' ];
        fprintf([save_string '\n']);
%------------------------------------------------------------------------
%add check points
%------------------------------------------------------------------------   

        if exist( ['my_output/' file_name '/checkpoint.mat'],'file' ) % If a checkpoint file exists, load it
                 fprintf('Checkpoint file found - Loading\n');
                 load(['my_output/' file_name '/' 'checkpoint.mat']);

        else %otherwise, start from the beginning
                 fprintf('No checkpoint file found - starting from beginning\n');


        %------------------------------------------------------------------------
        %construct material refractive index and nonlinear index
        %------------------------------------------------------------------------
                crystal_tilt= material_matrix_construct(x0,LN,my_input.sigma_crystal/cos(LN.alpha_crystal), my_setup_name);

        %------------------------------------------------------------------------
        %e_h (omega,x) define thz and ir in omega and kx space phi_ir(w,kx)
        %------------------------------------------------------------------------
                OP=initialize_OP;
                phi=initialization_input_field(e_h,k_vector,LN,OP);        

                m=1;                                                       %iteration number
                z_check_point=1500;                                        %save after 1500 steps
                RK=initialize_RK;                                          %initialize numerial method 
                data_out=initialize_data(z_mesh);

        end    
      
while z_mesh.z0<z_mesh.L_total
            
            [data_out,phi,z_mesh]=RK_slow_vary(crystal_tilt,data_out,phi,LN,k_vector,x0,z_mesh,RK,m,OP);
            my_plot_execute(z_check_point,save_string,z_mesh,data_out,m,x0,k_vector,phi,LN,OP);
            m=m+1;         

            if(rem(m,z_check_point)==0)
                   fprintf('Saving checkpoint\n');
                   save(['my_output/' file_name '/' 'checkpoint.mat'],'-v7.3');                   
               
            end
end
 
%------------------------------------------------------------------------
%save final, remove checkpoint 
%------------------------------------------------------------------------  
       
        celldata =[save_string '_' num2str(z_mesh.z0*1000,3) '_mm_'];
        plot_1d(OP,k_vector,phi,LN,data_out,celldata,z_mesh,m,x0);
        plot_2d_slow(OP,celldata,k_vector,phi,LN,x0,z_mesh.z0);
        close all
        delete(['my_output/' file_name '/' 'checkpoint.mat']);             %remove the checkpoint saving file
 toc


