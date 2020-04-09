function  x0=x_discritization(my_input,my_setup_name)
          x_scale=1;
         if(my_input.sigma_in>3e-3)
             x0=-(10.4*my_input.sigma_in/x_scale):(my_input.sigma_in/540/x_scale):(10.4*my_input.sigma_in/x_scale);
         elseif(strcmp(my_setup_name,'spatial_temporal')) 
             x0=-(15.4*my_input.sigma_in/x_scale):(my_input.sigma_in/320/x_scale):(15.4*my_input.sigma_in/x_scale);     
         else
             x0=-(15.4*my_input.sigma_in/x_scale):(my_input.sigma_in/220/x_scale):(15.4*my_input.sigma_in/x_scale);
         end
         
end 