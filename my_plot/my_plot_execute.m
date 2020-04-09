function my_plot_execute(z_check_point,save_string,z_mesh,data_out,m,x0,k_vector,phi,LN,OP)
   
                             
               if(m==1)
                celldata =[save_string '_' num2str(z_mesh.z0*1000,3) '_mm_'];
                out_couple_flag=0;
                plot_1d(OP,k_vector,phi,LN,data_out,celldata,z_mesh,m,x0);
                plot_2d_slow(OP,celldata,k_vector,phi,LN,x0,z_mesh.z0,out_couple_flag);
                close all
               elseif((data_out.eff(m-1)>data_out.eff(m))&&data_out.eff(m-1)>data_out.eff(m-2))
                    celldata =[save_string '_' num2str(z_mesh.z0*1000,3) '_mm_'];
                    plot_1d(OP,k_vector,phi,LN,data_out,celldata,z_mesh,m,x0);
                    out_couple_flag=1;
                    plot_2d_slow(OP,celldata,k_vector,phi,LN,x0,z_mesh.z0,out_couple_flag);
                    close all
                  
               end
                %save check points
               if(rem(m,z_check_point)==0)
                   fprintf('Saving checkpoint\n');
                               
                   celldata =[save_string '_' num2str(z_mesh.z0*1000,3) '_mm_'];
                   plot_1d(OP,k_vector,phi,LN,data_out,celldata,z_mesh,m,x0);
                   out_couple_flag=0;
                   plot_2d_slow(OP,celldata,k_vector,phi,LN,x0,z_mesh.z0,out_couple_flag);
                   close all
               end

                 fprintf('%d \n',m);
end