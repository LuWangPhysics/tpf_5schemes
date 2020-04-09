function plot_1d(OP,k_vector,phi,LN,data_out,celldata,z_mesh,m,x0)
                         figure
                         plot(LN.f_b./1e12,data_out.ir_spectrum);
                         xlabel('Frequency (Thz)');
                         title('IR spectrum');
                         ylabel('|I_{IR}|'); 
                        
                          savefig(gcf,[ celldata 'ir.fig'])

                        figure
                        plot(LN.f_t./1e12,data_out.thz_spectrum);
                        xlabel('Frequency (Thz)');
                        title('THz spectrum');
                        ylabel('|I_{THz}|'); 
                 
                        savefig(gcf,[ celldata 'thz.fig' ])
                        
                        figure
                        plot(z_mesh.dz.*(1:(m-1)).*1e3,data_out.eff(1:(m-1)))
                        xlabel('Distance (mm)')
                        title('Efficiency VS distance')
                        ylabel('Conversion efficiency \eta %')
                        savefig(gcf,[ celldata 'eff.fig' ])

     
                        figure
                        plot(x0.*1e3,data_out.spatial_thz(OP,phi,LN,k_vector,z_mesh.z0))
                        xlabel('x (mm)')
                        ylabel('thz fluence')
                        savefig(gcf,[ celldata '_thz_fluence.fig' ])
                        figure
                        plot(x0.*1e3,data_out.spatial_ir(OP,phi,LN,k_vector,z_mesh.z0))
                        xlabel('x (mm)')
                        ylabel('ir fluence')
                        savefig(gcf,[ celldata '_ir_fluence.fig' ])
                        
                     
end