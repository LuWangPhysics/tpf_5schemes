function plot_snapshot(LN,data_out,celldata,z_mesh,m)
                         figure
                         plot(LN.f_b./1e12,data_out.ir_spectrum);
                         xlabel('Frequency (Thz)');
                         title('IR spectrum');
                         ylabel('|E_{IR}|'); 
                        
                          savefig(gcf,[ celldata 'ir.fig'])

                        figure
                        plot(LN.f_t./1e12,data_out.thz_spectrum);
                        xlabel('Frequency (Thz)');
                        title('THz spectrum');
                        ylabel('|E_{THz}|'); 
                 
                        savefig(gcf,[ celldata 'thz.fig' ])
                        
                        figure
                        plot(z_mesh.dz.*(1:(m-1)).*1e3,data_out.eff(1:(m-1)))
                        xlabel('Distance (mm)')
                        title('Efficiency VS distance')
                        ylabel('Conversion efficiency \eta %')
                        savefig(gcf,[ celldata 'eff.fig' ])

  
                        
                        figure
                        plot(z_mesh.dz.*(1:length(data_out.energy_cons)).*1e3,data_out.energy_cons)
                        xlabel('z (mm)')
                        ylabel('energy_cons')
                        savefig(gcf,[ celldata '_energy_cons.fig' ])
end