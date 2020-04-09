function [data_out]=output_calculation(data_out,phi,k_vector,LN,m)

      data_out.ir_spectrum=data_out.spectrum(phi.ir,k_vector.kz_ir,(LN.omega_b./my_c.c),k_vector);
      data_out.thz_spectrum=data_out.spectrum(phi.thz,k_vector.kz_thz,(LN.omega_t./my_c.c),k_vector);
    

      
        thz_energy=sum(data_out.thz_spectrum).*LN.df;
        data_out.ir_energy(m)=sum(data_out.ir_spectrum).*LN.df;

        data_out.eff(m)=100*thz_energy/LN.F_pump;
        data_out.energy_cons(m)=thz_energy+data_out.ir_energy(m);


end
