classdef input_Lu
    properties 
        lambda_0; 
        f_c_t;
        tau_fwhm;                      
        d_beam_peak;
        sigma_in;                                                          %tpf dimension, sigma e-2 =sigma*sqrt(2);
        sigma_crystal;
        sigma_homo;                                                        %transvers dimension
        
        energy;
        w_high;                                                            %high and low frequency domian range
        w_low;
        deltaw; 
        
        GDD = 0;
        TOD=0;
        p;
        f1;
        f2;
        image_distance=1;
        name='LW'
        
    end
    
    methods
            function obj=input_Lu()
            end

            
            function obj=input_Lu_construct(obj,name_string,name_material)

                          L_alpha=1.5e-3;
                  switch name_string

                      case 'grating'
                             obj.deltaw =0.01;  
                              switch name_material
                                  case 'LN'

                                       if(obj.tau_fwhm>200e-15)    
                                              if(obj.sigma_in<L_alpha)
                                                  obj.d_beam_peak=L_alpha;
                                                  obj.w_low=32;
                                                  obj.w_high=18;
                                              else
                                                  obj.d_beam_peak=2.2*obj.sigma_in;
                                                  obj.w_low=44;
                                                  obj.w_high=21;
                                              end
                                      else
                                                  obj.d_beam_peak=L_alpha;
                                                  obj.w_low=42;
                                                  obj.w_high=32;  
                                      end

                                  case  'KTP'


                                     obj.d_beam_peak=1.5*L_alpha;
                                     obj.w_low=25;
                                     obj.w_high=20;
                              end



                              case 'tpf_NLES'
                                   obj.deltaw =0.01;  
                                  if(obj.tau_fwhm>200e-15)    
                                              obj.d_beam_peak=L_alpha;
                                              obj.w_low=25;
                                              obj.w_high=15;

                                  else
                                              obj.d_beam_peak=L_alpha;
                                              obj.w_low=42;
                                              obj.w_high=32;  
                                  end

                              case 'spatial_temporal'
                                  obj.deltaw =0.03;  
                                  obj.w_low=100;
                                  obj.w_high=75;
                                  obj.d_beam_peak=2*L_alpha;


                              case 'tpf_silicon'

                                  obj.deltaw =0.01;  
                                  if(obj.tau_fwhm>200e-15)    
                                          if(obj.sigma_in<L_alpha)
                                              obj.d_beam_peak=2*L_alpha;
                                              obj.w_low=18;
                                              obj.w_high=10;
                                          else
                                              obj.d_beam_peak=2.2*obj.sigma_in;
                                              obj.w_low=26;
                                              obj.w_high=15;
                                          end

                                  else
                                              obj.d_beam_peak=L_alpha;
                                              obj.w_low=32;
                                              obj.w_high=32;  
                                  end



                             %------------------------------
                             %Reflective enchelon
                             %--------------------------------
                              case 'echelon_reflect'
                                  obj.deltaw =0.01;  
                                  if(obj.tau_fwhm>200e-15)    
                                          if(obj.sigma_in<L_alpha)
                                              obj.d_beam_peak=2*L_alpha;
                                              obj.w_low=18;
                                              obj.w_high=10;
                                          else
                                              obj.d_beam_peak=2.2*obj.sigma_in;
                                              obj.w_low=23;
                                              obj.w_high=13;
                                          end

                                  else
                                              obj.d_beam_peak=L_alpha;
                                              obj.w_low=32;
                                              obj.w_high=32;  
                                  end


                              otherwise
                                  fprintf('scheme not found! check your spelling!!')
                  end


             end

    end
  
end
