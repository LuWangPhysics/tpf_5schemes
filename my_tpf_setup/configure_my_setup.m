function [k_vector,e_h,x0,LN,my_input]=configure_my_setup(my_name,my_input,LN,x0,e_f,alpha_tpf_aim)

switch my_name
    case 'grating'
        [k_vector,e_h,x0,LN,my_input]=tpf_grating(my_input,LN,x0,e_f,alpha_tpf_aim);
	LN.prism_comp=1;
    case 'echelon_reflect'
      	[k_vector,e_h,x0,LN,my_input]=tpf_echelon(my_input,LN,x0,e_f,alpha_tpf_aim);
	LN.prism_comp=1;
    case 'spatial_temporal'
        [k_vector,e_h,x0,LN,my_input]=tpf_spatial_temporal(my_input,LN,x0,e_f,alpha_tpf_aim);
        LN.prism_comp=1;
    case 'tpf_NLES'
        [k_vector,e_h,x0,LN,my_input]=tpf_NLES(my_input,LN,x0,e_f,alpha_tpf_aim);
	LN.prism_comp=0;
         case 'tpf_silicon'
        [k_vector,e_h,x0,LN,my_input]=tpf_silicon(my_input,LN,x0,e_f,alpha_tpf_aim);
       LN.prism_comp=1;
    otherwise
        print('no matching setup found');
        
end
    
end
