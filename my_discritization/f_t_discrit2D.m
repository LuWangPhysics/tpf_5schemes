function [s]=f_t_discrit2D(omega_0,my_input,deltaw,material_name)

%Central frequency in radians/second
%IR Frequencies 1.5-2.1
%Terahertz angular frequency 
T=my_c.T;
s=struct;  
s.omega_0=omega_0;                                                        
s.omega_b = (omega_0 + 2*pi*(-my_input.w_low:deltaw:my_input.w_high)*1e12)';

s.delta_omega=s.omega_b-omega_0;
                                               
s.f_b=(s.omega_b./(2*pi));    
s.omega_t=1e12*2*pi*linspace(deltaw,(my_input.w_high+my_input.w_low)/2,fix((length(s.omega_b))/2))';

s.f_t=s.omega_t./(2*pi);
s.df=s.f_b(2)-s.f_b(1);

% Vector of time points
s.t0  = -(length(s.omega_b)-1)/(2*length(s.omega_b)*s.df):1/(length(s.omega_b)*s.df):(length(s.omega_b)-1)/(2*length(s.omega_b)*s.df);
s.phase_modify=((my_input.w_low+my_input.w_high)/2-my_input.w_low)*1e12*2*pi.*s.t0';
         


switch material_name
    case 'LN'
  
        [s.n_ir_0,s.n_ir_b,s.n_thz_c,s.n_thz,s.alpha_thz,s.n_g,s.chi_2_0,s.n2_0]=LiNb03(omega_0,s.omega_b,s.omega_t,T,my_input.f_c_t);
   
    otherwise
        print('does not exisit, waiting to be implemented by the user!:)')
end

        
        
        d_k_dw = diff(s.omega_b.*s.n_ir_b)./diff(s.omega_b); %k in medium for all frequency ir
        s.n_g= d_k_dw(fix(my_input.w_low*1e12./s.df));  
 
end