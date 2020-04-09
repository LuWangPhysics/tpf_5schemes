function [n_g, v_g] =Gro_v(omega_0,omega,n_p)
c=3e8;
d_beta_ir = diff(omega.*n_p/c); %k in medium for all frequency ir
ng_w = c*d_beta_ir./diff(omega);
w_position=(omega_0-omega(1))/(omega(end)-omega(1));
n_g= ng_w(round((length(omega)*w_position)));  
v_g = c/n_g;
end