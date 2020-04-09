function [E_o] = field_homo(theta_tpf,n_g,m_out_o,sigma_in,x0,omega_b,e_f,omega_0)
%in the optical systerm coordinate. along the propagation direction
            %-------------------------------------------------------------
            %in optical system
            %e_w_b=e_f*exp(-(x0.^2./sigma_in^2./2).^n_gaussian);
            % my_sigma*sqrt(2)=sigma_gaussian; not the same as conventional
            % gaussian!!!!
            %--------------------------------------------------------------
            %q_in_b=ik\sigma^2/2
            q_in_b  =1i*sigma_in^2*omega_b/my_c.c;             

            %rotate to crystal coordinate
            %x z IR frame coordinate x0, z0 are tpf coordinate

            z=-sin(theta_tpf).*x0;
            x=x0.*cos(theta_tpf); 
            %z0=-0*tan(theta_tpf).*x;
            %x_1=0 theta_1=0,x_2=A_13,theta_2=A_23=m_out_o(2,3) add
            %coordinate chagne inside LN
            
            A_11=m_out_o(1,1)+z.*m_out_o(2,1);    
            A_12=m_out_o(1,2)+z.*m_out_o(2,2);
            A_13=m_out_o(1,3)+z.*m_out_o(2,3);
            k=omega_b/my_c.c;  
            
            %different material included in phi_x_1 and phi_x_2
            %divergence kx angle takes account for x direction shift
            phi_x_1=-k*sin(m_out_o(2,3))*A_13/2;
            %take phi tpf outside
            phi_x_2=0*k.*m_out_o(2,3);
            %phi_tpf=k.*m_out_o(2,3).*x0;
            phi_guoy=atan(z/abs(q_in_b));                      
            q_out_b = (A_11.*q_in_b+ A_12)./(m_out_o(2,1).*q_in_b+m_out_o(2,2));

            huygens=sqrt(abs(1./(A_11+A_12./q_in_b)).*exp(1i.*unwrap(angle(1./(A_11+A_12./q_in_b)))/2));
            %output field in IR coordinate.                        
 
            E_o = huygens.*exp((-1i*k./2./q_out_b).*(x-A_13).^2).*e_f.*exp(-1i*(phi_x_1+phi_x_2)).*exp(1i*phi_guoy);
  
 
               
end
