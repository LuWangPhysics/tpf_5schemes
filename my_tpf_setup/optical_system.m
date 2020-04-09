function [lens_M,M_optical,m23,s]=optical_system(my_input,LN)
%--------------------------------------------------------------------------
%Optical setup Matrices for imaging system
% grating-->s1-->lense m_l1--->s2 lens m_l2--->s3
%--------------------------------------------------------------------------
  
         s=my_input.image_distance.*my_input.sigma_in*tan(LN.alpha_crystal)*abs(LN.os_magnify);  %distance start inside the LN

         f1          = my_input.f1;
         f2          = my_input.f2;
         s1          = f1;                                                 %Grating to lens distance
         s2          = f1+f2;                                              %distance between two lens
         s3          = f2-s/LN.n_ir_0;                                     %Lens to crystal distance

        m_s1 = [1 s1 0;0 1 0;0 0 1];                                       %free propagation for s1 in vancumm
        m_l1  = [1 0 0;-1./f1 1 0;0 0 1];                                  %focusing
        m_s2 = [1 s2 0;0 1 0;0 0 1];                                       %free propagation for s2 in vancumm
        m_l2  = [1 0 0;-1./f2 1 0;0 0 1];                                  %focusing	
        m_s3 = [1 s3 0;0 1 0;0 0 1];                                       %free propagation for s2 in vancumm
        
        M_optical=m_s3*m_l2*m_s2*m_l1*m_s1;
        lens_M=f2/f1;
        m23=M_optical(2,2);
        %m23=(1-s2/f2+s1*(-(1/f2) - (1 - s2/f2)/f1));
       
end
