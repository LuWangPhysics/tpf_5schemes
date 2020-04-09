classdef material_matrix_construct

%construct the tilted crystal nonlinearity

   properties
        x0_;
        chi_2_0_;
        n2_0_;
        alpha_;
        sigma_;
        name_;
        echelon_w_; %beamlets translational size
   end

methods
    
    function obj= material_matrix_construct(x0,LN,sigma, setup_name)
        %constructor
        obj.x0_=x0;
        obj.chi_2_0_=LN.chi_2_0;
        obj.n2_0_=LN.n2_0;
        obj.alpha_=LN.alpha_tpf_aim;
        obj.sigma_=sigma;
        obj.name_=setup_name;

        if strcmp('tpf_NLES',setup_name)
            obj.echelon_w_=LN.echelon_w;
        end
        
    end

    
    
    function slice=chi_2(obj,z0)

      switch obj.name_
          case 'tpf_NLES'
              %creat a changing flattop to mimic the wedge effect . make
              %sure that the center of the wedge is the same as the center
              %of the beamlet
              N_total=fix((obj.x0_(end)-obj.x0_(1))/(obj.echelon_w_/cos(obj.alpha_)));
              slice=0;
              
              for n=1:N_total
                  a=z0*tan(obj.alpha_);
                  b=z0*tan(pi/2-obj.alpha_);
                  d=(a+b)/2;
                  flat_shift=(b-a)/2;
                  n_edge=(abs(obj.x0_(1))+obj.echelon_w_/2/cos(obj.alpha_))/(obj.echelon_w_/cos(obj.alpha_));
                  edge=(n_edge-floor(n_edge))*(obj.echelon_w_)/cos(obj.alpha_);
                  slice=slice+exp(-((obj.x0_-(obj.x0_(1)+edge+n*obj.echelon_w_/cos(obj.alpha_)-flat_shift)).^2./d.^2).^5);
              end
              slice(slice>1)=1;
              slice=slice.*obj.chi_2_0_;

          otherwise
               slice=obj.chi_2_0_.*(1-tanh(5e4*(obj.x0_+(obj.sigma_-z0./tan(obj.alpha_)))))/2;    
      end

       
    end
    
    
    

    function slice=n_2(obj,z0)
        switch obj.name_
          case 'tpf_NLES'
              N_total=fix((obj.x0_(end)-obj.x0_(1))/(obj.echelon_w_/cos(obj.alpha_)));
              slice=0;
              
              for n=1:N_total
                  a=z0*tan(obj.alpha_);
                  b=z0*tan(pi/2-obj.alpha_);
                  d=(a+b)/2;
                  flat_shift=(b-a)/2;
                  slice=slice+exp(-((obj.x0_-(obj.x0_(1)+n*obj.echelon_w_/cos(obj.alpha_)-flat_shift)).^2./d.^2).^5);
              end
              slice(slice>1)=1;
              slice=slice.*obj.n2_0_;

          otherwise
              slice=obj.n2_0_.*(1-tanh(5e4*(obj.x0_+(obj.sigma_-z0./tan(obj.alpha_)))))/2;
        end     
    end
  
    
    
end




end