function plot_2d_slow(OP,save_string,k_vector,phi,LN,x0,z0,flag)
        e_ir_wx=OP.propagation_phase(k_vector.kx_ir,x0,k_vector.kz_ir,z0).*OP.fft2_kxtox(phi.ir.*OP.diffraction_phase(k_vector.kx_ir,k_vector.kx,k_vector.kz_ir,z0),k_vector);
        e_thz_wx=OP.propagation_phase(k_vector.kx_thz,x0,k_vector.kz_thz,z0).*OP.fft2_kxtox(phi.thz.*OP.diffraction_phase(k_vector.kx_thz,k_vector.kx,k_vector.kz_thz,z0),k_vector);
       
        e_thz_wx_full=[conj(e_thz_wx(end:-1:1,:));zeros(1,length(x0));e_thz_wx];
        %reconstruct electric field  e(t,x)  
        e_ir_tx=exp(1i.*LN.phase_modify).*OP.fft2_ftot(e_ir_wx,LN);      
        e_thz_tx=OP.fft2_ftot(e_thz_wx_full,LN);
        
        save([ save_string '_thz_wx.mat'], 'e_thz_wx', '-v7.3'); 

        figure
        imagesc(x0,LN.f_t(1:fix(6e12/LN.df))./1e12,abs(e_thz_wx(1:fix(6e12/LN.df),:)));
        savefig(gcf,[ save_string '_thz_wx.fig' ],'compact')
         
        figure
        imagesc(k_vector.kx,LN.f_t./1e12,abs(phi.thz));
        savefig(gcf,[ save_string '_thz_wkx.fig' ],'compact')
        
        figure
        imagesc(x0,LN.f_b./1e12,abs(e_ir_wx));
        set(gca,'YDir','normal')
        savefig(gcf,[ save_string '_ir_wx.fig' ],'compact')
        


        figure
        imagesc(x0,LN.t0.*1e12,real(e_thz_tx));
        set(gca,'YDir','normal')
        savefig(gcf,[ save_string '_thz_tx.fig' ],'compact')
        save([ save_string '_thz_tx.mat'], 'e_thz_tx', '-v7.3'); 
        

        figure
        imagesc(x0,LN.t0.*1e12,abs(e_ir_tx));
        set(gca,'YDir','normal')
        savefig(gcf,[ save_string '_ir_tx.fig' ],'compact')
end
