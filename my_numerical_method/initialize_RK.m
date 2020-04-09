function RK=initialize_RK
        RK=struct;
        RK.array=[0,0.5,0.5,1];
        aa=cell(1,5);
        aa{1}=0;
        RK.ir_ds=aa;
        RK.ir_s=aa;
     
        RK.thz_ds=aa;
        RK.thz_s=aa;
end