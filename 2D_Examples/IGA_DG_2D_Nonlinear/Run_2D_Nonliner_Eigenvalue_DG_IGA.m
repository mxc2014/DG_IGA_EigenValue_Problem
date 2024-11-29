function Run_2D_Nonliner_Eigenvalue_DG_IGA()

%% DG-IGA method for the 2D nonlinear eigenvalue problem

 addpath('../../NURBS/')
 addpath('../IGA_2D_Grid/')
 addpath('../DG_Assemble_2D/')
   
 Refinement = 1:3;

 n_refine = length(Refinement);
 
 Lambda_h = zeros(n_refine,1);
 
 n_eles = zeros(n_refine,1);
 
 
 DOFs = zeros(n_refine,1);


 
pu = 1;

t =  0:2;



for j = 1:length(t)
    
    p = pu+t(j);
    for i=1:n_refine
    [n_eles(i),DOFs(i), Lambda_h(i,:)] = DG_IGA_2D_Nonliner_Eigenvalue(Refinement(i),t(j));
    end



save_EigenValue_Results(Lambda_h,n_eles,DOFs, p );

end


end