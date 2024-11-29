function Run_DG_IGA_2D_Linear_Eigenvalue()

% DG-IGA method for the 2D linear eigenvalue problem: Example 1 and Example
% 2

%%
addpath('../../NURBS/')
addpath('../DG_Assemble_2D/')
addpath('../IGA_2D_Grid/')



%%
 
 Example    =  'Example_1';     % 'Example_1' (Example 1) or 'Example_2' (Example 2)
 Refinement = 1:6;              % Uniformly refine each patch Refinement times

 n_refine = length(Refinement); % The number of meshes.
 
 Lambda_h = zeros(n_refine,4);  % For the 2D linear eigenvalue problem, we compute the first four eigenvalues.
 
 n_eles = zeros(n_refine,1);    % The number of elements in the whole domain
 
 
 DOFs = zeros(n_refine,1);      % The number of DOFs in the whole domain


 pu = 1;

 t = 0:2;



for j = 1:length(t)
    
    p = pu + t(j); % The degree of B-splines used for simulation
    
    for i=1:n_refine
    [n_eles(i),DOFs(i), Lambda_h(i,:)] = IGA_DG_2D_Linear_Eigenvalue(Refinement(i),Example, t(j));

    end


save_EigenValue_Results(Lambda_h,n_eles,DOFs, p, Example );

end


end