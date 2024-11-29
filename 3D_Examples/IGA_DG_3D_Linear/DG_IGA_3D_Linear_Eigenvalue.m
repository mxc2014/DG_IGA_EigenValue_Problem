function [lambda_h, n_eles, n_dofs] = DG_IGA_3D_Linear_Eigenvalue(nurbs_original, Refinement,t,mp)

datetime('now');


% DG-IGA method for the 3D linear eigenvalue problem


s = [0 0 0];

V_ext = @(x,y,z)-1.0/norm([x y z] - s);

% V_ext = @(x,y,z)(x^2+y^2+z^2)/2;

beta = 300;




n_subdomains = mp.n_subdomains;
 
 
%%

nurbs_refine_domains = cell(n_subdomains,1);

t_subdomains       = t*ones(n_subdomains,1);
Refine_subdomains  = Refinement*ones(n_subdomains,1);

%%
e = 14; Refine_subdomains(e) = Refinement + 1; % The sub-domain contains the singular point
%%

disp('******************************************************************************************')

disp('Begin to assemble the matrices for DG-IGA method...')

tic

for e = 1:n_subdomains
  
    [knotU,knotV,knotW] =  updateKnotVectors_3D(nurbs_original{e});%%  We want to have a regular mesh
 
    [knotU,knotV, knotW]=IGADegreeElevVolume(knotU,knotV,knotW,t_subdomains(e)); % Now we should consider the ratio of lengths in x- and y-directions
    pu = nurbs_original{e}.pu + t_subdomains(e);     pv = nurbs_original{e}.pv + t_subdomains(e); pw = nurbs_original{e}.pw + t_subdomains(e); 
    
    nurbs_refine_domains{e} = Iga_3d_grid(knotU,pu,knotV,pv,knotW,pw,Refine_subdomains(e));
    
    
end

n_dofs     =  nurbs_refine_domains{1}.n_dofs; % The total number of DOFs in all sub-domains

for e = 1:n_subdomains
    nurbs_refine_domains{e}.n_dofs_subDomains = nurbs_refine_domains{e}.n_dofs;
end
    

for e = 1:(n_subdomains-1)
    nurbs_refine_domains{e+1}.n_dofs_subDomains = n_dofs + nurbs_refine_domains{e+1}.n_dofs_subDomains;
    n_dofs = nurbs_refine_domains{e+1}.n_dofs_subDomains;
end

 

%%

A         = sparse(n_dofs,n_dofs);
M_ext     = A;
Mass_Mat  = sparse(n_dofs,n_dofs);
 
 

for e = 1:n_subdomains
    [stiffness_matrix, M_ext_domains, Mass_domains] = solve_laplace_A_M_3D(nurbs_original{e},nurbs_refine_domains{e},V_ext, e);
    
    if e==1
       row = 1:nurbs_refine_domains{1}.n_dofs_subDomains;
    else
      row = ( 1 + nurbs_refine_domains{e-1}.n_dofs_subDomains ): nurbs_refine_domains{e}.n_dofs_subDomains;
    end
    
    A(row,row)          =   stiffness_matrix;
    M_ext(row,row)      =   M_ext_domains;
    Mass_Mat(row,row)   =   Mass_domains;
%     disp('Assemble the matrices is OK')
end

 
uNoEs = nurbs_refine_domains{1}.uNoEs;


h = 1/uNoEs;

 n_eles = 0;
 
 for e = 1:n_subdomains
     n_eles = n_eles + nurbs_refine_domains{e}.NoEs;
 end



%% Update dof index from one sub-domain to its neighbour sub-domain
for e = 2:n_subdomains
%     disp([e, nurbs_refine_domains{e-1}.n_dofs])
    nurbs_refine_domains{e} = Update_Face_DoFs( nurbs_refine_domains{e-1}.n_dofs_subDomains, nurbs_refine_domains{e} );
end
%%  

nx = mp.nx;  % nx is the number of sub-intervals along the x-axis
ny = mp.ny;
nz = mp.nz;

ele_info = get3DElements_index(nx,ny,nz);

bottom2top_eles = ele_info.bottom2top_eles ;
back2front_eles = ele_info.back2front_eles ;
left2right_eles = ele_info.left2right_eles ;



%% From bottom to top

P = sparse(n_dofs,n_dofs);
S = P;

for k = 1:(nz-1)
%     disp('top-bottom case')
    bottom_ele_index = bottom2top_eles(k,:);
    top_ele_index    = bottom2top_eles(k+1,:);
    for i = 1:(nx*ny)
        [P_subdomains,S_subdomains] = IGA_DG_Top_Bottom_Face_Assemble(nurbs_original{bottom_ele_index(i)},nurbs_original{top_ele_index(i)},...
                                                       nurbs_refine_domains{bottom_ele_index(i)},nurbs_refine_domains{top_ele_index(i)}, n_dofs);
        
       P = P + P_subdomains;
       S = S + S_subdomains;
       
    end
end


 

%%  From left to right

for j = 1:(ny-1)
%     disp('left-right case')
    left_ele_index     = left2right_eles(j,:);
    right_ele_index    = left2right_eles(j+1,:);
    for i = 1:(nx*nz)
        [P_subdomains,S_subdomains] = IGA_DG_Left_Right_Face_Assemble(nurbs_original{left_ele_index(i)},nurbs_original{right_ele_index(i)},...
                                                       nurbs_refine_domains{left_ele_index(i)},nurbs_refine_domains{right_ele_index(i)}, n_dofs);
        
       P = P + P_subdomains;
       S = S + S_subdomains;
       
    end
end




%%  From back to front

for i = 1:(nx-1)
%     disp('back-front case')
    back_ele_index     = back2front_eles(i,:);
    front_ele_index    = back2front_eles(i+1,:);
    for d = 1:(ny*nz)
        [P_subdomains,S_subdomains] = IGA_DG_Back_Front_Face_Assemble(nurbs_original{back_ele_index(d)},nurbs_original{front_ele_index(d)},...
                                                       nurbs_refine_domains{back_ele_index(d)},nurbs_refine_domains{front_ele_index(d)}, n_dofs);
        
       P = P + P_subdomains;
       S = S + S_subdomains;
       
    end
end

%%
clear P_subdomains S_subdomains  stiffness_matrix  Mass_domains



Mat = A  + M_ext - S/2 - S'/2 + beta/h*P ;

clear S P


bnd_dof_index = getBoundaryDofs_index(nurbs_refine_domains,mp);




int_dof_index = 1:n_dofs;
int_dof_index(bnd_dof_index) = [];
Mat_int  = Mat(int_dof_index,int_dof_index);
Mass_int = Mass_Mat(int_dof_index,int_dof_index);

elapsedTime = toc;

fprintf('The running time for assembling matrices is %f seconds\n', elapsedTime);

disp('The matrices for DG-IGA method have been assembled!')

clear Mat Mass_Mat
%%
disp('===========================================================================================')

disp('Begin to compute the first eigenvalue of 3D linear eigenvalue problem ...')

n_eigenvalues = 1;


tic

[lambda_h,uh_int] = eigifp( Mat_int , Mass_int , n_eigenvalues); 
 
elapsedTime = toc;

fprintf('The running time for eigen solver is %f seconds \n', elapsedTime);

disp('The first eigenvalue of 3D linear eigenvalue problem has been coputed!')
  
 %%  Save uh and the NURBS information
  uh = zeros(n_dofs,n_eigenvalues); 
  uh(int_dof_index,:) = uh_int;
  
  nurbs = nurbs_refine_domains;
 
 
file_name_uh = strcat('./EigenFunction_Data/p_',num2str(pu));
file_name_uh = strcat(file_name_uh,'/uh_',num2str(Refinement));
file_name_uh = strcat(file_name_uh,'.mat');
save(file_name_uh,'uh')

file_name_nurbs = strcat('./EigenFunction_Data/p_',num2str(pu));
file_name_nurbs = strcat(file_name_nurbs,'/nurbs_',num2str(Refinement));
file_name_nurbs = strcat(file_name_nurbs,'.mat');
save(file_name_nurbs,'nurbs')
 
%%
fprintf('%s %d \n','The degree of B-splines is p = ', pu)
fprintf('No. of Elements = %d. No. of DOFs = %d, and the 1st Eigenvalue = %1.15f \n',n_eles, n_dofs, lambda_h);


end




