function [lambda_h,total_energy, n_eles, n_dofs] =  DG_IGA_3D_NonLinear_Eigenvalue(nurbs_original, Refinement,t,mp )

datetime('now')


s = [0 0  0]; % The singular point

V_ext = @(x,y,z) -2/norm([x y z] - s);% The external potential

rho = @(u) 2*u^2;

V_Har_bnd = @(x,y,z) 2/norm([x y z]);

% V_Har_bnd  = @(x,y,z)(x^2+y^2+z^2)/2;

beta = 100;



%%

n_subdomains = mp.n_subdomains;

nurbs_refine_domains = cell(n_subdomains,1);

t_subdomains       = t*ones(n_subdomains,1);
Refine_subdomains  = Refinement*ones(n_subdomains,1);

e = 14; Refine_subdomains(e) = Refinement + 1; % The sub-domain contains the singular point


%%

for e = 1:n_subdomains
  
    [knotU,knotV,knotW] =  updateKnotVectors_3D(nurbs_original{e});%%  We want to have a regular mesh
    [knotU,knotV, knotW]=IGADegreeElevVolume(knotU,knotV,knotW,t_subdomains(e)); % Now we should consider the ratio of lengths in x- and y-directions
    pu = nurbs_original{e}.pu + t_subdomains(e);     pv = nurbs_original{e}.pv + t_subdomains(e); pw = nurbs_original{e}.pw + t_subdomains(e); 
    
    nurbs_refine_domains{e} = Iga_3d_grid(knotU,pu,knotV,pv,knotW,pw,Refine_subdomains(e));
    
    
end


uNoEs = nurbs_refine_domains{1}.uNoEs;



h = 1/uNoEs;

 n_eles = 0;
 
 for e = 1:n_subdomains
     n_eles = n_eles + nurbs_refine_domains{e}.NoEs;
 end
 

n_dofs     =  nurbs_refine_domains{1}.n_dofs;
for e = 1:n_subdomains
    nurbs_refine_domains{e}.n_dofs_domains = nurbs_refine_domains{e}.n_dofs;
end

for e = 1:(n_subdomains-1) % Accumate the number of DOFs in each sub-domain
    nurbs_refine_domains{e+1}.n_dofs_domains = n_dofs + nurbs_refine_domains{e+1}.n_dofs_domains;
    n_dofs = nurbs_refine_domains{e+1}.n_dofs_domains;
end


%%
 
stiffness_matrix     = cell(n_subdomains,1);
V_ext_domains        = cell(n_subdomains,1);
Mass_domains         = cell(n_subdomains,1);
 
disp('******************************************************************************************')

disp('Begin to assemble the matrices for DG-IGA method...')

tic

for e = 1:n_subdomains
    [stiffness_matrix{e}, V_ext_domains{e}, Mass_domains{e}] = solve_laplace_A_M(nurbs_original{e},nurbs_refine_domains{e},V_ext);
end

 


A             = sparse(n_dofs,n_dofs);
V_ext_Mat     = A;
M             = A;
 

row = 1:nurbs_refine_domains{1}.n_dofs_domains;
A(row,row)         = stiffness_matrix{1};
V_ext_Mat(row,row) = V_ext_domains{1};
M(row,row)         = Mass_domains{1};

for e = 2:n_subdomains
    row = ( 1 + nurbs_refine_domains{e-1}.n_dofs_domains ): nurbs_refine_domains{e}.n_dofs_domains;
    A(row,row)         = stiffness_matrix{e}; 
    V_ext_Mat(row,row) = V_ext_domains{e};
    M(row,row)         = Mass_domains{e};
end

uh_idx = cell(n_subdomains,1);

for e = 1:n_subdomains
    
    if e==1
        uh_idx{e} = 1:nurbs_refine_domains{1}.n_dofs_domains;
    else
        uh_idx{e} = ( 1 + nurbs_refine_domains{e-1}.n_dofs_domains ): nurbs_refine_domains{e}.n_dofs_domains;
    end
        
end

%%






%% Update dof index from one sub-domain to its neighbour sub-domain

 

 

for e = 2:n_subdomains
    nurbs_refine_domains{e} = Update_Face_DoFs( nurbs_refine_domains{e-1}.n_dofs_domains, nurbs_refine_domains{e} );
end


%%  Handle the non-homogeneous Dirichlet boundary condition



 

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
    
    disp('top-bottom case')
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
    disp('left-right case')
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
    disp('back-front case')
    back_ele_index     = back2front_eles(i,:);
    front_ele_index    = back2front_eles(i+1,:);
    for d = 1:(ny*nz)
        [P_subdomains,S_subdomains] = IGA_DG_Back_Front_Face_Assemble(nurbs_original{back_ele_index(d)},nurbs_original{front_ele_index(d)},...
                                                       nurbs_refine_domains{back_ele_index(d)},nurbs_refine_domains{front_ele_index(d)}, n_dofs);
        
       P = P + P_subdomains;
       S = S + S_subdomains;
       
    end
end

clear P_subdomains S_subdomains

%%


A_DG = A/2   - S/2 - S'/2;

A_Poisson     = A - (S + S') + (beta/h)*P;
A_rhs         = A - (S + S');



Mat = A_DG +  V_ext_Mat + (beta/h)*P ;

clear S P

bnd_dof_index = getBoundaryDofs_index(nurbs_refine_domains,mp);
int_dof_index = 1:n_dofs;
int_dof_index(bnd_dof_index) = [];

A_Poisson(bnd_dof_index,:) = 0;
A_Poisson(:,bnd_dof_index) = 0;
A_Poisson(bnd_dof_index,bnd_dof_index) =speye(length(bnd_dof_index),length(bnd_dof_index));


Mass_Mat = M;
Mass_Mat(bnd_dof_index,:) = 0;
Mass_Mat(:,bnd_dof_index) = 0;
Mass_Mat(bnd_dof_index,bnd_dof_index) = speye(length(bnd_dof_index),length(bnd_dof_index));

Mass_int = M(int_dof_index,int_dof_index);

elapsedTime = toc;

fprintf('The running time for assembling matrices is %f seconds\n', elapsedTime);

disp('The matrices for DG-IGA method have been assembled!')

%% Handle the non-linear term

 u0_h = randn(n_dofs,1);
 u0_h(bnd_dof_index) = 0;

% u0_h = ones(n_dofs,1);
% Compute the L2 norm of uh in the subdomain 1





maxIt = 500;

tol = 1.0e-10;

alpha = 0.7;

lambda_old_h = 1000;

u1_h = zeros(n_dofs,1);

gh = L2_Projection_Boundary(nurbs_original,nurbs_refine_domains,mp,V_Har_bnd,n_dofs);

u0_h = normalized_uh(nurbs_original, nurbs_refine_domains, u0_h, n_subdomains);



disp('Begin the SCF iteration')

for k = 1:maxIt

    rhs = generate_Hartree_rhs(nurbs_original,nurbs_refine_domains, u0_h, rho, n_dofs, mp);

    rhs = rhs - A_rhs * gh;
    rhs(bnd_dof_index) = 0;
    Hartree_h = A_Poisson\rhs;
    Hartree_h = Hartree_h + gh;

    Hartree_Mat = sparse(n_dofs,n_dofs);
    
    row = 1:nurbs_refine_domains{1}.n_dofs_domains;
%     Hartree_subdomain = compute_Hartree(nurbs_original{1},nurbs_refine_domains{1},Hartree_h(row));

     Hartree_subdomain  = compute_Hartree_XC(nurbs_original{1},nurbs_refine_domains{1},Hartree_h(row), rho, u0_h(row));

     Hartree_Mat(row,row)  = Hartree_subdomain;
  
    
    
    for e = 2:n_subdomains
        row = ( 1 + nurbs_refine_domains{e-1}.n_dofs_domains ): nurbs_refine_domains{e}.n_dofs_domains;
 
        Hartree_subdomain     = compute_Hartree_XC(nurbs_original{e},nurbs_refine_domains{e},Hartree_h(row), rho, u0_h(row));
        Hartree_Mat(row,row)  = Hartree_subdomain;
    end
    
    Mat_Nonlinear = Mat + Hartree_Mat;
    
    Mat_int  = Mat_Nonlinear(int_dof_index,int_dof_index);
    
    
    [lambda_new_h,u1_h_int] = eigifp( Mat_int , Mass_int ,1); 
    u1_h(int_dof_index) = u1_h_int;
    
    format long
   
    
    
    if(abs(lambda_new_h - lambda_old_h)<=tol)
        break
    end
    
    
 
    
    u0_h = alpha*u1_h + (1-alpha)*u0_h;
    
    
    u0_h = normalized_uh(nurbs_original, nurbs_refine_domains, u0_h, n_subdomains);
    
    lambda_old_h = lambda_new_h;
    
   
    
    fprintf('The %d-th iteration, the first eigenvalue = %1.15f \n',k,lambda_old_h );
    
  
    
end


%%


u0_h  = normalized_uh(nurbs_original, nurbs_refine_domains, u0_h, n_subdomains);
uh    = u0_h;
nurbs = nurbs_refine_domains;


total_energy = get_total_energy(nurbs_original,nurbs_refine_domains,A,u0_h,rho,V_ext,Hartree_h,mp,uh_idx);

lambda_h = lambda_new_h;
 

fprintf('%s %d \n','The degree of B-splines is p = ', pu)
fprintf('No. of Elements = %d. No. of DOFs = %d, and 1st Eigenvalue = %1.15f. Total Energy = %1.15f \n',n_eles, n_dofs, lambda_h,total_energy);





degree_used = pu;

file_name_uh = strcat('./EigenFunction_Data/p_',num2str(degree_used));
file_name_uh = strcat(file_name_uh,'/uh_',num2str(Refinement));
file_name_uh = strcat(file_name_uh,'.mat');
save(file_name_uh,'uh')

file_name_nurbs = strcat('./EigenFunction_Data/p_',num2str(degree_used));
file_name_nurbs = strcat(file_name_nurbs,'/nurbs_',num2str(Refinement));
file_name_nurbs = strcat(file_name_nurbs,'.mat');
save(file_name_nurbs,'nurbs')

n_eles = 0;

end








