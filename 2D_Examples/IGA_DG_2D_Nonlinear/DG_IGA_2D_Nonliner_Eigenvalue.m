function [n_eles,n_dofs,lambda_h_new] = DG_IGA_2D_Nonliner_Eigenvalue(Refinement,t )


datetime('now')

DIM = 2;

n_test_eigenvalues = 1;
 
rho = @(u) u^2;

a = 1;   
c = 0.1;

s = [0,0];

V_ext = @(x,y)  -1/norm([x, y] - s);
 
%% 2D Quantum Harmonic Oscillator problem
% V_ext = @(x,y) (x^2 + y^2)/2;
%%


x = [-a, -c, c,a];
y = [-a, -c, c,a];

nx = length(x) - 1;
ny = length(y) - 1;
n_subdomains = nx*ny;

mp.nx = nx;
mp.ny = ny;

ConPts = zeros(n_subdomains,2,2,DIM);
% The first parameter is the number of elements
% The second parameter is the number of basis functions in x-direction
% The third parameter is the number of basis functions in y-direction
% The fourth parameter is the number dimensions

for j=1:ny
    for i=1:nx
        e = i + (j-1)*nx;
        ConPts(e,:,:,1) = [x(i) x(i)  ; x(i+1) x(i+1)];
        ConPts(e,:,:,2) = [y(j) y(j+1); y(j)   y(j+1)];
    end
end

pu = 1; pv = 1; % The ultimate degree of B splines basis functions used is (pu+t).
 
weights =[1 1;1 1];
knotU   =[0 0  1 1];
knotV   =[0 0  1 1];
nu = 2;
nv = 2;

nurbs_original = cell(n_subdomains,1);

for e = 1:n_subdomains
    
    nurbs_original{e}.ConPts  = zeros(nu,nv,DIM);
    nurbs_original{e}.weights = weights;
    nurbs_original{e}.pu      = pu;
    nurbs_original{e}.pv      = pv;
    nurbs_original{e}.knotU   = knotU;
    nurbs_original{e}.knotV   = knotV;
    
    for i = 1:nu
        for j=1:nv
            for d = 1:DIM
                nurbs_original{e}.ConPts(i,j,d) = ConPts(e,i,j,d);
            end
        end
    end
end


for j=1:ny
    for i=1:nx
        e = i+(j-1)*nx;
        nurbs_original{e}.a = x(i+1) - x(i);
    end
end

for j=1:ny
    for i=1:nx
        e = i+(j-1)*nx;
        nurbs_original{e}.b = y(j+1) - y(j);
    end
end





%%
                
nurbs_refine_domains = cell(n_subdomains,1);

t_subdomains       = t*ones(n_subdomains,1);
Refine_subdomains  = Refinement*ones(n_subdomains,1);

Refine_subdomains(5) = Refinement + 2; % For the sub-domain containning the singular point

for e = 1:n_subdomains
  
    [knotU,knotV] =  updateKnotVectors_2D(nurbs_original{e}.knotU,nurbs_original{e}.knotV,nurbs_original{e}.a,nurbs_original{e}.b);%%  We want to have a regular mesh
    [knotU,knotV] =  IGADegreeElevSurface(knotU,knotV,t_subdomains(e)); % Now we should consider the ratio of lengths in x- and y-directions
    pu = nurbs_original{e}.pu + t_subdomains(e);     pv = nurbs_original{e}.pv + t_subdomains(e); 
    
    nurbs_refine_domains{e} = IGA_2D_Grid(knotU,knotV,pu,pv, Refine_subdomains(e));
    
    
end




%%
 
n_dofs = nurbs_refine_domains{1}.n_dofs_domains;

for e = 1:(n_subdomains-1)
    nurbs_refine_domains{e+1}.n_dofs_domains  = nurbs_refine_domains{e+1}.n_dofs_domains + n_dofs;
    n_dofs =  nurbs_refine_domains{e+1}.n_dofs_domains;
end


subdomain_dofs_idx = cell(n_subdomains,1);
row = 1:nurbs_refine_domains{1}.n_dofs;
subdomain_dofs_idx{1} = row;

for e = 2:n_subdomains
    subdomain_dofs_idx{e} = (1  + nurbs_refine_domains{e-1}.n_dofs_domains ): nurbs_refine_domains{e}.n_dofs_domains;
end

 %%
 
 A      = sparse(n_dofs,n_dofs);
 M      = sparse(n_dofs,n_dofs);

disp('Begin to assemble the matrices for DG-IGA method...')

tic

for e = 1:n_subdomains
    [stiffness_Matrix,Mass_Mat] = solve_laplace_A_F(nurbs_original{e},nurbs_refine_domains{e},V_ext);
    row = subdomain_dofs_idx{e};
    A(row,row)      = stiffness_Matrix;
%     M_vert(row,row) = Mass_Mat_Vext;
    M(row,row)      = Mass_Mat;
    
end



%%
uNoEs = nurbs_refine_domains{1}.uNoEs;

h = 1.0/uNoEs;

beta  = 50 ; % The penalty parameter in DG scheme


%%

for e = 2:n_subdomains % Update the the dofs information that is across the boundary of two sub-domains
    nurbs_refine_domains{e} = Update_Edge_DoFs( nurbs_refine_domains{e-1}.n_dofs_domains,nurbs_refine_domains{e}   );
    
end

n_eles = 0;

for e = 1:n_subdomains
    
    n_eles = n_eles + nurbs_refine_domains{e}.NoEs;
end



%%



S = sparse(n_dofs,n_dofs);
P = S;


%%

dg_2D_domains_idx = get2DsubDomain_info(mp);

left2right_eles = dg_2D_domains_idx.dg_2D_domains_idx_left2right;
bottom2top_eles = dg_2D_domains_idx.dg_2D_domains_idx_bottom2top;

for i = 1:ny
    for j = 1:(nx-1)
        left_ele_idx  = left2right_eles(i,j);   % In the i-th row of elements for left to right, the index of left element
        right_ele_idx = left2right_eles(i,j+1); % In the i-th row of elements for left to right, the index of right element
    [P_subdomains,S_subdomains] = IGA_DG_Left_Right_Edge_Assemble(nurbs_original{left_ele_idx},nurbs_original{right_ele_idx},...
                                                                  nurbs_refine_domains{left_ele_idx},nurbs_refine_domains{right_ele_idx},n_dofs);
    P = P + P_subdomains;
    S = S + S_subdomains;
    end
end
 
 
%%


for i =1:nx
    for j = 1:(ny-1)
        bottom_ele_idx = bottom2top_eles(i,j);
        top_ele_idx    = bottom2top_eles(i,j+1);
    [P_subdomains,S_subdomains] = IGA_DG_Top_Bottom_Edge_Assemble(nurbs_original{bottom_ele_idx},nurbs_original{top_ele_idx},...
                                                                  nurbs_refine_domains{bottom_ele_idx},nurbs_refine_domains{top_ele_idx},n_dofs);
    P = P + P_subdomains;
    S = S + S_subdomains;
    end
end


clear P_subdomains S_subdomains 


clear stiffness_Matrix   Mass_Mat_Vext   Mass_Mat





%%
Mat = A - S/2 - S'/2  + (beta/h)*P; % + M_vert;

clear A S P M_vert


elapsedTime = toc;

fprintf('The running time for assembling matrices is %f seconds\n', elapsedTime);

disp('The matrices for DG-IGA method have been assembled!')

%%


% disp(size(Mat))



%%
 be = get2DsubDomain_bnd_dof(nurbs_refine_domains,mp);



int_dofs_idx = 1:n_dofs;
int_dofs_idx(be) = [];


uh0             = randn(n_dofs,1);
uh0(be)         = 0;

%%  Compute the L2 norm over all sub-domains, then do nomalization
 
uh0_L2_norm =  0;

for e = 1:n_subdomains
    row            =  subdomain_dofs_idx{e};
    L2_norm_e      =  compute_L2_uh(nurbs_original{e},nurbs_refine_domains{e},uh0(row));
    uh0_L2_norm    =  uh0_L2_norm + L2_norm_e*L2_norm_e;
end
   uh0_L2_norm     =  sqrt(uh0_L2_norm);
      
%    disp('L2 norm of eigenfunction uh is')
%    disp(uh0_L2_norm)
   
   uh0 = uh0/uh0_L2_norm;
   
   %%


   


   uh_1 = zeros(n_dofs,1);
%    
   alpha = 0.1;
 
   format long
 
 lambda_old = 100;
 
 max_it = 1000;
 
 err = 1.0e-12;
 
 disp('=============================================================================')
 
 disp('Begin the SCF iteration')

 
 for k = 1:max_it
     
  M_rho = sparse(n_dofs,n_dofs);
  
   for e = 1:n_subdomains
         row = subdomain_dofs_idx{e};
         M_rho_sub_domain_e    =  compute_Nonlinear_Mass_Mat(nurbs_original{e},nurbs_refine_domains{e},uh0(row),rho);
         M_rho(row,row)        =  M_rho_sub_domain_e ;
   end
%   
   Mat_new   =  Mat + M_rho;
%   
  [lambda_h_new,eigenvector_int_h] = eigifp(Mat_new(int_dofs_idx ,int_dofs_idx),M(int_dofs_idx ,int_dofs_idx ),n_test_eigenvalues);
%   
  uh_1(int_dofs_idx) = eigenvector_int_h(:,1);
  
  uh0 = alpha*uh0 + (1-alpha)*uh_1;
  
  uh0_L2_norm = 0;

 for e = 1:n_subdomains
     row            =  subdomain_dofs_idx{e};
     uh0_L2_norm    =  uh0_L2_norm +  compute_L2_uh(nurbs_original{e},nurbs_refine_domains{e},uh0(row));
      
 end
  
%    lambda_h_new = lambda_h_new';
 
   uh0_L2_norm     =  sqrt(uh0_L2_norm);
   uh0 = uh0/uh0_L2_norm;


if(abs(lambda_h_new(1) - lambda_old)<=err )
    break;
end


fprintf('The %d-th iteration, eigenvalue = %1.15f \n',k, lambda_h_new );
 

lambda_old = lambda_h_new(1);
%     
 end
 
disp('=============================================================================')
 

 uh = zeros(n_dofs,n_test_eigenvalues);
 uh(int_dofs_idx,:) = eigenvector_int_h;
 

 
 
%% Save the eigenfunction function and nurbs info
       
       uh_filename   =  strcat('./EigenFunction_Data/p_',num2str(pu));
       uh_filename   =  strcat(uh_filename,'/uh_');
       uh_filename   =  strcat(uh_filename,num2str(Refinement));
       uh_filename   =  strcat(uh_filename,'.mat');
       save(uh_filename,'uh')
       
       nurbs          = nurbs_refine_domains; 
       nurbs_filename = strcat('./EigenFunction_Data/p_',num2str(pu));
       nurbs_filename = strcat(nurbs_filename,'/nurbs_');
       nurbs_filename = strcat(nurbs_filename,num2str(Refinement));
       nurbs_filename = strcat(nurbs_filename,'.mat');
       save(nurbs_filename,'nurbs')
 
 
 %%

 
 fprintf('The degree of B-splines is %d. No. of Elements = %d. No. of DOFs = %d \n',pu, n_eles, n_dofs );
 fprintf('The first eigenvalue for 2D nonlinear eigenvalue problem (Example 3) is %1.15f \n',lambda_h_new);

 





end