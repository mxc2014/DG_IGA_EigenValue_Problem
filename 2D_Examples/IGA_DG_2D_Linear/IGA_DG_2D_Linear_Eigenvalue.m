function [n_eles,n_dofs,lambda_h]  = IGA_DG_2D_Linear_Eigenvalue(Refinement, Example, t)

DIM = 2;


fprintf('The test case is %s\n',Example)

%% For the domain with one singularity

if strcmp(Example, 'Example_1')

  V_ext = @(x,y) - 1.0/norm([x y]-[0 0]); 
 
  a = 1;    % The domain of Example 1 is Omega = [-1,1]^2
  c = 0.1;  % The patch that contains the origin is [-0.1,0.1]^2
 
  x = [-a,-c,c,a];
  y = x;
end



%%  For the domain with two singularities: Omega = [-1,1]^2

if strcmp(Example, 'Example_2')
    
a = 1;       c = 0.1;
x1 = -0.5;   x2 = 0.5;   
y1 =  0;     y2 = 0.;  
pnt_1 = [x1,y1];
pnt_2 = [x2,y2];
V_ext = @(x,y)  -1/norm([x, y] - pnt_1) -1/norm( [x, y] - pnt_2 ) ;

x = [-a, x1 - c, x1 + c, x2 - c, x2 + c, a];
y = [-a,y1-c,y1+c,a];


end


%% 2D Quantum Harmonic Oscillator problem

% V_ext = @(x,y) (x^2 + y^2)/2;
% a = 1; c = 1/2;
% x = [-a,-c,c,a];
% y = x;

%%


 

nx = length(x) - 1;
ny = length(y) - 1;
n_subdomains = nx*ny;



%%

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

%% The solution has one singularity

if strcmp(Example, 'Example_1')
    
 e = 5; Refine_subdomains(e) = Refinement + 2; % The patch that contains the singularity
 

end


%% The solution has two singularities

if strcmp(Example, 'Example_2')
    
e = 7;  Refine_subdomains(e) = Refinement + 2; % The patch that contains the singularity
e = 9;  Refine_subdomains(e) = Refinement + 2; % The patch that contains the singularity

end


%%

disp('Begin to assemble the matrices for DG-IGA method...')

tic

for e = 1:n_subdomains
  
    [knotU,knotV] =  updateKnotVectors_2D(nurbs_original{e}.knotU,nurbs_original{e}.knotV,nurbs_original{e}.a,nurbs_original{e}.b);%%  We want to have a regular mesh
    [knotU,knotV]  = IGADegreeElevSurface(knotU,knotV,t_subdomains(e)); % Now we should consider the ratio of lengths in x- and y-directions
    pu = nurbs_original{e}.pu + t_subdomains(e);     pv = nurbs_original{e}.pv + t_subdomains(e); 
    
    nurbs_refine_domains{e} = IGA_2D_Grid(knotU,knotV,pu,pv, Refine_subdomains(e));
    
    
end

%%

n_dofs = nurbs_refine_domains{1}.n_dofs_domains;

for e = 1:(n_subdomains-1)
    nurbs_refine_domains{e+1}.n_dofs_domains  = nurbs_refine_domains{e+1}.n_dofs_domains + n_dofs;
    n_dofs =  nurbs_refine_domains{e+1}.n_dofs_domains;
end


%%

for e = 2:n_subdomains % Update the the dofs information that is across the boundary of two sub-domains
    nurbs_refine_domains{e} = Update_Edge_DoFs( nurbs_refine_domains{e-1}.n_dofs_domains,nurbs_refine_domains{e}   );
    
end


n_eles = 0;

for e = 1:n_subdomains
    
    n_eles = n_eles + nurbs_refine_domains{e}.NoEs;
end
 

%%
% disp('Begin IGA+DG assemble')
tic

A      = sparse(n_dofs,n_dofs);
M      = sparse(n_dofs,n_dofs);

for e = 1:n_subdomains
    if e==1
      row = 1:nurbs_refine_domains{1}.n_dofs;
    else
      row = (1  + nurbs_refine_domains{e-1}.n_dofs_domains ): nurbs_refine_domains{e}.n_dofs_domains;
    end
    
    [stiffness_Matrix, Mass_Mat] = solve_laplace_A_F(nurbs_original{e},nurbs_refine_domains{e},V_ext);
     A(row,row)      = stiffness_Matrix;
     M(row,row)      = Mass_Mat;
end



%%
uNoEs = nurbs_refine_domains{1}.uNoEs;


h = 1.0/uNoEs;

beta  = 50;





%%



S = sparse(n_dofs,n_dofs);
P = S;


%%  The interface occurs at the left and right sub-domains

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
 
 
%%  The interface occurs at the lower and upper sub-domains


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



Mat = A - S/2 - S'/2  + (beta/h)*P ;


clear A S P M_vert




%%
 be = get2DsubDomain_bnd_dof(nurbs_refine_domains,mp);



int_dofs_idx = 1:n_dofs;
int_dofs_idx(be) = [];

n_eigenvalues = 4;
uh = zeros(n_dofs,n_eigenvalues);
 


Mat_int  = Mat(int_dofs_idx,int_dofs_idx);
M_int    = M(int_dofs_idx,int_dofs_idx);

elapsedTime = toc;

fprintf('The running time for assembling matrices is %f seconds\n', elapsedTime);

disp('The matrices for DG-IGA method have been assembled!')


clear Mat M

disp('===============================================================================')
disp('Begin to compute the eigenvalues ...')
tic
[lambda_h,uh_int] = eigifp( Mat_int , M_int ,n_eigenvalues); 
elapsedTime = toc;

fprintf('The running time for eigen solver is %f seconds \n', elapsedTime);



uh(int_dofs_idx,:) = uh_int;

fprintf('The degree of B-splines is p = %d. No. of Elements = %d. No. of DOFs = %d. \n',pu, n_eles, n_dofs);

fprintf('The first four eigenvalues are:\n');
format long
disp(lambda_h')
disp('===============================================================================')
 
%% Save the uh and the NURBS information

       uh_filename   =  strcat('./EigenFunction_Data/',Example,'/','p_',strcat(num2str(pu)));
       uh_filename   =  strcat( strcat(uh_filename,'/uh_'),num2str(Refinement));
       uh_filename   =  strcat(uh_filename,'.mat');
       
       
       nurbs_filename   =  strcat('./EigenFunction_Data/',Example,'/','p_',strcat(num2str(pu)));
       nurbs_filename   =  strcat( strcat(nurbs_filename,'/nurbs_'),num2str(Refinement));
       nurbs_filename   =  strcat(nurbs_filename,'.mat');
       
       nurbs = nurbs_refine_domains;
       save(uh_filename,'uh')
       save(nurbs_filename,'nurbs')
      

end


