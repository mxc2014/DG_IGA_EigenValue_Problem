
%% DG-IGA method for the 3D Kohn-Sham equation with a Helium atom (The Example 5 of our paper)


addpath('../../NURBS/')
addpath('../DG_Assemble_3D/')
addpath('../IGA_Grid_3D_Data/')
addpath('./DFT/')
addpath('./L2_Projection_Boundary/')
  
t = 0:2;    % The (pu+t) is the ultimate degree of B splines basis functions
DIM=3;

%% The physical domain is [-10,10]^3

a = 10;
c = 1;  % The patch that contains the singularity is [-1,1]^3
x = [-a, -c,c,  a];
 
y = x;
z = x;
  
nu=2; nv=2; nw=2;
pu=1; pv=1; pw=1;

knotU=[0 0  1 1]; knotV=[0 0  1 1];knotW=[0 0 1 1];
nx = length(x) - 1;
ny = length(y) - 1;
nz = length(z) - 1;
n_subdomains = nx*ny*nz;

mp.nx = nx;
mp.ny = ny;
mp.nz = nz;
mp.n_subdomains = n_subdomains;


sub_domains_lengths = zeros(n_subdomains,3);% The lengths of all sub-domains

ConPts = zeros(n_subdomains,nu,nv,nw,DIM);

for k=1:nz
    for j = 1:ny
        for i = 1:nx
            
            e = i + (j-1)*nx + ( k-1 )*nx*ny;
            
            ConPts(e,:,:,1,1) = [x(i) x(i);x(i+1)  x(i+1) ]; % the x-coordinates in the bottom face
            ConPts(e,:,:,2,1) = [x(i) x(i);x(i+1)  x(i+1) ]; % the x-coordinates in the top face
            
            ConPts(e,:,:,1,2) = [y(j) y(j+1);y(j)  y(j+1) ]; % the y-coordinates in the bottom face
            ConPts(e,:,:,2,2) = [y(j) y(j+1);y(j)  y(j+1) ]; % the y-coordinates in the top face
            
            ConPts(e,:,:,1,DIM) = [z(k)   z(k);  z(k)    z(k) ]; % the z-coordinates in the bottom face
            ConPts(e,:,:,2,DIM) = [z(k+1) z(k+1);z(k+1)  z(k+1) ]; % the z-coordinates in the top face
            
            sub_domains_lengths(e,:) = [x(i+1) - x(i), y(j+1) - y(j),   z(k+1) - z(k) ] ;
            
            
        end
    end
end

nurbs_original = cell(n_subdomains,1);

weights = ones(nu,nv,nw);

for e = 1:n_subdomains
    nurbs_original{e}.ConPts  = zeros(nu,nv,nw,DIM);
    nurbs_original{e}.weights = weights;
    nurbs_original{e}.pu      = pu;
    nurbs_original{e}.pv      = pv;
    nurbs_original{e}.pw      = pw;
    nurbs_original{e}.knotU   = knotU;
    nurbs_original{e}.knotV   = knotV;
    nurbs_original{e}.knotW   = knotW;
    nurbs_original{e}.Lengths = sub_domains_lengths(e,:);
    
    % We assume that the lengths of sub-domain in x-,y-, and z-directions are: a, b, and c 
    for i = 1:nu
        for j = 1:nv
            for k = 1:nw
                for d = 1:DIM
                nurbs_original{e}.ConPts(i,j,k,d) = ConPts(e,i,j,k,d);     
                end
            end
        end
    end
            
end


%%


Refinement= 1:4;

n_refine=length(Refinement);
 
n_dofs=zeros(n_refine,1);

Lambda_h = zeros(n_refine,1);
n_eles   = zeros(n_refine,1);
DOFs     = zeros(n_refine,1);
Energy   = zeros(n_refine,1);


for j = 1:length(t)
    p  = pu + t(j);
    for i=1:n_refine

[Lambda_h(i),Energy(i), n_eles(i), DOFs(i)] = DG_IGA_3D_NonLinear_Eigenvalue(nurbs_original,  Refinement(i),t(j),mp );
 
save_EigenValue_Results(Lambda_h,Energy, n_eles, DOFs, p );

end
end






