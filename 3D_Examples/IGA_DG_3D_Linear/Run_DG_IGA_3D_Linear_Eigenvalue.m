% DG-IGA method for the 3D linear eigenvalue problem with a singular potential (Example 4 of our paper)

addpath('../../NURBS/')
addpath('../DG_Assemble_3D/')
addpath('../IGA_Grid_3D_Data/')

  
t= 0:2;    % The (pu+t) is the ultimate degree of B splines basis functions

DIM=3;


  
%% Case I: rectangle domain
nu=2; nv=2; nw=2;
pu=1; pv=1; pw=1;

knotU=[0 0  1 1]; knotV=[0 0  1 1];knotW=[0 0 1 1];


 a = 2;    % The domain is [-2,2]^3
 c = 1/2;  % The patch that contains the singularity is [-0.5,0.5]^3


x = [-a, -c,c,   a];
y = [-a, -c,c,   a];
z = [-a, -c,c,   a];



nx = length(x) - 1;
ny = length(y) - 1;
nz = length(z) - 1;

mp.nx = nx;
mp.ny = ny;
mp.nz = nz;

n_subdomains = nx*ny*nz;

mp.n_subdomains = n_subdomains;

ConPts = zeros(n_subdomains,nu,nv,nw,DIM);

sub_domains_lengths = zeros(n_subdomains,3);% The lengths of all sub-domains

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


Refinement= 0:4;

n_refine=length(Refinement);
 
Lambda_h = zeros(n_refine,1);
n_eles   = zeros(n_refine,1);
DOFs     = zeros(n_refine,1);

for j = 1:length(t)
    p = pu + t(j);
    
    for i=1:size(Refinement,2)
        [Lambda_h(i), n_eles(i), DOFs(i)] = DG_IGA_3D_Linear_Eigenvalue(nurbs_original, Refinement(i),t(j),mp);
    end

    
save_EigenValue_Results(Lambda_h,n_eles, DOFs, p );

end

%%





