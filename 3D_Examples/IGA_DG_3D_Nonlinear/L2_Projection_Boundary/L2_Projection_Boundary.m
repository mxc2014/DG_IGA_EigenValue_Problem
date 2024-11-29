function gh = L2_Projection_Boundary(nurbs_original,nurbs_refine_domains,mp,g,n_dofs)

M     = sparse(n_dofs,n_dofs);
rhs   = zeros(n_dofs,1);

nx = mp.nx;
ny = mp.ny;
nz = mp.nz;



%% back and front faces

back_ele_index  = zeros(ny*nz,1);
front_ele_index = zeros(ny*nz,1);

i = 1;
local_idx = 1;
for k = 1:nz
    for j=1:ny
        back_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;
        local_idx = local_idx + 1;
    end
end


for e = 1:(ny*nz)
    ele_idx = back_ele_index(e);
    [M_back,rhs_back]= L2_Projection_Boundary_Back(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, n_dofs);
    M   = M   + M_back;
    rhs = rhs + rhs_back;
end


% back_ele_index

i = nx;
local_idx = 1;
for k = 1:nz
    for j=1:ny
        front_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;
        local_idx = local_idx + 1;
    end
end

for e = 1:(ny*nz)
    ele_idx = front_ele_index(e);
    [M_front,rhs_front]= L2_Projection_Boundary_Front(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, n_dofs);
    M   = M   + M_front;
    rhs = rhs + rhs_front;
end

% front_ele_index

%%  left and right faces


left_ele_index    = zeros(nx*nz,1);
right_ele_index   = zeros(nx*nz,1);

j = 1;
local_idx = 1;
for k=1:nz
    for i = 1:nx
        left_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;
        local_idx = local_idx + 1;
    end
end

for e = 1:(ny*nz)
    ele_idx = left_ele_index(e);
    [M_left,rhs_left]= L2_Projection_Boundary_Left(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, n_dofs);
    M   = M   + M_left;
    rhs = rhs + rhs_left;
end

% left_ele_index

j = ny;
local_idx = 1;
for k=1:nz
    for i = 1:nx
        right_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;
        local_idx = local_idx + 1;
    end
end

for e = 1:(ny*nz)
    ele_idx = right_ele_index(e);
 
    [M_right,rhs_right]= L2_Projection_Boundary_Right(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, n_dofs);
    M   = M   + M_right;
    rhs = rhs + rhs_right;
end

% right_ele_index

%% bottom and top faces        


bottom_ele_index  = zeros(nx*ny,1);
top_ele_index     = zeros(nx*ny,1);

k = 1;
local_idx = 1;
for j = 1:ny
    for i = 1:nx
        bottom_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;
        local_idx = local_idx + 1;
    end
end

for e = 1:(ny*nz)
    ele_idx = bottom_ele_index(e);
 
    [M_bottom,rhs_bottom]= L2_Projection_Boundary_Bottom(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, n_dofs);
    M   = M   + M_bottom;
    rhs = rhs + rhs_bottom;
end

% bottom_ele_index

k = nz;
local_idx = 1;
for j = 1:ny
    for i = 1:nx
        top_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;
        local_idx = local_idx + 1;
    end
end

for e = 1:(ny*nz)
    ele_idx = top_ele_index(e);
 
    [M_top,rhs_top]= L2_Projection_Boundary_Top(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, n_dofs);
    M   = M   + M_top;
    rhs = rhs + rhs_top;
end

% top_ele_index
%%
bnd_dof_index = getBoundaryDofs_index(nurbs_refine_domains,mp);
gh_bnd = M(bnd_dof_index,bnd_dof_index)\rhs(bnd_dof_index);

% gh_bnd

gh = zeros(n_dofs,1);
gh(bnd_dof_index) = gh_bnd;


end



