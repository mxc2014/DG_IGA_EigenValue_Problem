function L2_err = Compute_L2_Projection_Boundary_Error(nurbs_original,nurbs_refine_domains,mp,g,gh)

 

nx = mp.nx;
ny = mp.ny;
nz = mp.nz;

L2_err = 0;

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
    L2_err_ele = Compute_L2_Projection_Boundary_Bottom_Error(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, gh );
    L2_err = L2_err + L2_err_ele^2;
end



% top_ele_index

k = nz;
local_idx = 1;
for j = 1:ny
    for i = 1:nx
        top_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;  % top_ele_index
        local_idx = local_idx + 1;
    end
end

for e = 1:(ny*nz)
    ele_idx = top_ele_index(e);
    L2_err_ele = Compute_L2_Projection_Boundary_Top_Error(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, gh );
    L2_err = L2_err + L2_err_ele^2;
end


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
    L2_err_ele = Compute_L2_Projection_Boundary_Left_Error(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, gh );
    L2_err = L2_err + L2_err_ele^2;
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
    L2_err_ele = Compute_L2_Projection_Boundary_Right_Error(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, gh );
    L2_err = L2_err + L2_err_ele^2;
end




% right_ele_index



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
      L2_err_ele = Compute_L2_Projection_Boundary_Back_Error(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, gh );
      L2_err = L2_err + L2_err_ele^2;
end




% back_ele_index

i = nx;
local_idx = 1;
for k = 1:nz
    for j=1:ny
        front_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;  % front_ele_index
        local_idx = local_idx + 1;
    end
end

for e = 1:(ny*nz)
     ele_idx = front_ele_index(e);
     L2_err_ele = Compute_L2_Projection_Boundary_Front_Error(nurbs_original{ele_idx}, nurbs_refine_domains{ele_idx},g, gh );
     L2_err = L2_err + L2_err_ele^2;
end




%%

 L2_err = sqrt(L2_err);
 
%  disp(L2_err)

end



