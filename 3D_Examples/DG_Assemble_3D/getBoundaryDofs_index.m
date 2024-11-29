function bnd_dof_index = getBoundaryDofs_index(nurbs_refine_domains,mp)
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

% back_ele_index

i = nx;
local_idx = 1;
for k = 1:nz
    for j=1:ny
        front_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;
        local_idx = local_idx + 1;
    end
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

% left_ele_index

j = ny;
local_idx = 1;
for k=1:nz
    for i = 1:nx
        right_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;
        local_idx = local_idx + 1;
    end
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

% bottom_ele_index

k = nz;
local_idx = 1;
for j = 1:ny
    for i = 1:nx
        top_ele_index(local_idx) = i + (j-1)*nx + (k-1)*nx*ny;
        local_idx = local_idx + 1;
    end
end

% top_ele_index

%% back and front faces
n_dofs_back = 0;
n_dofs_front = 0;

for i = 1:ny*nz
      n_dofs_back   = n_dofs_back   + nurbs_refine_domains{back_ele_index(i)}.n_dofs_back_face ;
      n_dofs_front  = n_dofs_front  + nurbs_refine_domains{front_ele_index(i)}.n_dofs_front_face;
end


back_face_dofs  =  zeros(1,n_dofs_back);
front_face_dofs =  zeros(1,n_dofs_front);

local_idx = 1;

for e = 1:ny*nz
    tmp_back  = nurbs_refine_domains{back_ele_index(e)}.back_face_dof;
    tmp_front = nurbs_refine_domains{front_ele_index(e)}.front_face_dof;
    for i = 1:length(tmp_back)
        back_face_dofs(local_idx)  = tmp_back(i);
        front_face_dofs(local_idx) = tmp_front(i);
        local_idx = local_idx + 1;
    end
end



%% left and right faces

n_dofs_left  = 0;
n_dofs_right = 0;

for i = 1:nx*nz
    n_dofs_left   = n_dofs_left   + nurbs_refine_domains{left_ele_index(i)}.n_dofs_left_face;
    n_dofs_right  = n_dofs_right  + nurbs_refine_domains{right_ele_index(i)}.n_dofs_right_face;
end


left_face_dofs  =  zeros(1,n_dofs_left);
right_face_dofs =  zeros(1,n_dofs_right);

local_idx = 1;

for e = 1:nx*nz
    tmp_left  = nurbs_refine_domains{left_ele_index(e)}.left_face_dof;
    tmp_right = nurbs_refine_domains{right_ele_index(e)}.right_face_dof;
    for i = 1:length(tmp_left)
        left_face_dofs(local_idx)  = tmp_left(i);
        right_face_dofs(local_idx) = tmp_right(i);
        local_idx = local_idx + 1;
    end
end



%% bottom and top faces

n_dofs_bottom = 0;
n_dofs_top = 0;

for i = 1:nx*ny
    n_dofs_bottom   = n_dofs_bottom   + nurbs_refine_domains{bottom_ele_index(i)}.n_dofs_bottom_face;
    n_dofs_top      = n_dofs_top      + nurbs_refine_domains{top_ele_index(i)}.n_dofs_top_face;
end

bottom_face_dofs  =  zeros(1,n_dofs_bottom);
top_face_dofs     =  zeros(1,n_dofs_top);

local_idx = 1;

for e = 1:nx*ny
    tmp_bottom  = nurbs_refine_domains{bottom_ele_index(e)}.bottom_face_dof;
    tmp_top     = nurbs_refine_domains{top_ele_index(e)}.top_face_dof;
    for i = 1:length(tmp_bottom)
        bottom_face_dofs(local_idx)  = tmp_bottom(i);
        top_face_dofs(local_idx) = tmp_top(i);
        local_idx = local_idx + 1;
    end
end
%%

bnd_dof_index  = [bottom_face_dofs,top_face_dofs,left_face_dofs,right_face_dofs,back_face_dofs,front_face_dofs];

bnd_dof_index  = unique(bnd_dof_index); 

end