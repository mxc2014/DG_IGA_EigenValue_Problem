function ele_info = get3DElements_index(nx,ny,nz)

bottom2top_eles = zeros(nz,nx*ny);
back2front_eles = zeros(nx,ny*nz);
left2right_eles = zeros(ny,nx*nz);


for k = 1:nz
    local_idx = 1;
    for j = 1:ny
        for i = 1:nx
          e = i + (j-1)*nx + (k-1)*nx*ny;
          bottom2top_eles(k,local_idx) = e;
          local_idx = local_idx + 1;
        end
    end
end

for i = 1:nx
    local_idx = 1;
    for k = 1:nz
        for j = 1:ny
          e = i + (j-1)*nx + (k-1)*nx*ny;
          back2front_eles(i,local_idx) = e;
          local_idx = local_idx + 1;
        end
    end
end


for j = 1:ny
    local_idx = 1;
    for k = 1:nz
        for i = 1:nx
          e = i + (j-1)*nx + (k-1)*nx*ny;
          left2right_eles(j,local_idx) = e;
          local_idx = local_idx + 1;
        end
    end
end

ele_info.bottom2top_eles = bottom2top_eles;
ele_info.back2front_eles = back2front_eles;
ele_info.left2right_eles = left2right_eles;



    


end

