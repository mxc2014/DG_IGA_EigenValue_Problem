function nurbsInfo=Iga_3d_grid(knotU,pu,knotV,pv,knotW,pw,Refinement)

DIM=3;

% addpath('../NURBS/')

[Ubar,Vbar,Wbar,n_dofs]=IGAknotRefineVolume(knotU,pu,knotV,pv,knotW,pw,Refinement);

nurbsInfo.Ubar=Ubar;     % The knot vector in the u-direction.
nurbsInfo.Vbar=Vbar;     % The knot vector in the v-direction.
nurbsInfo.Wbar=Wbar;   % The knot vector in the w-direction.
nurbsInfo.n_dofs=n_dofs; % The number of DOFs in the refined grid.
nurbsInfo.pu = pu;
nurbsInfo.pv = pv;
nurbsInfo.pw = pw;

UBreaks=unique(Ubar);   % u 方向上节点向量中的断点.
VBreaks=unique(Vbar);   % v 方向上节点向量中的断点.
WBreaks=unique(Wbar); % w 方向上节点向量中的断点.

nurbsInfo.UBreaks=UBreaks;
nurbsInfo.VBreaks=VBreaks;
nurbsInfo.WBreaks=WBreaks;


N1=length(Ubar)-pu-1;  %  u 方向上基函数的个数。
N2=length(Vbar)-pv-1;  %  v 方向上基函数的个数。
N3=length(Wbar)-pw-1;%  w 方向上基函数的个数。

nurbsInfo.N1=N1;
nurbsInfo.N2=N2;
nurbsInfo.N3=N3;

uNoEs=length(UBreaks)-1;     %   u 方向上的区间数。
vNoEs=length(VBreaks)-1;      %  v 方向上的区间数。
wNoEs=length(WBreaks)-1;    %  w 方向上的区间数。
NoEs=uNoEs*vNoEs*wNoEs;   %  计算区域上的区间总数。

nurbsInfo.uNoEs=uNoEs;
nurbsInfo.vNoEs=vNoEs;
nurbsInfo.wNoEs=wNoEs;
nurbsInfo.NoEs=NoEs;

Eledof=(pu+1)*(pv+1)*(pw+1);     % 一个单元上的自由度总个数。
Element=zeros(NoEs,Eledof);         % 存储网格里每个单元上的自由度的编号.
knotSpanIndex=zeros(NoEs,DIM);  % 存储每个单元上的参数开始坐标的 knot span index. 
Coordinate=zeros(NoEs,2*DIM);     % 存储每个参数单元上的坐标起止值: $[u_i, u_{i+1}] *  [v_j, v_{j+1}]*  [w_k, w_{k+1}]$.

UKnotSpan = zeros(1,uNoEs);
VKnotSpan = zeros(1,vNoEs);
WKnotSpan = zeros(1,wNoEs);

 for k1=1:wNoEs   % 循环w方向上的全部单元.
     w_span=findspan(Wbar,pw,WBreaks(k1)); %当前单元e上的w方向上的节点张成区间的index,即$[w_k, w_{k+1}]$.
     WKnotSpan(k1) = w_span;
     for j1=1:vNoEs % 循环v方向上的全部单元.
     v_span=findspan(Vbar,pv,VBreaks(j1));      %当前单元e上的v方向上的节点张成区间的index,即$[v_j, v_{j+1}]$.
     VKnotSpan(j1) = v_span;
         for i1=1:uNoEs %循环 u方向上的全部单元.
     u_span=findspan(Ubar,pu,UBreaks(i1));      %当前单元e上的u方向上的节点张成区间的index,即$[u_i, u_{i+1}]$.
     UKnotSpan(i1) = u_span;
	
	e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
    % 注意，参数区域的网格单元的编号顺序是：　是先把最底部那一行单元从左到右，再从下到上来排列的.
	Coordinate(e,:)=[UBreaks(i1:i1+1),VBreaks(j1:j1+1),WBreaks(k1:k1+1)];% 存储当前单元的四个参数坐标.
    
    
    
	knotSpanIndex(e,:)=[u_span,v_span,w_span];
	
    local_index = 1;
    
    for k2=(w_span-pw):w_span
        for j2=(v_span-pv):v_span
            for i2=(u_span-pu):u_span
     global_index = i2+(j2-1)*N1+(k2-1)*N1*N2;
     Element(e,local_index) = global_index;
     local_index = local_index+1;
            end
        end
    end
     % 注意，这里全局自由度的编号是先让v和w方向的index固定，把u方向上的index变化，再让w方向上的index固定，让v方向的index变化，也就是:
     % N_{1,1,1}, N_{2,1,1}, ... , N_{m,1,1}; N_{1,2,1},N_{2,2,1},..., N_{m,2,1}; ...;
     % N_{1,n,1},N_{2,n,1},..., N_{m,n,1}.
end
    end
end




nurbsInfo.Element=Element;
nurbsInfo.Coordinate=Coordinate;
nurbsInfo.knotSpanIndex=knotSpanIndex;


Element_w_0 = zeros(uNoEs*vNoEs,(pu+1)*(pv+1)); % The dofs index in the elements lying on the surface w=0;
% These elements are the 2D elements, not the elements of the 3D mesh.


for j1=1:vNoEs
    for i1=1:uNoEs
        e = i1 + (j1-1)*uNoEs; % (i,j)处的单元的全局编号为e.
        u_span=findspan(Ubar,pu,UBreaks(i1));  
        v_span=findspan(Vbar,pv,VBreaks(j1));
        local_index = 1;
   for j2=(v_span - pv):v_span
       for i2=(u_span - pu):u_span
           global_index = i2 + (j2-1)*N1;
           Element_w_0(e,local_index)=global_index;
           local_index = local_index + 1;
       end
   end
   
    end
end

nurbsInfo.Element_w_0 = Element_w_0;



Element_w_1 = zeros(uNoEs*vNoEs,(pu+1)*(pv+1)); % The dofs index in the elements lying on the surface w=0;
% These elements are the 2D elements, not the elements of the 3D mesh.


for j1=1:vNoEs
    for i1=1:uNoEs
        e = i1 + (j1-1)*uNoEs; % (i,j)处的单元的全局编号为e.
        u_span=findspan(Ubar,pu,UBreaks(i1));  
        v_span=findspan(Vbar,pv,VBreaks(j1));
        local_index = 1;
   for j2=(v_span - pv):v_span
       for i2=(u_span - pu):u_span
           global_index = i2 + (j2-1)*N1;
           Element_w_0(e,local_index)=global_index;
           local_index = local_index + 1;
       end
   end
   
    end
end

nurbsInfo.Element_w_0 = Element_w_0;


u_n_ele_dofs = pu + 1;
v_n_ele_dofs = pv + 1;
w_n_ele_dofs = pw + 1;

%% Bottom face
bottom_face_dof =1:N1*N2;                              %  The face  w=0，即底面。
nurbsInfo.bottom_face_dof  = bottom_face_dof ;
nurbsInfo.n_dofs_bottom_face = length(bottom_face_dof); 



bottom_face_eles_2layers_dofs   = zeros(uNoEs*vNoEs, 2*u_n_ele_dofs*v_n_ele_dofs);


for j = 1:vNoEs
for i = 1:uNoEs
   ele_idx   = i + (j-1)*uNoEs;
   local_idx = 1;

   for k1 = 1:2 % The 2 layers dofs near bottom face
      for j1 = VKnotSpan(j)-pv:VKnotSpan(j)
         for i1 = UKnotSpan(i)-pu:UKnotSpan(i)
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        bottom_face_eles_2layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
   end
end
end

nurbsInfo.bottom_face_eles_2layers_dofs = bottom_face_eles_2layers_dofs;
   
   
bottom_face_eles_1_layers_dofs   = zeros(uNoEs*vNoEs, u_n_ele_dofs*v_n_ele_dofs);


for j = 1:vNoEs
for i = 1:uNoEs
   ele_idx   = i + (j-1)*uNoEs;
   local_idx = 1;

   for k1 = 1:1 % The 1 layer dofs on bottom face
      for j1 = VKnotSpan(j)-pv:VKnotSpan(j)
         for i1 = UKnotSpan(i)-pu:UKnotSpan(i)
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        bottom_face_eles_1_layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
   end
end
end

nurbsInfo.bottom_face_eles_1_layers_dofs = bottom_face_eles_1_layers_dofs;

   
%% Top face case
top_face_dof        = (1:N1*N2) + (N3-1)*N1*N2; % The face  w=1，即顶面。
nurbsInfo.top_face_dof  = top_face_dof;
nurbsInfo.n_dofs_top_face = length(top_face_dof); 

top_face_eles_2layers_dofs   = zeros(uNoEs*vNoEs,2*u_n_ele_dofs*v_n_ele_dofs);


for j = 1:vNoEs
for i = 1:uNoEs
   ele_idx   = i + (j-1)*uNoEs;
   local_idx = 1;

   for k1 = (N3-1):N3  % The 2 layers dofs near the top face
      for j1 = VKnotSpan(j)-pv:VKnotSpan(j)
         for i1 = UKnotSpan(i)-pu:UKnotSpan(i)
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        top_face_eles_2layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
end
   
  
end
end

nurbsInfo.top_face_eles_2layers_dofs = top_face_eles_2layers_dofs;



top_face_eles_1_layers_dofs   = zeros(uNoEs*vNoEs,u_n_ele_dofs*v_n_ele_dofs);


for j = 1:vNoEs
for i = 1:uNoEs
   ele_idx   = i + (j-1)*uNoEs;
   local_idx = 1;

   for k1 = N3:N3  % The 1 layer dofs on the top face
      for j1 = VKnotSpan(j)-pv:VKnotSpan(j)
         for i1 = UKnotSpan(i)-pu:UKnotSpan(i)
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        top_face_eles_1_layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
   end 
end
end

nurbsInfo.top_face_eles_1_layers_dofs = top_face_eles_1_layers_dofs;


%%  Front face: u = 1, i = N1 

front_face_dof    =  zeros(1,N2*N3);                  %  The face  u=1，即正面。
local_index = 1;
for k=1:N3
    for j=1:N2
        global_index = N1 + (j-1)*N1 +  (k-1)*N1*N2;
        front_face_dof(local_index) = global_index;
        local_index = local_index + 1;
    end
end

nurbsInfo.front_face_dof  =  front_face_dof ;

nurbsInfo.n_dofs_front_face = length(front_face_dof); 



front_face_eles_2layers_dofs   = zeros(vNoEs*wNoEs, 2*v_n_ele_dofs*w_n_ele_dofs);

for k = 1:wNoEs
for j = 1:vNoEs
   ele_idx   = j + (k-1)*vNoEs;
   local_idx = 1;

   for k1 = WKnotSpan(k)-pw:WKnotSpan(k)
      for j1 = VKnotSpan(j)-pv:VKnotSpan(j)
         for i1 = (N1-1):N1
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        front_face_eles_2layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
end
   
  
end
end

nurbsInfo.front_face_eles_2layers_dofs = front_face_eles_2layers_dofs;



front_face_eles_1_layers_dofs   = zeros(vNoEs*wNoEs, v_n_ele_dofs*w_n_ele_dofs);

for k = 1:wNoEs
for j = 1:vNoEs
   ele_idx   = j + (k-1)*vNoEs;
   local_idx = 1;

   for k1 = WKnotSpan(k)-pw:WKnotSpan(k)
      for j1 = VKnotSpan(j)-pv:VKnotSpan(j)
         for i1 = N1:N1
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        front_face_eles_1_layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
   end
end
end

nurbsInfo.front_face_eles_1_layers_dofs = front_face_eles_1_layers_dofs;




%% Back face case: u = 0, and i = 1

back_face_dof    =  zeros(1,N2*N3);                  %  The face  u=0，即背面。
local_index = 1;
for k=1:N3
    for j=1:N2
        global_index = 1 + (j-1)*N1 +  (k-1)*N1*N2;
        back_face_dof(local_index) = global_index;
        local_index = local_index + 1;
    end
end

nurbsInfo.back_face_dof = back_face_dof;
nurbsInfo.n_dofs_back_face = length(back_face_dof); 

back_face_eles_2layers_dofs   = zeros(vNoEs*wNoEs,2*v_n_ele_dofs*w_n_ele_dofs);

for k = 1:wNoEs
for j = 1:vNoEs
   ele_idx   = j + (k-1)*vNoEs;
   local_idx = 1;

   for k1 = WKnotSpan(k)-pw:WKnotSpan(k)
      for j1 = VKnotSpan(j)-pv:VKnotSpan(j)
         for i1 = 1:2
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        back_face_eles_2layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
   end
end
end


nurbsInfo.back_face_eles_2layers_dofs = back_face_eles_2layers_dofs;


back_face_eles_1_layers_dofs   = zeros(vNoEs*wNoEs,v_n_ele_dofs*w_n_ele_dofs);

for k = 1:wNoEs
for j = 1:vNoEs
   ele_idx   = j + (k-1)*vNoEs;
   local_idx = 1;

   for k1 = WKnotSpan(k)-pw:WKnotSpan(k)
      for j1 = VKnotSpan(j)-pv:VKnotSpan(j)
         for i1 = 1:1
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        back_face_eles_1_layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
   end
end
end


nurbsInfo.back_face_eles_1_layers_dofs = back_face_eles_1_layers_dofs;


 
%%  Left face case: v = 0, j = 1

        
left_face_dof = zeros(1,N1*N3);  %  The face  v=0，即左面。
local_index = 1;
for k=1:N3
    for i=1:N1
        global_index = i +  (k-1)*N1*N2;
        left_face_dof(local_index) = global_index;
        local_index = local_index + 1;
    end
end

nurbsInfo.left_face_dof = left_face_dof;
nurbsInfo.n_dofs_left_face = length(left_face_dof); 


left_face_eles_2layers_dofs   = zeros(uNoEs*wNoEs, 2*u_n_ele_dofs*w_n_ele_dofs);


for k = 1:wNoEs
for i = 1:uNoEs
   ele_idx   = i + (k-1)*uNoEs;
   local_idx = 1;

   for k1 = WKnotSpan(k)-pw:WKnotSpan(k)
      for j1 = 1:2
         for i1 = UKnotSpan(i)-pu:UKnotSpan(i)
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        left_face_eles_2layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
   end
end
end

nurbsInfo.left_face_eles_2layers_dofs = left_face_eles_2layers_dofs;


left_face_eles_1_layers_dofs   = zeros(uNoEs*wNoEs, u_n_ele_dofs*w_n_ele_dofs);


for k = 1:wNoEs
for i = 1:uNoEs
   ele_idx   = i + (k-1)*uNoEs;
   local_idx = 1;

   for k1 = WKnotSpan(k)-pw:WKnotSpan(k)
      for j1 = 1:1
         for i1 = UKnotSpan(i)-pu:UKnotSpan(i)
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        left_face_eles_1_layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
   end
end
end

nurbsInfo.left_face_eles_1_layers_dofs = left_face_eles_1_layers_dofs;



%% Right face case: v = 1, j = N2

right_face_dof = zeros(1,N1*N3);  %  The face  v=1，即右面。
local_index = 1;
for k=1:N3
    for i=1:N1
        global_index = i +(N2-1)*N1 +  (k-1)*N1*N2;
        right_face_dof(local_index) = global_index;
        local_index = local_index + 1;
    end
end


nurbsInfo.right_face_dof = right_face_dof;
nurbsInfo.n_dofs_right_face = length(right_face_dof); 


right_face_eles_2layers_dofs   = zeros(uNoEs*wNoEs, 2*u_n_ele_dofs*w_n_ele_dofs);


for k = 1:wNoEs
for i = 1:uNoEs
   ele_idx   = i + (k-1)*uNoEs;
   local_idx = 1;

   for k1 = WKnotSpan(k)-pw:WKnotSpan(k)
      for j1 = (N2-1):N2
         for i1 = UKnotSpan(i)-pu:UKnotSpan(i)
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        right_face_eles_2layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
   end
end
end

nurbsInfo.right_face_eles_2layers_dofs = right_face_eles_2layers_dofs;


right_face_eles_1_layers_dofs   = zeros(uNoEs*wNoEs, u_n_ele_dofs*w_n_ele_dofs);


for k = 1:wNoEs
for i = 1:uNoEs
   ele_idx   = i + (k-1)*uNoEs;
   local_idx = 1;

   for k1 = WKnotSpan(k)-pw:WKnotSpan(k)
      for j1 = N2:N2
         for i1 = UKnotSpan(i)-pu:UKnotSpan(i)
        global_idx = i1 + (j1-1)*N1 + (k1-1)*N1*N2;
        right_face_eles_1_layers_dofs(ele_idx, local_idx) = global_idx;
        local_idx = local_idx + 1;
end
end
   end
end
end

nurbsInfo.right_face_eles_1_layers_dofs = right_face_eles_1_layers_dofs;



%%





bnd_ele = zeros(2*uNoEs*vNoEs+2*uNoEs*wNoEs+2*vNoEs*wNoEs,1);

index  = 1;

     k1=1;
     for j1=1:vNoEs     % 循环v方向上的全部单元.
         for i1=1:uNoEs % 循环u方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
     k1=wNoEs;
     for j1=1:vNoEs     % 循环v方向上的全部单元.
         for i1=1:uNoEs % 循环u方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
     
     j1=1;
     for k1=1:wNoEs   % 循环w方向上的全部单元.
         for i1=1:uNoEs % 循环u方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
     j1=vNoEs;
     for k1=1:wNoEs   % 循环w方向上的全部单元.
         for i1=1:uNoEs % 循环u方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
     
     i1=1;
     for k1=1:wNoEs   % 循环w方向上的全部单元.
         for j1=1:vNoEs % 循环v方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
     i1=uNoEs;
     for k1=1:wNoEs   % 循环w方向上的全部单元.
         for j1=1:vNoEs % 循环v方向上的全部单元.
	     e=i1+(j1-1)*uNoEs + (k1-1)*uNoEs *vNoEs ;%=== 第    (i1,j1,k1)　号单元的全局编号。
         bnd_ele(index) = e;
         index = index + 1;
         end
     end
     
   bnd_ele = unique(bnd_ele);
   
   all_ele = [1:uNoEs*vNoEs*wNoEs]';
       
   all_ele(bnd_ele) = [];
   
   interior_ele = all_ele;
   
   
   nurbsInfo.bnd_ele = bnd_ele;
   
   nurbsInfo.interior_ele = interior_ele;
   
   
   n_dofs_w_0_u_0   =  N1*N2 + N2*N3 - N2; 
   dofs_w_0_u_0 = cell(2,1);
   dofs_w_0_u_0{1} = 1:N1*N2;
   dofs_w_0_u_0{2} = [1+( (1:N2) - 1)*N1,  (N1*N2+1):n_dofs_w_0_u_0];
   
   
  nurbsInfo.dofs_w_0_u_0 = dofs_w_0_u_0;
  nurbsInfo.n_dofs_w_0_u_0 = n_dofs_w_0_u_0;
     

end
