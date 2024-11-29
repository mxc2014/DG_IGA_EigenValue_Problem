function [P,S]= IGA_DG_Left_Right_Face_Assemble(nurbs_original_1,nurbs_original_2,nurbs_refine_1,nurbs_refine_2, n_dofs)
%%
ConPts_o_1  = nurbs_original_1.ConPts;
weights_o_1 = nurbs_original_1.weights;
knotU_o_1   = nurbs_original_1.knotU;
knotV_o_1   = nurbs_original_1.knotV;
knotW_o_1   = nurbs_original_1.knotW;
pu_o_1      = nurbs_original_1.pu;
pv_o_1      = nurbs_original_1.pv;
pw_o_1      = nurbs_original_1.pw;

knotU_1      = nurbs_refine_1.Ubar; 
knotV_1      = nurbs_refine_1.Vbar; 
knotW_1      = nurbs_refine_1.Wbar;
pu_1         = nurbs_refine_1.pu;
pv_1         = nurbs_refine_1.pv;
pw_1         = nurbs_refine_1.pw;
uNoEs_1      = nurbs_refine_1.uNoEs;
vNoEs_1      = nurbs_refine_1.vNoEs;
wNoEs_1      = nurbs_refine_1.wNoEs;
%%


%%
ConPts_o_2  = nurbs_original_2.ConPts;
weights_o_2 = nurbs_original_2.weights;
knotU_o_2   = nurbs_original_2.knotU;
knotV_o_2   = nurbs_original_2.knotV;
knotW_o_2   = nurbs_original_2.knotW;
pu_o_2      = nurbs_original_2.pu;
pv_o_2      = nurbs_original_2.pv;
pw_o_2      = nurbs_original_2.pw;

knotU_2      = nurbs_refine_2.Ubar; 
knotV_2      = nurbs_refine_2.Vbar; 
knotW_2      = nurbs_refine_2.Wbar;
pu_2         = nurbs_refine_2.pu;
pv_2         = nurbs_refine_2.pv;
pw_2         = nurbs_refine_2.pw;
uNoEs_2      = nurbs_refine_2.uNoEs;
vNoEs_2      = nurbs_refine_2.vNoEs;
wNoEs_2      = nurbs_refine_2.wNoEs;
%%

UBreaks_1    = nurbs_refine_1.UBreaks; 
WBreaks_1    = nurbs_refine_1.WBreaks; 

UBreaks_2    = nurbs_refine_2.UBreaks; 
WBreaks_2    = nurbs_refine_2.WBreaks; 

Ubreaks      = [UBreaks_1,UBreaks_2];
Ubreaks      = unique(Ubreaks);

u_edge_element_idx_1 =  Merge_edge_elements_idx(UBreaks_1,Ubreaks);
u_edge_element_idx_2 =  Merge_edge_elements_idx(UBreaks_2,Ubreaks);



NoEs_u       = length(Ubreaks)-1;

Wbreaks      = [WBreaks_1,WBreaks_2];
Wbreaks      = unique(Wbreaks);
w_edge_element_idx_1 =  Merge_edge_elements_idx(WBreaks_1,Wbreaks);
w_edge_element_idx_2 =  Merge_edge_elements_idx(WBreaks_2,Wbreaks);

NoEs_w       = length(Wbreaks)-1;

NoEs         = NoEs_u*NoEs_w;


ele_idx_on_edge_1 = zeros(NoEs,1);
ele_idx_on_edge_2 = zeros(NoEs,1);

merge_coordinate  = zeros(NoEs,4);


for j = 1:NoEs_w
for i = 1:NoEs_u
 merge_ele_idx = i + (j-1)*NoEs_u;  % The element index in the finner mesh
 ele_idx_1     = u_edge_element_idx_1(i) + ( w_edge_element_idx_1(j)  -1 )*uNoEs_1; % The element index in the coarser mesh from patch 1
 ele_idx_2     = u_edge_element_idx_2(i) + ( w_edge_element_idx_2(j)  -1 )*uNoEs_2; % The element index in the coarser mesh from patch 2
 ele_idx_on_edge_1( merge_ele_idx ) = ele_idx_1 ;
 ele_idx_on_edge_2( merge_ele_idx ) = ele_idx_2 ;
 merge_coordinate(merge_ele_idx,:)  = [Ubreaks(i:i+1), Wbreaks(j:j+1)];
end
end


% ele_idx_on_edge_1
% ele_idx_on_edge_2

% merge_coordinate


right_face_eles_2layers_dofs_plus       = nurbs_refine_1.right_face_eles_2layers_dofs;
left_face_eles_2layers_dofs_minus       = nurbs_refine_2.left_face_eles_2layers_dofs;

% nurbs_refine_1.Element
% right_face_eles_2layers_dofs_plus
% disp('=====================')
% left_face_eles_2layers_dofs_minus


np      = max(pu_1,pu_2) +1 ;
[gp,gw] = grule(np);

%%  The right face of the left volume
v_right = 1.0;
Vders_plus  = bspbasisDers(knotV_1,pv_1,v_right,1); % For lower volume, the top face:    w = 1
Nv_plus     = Vders_plus(1,end-1:end);  %%  This needs special attention
DNv_plus    = Vders_plus(2,end-1:end);
%%
v_left = 0.0;
Vders_minus = bspbasisDers(knotV_2,pv_2,v_left,1); % For upper volume, the bottom face: w = 0
Nv_minus    = Vders_minus(1,1:2);
DNv_minus   = Vders_minus(2,1:2);
%%


P = sparse(n_dofs,n_dofs);
S = P;


normal = [0;1;0];  % The unit outward normal is rightward

for e = 1:NoEs
    ue = merge_coordinate(e,1:2);  we = merge_coordinate(e,3:4);
    
%     disp('ue = ')
%     disp(ue)
%     disp('ve = ')
%     disp(ve)
    
    ua = ue(1);  ub = ue(2);  uJ = (ub - ua)/2;
    wa = we(1);  wb = we(2);  wJ = (wb - wa)/2;
    face_idx_plus  = ele_idx_on_edge_1(e);
    face_idx_minus = ele_idx_on_edge_2(e);
    
    face_dofs_plus  = right_face_eles_2layers_dofs_plus(face_idx_plus,:);
    face_dofs_minus = left_face_eles_2layers_dofs_minus(face_idx_minus,:);
    face_dofs       = [face_dofs_plus,face_dofs_minus];
    n_face_dofs     = length(face_dofs);
    
    face_jump_Ae    = zeros(n_face_dofs,n_face_dofs);
    face_average_Ae = zeros(n_face_dofs,n_face_dofs);


   for i = 1:np
  u = ( (ub - ua )*gp(i) + ua + ub  )/2;
  Uders_plus = bspbasisDers(knotU_1,pu_1,u,1);
  Nu_plus    = Uders_plus(1,:);   DNu_plus = Uders_plus(2,:);

  Uders_minus = bspbasisDers(knotU_2,pu_2,u,1);
  Nu_minus    = Uders_minus(1,:);   DNu_minus = Uders_minus(2,:);
  
for j = 1:np
  w = ( (wb - wa )*gp(j) + wa + wb  )/2;

[~,DF_plus]  = NurbsVolume(ConPts_o_1,weights_o_1,knotU_o_1,pu_o_1,u,knotV_o_1,pv_o_1,v_right,knotW_o_1,pw_o_1,w);
[~,DF_minus] = NurbsVolume(ConPts_o_2,weights_o_2,knotU_o_2,pu_o_2,u,knotV_o_2,pv_o_2,v_left ,knotW_o_2,pw_o_2,w);

A = DF_plus(2,1)*DF_plus(3,2) - DF_plus(2,2)*DF_plus(3,1);
B = DF_plus(3,1)*DF_plus(1,2) - DF_plus(1,1)*DF_plus(3,2);
C = DF_plus(1,1)*DF_plus(2,2) - DF_plus(1,2)*DF_plus(2,1);
J_plus = sqrt(A*A+B*B+C*C)*uJ*wJ*gw(i)*gw(j);





  Wders_plus = bspbasisDers(knotW_1,pw_1,w,1);
  Nw_plus    = Wders_plus(1,:);   DNw_plus = Wders_plus(2,:);

  Wders_minus = bspbasisDers(knotW_2,pw_2,w,1);
  Nw_minus    = Wders_minus(1,:);   DNw_minus = Wders_minus(2,:);

  vw_basis_plus = Nv_plus'* Nw_plus ; vw_basis_plus = vw_basis_plus(:)';
  basis_plus    = Nu_plus'*vw_basis_plus;    basis_plus = basis_plus(:);
  vw_basis_minus= Nv_minus'* Nw_minus; vw_basis_minus = vw_basis_minus(:)';
  basis_minus   = Nu_minus'*vw_basis_minus;  basis_minus = basis_minus(:);

  face_jump_basis = [basis_plus;-basis_minus];
  
  face_jump_Ae =  face_jump_Ae  +  face_jump_basis*face_jump_basis'*J_plus;

  Du_vw_basis_plus = DNu_plus'*vw_basis_plus;   Du_vw_basis_plus = Du_vw_basis_plus(:); 
  Du_vw_basis_minus= DNu_minus'*vw_basis_minus; Du_vw_basis_minus = Du_vw_basis_minus(:); 
  
  

  Dv_w_basis_plus  = DNv_plus'* Nw_plus;  Dv_w_basis_plus = Dv_w_basis_plus(:)'; 
  Dv_uw_basis_plus = Nu_plus'*Dv_w_basis_plus;  Dv_uw_basis_plus = Dv_uw_basis_plus(:);
  Dv_w_basis_minus  = DNv_minus'* Nw_minus;     Dv_w_basis_minus = Dv_w_basis_minus(:)'; 
  Dv_uw_basis_minus = Nu_minus'*Dv_w_basis_minus;  Dv_uw_basis_minus = Dv_uw_basis_minus(:);

  Dw_v_basis_plus   = Nv_plus'* DNw_plus; Dw_v_basis_plus = Dw_v_basis_plus(:)'; 
  Dw_uv_basis_plus  = Nu_plus'*Dw_v_basis_plus;  Dw_uv_basis_plus = Dw_uv_basis_plus(:);
  Dw_v_basis_minus  = Nv_minus'* DNw_minus; Dw_v_basis_minus = Dw_v_basis_minus(:)'; 
  Dw_uv_basis_minus = Nu_minus'*Dw_v_basis_minus;  Dw_uv_basis_minus = Dw_uv_basis_minus(:);

  basis_grad_plus  = [Du_vw_basis_plus,Dv_uw_basis_plus,Dw_uv_basis_plus]/DF_plus;
  basis_grad_minus = [Du_vw_basis_minus,Dv_uw_basis_minus,Dw_uv_basis_minus]/DF_minus;
  face_average_basis = [basis_grad_plus;basis_grad_minus]*normal/2;
  face_average_basis = face_average_basis'; % Now it is a row vector.
  
  face_average_Ae  = face_average_Ae + face_jump_basis* face_average_basis*J_plus;
end
   end

   P(face_dofs,face_dofs) = P(face_dofs,face_dofs) +  face_jump_Ae;
   S(face_dofs,face_dofs) = S(face_dofs,face_dofs) +  face_average_Ae;
end


end