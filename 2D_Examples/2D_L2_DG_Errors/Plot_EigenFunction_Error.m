function  Plot_EigenFunction_Error(nurbs_original,nurbs_refine,uh,ref_nurbs,ref_uh)

% This script is used to compute the $L^2$ norm, $H^1$ semi-norm and $H^2$
% semi-norm  error between $u$ and the IGA solution $u_h$.



ConPts_o  = nurbs_original.ConPts ;
weights_o = nurbs_original.weights;
knotU_o    = nurbs_original.knotU; 
knotV_o    = nurbs_original.knotV; 
 
pu_o         =  nurbs_original.pu;
pv_o         =  nurbs_original.pv;


pu   =  nurbs_refine.pu;
pv   =  nurbs_refine.pv;

m = nurbs_refine.m;




 


knotU = nurbs_refine.Ubar;
knotV = nurbs_refine.Vbar;
 

 

 







m_pnts  = 200;
n_pnts  = m_pnts;



UBreaks  = nurbs_refine.UBreaks;  
VBreaks  = nurbs_refine.VBreaks;    

u_pnts = linspace(UBreaks(1),UBreaks(end),m_pnts+1);
v_pnts = linspace(VBreaks(1),VBreaks(end),n_pnts+1);

X = zeros(m_pnts+1,n_pnts+1);
Y = X;

 for j=1:(n_pnts+1)
    for i=1:(m_pnts+1)
     S= PointOnNurbsSurface(ConPts_o,weights_o,knotU_o,pu_o,u_pnts(i),knotV_o,pv_o,v_pnts(j));
     X(i,j) = S(1);
     Y(i,j) = S(2);
    end
 end

X = X';
Y = Y';



Z = zeros(m_pnts+1,n_pnts+1);

for j=1:(n_pnts+1)
     v_span = findspan(knotV,pv,v_pnts(j));
     Nv = bsplinebasis(knotV,pv,v_pnts(j));
     v_span_index = (v_span - pv): v_span;
     
    for i=1:(m_pnts+1)
        u_span = findspan(knotU,pu,u_pnts(i));
        Nu = bsplinebasis(knotU,pu,u_pnts(i));
       
        u_span_index = (u_span - pu): u_span;
        
        
        for j1=1:(pv+1)
            for i1=1:(pu+1)
                index = u_span_index(i1) + (v_span_index(j1)-1)*m ; 
                Z(i,j) = Z(i,j) + uh(index)*Nu(i1)*Nv(j1);
            end
        end
        
        
        ref_uh_at_uv   = reference_uh( ref_nurbs,ref_uh,u_pnts(i),v_pnts(j) );
        
        Z(i,j) = abs(   Z(i,j) -  ref_uh_at_uv );
        
        
    end
end

Z = Z';

% pcolor(X,Y,Z)
s = mesh(X,Y,Z);
% title('Contour plot of error function')
set(gca,'FontName','Times New Roman','FontSize',35,'LineWidth',1);

% colormap(jet)

 s.FaceColor = 'interp'; 

 s.MeshStyle = 'none';

% axis equal
xlim([-5,5])
ylim([-5,5])
colorbar



 hold on




     
end