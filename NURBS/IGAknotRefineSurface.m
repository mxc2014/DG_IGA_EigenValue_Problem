function [Ubar,Vbar,dof]=IGAknotRefineSurface(knotU,knotV,pu,pv,Refinement)

% Ubar=knotU;Vbar=knotV;

%% To unequally distribute the knots

% The interval is [-a,a];

a = 2;

xL = a;

x(1) = 0;
x(2) = a/2;
i=2;







% Ubar = [knotU(1:pu+1), inner_knot, knotU(end-pu:end) ];

%%


Ubar = knotU;

% Ubar = [0 0 0.1 0.2 1 1];

Vbar = knotV;

for i=1:Refinement
	 UBreks=unique(Ubar);VBreks=unique(Vbar);
     Xu=(UBreks(1:end-1)+UBreks(2:end))/2;
     Xv=(VBreks(1:end-1)+VBreks(2:end))/2;
     temp=[Ubar,Xu];Ubar=temp;Ubar=sort(Ubar); % Ubar=[Ubar,Xu];Ubar=sort(Ubar);
     temp=[Vbar,Xv];Vbar=temp;Vbar=sort(Vbar);% Vbar=[Vbar,Xv]; Vbar=sort(Vbar);
end

mu=length(Ubar)-pu-1;nv=length(Vbar)-pv-1;
dof=mu*nv;
