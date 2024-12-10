function [tv,yv]=EulerSyst_4(f,tspan,y0,n,theta, A, h, g)
%Solve y'=f(t,y) with inital condition y(tspan(1))=y0
%up to time tspan(2).
%y0 should be a column vector, 
%f is assumed to return a column vector of the same length.

Dt=(tspan(2)-tspan(1))/n; 
tv=linspace(tspan(1), tspan(2), n + 1); %tv a row vector.
p=length(y0); 
yv=zeros(p,n+1);
yv(:,1)=y0;

M = eye(size(A)) - Dt * theta * A;

for ii=1:n
    t_next = tv(ii+1);
    
    g_next = g(t_next);

    rhs = yv(:,ii) + Dt * ((1 - theta)*f(tv(ii),yv(:,ii)) + theta * g_next);
    
    yv(:, ii+1) = M \ rhs; 
end
