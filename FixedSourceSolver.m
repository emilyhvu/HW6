function [anal_phi ,phi, mesh] = FixedSourceSolver(a,D,Sigma_A,S,h)
%Question 2

mesh=[-a:h:a];

%fake_mesh=[]; %unecessary for this code but would be used for non-uniform spacing
%for i=1:length(mesh)-1 
%    fake_mesh=[fake_mesh mesh(i)+((mesh(i+1)-mesh(i)))/2];
%end
%fake_mesh=[mesh(1) fake_mesh mesh(length(mesh))];

A=zeros(length(mesh)-2,length(mesh)); 
for i=1:length(mesh)-2
    A(i,i:i+2)=[-D/h^2 2*D/h^2+Sigma_A -D/h^2]; 
end
A=A(:,2:length(mesh)-1);
b=(ones(1,length(mesh)-2)*S)';

phi=A\b; %NUMERICAL (FINITE DIFF METHOD)
phi=[0;phi;0];

plot(mesh,phi,'r+')
title('Fixed-Source Diffusion Equation')
xlabel('x')
ylabel('phi(x)')
hold all

f=@computeAnalytical; %ANALYTICAL (SOLVED ODE BY HAND)
anal_phi=[];
for i=1:length(mesh)
    anal_phi=[anal_phi f(mesh(i))];
end
anal_phi=anal_phi';

plot(mesh,anal_phi, 'bl-')
legend('Numerical', 'Analytical')

end