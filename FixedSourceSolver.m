function [anal_phi ,phi, mesh] = FixedSourceSolver(a,D,Sigma_A,S,h)

mesh=[-a:h:a];
BD=0;

%fake_mesh=[]; %unecessary for this code but would be used for non-uniform spacing
%for i=1:length(mesh)-1 
%    fake_mesh=[fake_mesh mesh(i)+((mesh(i+1)-mesh(i)))/2];
%end
%fake_mesh=[mesh(1) fake_mesh mesh(length(mesh))];

A=zeros(length(mesh),length(mesh)+2);
for i=1:length(mesh)
    A(i,i:i+2)=[-D/h^2 2*D/h^2+Sigma_A -D/h^2]; 
end
A=A(:,2:length(mesh)+1);
b=(ones(1,length(mesh))*S)';

phi=A\b;

plot(mesh,phi,'bl')
title('Fixed-Source Diffusion Equation')
xlabel('x')
ylabel('phi(x)')
hold all

f=@computeAnalytical;
anal_phi=[];
for i=1:length(mesh)
    anal_phi=[anal_phi f(mesh(i))];
end
anal_phi=anal_phi';

plot(mesh,anal_phi, 'r+')
legend('Numerical', 'Analytical')

end