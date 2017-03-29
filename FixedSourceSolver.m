function [anal_phi ,phi, mesh] = FixedSourceSolver(a,D,Sigma_A,S,h)
%Question 2

mesh=[-a:h:a];

A=zeros(length(mesh)-2,length(mesh)); 
for i=1:length(mesh)-2
    A(i,i:i+2)=[-D/h^2 2*D/h^2+Sigma_A -D/h^2]; 
end
A=A(:,2:length(mesh)-1);
b=(ones(1,length(mesh)-2)*S)';

dim=size(A,1); %Gauss Elem
Ab=[A b];
for i=2:dim %Forward Sub
   Ab(i-1,:)=Ab(i-1,:)./Ab(i-1,i-1); %Normalize
   Ab(i,:)=Ab(i,:)-Ab(i-1,:).*Ab(i,i-1);
end
A=Ab(1:dim,1:dim);
b=Ab(:,dim+1);

phi=zeros(dim,1);
phi(dim)=b(dim)/A(dim,dim);
for i=dim-1:-1:1 %Backward Sub
    phi(i)=b(i)-A(i,i+1)*phi(i+1);
end
phi=[0;phi;0]; %NUMERICAL (FINITE DIFF METHOD)

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