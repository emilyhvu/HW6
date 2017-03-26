function [eigenvalue,k,iter] = EigenSolver(a,D,Sigma_A,vSigma_F,h,e1,e2)
%Question 4 

mesh=[-a:h:a];
A=zeros(length(mesh)-2,length(mesh)); 
for i=1:length(mesh)-2
    A(i,i:i+2)=[-D/h^2 2*D/h^2+Sigma_A -D/h^2]; 
end
A=A(:,2:length(mesh)-1);

k_guess=0.755; %initial k value
phi_guess=ones(length(mesh)-2,1); %initial phi vector

%0th iteration
k=k_guess; %k
phi=phi_guess/norm(phi_guess,2); %step1: normalized phi vector
Q=vSigma_F*phi; %step3: initial fission source

%itital iteration
iter=1;
%
phi_next=myGaussSeidel(A,length(mesh)-2,e2,1/k*Q,phi);%phi vector
%
Q_next=vSigma_F*phi_next;%next fission source
k_next=k*((h/2*Q_next(1)+h*(sum(Q_next)-Q_next(1)))/((h/2*Q(1)+h*(sum(Q)-Q(1)))));
phi=phi_next;

e1_check=abs((k_next-k)/k_next);
%e2_check=abs((norm(phi,2)-norm(phi_prev,2))/norm(phi,2));

while e1_check>e1 %||e2_check>e2
    iter=iter+1;
    %
    phi_next=myGaussSeidel(A,length(mesh)-2,e2,1/k_next*Q_next,phi);
    %
    k=k_next;
    Q=Q_next;
    Q_next=vSigma_F*phi_next;%next fission source
    k_next=k*((h/2*Q_next(1)+h*(sum(Q_next)-Q_next(1)))/((h/2*Q(1)+h*(sum(Q)-Q(1)))));
    phi=phi_next;
    
    e1_check=abs((k_next-k)/k_next);
end
phi=[0;phi;0];
eigenvalue=1/k;

plot(mesh,phi,'b')
hold all
title('Eigenvector Solution')
ylabel('phi(x)')
xlabel('x')

end

