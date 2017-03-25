function [iter, eigenvalue, k] = EigenSolver(a,D,Sigma_A,vSigma_F,h,k_guess,phi_guess,e1,e2)

mesh=[-a:h:a];
A=zeros(length(mesh),length(mesh)+2);
for i=1:length(mesh) %step2: compute elements of A
    A(i,i:i+2)=[-D/h^2 2*D/h^2+Sigma_A -D/h^2];
end
A=A(:,2:length(mesh)+1);
F=eye(length(mesh))*vSigma_F; %unecessary

%initial guesses
k_prev=k_guess; %k
phi_prev=phi_guess/norm(phi_guess,2); %step1: normalized phi vector
Q_prev=vSigma_F*phi_prev; %step3: initial fission source

%itital iteration
iter=1;
%
phi=A\(1/k_prev*Q);%phi vector
%
Q=vSigma_F*phi;%next fission source
k=k_prev*((h/2*Q(1)+h*(sum(Q)-Q(1)))/((h/2*Q_prev(1)+h*(sum(Q_prev)-Q_prev(1)))));

e1_check=abs((k-k_prev)/k);
e2_check=abs((norm(phi,2)-norm(phi_prev,2))/norm(phi,2));

while e1_check>e1||e2_check>e2
    iter=iter+1;
    %
    phi=A\(1/k*Q);
    %
    k_prev=k;
    Q_prev=Q;
    Q=vSigma_F*phi;%next fission source
    k=k_prev*((h/2*Q(1)+h*(sum(Q)-Q(1)))/((h/2*Q_prev(1)+h*(sum(Q_prev)-Q_prev(1)))));

    e1_check=abs((k-k_prev)/k);
    e2_check=abs((norm(phi,2)-norm(phi_prev,2))/norm(phi,2)); 
end

eigenvalue=1/k;

end

