function [GS, GS_iter] = myGaussSeidel(n,error,tol)

A=zeros(n,n+2); %A matrix made of zeros with n rows and n+2 columns
for i=1:n %for every row in A
    A(i,i:i+2)=[-1 3 -1]; %replace corresponding indices with recurring coefficients
end
A=A(:,2:n+1); %get rid of first and last columns
b=100*ones(n,1); %b vector with n 100s
x_0=zeros(n,1); %initial guess of all zeros
    
%gauss seidel method
LplusD=A;
for i=1:n-1 %creates sum of lower and diagonal matrix by turning all upper elements to 0
        LplusD(i,i+1)=0;
end

U=A-LplusD; %creates upper matrix by subtracting (L+D) from A
P=LplusD\(-U); %multiply -U and L+D

GS=P*x_0+LplusD\b; %first iteration with backslash to take inverse of L+D
GS_iter=1;
x=A\b; %actual solution
e_vector=GS-x; 
e=norm(e_vector,2); %normalizing error vector

while e>error %tolerance
    GS_prev=GS; %store previous vector
    GS=P*GS_prev+LplusD\b; %uses last guess
    GS_iter=GS_iter+1; %counts iterations
    
    e_vector=GS-GS_prev; %error vector
    if strcmpi(tol,'absolute')
        e=norm(e_vector,2); %normalize to get absolute error 
    elseif strcmpi(tol,'relative')
        e=norm(e_vector,2)/norm(GS,2); %finds relative error
    end
end
end