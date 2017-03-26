function [GS] = myGaussSeidel(A,n,e2,b,x_0)
%Question 4

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
e=norm(e_vector,2)/norm(GS,2); %normalizing error vector

while e>e2 %tolerance
    GS_prev=GS; %store previous vector
    GS=P*GS_prev+LplusD\b; %uses last guess
    GS_iter=GS_iter+1; %counts iterations
    
    e_vector=GS-GS_prev; %error vector
    e=norm(e_vector,2)/norm(GS,2); %finds relative error
end
end