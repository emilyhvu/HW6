function [Phi] = Thomas(A,b)

dim=size(A,1);
for i=2:dim
    L=A(i,i-1);
    D_Current=A(i,i);
    D_Previous=A(i-1,i-1);
    U=A(i-1,i);
    A(i,i)=D_Current-(L/D_Previous)*U;
    b(i)=b(i)-(L/D_Previous)*b(i-1);
    A(i,i-1)=0;
end

Phi=zeros(dim+2,1);
Phi(dim)=b(dim)/A(dim,dim);
for i=dim-1:-1:2
    U=A(i-1,i);
    D=A(i,i);
    Phi(i)=(b(i)-U*Phi(i+1))/D;
end

end

