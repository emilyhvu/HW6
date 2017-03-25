function [phi_anal] = computeAnalytical(x)

a=4;
D=1;
Sigma_A=0.2;
S=8;
L=sqrt(D/Sigma_A);
C1=-S/Sigma_A*exp(a/L)/(exp(2*a/L)+1);

phi_anal=C1*(exp(x/L)+exp(-x/L))+S/Sigma_A;

end