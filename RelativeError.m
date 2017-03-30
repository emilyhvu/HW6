function [x_axis, y_axis] = RelativeError(a,D,Sigma_A,S,h_values)
%Question 3

x_axis=[];
y_axis=[];

for i=1:length(h_values)
    h=h_values(i); %pass in different h values
    [anal_phi ,phi, mesh] = FixedSourceSolver(a,D,Sigma_A,S,h); %outputs numerical and analytical solutions
    mesh_length=length(mesh)-1;
    error=(anal_phi-phi)./anal_phi;  
    error=error(2:end-1); %remove first/last elements
    
    x_axis=[x_axis mesh_length];
    y_axis=[y_axis max(error)];
end
close all

loglog(x_axis,y_axis,'bo-')
hold all
%plot(h_values,y_axis)

title('Maximum Relative Error vs. Total Number of Meshes')
xlabel('Total Number of Meshes')
ylabel('Maximum Relative Error')
grid on

end
