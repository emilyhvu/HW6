function [x_axis, y_axis] = RelativeError(a,D,Sigma_A,S,h_values)

x_axis=[];
y_axis=[];

for i=1:length(h_values)
    h=h_values(i);
    [anal_phi ,phi, mesh] = FixedSourceSolver(a,D,Sigma_A,S,h);
    error=abs(norm(phi,2)-norm(anal_phi,2))/norm(phi,2);
    mesh_length=length(mesh)-1;
    
    x_axis=[x_axis mesh_length];
    y_axis=[y_axis error];
end
close all

plot(x_axis,y_axis,'bo-')
title('Maximum Relative Error vs. Total Number of Meshes')
xlabel('Total Number of Meshes')
ylabel('Maximum Relative Error')

end
