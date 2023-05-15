% create mesh, receiving gd, sf
pderect([0 1])
model = createpde;
dl = decsg(gd);
geometryFromEdges(model,dl);
%pdegplot(model,"EdgeLabels","on")

applyBoundaryCondition(model,"dirichlet","Edge",1:4,"u",0);
specifyCoefficients(model,"c",1,"a",0,"d",1, "m",0,"f",0);
mesh_Hmax = generateMesh(model);
figure
pdemesh(model)
results = solvepdeeig(model,[-Inf,20]);

%Pick one of the simple eigenvalues
eig_1 = results.Eigenvectors(:);
max(abs(eig_1))
eig_1_val = results.Eigenvalues(1);
figure
pdeplot(model,'XYData',eig_1)

nodes_pde = results.Mesh.Nodes;

% calculate analytical solution
analytical_sol = zeros(length(nodes_pde),1);
%errors = zeros(length(nodes_pde),1);
for index = 1:length(nodes_pde)
    analytical_sol(index,1) = 2*sin(pi*nodes_pde(1,index))*sin(pi*nodes_pde(2,index));
end

analytical_sol
%xlin = linspace(min(nodes_pde(1,:)),max(nodes_pde(1,:)),100);
%ylin = linspace(min(nodes_pde(2,:)),max(nodes_pde(2,:)),100);
%[X,Y] = meshgrid(xlin,ylin);

% visualize analytical solution on gridpoints
%Z = griddata(nodes_pde(1,:),nodes_pde(2,:),eig_1,X,Y,'v4');
%figure
%mesh(X,Y,Z)
%axis tight; hold on
%plot3(nodes_pde(1,:),nodes_pde(2,:),eig_1,'.','MarkerSize',15)

% calculate errors
%errors(:,1) = analytical_sol(:,1) - eig_1(:,1)
%Z_err = griddata(nodes_pde(1,:),nodes_pde(2,:),errors,X,Y,'v4');
%figure
%mesh(X,Y,Z_err)
%axis tight; hold on
%plot3(nodes_pde(1,:),nodes_pde(2,:),errors,'.','MarkerSize',15)