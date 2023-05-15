%% Calculate analytical solution
% meta information
eig = 1;

% Setup the pde
model = createpde;
dl = decsg(gd_ex_2,sf_ex_2,ns_ex_2);
geometryFromEdges(model,dl);
pdegplot(model,"EdgeLabels","on")
applyBoundaryCondition(model,"dirichlet","Edge",1:4,"u",0);
specifyCoefficients(model,"c",1,"a",0,"d",1, "m",0,"f",0);
mesh_pde = generateMesh(model);
p = mesh_pde.Nodes;
% visualize mesh
pdemesh(model);

%% Calculate analytical solution for all mesh points
ana_sol = zeros(length(p),1);
if eig == 1
    for index = 1:length(p)
        ana_sol(index,1) = 2*sin(pi*p(1,index))*sin(pi*p(2,index));
    end
elseif eig == 2
    for index = 1:length(p)
        ana_sol(index,1) = sin(2*pi*p(1,index))*sin(pi*p(2,index));
    end
elseif eig == 3
    for index = 1:length(p)
        ana_sol(index,1) = sin(pi*p(1,index))*sin(2*pi*p(2,index));
    end
elseif eig == 4
    for index = 1:length(p)
        ana_sol(index,1) = sin(2*pi*p(1,index))*sin(2*pi*p(2,index));
    end
end
%% Calculate the (scaled) numerical solution
u = solvepdeeig(model,[-Inf,100]);

%Pick one of the eigenvalues
num_sol = u.Eigenvectors(:,eig);
num_eig = u.Eigenvalues

% scale eigenfunction to analytical solution using middle point of the mesh
scaling_fac = ana_sol(round(length(ana_sol)/2))/num_sol(round(length(ana_sol)/2));
num_sol = scaling_fac*num_sol;

%% Calculate error
errs = zeros(length(p),1);


for index = 1:length(p)
    errs(index,1) = abs(ana_sol(index,1)-num_sol(index,1));
end
max(errs)
%% Visualization of numerical solution and analytical solution
xlin = linspace(min(p(1,:)),max(p(1,:)),100);
ylin = linspace(min(p(2,:)),max(p(2,:)),100);
[X,Y] = meshgrid(xlin,ylin);

% visualize analytical solution on gridpoints
Z_ana = griddata(p(1,:),p(2,:),ana_sol,X,Y,'v4');
figure('NumberTitle', 'off', 'Name', 'Analytical solution');
mesh(X,Y,Z_ana);
axis tight; hold on
plot3(p(1,:),p(2,:),ana_sol,'.','MarkerSize',15);
title('Analytical solution');

% visualize numerical solution on gridpoints
Z_num = griddata(p(1,:),p(2,:),num_sol,X,Y,'v4');
figure('NumberTitle', 'off', 'Name', 'Numerical solution');
mesh(X,Y,Z_num);
axis tight; hold on
plot3(p(1,:),p(2,:),num_sol,'.','MarkerSize',15);
title('Numerical solution');
hold off
% visualize error
Z_error = griddata(p(1,:),p(2,:),errs,X,Y,'v4');
figure('NumberTitle', 'off', 'Name', 'Error')
mesh(X,Y,Z_error);
axis tight; hold on
plot3(p(1,:),p(2,:),errs,'.','MarkerSize',15);
title('Error');

%% Task 3: Error Convergence

x_3 = linspace(0,4,5);
y_3_num = [l_init(1), l_1(1), l_2(1), l_3(1), l_4(1)];
y_3_ana = ones(5)*2*pi^2;
figure
plot(x_3, y_3_num);
hold on
plot(x_3, y_3_ana);
axis tight; hold on
title("Error Convergence of Eigenvalues, First Eigenvalue");

x_3 = linspace(0,4,5);
y_3_num = [l_init(2), l_1(2), l_2(2), l_3(2), l_4(2)];
y_3_ana = ones(5)*5*pi^2;
figure
plot(x_3, y_3_num);
hold on
plot(x_3, y_3_ana);
axis tight; hold on
legend("Numerical Eigenvalue", "Analytical Eigenvalue");
title("Error Convergence of Eigenvalues, Second Eigenvalue");