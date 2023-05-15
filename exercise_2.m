%% Calculate analytical solution
% meta information
eig = 1;
hmaxs = linspace(0.2,0.002,5);
% calculate errors for first 5 eigenvalues
errors = zeros(5,length(hmaxs));
n_nodes = zeros(length(hmaxs),1);
% Setup the pde
model = createpde;
dl = decsg(gd_ex_2,sf_ex_2,ns_ex_2);
geometryFromEdges(model,dl);
pdegplot(model,"EdgeLabels","on")
applyBoundaryCondition(model,"dirichlet","Edge",1:4,"u",0);
specifyCoefficients(model,"c",1,"a",0,"d",1, "m",0,"f",0);

for i = 1:length(hmaxs)
    mesh_pde = generateMesh(model,"Hmax",hmaxs(i),"GeometricOrder","linear");
    p = mesh_pde.Nodes;
    pdemesh(model);
    n_nodes(i) = length(p);
    % Calculate the (scaled) numerical solution
    u = solvepdeeig(model,[-Inf,300]);
    
    % Pick one of the eigenvalues
    num_eigs = u.Eigenvalues(1:5).';
    % Analytical eigenvalues
    ana_eigs = [2*pi^2, 5*pi^2, 5*pi^2, 8*pi^2, 10*pi^2];
    
    errs = abs(num_eigs - ana_eigs);
    for j=1:5
        errors(j,i) = errs(j);
    end
end
%%
loglog(n_nodes,errors(1,:),'-x','DisplayName','\lambda_{11}');
hold on
loglog(n_nodes,errors(2,:),'-x','DisplayName','\lambda_{12}');
loglog(n_nodes,errors(3,:),'-x','DisplayName','\lambda_{21}');
loglog(n_nodes,errors(4,:),'-x','DisplayName','\lambda_{22}');
loglog(n_nodes,errors(5,:),'-x','DisplayName','\lambda_{13}');
hold off
lgd = legend;
%set(gca, 'XDir', 'reverse');
xlabel('number of nodes')
ylabel('absolute error to anal. sol.')
title('Convergence analysis of different eigenvalues')
hold off
%%
%% Determine order
rates = zeros(5,1);
constants = zeros(5,1);
for elem = 1:5
    log_node_points = log(n_nodes);
    log_errs = log(errors(elem,:));
    a = polyfit(log_node_points, log_errs, 1);
    rates(elem) = a(1);
    constants(elem) = a(2);
end
