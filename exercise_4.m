%%
%hmaxs = [0.2,0.12,0.04];
hmaxs = logspace(-1,-2,10);
hmaxs(end) = [];
errors = zeros(8,length(hmaxs));
node_points = zeros(8,length(hmaxs));
for alpha = 4
    if alpha == 1
        gd = gd_0;
        ns = ns_0;
        sf = sf_0;
    elseif alpha == 2
        gd = gd_1;
        ns = ns_1;
        sf = sf_1;
    elseif alpha == 3
        gd = gd_2;
        ns = ns_2;
        sf = sf_2;
    elseif alpha == 4
        gd = gd_3;
        ns = ns_3;
        sf = sf_3;
    elseif alpha == 5
        gd = gd_4;
        ns = ns_4;
        sf = sf_4;
    elseif alpha == 6
        gd = gd_5;
        ns = ns_5;
        sf = sf_5;
    elseif alpha == 7
        gd = gd_6;
        ns = ns_6;
        sf = sf_6;
    elseif alpha == 8
        gd = gd_7;
        ns = ns_7;
        sf = sf_7;
    end
    model = createpde;
    dl = decsg(gd,sf,ns);
    geometryFromEdges(model,dl);
    pdegplot(model,"EdgeLabels","on")
    if alpha == 1
        applyBoundaryCondition(model,"dirichlet","Edge",1:3,"u",0);
    else
        applyBoundaryCondition(model,"dirichlet","Edge",1:5,"u",0);
    end
    specifyCoefficients(model,"c",1,"a",0,"d",1, "m",0,"f",0);

    ana_hmax = 0.01;
    mesh_Hmax = generateMesh(model,"Hmax",0.01,"GeometricOrder","linear");
    pdemesh(model);
    ana_results = solvepdeeig(model,[-Inf,20]);
    
    for i = 1:length(hmaxs)
        mesh_Hmax = generateMesh(model,"Hmax",hmaxs(i), "GeometricOrder","linear");
        
        node_size = size(mesh_Hmax.Nodes);
        node_points(alpha,i) = node_size(2);
        %node_points(alpha,i) = hmaxs(i);
        
 
        pdemesh(model);
        result_num = solvepdeeig(model,[-Inf,20]);
        num_eig = result_num.Eigenvalues(1);
        errors(alpha,i) = abs(num_eig-ana_results.Eigenvalues(1));
    end
    
end
disp('job done')
%%
loglog(node_points(1,:),errors(1,:),'-x','DisplayName',append('\alpha_', sprintf('%0.0f',0)));
hold on
loglog(node_points(2,:),errors(2,:),'-x','DisplayName',append('\alpha_', sprintf('%0.0f',1)));
loglog(node_points(3,:),errors(3,:),'-x','DisplayName',append('\alpha_', sprintf('%0.0f',2)));
loglog(node_points(4,:),errors(4,:),'-x','DisplayName',append('\alpha_', sprintf('%0.0f',3)));
loglog(node_points(5,:),errors(5,:),'-x','DisplayName',append('\alpha_', sprintf('%0.0f',4)));
loglog(node_points(6,:),errors(6,:),'-x','DisplayName',append('\alpha_', sprintf('%0.0f',5)));
loglog(node_points(7,:),errors(7,:),'-x','DisplayName',append('\alpha_', sprintf('%0.0f',6)));
loglog(node_points(8,:),errors(8,:),'-x','DisplayName',append('\alpha_', sprintf('%0.0f',7)));
hold off
lgd = legend;
%set(gca, 'XDir', 'reverse');
xlabel('number of nodes')
ylabel('absolute error to anal. sol.')
title('Convergence analysis on non-convex domain')
hold off
%% Determine order
rates = zeros(1,8);
constants = zeros(1,8);
for alpha = 1:8
    log_node_points = log(node_points(alpha,:));
    log_errs = log(errors(alpha,:));
    a = polyfit(log_node_points, log_errs, 1);
    rates(alpha) = a(1);
    constants(alpha) = a(2);
end
%%


    
