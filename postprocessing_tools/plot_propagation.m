%%
% BEGIN: User-modified variables
% format of the simulation file ('simfile'):
% <par. 1> ... <par. n> <sigma> <log-lik.> <sim.res. x 1> <sim.res. y 1> ... <sim.res. x N> <sim.res. y N>
% format of the data file ('datfile'):
% <data x i> <data y i> <data sigma i>

% data
simfile = 'data/propagation_sim.txt';
datfile = 'data/propagation_dat.txt';
dim = 1; % # model params

% plotting
isplot = 1; % shall we plot?
fignum = 1; % number of the figure
name_x = 'shear rate'; name_y = 'TTF'; % axes labels
lim_x = [0 13250]; lim_y = [0.0 2600]; % axes limits
nmesh_y = 1000; % number of mesh points in y direction
c1 = [44 157 255]/255; % colors
c2 = [204 0 102]/255;
% END: User-modified variables

%%
% Load data and simulation results
dat    = load(datfile);
dat_x  = dat(:, 1)';
dat_y  = dat(:, 2)';
dat_s  = dat(:, 3)'; % data sigma

sim    = load(simfile);
sim_n  = size(sim, 1);
sim_s  = sim(:, dim+1); % simulation sigma
sim_ll = sim(:, dim+2); % log-likelihood
sim_x  = sim(:, dim+3:2:end);
sim_y  = sim(:, dim+4:2:end);
sim_d  = size(sim_x, 2); % dimension of the simulation result

%%
% Construct a mesh
min_val_y = min(sim_y) - 0.3*abs(min(sim_y));
max_val_y = max(sim_y) + 0.3*abs(min(sim_y));
mesh_step_y = (max_val_y-min_val_y)/(nmesh_y-1);
mesh_y = zeros(nmesh_y, sim_d);
for i=1:nmesh_y; for j=1:sim_d
    mesh_y(i,j) = min_val_y(j) + (i-1)*mesh_step_y(j);
end; end

%%
% Compute PDF
fprintf('Computing PDF...\n');
pdf_y = zeros(nmesh_y, sim_d);
if prod(sim_s) > 0
    for k=1:sim_n
        sigma = sim_s(k);
        mu    = sim_y(k,:);
        for j=1:nmesh_y
            pdf_y(j,:) = pdf_y(j,:) + normpdf(mesh_y(j,:), mu, sqrt(sigma));
        end
    end
else
    for j=1:sim_d
        [pdf_y(:,j),edges] = histcounts(sim_y(:,j), nmesh_y, 'Normalization', 'probability');
        mesh_y(:,j) = 0.5*(edges(1:end-1) + edges(2:end));
        mesh_step_y(j) = mesh_y(2,j)-mesh_y(1,j);
    end
end
pdf_y = pdf_y ./ repmat(trapz(pdf_y).*mesh_step_y, nmesh_y, 1);

%%
% Compute CDF
fprintf('Computing CDF...\n');
cdf_y = zeros(nmesh_y, sim_d);
for j=1:sim_d
    cdf_y(1,j) = mesh_step_y(j)*pdf_y(1,j);
    for i=2:nmesh_y
        cdf_y(i,j) = cdf_y(i-1,j) + mesh_step_y(j)*pdf_y(i,j);
    end
end

%%
% Compute quantiles
fprintf('Computing quantiles...\n');
myqq = [0.05 0.95];
nqq = length(myqq);
q_y = zeros(nqq, sim_d);
for k=1:nqq
    for j=1:sim_d
        if myqq(k) <= cdf_y(1,j)
            q_y(k,j) = mesh_y(1,j);
        end
        for i=1:nmesh_y-1
            if cdf_y(i,j) < myqq(k) && myqq(k) <= cdf_y(i+1,j)
                q_y(k,j) = mesh_y(i+1,j);
            end
        end
    end
end

%%
% Compute best values
[~, best_id] = max(sim_ll);
sim_best_x = sim_x(best_id, :);
sim_best_y = sim_y(best_id, :);

%%
% Plot
if isplot
    set(0, 'DefaultTextFontSize', 30)
    set(0, 'DefaultAxesFontSize', 30)
    set(0, 'DefaultAxesFontName', 'Times')

    % create continuous x value array for plotting
    X = [sim_x(1,:), fliplr(sim_x(1,:))];

    % create continuous y values
    Y = zeros(nqq, sim_d);
    for i=1:2:nqq
        Y(i,   :) =        q_y(i,   :);
        Y(i+1, :) = fliplr(q_y(i+1, :));
    end

    % plot
    fig = figure(fignum); fignum = fignum + 1;
    grid off; box on
    xlabel(name_x); ylabel(name_y); xlim(lim_x); ylim(lim_y)
    hold on
    for i=1:2:length(myqq)
        h = fill(X(:), [Y(i,:) Y(i+1,:)], c1);
        alpha(0.25)
        set(h, 'EdgeColor', 'None');
    end
    plot(sim_best_x, sim_best_y, 'Color', c1, 'LineWidth', 3)
    errorbar(dat_x, dat_y, dat_s, 'Marker' , '.', 'Color', c2, 'LineWidth', 3, 'MarkerSize', 40)
    hold off
end
