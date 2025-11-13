%==========================================================================
% LigeiaMare_FlowVisualization
%
% Author: Eloy Solanes Gracia
% Year: 2025
%
% ------------
% Description:
% ------------
%
% This script computes and visualizes the velocity fields (u and w)
% in a semi-enclosed basin of Ligeia Mare (Titan) for different
% rheological indices n (e.g., n = 0.6, 1.0, 1.2).
% It computes the fluid column for each discretization point in the
% x-domain.
%
% The bathymetric profile of the basin is incorporated into the model.
% Depending on the selected mode:
%   Mode 1 -> Plots the horizontal velocity component u(x, z)
%   Mode 2 -> Computes and plots the vertical velocity component w(x, z)
%
% The script loads precomputed matrices of velocity (U_matriz_nXX.mat)
% and normalizes them for visualization.
%
%
%==========================================================================
clear all
close all
clc

% mode 1 --> Plot u(x, z)
% mode 2 --> Plot w(x, z)
mode = 1;

% Bathymetry profile (Mastrogiuseppe et al. (2014))
x_dom = linspace(0, 300, 300);
bathy_profile = zeros(1, length(x_dom));
for i = 1:length(x_dom)
    x = x_dom(i);
    bathy_profile(i) = -0.000018 * x.^3 + 0.012322 * x.^2 - 1.955119 * x - 60.654995;
end

% Matrix to store vertical coordinates
Z_matriz = zeros(400, length(x_dom));

% Generate depth levels for each x
for j = 1:length(x_dom)
    H = -bathy_profile(j);
    z = linspace( 0, -H, 400)';
    Z_matriz(:,j) = z;
end

% Load precomputed velocity matrices (400x300)
U_n06_struct = load('U_matriz_n06.mat');
U_n10_struct = load('U_matriz_n10.mat');
U_n12_struct = load('U_matriz_n12.mat');

U_n06 = U_n06_struct.U_matriz;
U_n10 = U_n10_struct.U_matriz;
U_n12 = U_n12_struct.U_matriz;

if mode == 1
    % Variables
    U_max_06 = 0;
    U_max_10 = 0;
    U_max_12 = 0;
    
    % Iterate throgh the three matrices to find their maximum values
    for i = 1:3
        U_max = 0;
        if i == 1
            U_matriz = U_n06;
        elseif i == 2
            U_matriz = U_n10;
        elseif i == 3
            U_matriz = U_n12;
        end
    
        % Loop through each element to find the maximum absolute velocity value
        for j = 1:400
            for k = 1:300
                U = abs(U_matriz(j,k));
                
                if U > U_max
                    U_max = U;
                end
            end
        end

        % Store the maximum velocity for each U matrix
        if i == 1
            U_max_06 = U_max;
        elseif i == 2
            U_max_10 = U_max;
        elseif i == 3
            U_max_12 = U_max;
        end
    end
      
    % Normalize the velocity matrices by the largest maximum among the three
    U_n06_Norm = U_n06 / U_max_12;
    U_n10_Norm = U_n10 / U_max_12;
    U_n12_Norm = U_n12 / U_max_12;

    % Bathymetry mesh (X–Z grid)
    [X, Z] = meshgrid(x_dom, linspace(0, max(bathy_profile), 400));
    for i = 1:length(x_dom)
         Z(:, i) = linspace(0, bathy_profile(i), 400);
    end
    
    % Plot results for n = 0.6
    figure
    surf(X, Z, U_n06_Norm, 'EdgeColor', 'none')
    colormap(turbo(18))
    colorbar
    clim([-0.8 1])
    shading flat
    xlabel('Horizontal distance [km]', 'Interpreter', 'latex')
    ylabel('Depth [m]', 'Interpreter', 'latex')
    zlabel('Normalized velocity', 'Interpreter', 'latex')
    title('u(x,z) for n = 0,6', 'Interpreter', 'latex')
    set(gcf, 'Color', 'w')
    view([0, 90])

    % Plot results for n = 1
    figure
    surf(X, Z, U_n10_Norm, 'EdgeColor', 'none')
    colormap(turbo(18))
    colorbar
    clim([-0.8 1])
    shading flat
    xlabel('Horizontal distance [km]', 'Interpreter', 'latex')
    ylabel('Depth [m]', 'Interpreter', 'latex')
    zlabel('Normalized velocity', 'Interpreter', 'latex')
    title('u(x,z) for n = 1', 'Interpreter', 'latex')
    set(gcf, 'Color', 'w')
    view([0, 90])
    
    % Plot results for n = 1.2
    figure
    surf(X, Z, U_n12_Norm, 'EdgeColor', 'none')
    %shading interp
    %colormap turbo
    colormap(turbo(18))
    colorbar
    clim([-0.8 1])
    shading flat
    xlabel('Horizontal distance [km]', 'Interpreter', 'latex')
    ylabel('Depth [m]', 'Interpreter', 'latex')
    zlabel('Normalized velocity', 'Interpreter', 'latex')
    title('u(x,z) for n = 1,2', 'Interpreter', 'latex')
    set(gcf, 'Color', 'w')
    view([0, 90])
end


if mode == 2
    
    % Create the width profile b(x)
    b_min = 40;          % Minimum width [km]
    b_max = 50;          % Maximum width [km]
    %b_x = b_min + (b_max - b_min) * (1 - ((2 * x_dom / 300) - 1).^2);
    
    % Shifted parabolic profile
    b_x = b_min + (b_max - b_min) .* ( 1 - ((x_dom - 230)./(300/2)).^2 );
    [b, ~] = meshgrid(b_x, linspace(0, 1, 400));
   
    % Variables
    W_max_06 = 0;
    W_max_10 = 0;
    W_max_12 = 0;
       
    % Loop through the three velocity matrices
    for i = 1:3
        W_max = 0;
        if i == 1
            U_matriz = U_n06;
            bu = b .* U_matriz;
        elseif i == 2
            U_matriz = U_n10;
            bu = b .* U_matriz;
        elseif i == 3
            U_matriz = U_n12;
            bu = b .* U_matriz;
        end
        
        % Compute dw/dz and integrate to obtain w(x,z)
        dx = x_dom(2) - x_dom(1);
        dbu_dx = gradient(bu, dx, 2);
        dw_dz = - dbu_dx ./ b;
        
        W_matriz = zeros(size(U_matriz));
        for j = 1:length(x_dom)
            W_matriz(:,j) = cumtrapz(Z_matriz(:,j), dw_dz(:,j));
        end
        
        % Loop through each element to find the maximum absolute value of w
        for j = 1:400
            for k = 1:300
                W = abs(W_matriz(j,k));
                
                if W > W_max
                    W_max = W;
                end
            end
        end
        
        % Store results for each rheological index
        if i == 1
            W_n06 = W_matriz;
            W_max_06 = W_max;
        elseif i == 2
            W_n10 = W_matriz;
            W_max_10 = W_max;
        elseif i == 3
            W_n12 = W_matriz;
            W_max_12 = W_max;
        end
    end
    
    % Normalize the vertical velocity matrices using the largest maximum
    W_n06_Norm = W_n06 / W_max_12;
    W_n10_Norm = W_n10 / W_max_12;
    W_n12_Norm = W_n12 / W_max_12;

    % Bathymetry mesh (X–Z grid)
    [X, Z] = meshgrid(x_dom, linspace(0, max(bathy_profile), 400));
    for i = 1:length(x_dom)
         Z(:, i) = linspace(0, bathy_profile(i), 400);
    end

     % Plot results for n = 0.6
    figure
    surf(X, Z, W_n06_Norm, 'EdgeColor', 'none')
    colormap(turbo(30))
    colorbar
    clim([-0.5 1])
    xlabel('Horizontal distance [km]', 'Interpreter', 'latex')
    ylabel('Depth [m]', 'Interpreter', 'latex')
    zlabel('Normalized velocity', 'Interpreter', 'latex')
    title('w(x,z) for n = 0,6', 'Interpreter', 'latex')
    set(gcf, 'Color', 'w')
    view([0, 90])
    
     % Plot results for n = 1
    figure
    surf(X, Z, W_n10_Norm, 'EdgeColor', 'none')
    colormap(turbo(30))
    colorbar
    clim([-0.5 1])
    xlabel('Horizontal distance [km]', 'Interpreter', 'latex')
    ylabel('Depth [m]', 'Interpreter', 'latex')
    zlabel('Normalized velocity', 'Interpreter', 'latex')
    title('w(x,z) for n = 1', 'Interpreter', 'latex')
    set(gcf, 'Color', 'w')
    view([0, 90])
    
     % Plot results for n = 1.2
    figure
    surf(X, Z, W_n12_Norm, 'EdgeColor', 'none')
    colormap(turbo(30))
    colorbar
    clim([-0.5 1])
    xlabel('Horizontal distance [km]', 'Interpreter', 'latex')
    ylabel('Depth [m]', 'Interpreter', 'latex')
    zlabel('Normalized velocity', 'Interpreter', 'latex')
    title('w(x,z) for n = 1,2', 'Interpreter', 'latex')
    set(gcf, 'Color', 'w')
    view([0, 90])
end