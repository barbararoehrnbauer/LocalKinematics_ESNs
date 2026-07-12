%% Create field plots in the reference configuration

%
% The plots are generated in the reference configuration, as in the original
% plotting workflow. Therefore, all contour plots use the undeformed facet-center
% coordinates:
%
%   plotResult.facetCenters.referenceX
%   plotResult.facetCenters.referenceY

% Required workspace input:
%   Either currentResult must exist, or resultsBySetting must contain at least
%   one result structure.

%% Select result to plot

if exist('currentResult', 'var')
    plotResult = currentResult;
% elseif exist('resultsBySetting', 'var')
%     plotResult = resultsBySetting{end}{end};
else
    error('No result structure found. Run the Q4 kinematic analysis first.')
end

%% Plot settings

res = 100;
halfNetworkWidth = plotResult.networkWidth / 2;

% Reference facet-center coordinates. All plots below are evaluated and shown
% at these undeformed positions.
referenceX = plotResult.facetCenters.referenceX;
referenceY = plotResult.facetCenters.referenceY;

% Displacement fields at reference facet-center positions.
displacementX = plotResult.displacement.X;
displacementY = plotResult.displacement.Y;
displacementMagnitude = plotResult.displacement.magnitude;

% Averaged kinematic fields at reference facet-center positions.
% These averaged fields are intended for visualization. The statistical
% kinematic fingerprint should be based on the non-averaged local fields.
E11 = plotResult.averaged.E11;
E22 = plotResult.averaged.E22;
E12 = plotResult.averaged.E12;
deformationTypeExponent = plotResult.averaged.m;
rotationDeg = plotResult.averaged.R_deg;

%% x-displacement field in the reference configuration

figure
contourf(mean(referenceX), mean(referenceY'), displacementX, res, 'LineStyle', 'none')
colormap(jet)
colorbar
title('Displacement field U_x, FEM-MATLAB')
xlabel('X coordinate [\mum]')
ylabel('Y coordinate [\mum]')
axis equal
caxis([-3.5 4.5])
xlim([-halfNetworkWidth, halfNetworkWidth])
ylim([-halfNetworkWidth, halfNetworkWidth])
hold off

%% y-displacement field in the reference configuration

figure
contourf(mean(referenceX), mean(referenceY'), displacementY, res, 'LineStyle', 'none')
colormap(jet)
colorbar
title('Displacement field U_y, FEM-MATLAB')
xlabel('X coordinate [\mum]')
ylabel('Y coordinate [\mum]')
axis equal
caxis([-2.6 1.6])
xlim([-halfNetworkWidth, halfNetworkWidth])
ylim([-halfNetworkWidth, halfNetworkWidth])
hold off

%% Displacement magnitude in the reference configuration

figure
contourf(mean(referenceX), mean(referenceY'), displacementMagnitude, res, 'LineStyle', 'none')
colormap(jet)
colorbar
title('Displacement magnitude |U|')
xlabel('X coordinate [\mum]')
ylabel('Y coordinate [\mum]')
axis equal
xlim([-halfNetworkWidth, halfNetworkWidth])
ylim([-halfNetworkWidth, halfNetworkWidth])
hold off

%% Green-Lagrange strain E11 in the reference configuration

figure
contourf(mean(referenceX), mean(referenceY'), E11, res, 'LineStyle', 'none')
colormap(jet)
colorbar
title('Green-Lagrange strain E_{11}, FEM-MATLAB')
xlabel('X coordinate [\mum]')
ylabel('Y coordinate [\mum]')
axis equal
xlim([-halfNetworkWidth, halfNetworkWidth])
ylim([-halfNetworkWidth, halfNetworkWidth])
caxis([0.022 0.068])
% caxis([0.022 0.068])
hold off

%% Green-Lagrange strain E22 in the reference configuration

figure
contourf(mean(referenceX), mean(referenceY'), E22, res, 'LineStyle', 'none')
colormap(jet)
colorbar
title('Green-Lagrange strain E_{22}, FEM-MATLAB')
xlabel('X coordinate [\mum]')
ylabel('Y coordinate [\mum]')
axis equal
xlim([-halfNetworkWidth, halfNetworkWidth])
ylim([-halfNetworkWidth, halfNetworkWidth])
caxis([-0.05 0.002])
hold off

%% Green-Lagrange strain E12 in the reference configuration

figure
contourf(mean(referenceX), mean(referenceY'), E12, res, 'LineStyle', 'none')
colormap(jet)
colorbar
title('Green-Lagrange strain E_{12}, FEM-MATLAB')
xlabel('X coordinate [\mum]')
ylabel('Y coordinate [\mum]')
axis equal
xlim([-halfNetworkWidth, halfNetworkWidth])
ylim([-halfNetworkWidth, halfNetworkWidth])
caxis([-0.012 0.016]);
hold off

%% Deformation type exponent m in the reference configuration

figure
contourf(mean(referenceX), mean(referenceY'), deformationTypeExponent, res, 'LineStyle', 'none')
colormap(jet)
colorbar
title('Deformation type exponent m, FEM-MATLAB')
xlabel('X coordinate [\mum]')
ylabel('Y coordinate [\mum]')
axis equal
xlim([-halfNetworkWidth, halfNetworkWidth])
ylim([-halfNetworkWidth, halfNetworkWidth])
caxis([-1 1])
hold off

%% Rotation R in the reference configuration

figure
contourf(mean(referenceX), mean(referenceY'), rotationDeg, res, 'LineStyle', 'none')
colormap(jet)
colorbar
title('Rotation R, FEM-MATLAB')
xlabel('X coordinate [\mum]')
ylabel('Y coordinate [\mum]')
axis equal
caxis([-1 1])
xlim([-halfNetworkWidth, halfNetworkWidth])
ylim([-halfNetworkWidth, halfNetworkWidth])
hold off

%% Example: E11 contour plot interpolated in the reference configuration
% This section replaces the former current-configuration example. It is only
% needed if an interpolated contour plot is desired. The coordinates remain the
% undeformed reference coordinates.

xRow = referenceX(:);
yRow = referenceY(:);
zRow = E11(:);

validInterpolationPoint = ~isnan(xRow) & ~isnan(yRow) & ~isnan(zRow);
xRow = xRow(validInterpolationPoint);
yRow = yRow(validInterpolationPoint);
zRow = zRow(validInterpolationPoint);

if numel(xRow) >= 3

    gridStepX = plotResult.facetWidth / 10;
    gridStepY = gridStepX;

    [xGridReference, yGridReference] = meshgrid(min(xRow):gridStepX:max(xRow), ...
                                               min(yRow):gridStepY:max(yRow));

    zGridReference = griddata(xRow, yRow, zRow, xGridReference, yGridReference);

    % Mask the interpolated field outside the measured reference domain.
    boundaryIndex = boundary(xRow, yRow, 0.1);
    measuredReferenceDomain = polyshape(xRow(boundaryIndex), yRow(boundaryIndex), ...
                                        'Simplify', false);

    isInsideReferenceDomain = isinterior(measuredReferenceDomain, ...
                                         xGridReference(:), yGridReference(:));
    isInsideReferenceDomain = reshape(isInsideReferenceDomain, size(xGridReference));
    zGridReference(~isInsideReferenceDomain) = NaN;

%     figure
%     contourf(xGridReference, yGridReference, zGridReference, res, 'LineStyle', 'none')
%     colormap(jet)
%     axis equal
%     xlabel('X coordinate [\mum]')
%     ylabel('Y coordinate [\mum]')
%     title('Green-Lagrange strain E_{11}, reference configuration')
%     hold on
% 
%     plot([min(referenceX(:)), max(referenceX(:)), max(referenceX(:)), min(referenceX(:)), min(referenceX(:))], ...
%          [min(referenceY(:)), min(referenceY(:)), max(referenceY(:)), max(referenceY(:)), min(referenceY(:))], ...
%          'black')
% 
%     colorbar
%     xlim([-halfNetworkWidth, halfNetworkWidth])
%     ylim([-halfNetworkWidth, halfNetworkWidth])
%     hold off
else
    warning('The interpolated reference-configuration plot was skipped because fewer than three valid points are available.')
end
