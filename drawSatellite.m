function drawSatellite(DCM, satellite_patch)
    % Function to update satellite orientation in the animation
    
    % Define satellite vertices (for visualization)
    satellite_vertices = [1 0 0; -0.5 0.5 0; -0.5 -0.5 0; 0 0 1];
    
    % Rotate satellite vertices using DCM
    rotated_satellite_vertices = (DCM * satellite_vertices')';
    
    % Update satellite position in the animation
    set(satellite_patch, 'XData', rotated_satellite_vertices(:, 1), ...
                         'YData', rotated_satellite_vertices(:, 2), ...
                         'ZData', rotated_satellite_vertices(:, 3));
end