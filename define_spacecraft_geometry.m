function define_spacecraft_geometry(length, width, height)

% Example: Cube-shaped spacecraft

    % Visualize spacecraft geometry
    figure;
    hold on;
    axis equal;
    xlabel('X-axis (m)');
    ylabel('Y-axis (m)');
    zlabel('Z-axis (m)');
    title('Spacecraft Geometry');

    % Plot spacecraft components
    % Example: Cube-shaped spacecraft
    drawCube([0, 0, 0], length, width, height);

    % Add labels for spacecraft components
    text(length/2, 0, 0, 'Length', 'HorizontalAlignment', 'center');
    text(0, width/2, 0, 'Width', 'HorizontalAlignment', 'center');
    text(0, 0, height/2, 'Height', 'HorizontalAlignment', 'center');

    % Generate mesh using PDE Toolbox
    g = decsg([3 4 0 10 10 0 0 5 5], 'R1+R2-R3'); % Define a rectangle with corners at (0,0) and (10,5)
    msh = generateMesh(g); % Generate mesh
    nodes = msh.Nodes; % Get mesh nodes
    elements = msh.Elements; % Get mesh elements

    % Visualize mesh
    trisurf(elements, nodes(:,1), nodes(:,2), 'FaceColor', 'cyan', 'EdgeColor', 'none');
end


function drawCube(center, length, width, height)
    % Draw a cube with specified dimensions and center
    vertices = [
        center + [-length/2, -width/2, -height/2];
        center + [length/2, -width/2, -height/2];
        center + [length/2, width/2, -height/2];
        center + [-length/2, width/2, -height/2];
        center + [-length/2, -width/2, height/2];
        center + [length/2, -width/2, height/2];
        center + [length/2, width/2, height/2];
        center + [-length/2, width/2, height/2];
    ];
    faces = [
        1, 2, 3, 4;
        5, 6, 7, 8;
        1, 2, 6, 5;
        2, 3, 7, 6;
        3, 4, 8, 7;
        4, 1, 5, 8;
    ];
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'blue', 'FaceAlpha', 0.5);

end

