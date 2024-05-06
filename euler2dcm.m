function DCM = euler2dcm(phi, theta, psi)
    % Euler angles to direction cosine matrix (DCM) conversion
    
    % Check if input arguments are provided
    if nargin < 3
        error('Insufficient input arguments. Provide phi, theta, and psi.');
    end
    
    % Calculate direction cosine matrix elements
    Rz = [cos(psi), -sin(psi), 0;
          sin(psi), cos(psi), 0;
          0, 0, 1];
    Ry = [cos(theta), 0, sin(theta);
          0, 1, 0;
          -sin(theta), 0, cos(theta)];
    Rx = [1, 0, 0;
          0, cos(phi), -sin(phi);
          0, sin(phi), cos(phi)];
    
    % Compute DCM
    DCM = Rz * Ry * Rx;
end
