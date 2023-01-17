clear all; close all; clc;

global K_global b_global M N

% Define constants
epsilon_naught = 8.854187817E-12; % Permittivity of free space (F/m)
mu_naught = 4*pi*10e-7; % Permeability of vaccum in magnetic field (Hm)
freq = 1e6; % 1 GHz
k = 2*pi*freq/3e8; % Check definition

% TODO: complete
M = 0;
N = 0

% Initialize the global matrices
K_global = zeros(num_unknowns);
b_global = zeros(num_unknowns,1);

% loop over all i and k elements
for i = 1:num_elements
   
    % Initialize matrices
    Kji = zeros(3);
    b = zeros(num_unknowns,1);

    % Get the node coordinates
    Xi = [x1_i, x2_i, x3_i];
    Yi = [y1_i, y2_i, y3_i];


    for j = 1:num_elements
        % Get the coordinates of this element's nodes
        Xj = [x1_j, x2_j, x3_j];
        Yj = [y1_j, y2_j, y3_j];


        % Initialize the summations
        A_sum = 0;
        B_sum = 0;
        y_sum = 0; % RHS

        for p = 1:3

            % Calculate xp
            x_p = zeta_i(p,:) * ki_node_coord;

            % Calculate lambda p 
            lambda_p =  1/(2*Area_i)*[[zeta_i(p,3)*T1_i - zeta_i(p,2)*T2_i]; ...
                        [zeta_i(p,1)*T2_i - zeta_i(p,3)*T0_i]; ...
                        [zeta_i(p,2)*T0_i - zeta_i(p,1)*T1_i]];
    

            % Calculate local incident field vector (right-hand side)
            z = 0;
            x_hat = x_p/norm(x_p);
            e_inc = x_hat * exp(-1i*k*z); % = x_hat
            e_vec = [e_inc; e_inc; e_inc];
            y_sum = y_sum + dot(lambda_p,e_vec)' * 1/6;

            for q = 1:3
                % Calculate yq
                y_q = zeta_j(q,:) * kj_node_coord;

                % Calculate R
                Rpq = sqrt((x_p(1)-y_q(1))^2 + (x_p(2)-y_q(2))^2 + (x_p(3)-y_q(3))^2); % Scalar

                % Calculate Fpq
                Fpq = exp(-1i * k * Rpq) / (4*pi * Rpq); % Complex number
    
                % Calculate Apq
                A_val = Fpq * lambda_p' * lambda_q;
                A_sum = A_sum + A_val;

                % Calculate Bqp
                B_sum = B_sum + Fpq * ones(3,3);
            
            end %q            
        end %p

        % Calculate local K matrix
        Kij = []

        % Calculate the total electric field (b)
        b = []
        
        % Add any values that contribute to the unknowns to the globals
        map_to_global(i,j,Kij,b);

    end %j
end %i
   
% Solve for the permittivity coefficients
a = K \ b_global;


% Plot recovvered permittivity vs original permittivity


% Function for mapping from local matrices to global matrix
function map_to_global(i,j,K,b)
    global K_global b_global M N

    % Map to global matrix
    % for row = 1:3
    %     for col = 1:3
    %         % If either element holds a known value here, skip.
    %         if(state_rows(row) == -1 || state_cols(col) == -1)
    %             continue;
    %         end
    %         % Multiply by +1 or -1 depending on sign, and add to the 
    %         %   global matrices
    %         K_global(state_rows(row),state_cols(col)) = K_global(state_rows(row),state_cols(col)) + (dir_rows(row) * dir_cols(col)*Z(row,col));
    %         b_global(state_rows(row),:) = b_global(state_rows(row),:) + (dir_rows(row) * y(row,:));

    %     end %col
    % end %row
end % map ftn