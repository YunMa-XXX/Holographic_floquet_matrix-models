close all
clear all
clc

set(groot,'defaulttextinterpreter','latex','defaulttextFontWeight','bold')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultColorbarTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultAxesFontSize',20)
set(groot,'defaultLineLineWidth',1.5)


%% Build Clifford gates: Hadamard H, Phase S, CNOT
I = eye(2);
CNOT = [1 1 0 0; 0 1 0 0 ; 0 0 1 0; 0 0 1 1]; % binary form
H = [0 1;1 0]; % binary form
S = [1 0;1i 1]; % binary form

%% Build the four-qubit Clifford gate
gate = phase(4) * controlled(4,1) * controlled(1,2) * controlled(2,3) * controlled(3,4) * hadamard(1); %% the one we originally use
gate = mod(gate,2);

% Q = [zeros(length(gate)/2), eye(length(gate)/2);eye(length(gate)/2), zeros(length(gate)/2)]; % check if gate is Clifford
% symplectic = mod(gate*Q*gate',2);
% check1 = sum(sum(symplectic ~= Q));

%% Main body
N_and_ts_vec = [];
N_values = [8, 12, 16, 20, 24, 28, 32, 36, 40, 44, 48, 52, 56, 60, 64];
% N_values = 2*(4:32);

for runs_ind = 1:30
    colors = lines(length(N_values)+1); 
    % colors(3,:) = [];
    colors = colors(1:length(N_values),:);

    t_final = 40; % final time
     
    p_log = [];
    p_OSG = [];
    p_sat = [];
    legend_OSG_vec = [];
    ts_vec = [];
    
    figure(1)
    for N_ind = 1:length(N_values) 
        N = N_values(N_ind);
        n = N^2; 
        U_permute = zeros(2*n);
        U_permute_one_kind = zeros(n);
        sigma = zeros(N);
        P = zeros(2*n);
        pauli_string = zeros(2*n,1);
        initial = randi([1,2*n]);
        % initial = 23;
        pauli_string(initial,1) = 1;
        pairs_vec = [1];
        time_vec = [0];
    
        %% Build the four-qubit Clifford gate
        D = kron(eye(2*n/8), gate);
        
        %% Find U_permute (and thus also need to find sigma first):
        v1 = (1:N)';
        for k = 1:N
            if mod(k,2) == 1
                v1(k,2) = (k+1)/2;
            else
                v1(k,2) = (N+v1(k,1))/2; % first column in v1 records [1;2;...;n], second column records result after permuting v1 according the the predetermined permutation rule
            end
        end
        for k = 1:N
            sigma(v1(k,2),v1(k,1)) = 1; % We've built the permutation matrix sigma
        end
        for kkk = 1:n
            kkk_row_ind = floor(kkk/N) + 1;
            kkk_col_ind = mod(kkk,N);
            if kkk_col_ind == 0
                kkk_row_ind = kkk_row_ind - 1;
                kkk_col_ind = N;
            end
            kkk_row_perm_ind = find(sigma(kkk_row_ind,:) == 1);
            kkk_col_perm_ind = find(sigma(kkk_col_ind,:) == 1);
            U_permute_one_kind(N*(kkk_row_ind-1)+kkk_col_ind, N*(kkk_row_perm_ind-1)+kkk_col_perm_ind) = 1;
        end
        U_permute = kron(eye(2), U_permute_one_kind); % We've built U_permute
    
        %% Find index change permutation matrix P:
        [i, j] = meshgrid(1:2*N, 1:N); 
        v2 = [i(:), j(:)]; 
        v3 = [];
        a1 = [1 1;1 4;4 4;4 1]; % [X11; X14; X44; X41] % Uncomment if want Rule 1 (STARTS HERE!)
        a1 = [a1; a1(:,1)+N a1(:,2)]; % [X11; X14; X44; X41; Z11; Z14; Z44; Z41]
        a2 = [2 2;2 3;3 3;3 2]; % [X22; X23; X33; X32]
        a2 = [a2; a2(:,1)+N a2(:,2)]; % [X22; X23; X33; X32; Z22; Z23; Z33; Z32]
        a3 = [1 3;3 4;4 2;2 1]; % [X13; X34; X42; X21]
        a3 = [a3; a3(:,1)+N a3(:,2)]; % [X13; X34; X42; X21; Z13; Z34; Z42; Z21]
        a4 = [1 2;2 4;4 3;3 1]; % [X12; X24; X43; X31]
        a4 = [a4; a4(:,1)+N a4(:,2)]; % [X12; X24; X43; X31; Z12; Z24; Z43; Z31]
        a = [a1;a2;a3;a4];
        for i = 1:N/4;
            v3(4^2*2*(i-1)+1:4^2*2*i,:) = a + 4*(i-1); % all 8x8 diagonal blocks containing both X and Z gates
        end   % Uncomment if want Rule 1 (ENDS HERE!)
        % v3 = [0,0]; % Uncomment if want Rule 2
        v4 = setdiff(v2, v3,'rows');
        for i = 1:length(v4)
            if ~ismember(v4(i,:), v3, 'rows') % for each row in v4 this should always be satisfied
                if mod(v4(i,1),2) == 0
                    p_i = 0;
                else
                    p_i = 1;
                end
                if mod(v4(i,2),2) == 0
                    p_j = 0;
                else
                    p_j = 1;
                end
                b = [v4(i,1) v4(i,2); v4(i,2) v4(i,1)+(-1)^(p_i+1); v4(i,1)+(-1)^(p_i+1) v4(i,2)+(-1)^(p_j+1); v4(i,2)+(-1)^(p_j+1) v4(i,1)]; % rule of rearrangement for interaction (i.e. change from normal set to predetermined set) for any 8x8 off-diagonal block (i.e.containing both X and Z gates)
                b = [b; b(:,1)+N b(:,2)];
                v3 = [v3; b];
            end 
        end
        % v3(1,:) = []; % Uncomment if want Rule 2
        v2(:,3) = (v2(:,1)-1)*N + v2(:,2);
        v3(:,3) = (v3(:,1)-1)*N + v3(:,2);
        for k = 1:length(v2)
            P(v2(k,3),v3(k,3)) = 1; % converts the normal index set to the cyclic connectivity index set
        end
        
        %% Now we begin the time evolution of the system
        t = 0;
        while t ~= t_final
            t = t + 1;
            time_vec = [time_vec, t];
    
            last_size = sum(sum(pauli_string));
    
            pauli_string = U_permute * pauli_string;
            pauli_string = P * pauli_string;
            pauli_string = pauliprod(D,pauli_string);
            pauli_string = P' * pauli_string;
            
            pairs = sum(pauli_string(1:length(pauli_string)/2) | pauli_string((length(pauli_string)/2+1): length(pauli_string))); % number of single-qubit Pauli operators in the Pauli string
            pairs_vec = [pairs_vec, pairs];
        end
    
        c = colors(N_ind,:);
        p_OSG = [p_OSG, semilogy(time_vec, pairs_vec,'Color', c)]; 
        a = log(N^2)/log(4);
        b = interp1(time_vec, pairs_vec, a, 'linear', 'extrap');   
        hold on
        plot(a, b, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', c, 'LineWidth', 1.5); % when t=log_4 (N^2), what's the corresponding operator size
        p_log = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
        hold on
        plot(time_vec, 3*n/4 * ones(1,length(time_vec)),'--','Color', c,'HandleVisibility','off') % the typical operator size at late times 
        set(gca,'yscale','log')
        legend_OSG_vec = [legend_OSG_vec, strcat('$N^2$=', string(n))];
        hold on
    
        y_target = 1/2 * N^2*(1-0.00);
        % 1) first index where y reaches/exceeds the threshold
        y_index_max = find(pairs_vec >= y_target, 1, 'first');
        
        if isempty(y_index_max)
            error('y never reaches the target value.');
        end
        
        new_x = time_vec(y_index_max-1:y_index_max);
        new_y = pairs_vec(y_index_max-1:y_index_max);
        
        % 2) optional sanity check: monotone nondecreasing
        if any(diff(new_y) < 0)
            warning('new_y is not monotone nondecreasing; 2a may be unreliable. Consider 2b on the full data.');
        end
        
        x_s = interp1(new_y, new_x, y_target, 'linear'); 
        y_s = y_target;
        plot(x_s, y_s, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'LineWidth', 1.5);
        p_sat = plot(NaN, NaN, 'o', 'MarkerEdgeColor', [0.5 0.5 0.5], 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
        ts_vec = [ts_vec,x_s];
    end
    N_and_ts_vec = [N_and_ts_vec;ts_vec];
end

%% Plotting
figure(2)
X = log(1/2*N_values.^2);
Y = N_and_ts_vec;
mean_Y = mean(Y, 1);        
std_Y  = std(Y, 0, 1);      
errorbar(X, mean_Y, std_Y, 'o-')
xlabel('$\ln(N^2/2)$')
ylabel('$t_s$')
title('Rule 1')
% title('Rule 2')
hold on
W = 1 ./ std_Y.^2;
Xmat = [X(:), ones(length(X),1)];   
Wmat = diag(W);
beta = (Xmat' * Wmat * Xmat) \ (Xmat' * Wmat * mean_Y(:));
a = beta(1);
b = beta(2);
Y_fit = a*X + b;
plot(X, mean_Y, 'o'); 
hold on
fitting_curve = plot(X, Y_fit, 'r', 'LineWidth', 1.5)
Cov = inv(Xmat' * Wmat * Xmat);   
sigma_a = sqrt(Cov(1,1));         
sigma_b = sqrt(Cov(2,2));         
legend([fitting_curve], strcat('slope = ', num2str(a), '$\pm$', num2str(sigma_a)))

%% Define functions
function vf = hadamard(i)
    I = eye(8);
    vf = I;
    vf(:,[i,i+4]) = vf(:,[i+4,i]);
end

function vf = phase(i)
    I = eye(8);
    vf = I;
    vf(i+4,:) = I(i,:)+I(i+4,:);
end

function vf = controlled(i,j) % i: controlled qubit, j: target qubit
    I = eye(8);
    vf = [];
    for k = 1:8
        vf(:,k) = I(:,k);
        if k == i
            vf(j,k) = 1;
        elseif k == j+4
            vf(i+4,k) = 1;
        end
    end
end

function sf = pauliprod(A,si)
    pos = find(si' == 1);
    sf = 0;
    for i = pos
        sf = sf + A(:,i);
    end
    sf = mod(sf,2);
end
