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
gate = phase(4) * controlled(4,1) * controlled(1,2) * controlled(2,3) * controlled(3,4) * hadamard(1);
gate = mod(gate,2);

% Q = [zeros(length(gate)/2), eye(length(gate)/2);eye(length(gate)/2), zeros(length(gate)/2)]; % check if gate is Clifford
% symplectic = mod(gate*Q*gate',2);
% check1 = sum(sum(symplectic ~= Q));

%% Main body
legend_vec = [];
rate_vec = [];
N_values = [2^(3), 2^(4),2^(5)]; % N is a power of two
% N_values = [12, 20, 36,52]; % N contains some odd factor

colors = lines(length(N_values)+1); 
colors(3,:) = [];
colors = colors(1:length(N_values),:);
t_final = 40; % final time

p_log = [];
p_OSG = [];
legend_OSG_vec = [];

figure(1)
for N_ind = 1:length(N_values) % loop wrt N
    N = N_values(N_ind);
    n = N^2;
    L = zeros(n);
    U_permute = zeros(2*2*n);
    U_permute_one_kind = zeros(n);
    sigma = zeros(N);
    P = zeros(2*2*n);
    P_one_layer = zeros(2*n);
    D = eye(2*n);
    pauli_string = zeros(2*2*n,1);
    initial = randi([1,2*2*n]); % random initial single-qubit operator
    % initial = 1217; % if you want an assigned initial single-qubit operator
    pauli_string(initial,1) = 1;
    pairs_vec = [1];
    time_vec = [0];

    %% We want to apply the four-qubit Clifford gate U_int on each cyclically connected subset
    D = kron(eye(2*2*n/8), gate); % gate D works as the U_int

    %% Build matrix L to take the transpose of the second layer
    first_vec = 1:n;
    first_matrix = reshape(first_vec, N, N)';
    second_vec = reshape(first_matrix, 1, n);
    for k = 1:n
        L(first_vec(k),second_vec(k)) = 1;
    end
    L = blkdiag(eye(n), L, eye(n), L);


    %% Find U_permute (and thus also need to find sigma first)
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
        U_permute_one_kind(N*(kkk_row_ind-1)+kkk_col_ind, N*(kkk_row_perm_ind-1)+kkk_col_perm_ind) = 1; % we build U_permute
    end
    U_permute = kron(eye(4), U_permute_one_kind);

    %% Find index change permutation matrix P:
    [i, j] = meshgrid(1:4*N, 1:N); 
    v2 = [i(:), j(:)]; 
    v3 = [];
    % a1 = [1 1;1+N 4;4 4;4+N 1]; % [X11^1; X14^2; X44^1; X41^2] % Uncomment if want Rule 1 (STARTS HERE!)
    % a1 = [a1; a1(:,1)+2*N a1(:,2)]; % [X11^1; X14^2; X44^1; X41^2; Z11^1; Z14^2; Z44^1; Z41^2]
    % a2 = [2 2;2+N 3;3 3;3+N 2]; % [X22; X23; X33; X32]
    % a2 = [a2; a2(:,1)+2*N a2(:,2)]; % [X22; X23; X33; X32; Z22; Z23; Z33; Z32]
    % a3 = [1 3;3+N 4;4 2;2+N 1]; % [X13; X34; X42; X21]
    % a3 = [a3; a3(:,1)+2*N a3(:,2)]; % [X13; X34; X42; X21; Z13; Z34; Z42; Z21]
    % a4 = [1 2;2+N 4;4 3;3+N 1]; % [X12; X24; X43; X31]
    % a4 = [a4; a4(:,1)+2*N a4(:,2)]; % [X12; X24; X43; X31; Z12; Z24; Z43; Z31]
    % a5 = [1+N 1;1 4;4+N 4;4 1]; % [X11; X14; X44; X41]
    % a5 = [a5; a5(:,1)+2*N a5(:,2)]; % [X11; X14; X44; X41; Z11; Z14; Z44; Z41]
    % a6 = [2+N 2;2 3;3+N 3;3 2]; % [X22; X23; X33; X32]
    % a6 = [a6; a6(:,1)+2*N a6(:,2)]; % [X22; X23; X33; X32; Z22; Z23; Z33; Z32]
    % a7 = [1+N 3;3 4;4+N 2;2 1]; % [X13; X34; X42; X21]
    % a7 = [a7; a7(:,1)+2*N a7(:,2)]; % [X13; X34; X42; X21; Z13; Z34; Z42; Z21]
    % a8 = [1+N 2;2 4;4+N 3;3 1]; % [X12; X24; X43; X31]
    % a8 = [a8; a8(:,1)+2*N a8(:,2)]; % [X12; X24; X43; X31; Z12; Z24; Z43; Z31]
    % a = [a1;a5;a2;a6;a3;a7;a4;a8];
    % for i = 1:N/4;
    %     v3(4^2*2*2*(i-1)+1:4^2*2*2*i,:) = a + 4*(i-1); 
    % end % Uncomment if want Rule 1 (ENDS HERE!)
    v3 = [0,0]; % Uncomment if want Rule 2
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
            if v4(i,1) <= N
                b = [v4(i,1) v4(i,2); v4(i,2)+N v4(i,1)+(-1)^(p_i+1); v4(i,1)+(-1)^(p_i+1) v4(i,2)+(-1)^(p_j+1); v4(i,2)+(-1)^(p_j+1)+N v4(i,1)]; 
            end
            b = [b; b(:,1)+2*N b(:,2)];
            v3 = [v3; b];
        end 
    end
    v3(1,:) = []; % Uncomment if want Rule 2
    [~, idx_v2_in_v3] = ismember(v2, v3, 'rows'); % check if every row in v2 appears in v3
    rows_not_in_v3 = v2(idx_v2_in_v3 == 0, :);
    [uniqueRows, ~, idx] = unique(v3, 'rows'); % check if every row in v3 is different
    rowCounts = histcounts(idx, 1:max(idx)+1);
    numRepeatedRows = sum(rowCounts > 1);

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

        pauli_string = L * pauli_string;
        pauli_string = U_permute * pauli_string;
        pauli_string = L' * pauli_string;
        pauli_string = P * pauli_string;
        pauli_string = pauliprod(D,pauli_string);
        pauli_string = P' * pauli_string;
        
        thissss_pauli = find(pauli_string == 1);
        
        pairs = sum(pauli_string(1:length(pauli_string)/2) | pauli_string((length(pauli_string)/2+1): length(pauli_string))); % number of single-qubit Pauli operators in the Pauli string
        pairs_vec = [pairs_vec, pairs];
    end

    set(groot,'defaulttextinterpreter','latex','defaulttextFontWeight','bold')
    set(groot,'defaultAxesTickLabelInterpreter','latex')
    set(groot,'defaultColorbarTickLabelInterpreter','latex')
    set(groot,'defaultLegendInterpreter','latex')

    c = colors(N_ind,:);
    p_OSG = [p_OSG, semilogy(time_vec, pairs_vec,'Color', c)]; 
    a = log(2*N^2)/log(4);
    b = interp1(time_vec, pairs_vec, a, 'linear', 'extrap');  
    hold on
    plot(a, b, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', c, 'LineWidth', 1.5); % when t=log_4 (N^2), what's the corresponding operator size
    p_log = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
    hold on
    plot(time_vec, 3*2*n/4 * ones(1,length(time_vec)),'--','Color', c,'HandleVisibility','off') %% the typical operator size at late times 
    set(gca,'yscale','log')
    legend_OSG_vec = [legend_OSG_vec, strcat('$2N^2$=', string(2*n))];
    hold on
end

%% Plotting
xlabel('$t$')
ylabel('$n(t)$')
legend([p_log, p_OSG],['$\log_4{2N^2}$',legend_OSG_vec],'Location','best')
% title('Rule 1')
title('Rule 2')
% exportgraphics(gcf,'/Users/yun/Desktop/research/fast scrambler/official_figures/double-layer/double_diff_OSG.png','Resolution',300)


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
