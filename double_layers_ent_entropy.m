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
legend_entropy_vec = [];
N_values = [2^3, 2^4]; % N is a power of two
% N_values = [12,24]; % N contains some odd factor

t_final = 8; % final time

colors = lines(4+1);
colors(3,:) = [];
colors = colors(1:4,:);
color_ind = 0;

p_entropy = []; % this vector stores the curves for entanglement entropy growth
p_dashed = [];

figure(1)
for N_ind = 1:length(N_values) 
    N = N_values(N_ind);
    n = N^2; 
    L = zeros(n);
    U_permute = zeros(2*2*n);
    U_permute_one_kind = zeros(n);
    sigma = zeros(N);
    P = zeros(2*2*n);
    P_one_layer = zeros(2*n);

    %% We want to apply the four-qubit Clifford gate U_int on each cyclically connected subset
    D = kron(eye(2*2*n/8), gate); % gate D works as the U_int

    %% Build matrix L to take the transpose of the second layer
    first_vec = 1:n;
    first_matrix = reshape(first_vec, N, N)';
    second_vec = reshape(first_matrix, 1, n);
    for k = 1:n
        L(first_vec(k),second_vec(k)) = 1; 
    end
    L_RUC = blkdiag(eye(n), L);
    L = blkdiag(eye(n), L, eye(n), L);


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
    U_permute = kron(eye(4), U_permute_one_kind); % we've built U_permute

    %% Find index change permutation matrix P:
    [i, j] = meshgrid(1:4*N, 1:N); 
    v2 = [i(:), j(:)]; 
    v3 = [];
    % a1 = [1 1;1+N 4;4 4;4+N 1]; % [X11^1; X14^2; X44^1; X41^2]  % Uncomment if want Rule 1 (STARTS HERE!)
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
    %     v3(4^2*2*2*(i-1)+1:4^2*2*2*i,:) = a + 4*(i-1); % cyclic connectivity rule for all 8x8 diagonal blocks (i.e.containing both X and Z gates)
    % end         % Uncomment if want Rule 1 (ENDS HERE!)
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
    

    %% Define subsystem A
    choice_of_qubit_vec = [randi([1,2*N^2]),randi([1,2*N^2])]; % we pick one reference qubit here
    % A_size_vec = sort([randi([1,floor(2*N^2/6)]), randi([2*N^2/2-1,2*N^2/2])]); % We will set A as the subset of qubits that have interacted with the reference qubit during the time evoluton between t=0 and t=t_final
    A_size_vec = sort([2*N^2/4, randi([2*N^2/4+1,2*N^2/2])]); % set the sizes of A when only HALF of the qubits are involved in the interactions 

    gate_RUC = ones(4,4);
    D_RUC = kron(eye(2*n/4), gate_RUC);
    P_RUC = zeros(n);
    [i, j] = meshgrid(1:2*N, 1:N); 
    v2_RUC = [i(:), j(:)]; 
    v3_RUC = [];
    % a1_RUC = [1 1;1+N 4;4 4;4+N 1]; % [X11^1; X14^2; X44^1; X41^2] % Uncomment if want Rule 1 (STARTS HERE!)
    % a2_RUC = [2 2;2+N 3;3 3;3+N 2]; % [X22; X23; X33; X32]
    % a3_RUC = [1 3;3+N 4;4 2;2+N 1]; % [X13; X34; X42; X21]
    % a4_RUC = [1 2;2+N 4;4 3;3+N 1]; % [X12; X24; X43; X31]
    % a5_RUC = [1+N 1;1 4;4+N 4;4 1]; % [X11; X14; X44; X41]
    % a6_RUC = [2+N 2;2 3;3+N 3;3 2]; % [X22; X23; X33; X32]
    % a7_RUC = [1+N 3;3 4;4+N 2;2 1]; % [X13; X34; X42; X21]
    % a8_RUC = [1+N 2;2 4;4+N 3;3 1]; % [X12; X24; X43; X31]
    % a_RUC = [a1_RUC;a5_RUC;a2_RUC;a6_RUC;a3_RUC;a7_RUC;a4_RUC;a8_RUC];
    % for i = 1:N/4;
    %     v3_RUC(4^2*2*(i-1)+1:4^2*2*i,:) = a_RUC + 4*(i-1); % rule of rearrangement for interaction (i.e. change from normal set to predetermined set) for all 8x8 diagonal blocks (i.e.containing both X and Z gates)
    % end   % Uncomment if want Rule 1 (ENDS HERE!)
    v3_RUC = [0,0]; % Uncomment if want Rule 2
    v4_RUC = setdiff(v2_RUC, v3_RUC,'rows');
    for i = 1:length(v4_RUC)
        if ~ismember(v4_RUC(i,:), v3_RUC, 'rows') % for each row in v4 this should always be satisfied
            if mod(v4_RUC(i,1),2) == 0
                p_i = 0;
            else
                p_i = 1;
            end
            if mod(v4_RUC(i,2),2) == 0
                p_j = 0;
            else
                p_j = 1;
            end

            if v4_RUC(i,1) <= N
                b_RUC = [v4_RUC(i,1) v4_RUC(i,2); v4_RUC(i,2)+N v4_RUC(i,1)+(-1)^(p_i+1); v4_RUC(i,1)+(-1)^(p_i+1) v4_RUC(i,2)+(-1)^(p_j+1); v4_RUC(i,2)+(-1)^(p_j+1)+N v4_RUC(i,1)]; 
            end
            v3_RUC = [v3_RUC; b_RUC];
        end 
    end
    v3_RUC(1,:) = []; % Uncomment if want Rule 2
    [~, idx_v2_in_v3] = ismember(v2_RUC, v3_RUC, 'rows'); % check if every row in v2 appears in v3
    rows_not_in_v3 = v2(idx_v2_in_v3 == 0, :);
    [uniqueRows, ~, idx] = unique(v3, 'rows'); % check if every row in v3 is different
    rowCounts = histcounts(idx, 1:max(idx)+1);
    numRepeatedRows = sum(rowCounts > 1);

    v2_RUC(:,3) = (v2_RUC(:,1)-1)*N + v2_RUC(:,2);
    v3_RUC(:,3) = (v3_RUC(:,1)-1)*N + v3_RUC(:,2);
    
    for k = 1:length(v2_RUC)
        P_RUC(v2_RUC(k,3),v3_RUC(k,3)) = 1; % converts the normal index set to the cyclic connectivity index set
    end

    for A_size_ind = 1:length(A_size_vec)
        t = 0;
        color_ind = color_ind + 1;
        c = colors(color_ind,:);
        ind_qubit_vec = zeros(2*n,1); 
        new_qubit_vec = [];
        A_size = A_size_vec(A_size_ind);
        choice_of_qubit = choice_of_qubit_vec(A_size_ind);
        A = [choice_of_qubit];
        ind_qubit_vec(choice_of_qubit) = 1;
        while t ~= t_final
            new_qubit_vec = [];
            new_qubit_vec = ind_qubit_vec;
            new_qubit_vec = L_RUC * new_qubit_vec;
            new_qubit_vec = kron(eye(2), U_permute_one_kind) * new_qubit_vec;
            new_qubit_vec = L_RUC' * new_qubit_vec;
            new_qubit_vec = P_RUC * new_qubit_vec;
            new_qubit_vec = matprod_RUC(D_RUC,new_qubit_vec);
            new_qubit_vec = P_RUC' * new_qubit_vec;
            ind_qubit_vec = new_qubit_vec;
            t = t + 1;

            [qubit_row, qubit_col] = find(ind_qubit_vec ~= 0);
            A = [A;qubit_row];
            A = unique(A,'stable');
            if length(A) == 2*N^2
                break
            end
        end
        A = unique(A,'stable'); % Now we find the subset of qubits that have interacted with the reference qubit during the time evoluton between t=0 and t=t_final
        A = A(1:A_size); % Set A to contain the first number |A| of qubits
        Ac = setdiff(1:2*n,A); % complement set A^c

        %% Now we begin the time evolution of the system
        t = 0;
        time_vec = [0];
        ind_stabilizer_vec = [zeros(2*n);eye(2*n)]; % for ground state, a set of independent stabilizers is just all the Z operators
        new_stabilizer_vec = [];
        entropy_vec = []; 
        projA = ind_stabilizer_vec;
        projA([Ac, Ac+2*n],:) = 0; % projection of each independent stabilizer onto A
        rank_projA = rank(projA');
        SA = rank_projA - length(A);
        entropy_vec = [entropy_vec, SA]; % entanglement entropy
        while t ~= t_final
            new_stabilizer_vec = [];
            new_stabilizer_vec = ind_stabilizer_vec;
            new_stabilizer_vec = L * new_stabilizer_vec;
            new_stabilizer_vec = U_permute * new_stabilizer_vec;
            new_stabilizer_vec = L' * new_stabilizer_vec;
            new_stabilizer_vec = P * new_stabilizer_vec;
            new_stabilizer_vec = paulimatprod(D,new_stabilizer_vec);
            new_stabilizer_vec = P' * new_stabilizer_vec;
            ind_stabilizer_vec = new_stabilizer_vec;

            projA = ind_stabilizer_vec;
            projA([Ac, Ac+2*n],:) = 0;
            rank_projA = rank(projA');
            SA = rank_projA - length(A);
            entropy_vec = [entropy_vec, SA];
            t = t + 1;
            time_vec = [time_vec, t];
        end
        if length(A) > (2*N^2)/2 %% if ALL of the qubits are involved in the interaction
        % if length(A) > (2*N^2)/4 %% if HALF of the qubits are involved in the interaction
            a = log(A_size)/log(4);
            b = interp1(time_vec, entropy_vec/A_size, a, 'linear', 'extrap');  
            p_dashed = [p_dashed, plot(time_vec,((2*N^2/2-length(A))/length(A))*ones(1,length(time_vec)),'Color', c,'LineStyle', '--')];
            hold on
            l1 = plot(nan, nan, 'k--');
            hold on

            a = log(N^2 - A_size)/log(4);
            b = interp1(time_vec, entropy_vec/A_size, a, 'linear', 'extrap'); 
            plot(a, b, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', c, 'LineWidth', 1.5);
            hold on
            p_log = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
            hold on
            p_entropy = [p_entropy, plot(time_vec,entropy_vec/A_size,'Color', c)];
            hold on
            legend_entropy_vec = [legend_entropy_vec,strcat(strcat('$2N^2=$', string(2*N^2)), ',',' ', strcat('$|A|=$', string(length(A))))];
            hold on
        else
            a = log(A_size)/log(4);
            b = interp1(time_vec, entropy_vec/A_size, a, 'linear', 'extrap');   
            plot(a, b, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', c, 'LineWidth', 1.5); % when t=log_4 (N^2), what's the corresponding entanglement entropy
            hold on
            p_log = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
            hold on
            p_entropy = [p_entropy, plot(time_vec,entropy_vec/length(A),'Color', c)];
            hold on
            legend_entropy_vec = [legend_entropy_vec,strcat(strcat('$2N^2=$', string(2*N^2)), ',',' ', strcat('$|A|=$', string(length(A))))];
            hold on
        end
    end
end
set(groot,'defaulttextinterpreter','latex','defaulttextFontWeight','bold')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultColorbarTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
legend([p_log, p_entropy],['$\log_4{|A|}$',legend_entropy_vec],'Location','best'); %% if ALL of the qubits are involved in the interaction
% legend([p_log, l1, p_entropy],['$\log_4{|A|}$','$\frac{|A^c|}{|A|}$',legend_entropy_vec],'Location','best'); %% if HALF of the qubits are involved in the interaction
% title('Rule 1')
title('Rule 2')
xlabel('$t$')
ylabel('$\frac{S_A(t)}{|A|}$')
% exportgraphics(gcf,'/Users/yun/Desktop/research/fast scrambler/official_figures/double-layer/double_diff_EE_odd.png','Resolution',300)

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

function mf = paulimatprod(A,M)
    mf = [];
    [numrow, numcol] = size(M);
    for i = 1:numcol
        mf = [mf, pauliprod(A,M(:,i))];
    end
end

function mf_RUC = matprod_RUC(B,K)
    mf_RUC = [];
    [numrow_RUC, numcol_RUC] = size(K);
    for i = 1:numcol_RUC
        mf_RUC = [mf_RUC, B*K(:,i)];
    end
end