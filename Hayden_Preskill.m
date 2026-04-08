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

Q = [zeros(length(gate)/2), eye(length(gate)/2);eye(length(gate)/2), zeros(length(gate)/2)]; % check if gate is Clifford
symplectic = mod(gate*Q*gate',2);
check1 = sum(sum(symplectic ~= Q));

%% Main body
N_values = 8; % We need an N by N matrix of qubits

tf = 10; % final time

for N = N_values
    n = N^2; % in total we have N^2 qubits
    U_permute = zeros(2*n);
    U_permute_one_kind = zeros(n);
    sigma = zeros(N);
    P = zeros(2*n);

    %% We want to apply the four-qubit Clifford gate U_int on each cyclically connected subset
    D = kron(eye(2*n/8), gate); % gate D works as the U_int

    %% Find U_permute (and thus also need to find sigma first)
    v1 = (1:N)';
    for k = 1:N
        if mod(k,2) == 1
            v1(k,2) = (k+1)/2;
        else
            v1(k,2) = (N+v1(k,1))/2; % first column in v1 records [1;2;...;n], second column records result after permuting v1 according to the cyclic connectivity rule
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
        U_permute_one_kind(N*(kkk_row_ind - 1)+kkk_col_ind, N*(kkk_row_perm_ind - 1)+kkk_col_perm_ind) = 1;
    end
    U_permute = kron(eye(2), U_permute_one_kind); % we've built U_permute
    
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
    end    % Uncomment if want Rule 1 (ENDS HERE!)
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
            b = [v4(i,1) v4(i,2); v4(i,2) v4(i,1)+(-1)^(p_i+1); v4(i,1)+(-1)^(p_i+1) v4(i,2)+(-1)^(p_j+1); v4(i,2)+(-1)^(p_j+1) v4(i,1)];
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


    %% Define subsystem A
     % choice_of_qubit = randi([1,n]); % choose the reference qubit randomly
     choice_of_qubit = 39; % if you wanna choose a fixed reference qubit
     q_pool = [choice_of_qubit]; % q_pool will record all the qubits that the reference qubit could have spread to under random unitary dynamics
     ind_qubit_vec = zeros(n,1); 
     ind_qubit_vec(choice_of_qubit) = 1; 
     new_qubit_vec = [];

     gate_RUC = ones(4,4); 
     D_RUC = kron(eye(n/4), gate_RUC);
     P_RUC = zeros(n);
    [i, j] = meshgrid(1:N, 1:N); 
    v2_RUC = [i(:), j(:)]; 
    v3_RUC = [];
    a1_RUC = [1 1;1 4;4 4;4 1]; % [11; 14; 44; 41] % Uncomment if want Rule 1 (STARTS HERE!)
    a2_RUC = [2 2;2 3;3 3;3 2]; % [22; 23; 33; 32] 
    a3_RUC = [1 3;3 4;4 2;2 1]; % [13; 34; 42; 21]
    a4_RUC = [1 2;2 4;4 3;3 1]; % [12; 24; 43; 31]
    a_RUC = [a1_RUC;a2_RUC;a3_RUC;a4_RUC];
    for i = 1:N/4;
        v3_RUC(4^2*(i-1)+1:4^2*i,:) = a_RUC + 4*(i-1); % all 8x8 diagonal blocks containing both X and Z gates
    end % Uncomment if want Rule 1 (ENDS HERE!)
    % v3_RUC = [0,0]; % Uncomment if want Rule 2
    v4_RUC = setdiff(v2_RUC, v3_RUC,'rows');
    for i = 1:length(v4_RUC)
        if ~ismember(v4_RUC(i,:), v3_RUC, 'rows') % for each row in v4 this should always satisfy
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
            b = [v4_RUC(i,1) v4_RUC(i,2); v4_RUC(i,2) v4_RUC(i,1)+(-1)^(p_i+1); v4_RUC(i,1)+(-1)^(p_i+1) v4_RUC(i,2)+(-1)^(p_j+1); v4_RUC(i,2)+(-1)^(p_j+1) v4_RUC(i,1)]; 
            v3_RUC = [v3_RUC; b];
        end 
    end
    % v3_RUC(1,:) = []; % Uncomment if want Rule 2
    v2_RUC(:,3) = (v2_RUC(:,1)-1)*N + v2_RUC(:,2);
    v3_RUC(:,3) = (v3_RUC(:,1)-1)*N + v3_RUC(:,2);
    for k = 1:length(v2_RUC)
        P_RUC(v2_RUC(k,3),v3_RUC(k,3)) = 1; % converts the normal index set to the cyclic connectivity index set
    end

     t = 0;
     while t ~= tf
            new_qubit_vec = [];
            new_qubit_vec = ind_qubit_vec;
            new_qubit_vec = U_permute_one_kind * new_qubit_vec;
            new_qubit_vec = P_RUC * new_qubit_vec;
            new_qubit_vec = matprod_RUC(D_RUC, new_qubit_vec);
            new_qubit_vec = P_RUC' * new_qubit_vec;
            ind_qubit_vec = new_qubit_vec;
            t = t + 1;

            [qubit_row, qubit_col] = find(ind_qubit_vec ~= 0); 
            q_pool = [q_pool;qubit_row];
            q_pool = unique(q_pool,'stable');
            if length(q_pool) == N^2
                break
            end
     end
     q_pool = unique(q_pool,'stable'); % Now we get the q_pool as the set of qubits that the reference qubit could have spread to during the time evolution (recorded according to the first time it has appeared in the set)

    %% In the following we begin our main calculation/measurement for the H-P protocol
    saturated = [];
    non_saturated = [];

     for qqq = 1:length(q_pool) % set the size of q_pool to range from 1 to all the contained qubits
         q_vec = q_pool(1:qqq);
         interacting_qubits_A = [q_vec]; % define A to be the first qqq qubits in the subset of qubits that could have been spreaded by the refernece qubit during the time evolution
         ind_stabilizer_vec = eye(2*n); % Define initial set of independent stabilizers: For an arbitrary Clifford state, the stabilizers are composed of X and Z operators, so we essentially need both kinds
         new_stabilizer_vec = [];
         thisrank = [];
         qubit_row = [];
         qubit_col = [];

         %% Now we begin the time evolution of the system; recall that previously we loop wrt time just to get the "possibly infected pool"
         t = 0; % dynamics starts from t=0
         while t ~= tf
            new_stabilizer_vec = [];
            new_stabilizer_vec = ind_stabilizer_vec;
            new_stabilizer_vec = U_permute * new_stabilizer_vec;
            new_stabilizer_vec = P * new_stabilizer_vec;
            new_stabilizer_vec = paulimatprod(D,new_stabilizer_vec); % well we don't need it here, we can just apply mod(D * new_stabilizer_vec, 2) which is essentially the same. Somehow I defined a bunch of functions I forgot why
            new_stabilizer_vec = P' * new_stabilizer_vec; 
            ind_stabilizer_vec = new_stabilizer_vec;
            t = t + 1;

             k = qqq;
             r_vec = interacting_qubits_A(1:k);
             r = length(r_vec);
             q = length(q_vec); % Here one can realize that essentially r=q and r_vec = q_vec
    
             mat_K = ind_stabilizer_vec([r_vec,n+r_vec],setdiff(n+1:2*n,q_vec+n)); % operators in A^c
             mat_G = ind_stabilizer_vec([r_vec,n+r_vec],[q_vec,q_vec+n]); % operators in A
             R = rref([mat_K mat_G]); % maybe rref is unnecessary here
             if rank(R) == rank(mat_K) 
                 saturated = [saturated, [r;t;1]];
             else
                 non_saturated = [non_saturated, [r;t;0]];
             end  
         end
     end
end

if length(saturated) ~= 0
    p1 = scatter(saturated(1,:), saturated(2,:),10,'r','filled');
end
hold on
if length(non_saturated) ~= 0
    p2 = scatter(non_saturated(1,:), non_saturated(2,:),10,'b','filled');
end
hold on
p3 = plot(floor(N^2/3)*ones(1,tf),1:tf,'k-.','Linewidth', 1);
legend([p1, p2, p3], {'saturated', 'not saturated', strcat('$r=$', string(floor(N^2/3)))},'Location','best')
xlabel('$|B|$')
ylabel('$t$')
title(strcat('Rule 1, $N^2$=', string(n)))
% title(strcat('Rule 2, $N^2$=', string(n)))
pbaspect([1 0.8 1])
% exportgraphics(gcf,'/Users/yun/Desktop/research/fast scrambler/official_figures/single-layer/single_same_HP.png','Resolution',300)

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