% 1. choose N divisible by 4.  initialize N x N matrix X with one random 1 and the rest of the entries zero.
% 2. apply a random permutation matrix Sigma to X:   X -> Sigma^T X Sigma.   (the way to generate random permutations is discussed here: https://www.mathworks.com/matlabcentral/answers/347033-generate-random-permutation-matrix, but i am suspicious whether anyone s actual code is working.  you could check if they do)
% 3. build a subroutine that checks every cluster of the matrix in the same way that we did on the board today for the diagonal/off-diagonal blocks.  for each cluster with matrix elements {c1,c2,c3,c4} , send X(c1) -> X(c1)+X(c2)+X(c3)+X(c4).
% 4. for each matrix element which is larger than 0 in X, set it to 1.   evaluate and record size = \sum_{i,j} X(i,j)
% 5. repeat steps 2-4 for the desired number of time steps
% 6. plot size as a function of time, compare to exponential growth (easiest by plotting log(size) vs. t)
%%
close all
clear all
clc

set(groot,'defaulttextinterpreter','latex','defaulttextFontWeight','bold')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultColorbarTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
set(groot,'defaultAxesFontSize',20)
set(groot,'defaultLineLineWidth',1.5)

legend_vec = [];
saturation_time = [];
N_values = [2^(3), 12, 2^(4),40,2^(6)]; % We need an NxN matrix of qubits

colors = lines(length(N_values)+1); 
% colors(3,:) = [];
colors = colors(1:length(N_values),:);

p_RUC = [];
p_log = [];
figure(1)
for N_ind = 1:length(N_values)
    N = N_values(N_ind);
    initial_i = randi([1,N]);  %% row index of the initially infected qubit
    initial_j = randi([1,N]); %% column index of the initially infected qubit
    % initial_i = 22;   
    % initial_j = 29; 
    A = zeros(N);
    A(initial_i,initial_j) = 1;
    perm_A = A;
    size_vec = [1]; 
    time_vec = [0];
    t = 0; % initial time of the evolution
    t_final = 20; % final time of the evolution
    diag_block = [];
    for b = 1:(N/4)
        diag_block(b,:) = ((b-1)*4+1):((b-1)*4+4); % define the 4x4 diagonal blocks
    end
    
    while t ~= t_final
        t = t + 1;
        time_vec = [time_vec, t];
            [ones_x, ones_y] = find(perm_A == 1);
            for k = 1:length(ones_x)
                if mod(ones_x(k),2) == 1
                    ones_x(k) = (ones_x(k)+1)/2;
                else
                    ones_x(k) = (N+ones_x(k))/2;
                end
                if mod(ones_y(k),2) == 1
                    ones_y(k) = (ones_y(k)+1)/2;
                else
                    ones_y(k) = (N+ones_y(k))/2;
                end
            end
            perm_A(find(perm_A == 1)) = 0;
            perm_A((ones_y-1)*N+ones_x) = 1;
        [ones_i, ones_j] = find(perm_A == 1);
        last_size = sum(sum(perm_A));
        R1 = {'(1,1)','(1,4)','(4,4)','(4,1)'}; % define the cyclic connectivity rule for 4x4 "digaonal" blocks
        R2 = {'(2,2)','(2,3)','(3,3)','(3,2)'};
        R3 = {'(1,3)','(3,4)','(4,2)','(2,1)'};
        R4 = {'(1,2)','(2,4)','(4,3)','(3,1)'};
        R = [R1; R2; R3; R4];
        R = string(R);
        for k = 1:last_size
            [r1,c1] = find(diag_block == ones_i(k));
            [r2,c2] = find(diag_block == ones_j(k));
            % if r1 == r2 % Uncomment if want Rule 1 (STARTS HERE!). This line is saying, if (i,j) is in some 4x4 "diagonal" block,
            %     I = mod(ones_i(k),4);
            %     J = mod(ones_j(k),4);
            %     if I == 0 
            %         I = 4;
            %     end
            %     if J == 0
            %         J = 4;
            %     end
            %     [row,col] = find(R == ['(',num2str(I),',',num2str(J),')']);
            %     for l = 1:length(R(row,:))
            %         this_pair = convertStringsToChars(R(row,l));
            %         perm_A((ceil(min(ones_i(k),ones_j(k)) / 4) - 1) * 4 + str2double(this_pair(2)),(ceil(min(ones_i(k),ones_j(k)) / 4) - 1) * 4 + str2double(this_pair(4))) = 1;
            %     end
            % else % Uncomment if want Rule 1 (ENDS HERE!)
                if mod(ones_i(k),2) == 0 % define the cyclic connectivity rule for 2x2 "off-digaonal" blocks
                    p_i = 0;
                else
                    p_i = 1;
                end
                if mod(ones_j(k),2) == 0
                    p_j = 0;
                else
                    p_j = 1;
                end
                perm_A(ones_j(k), ones_i(k) + (-1) ^ (p_i + 1)) = 1;
                perm_A(ones_i(k) + (-1) ^ (p_i + 1), ones_j(k) + (-1) ^ (p_j + 1)) = 1;
                perm_A(ones_j(k) + (-1) ^ (p_j + 1), ones_i(k)) = 1;
            % end % Uncomment if want Rule 1
        end
        size = sum(sum(perm_A)); % number of "ones" / infected qubits of the system
        size_vec = [size_vec, size];
    end
    c = colors(N_ind,:);
    saturation_time = [saturation_time, t]; % time that the size saturates at N^2
    p_RUC = [p_RUC, semilogy(time_vec, size_vec, 'Color', c)]; 
    hold on
    plot(time_vec, N^2 * ones(1,length(time_vec)),'--','Color', c,'HandleVisibility','off') % dashed upper bound of the size
    hold on
    a = log(N^2)/log(4);
    b = interp1(time_vec, size_vec, a, 'linear', 'extrap');  
    p2 = semilogy(a, b, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', c, 'LineWidth', 1.5);
    p_log = plot(NaN, NaN, 'o', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'none', 'LineWidth', 1.5);
    hold on
    legend_vec = [legend_vec, strcat('N=', string(N))];
    hold on
end
legend([p_log, p_RUC],['$\log_4{N^2}$',legend_vec],'Location','best')
% title('Rule 1')
title('Rule 2')
set(groot,'defaulttextinterpreter','latex','defaulttextFontWeight','bold')
set(groot,'defaultAxesTickLabelInterpreter','latex')
set(groot,'defaultColorbarTickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex')
xlabel('$t$')
ylabel('$n(t)$')
% lgd.Color = 'none';   % no background
% lgd.Box   = 'off'; 
% exportgraphics(gcf,'size_time_2^10.png','Resolution',300)

