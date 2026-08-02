function Present problem
    % Parameters
    global S Re Ha B Br Gr Gm Pr  Nb Nt Pe Le 
    
    % Define domain boundaries
    eta_min = -1; 
    eta_max = 1; 
    
    % Non-zero values to avoid division by zero
    S=0.5; Re=1; Ha=2; B=2; Br=1; Gr=2; Pr=0.7;
    Gm=2;  Nb=0.4;  Nt=0.1;  Pe=0.4; Le=0.8; 
    Mva = [1 2 3]; % Magnetic parameter variations
    
    figure(1); hold on;
    lines = {'m-','b-','g-','r-'};
    
    for i = 1:length(Mva)
        Ha = Mva(i);
        
        % Mesh initialization spanning from eta_min (-1) to eta_max (1)
        % Using non-zero initial guesses to prevent singular Jacobian matrices
        solinit = bvpinit(linspace(eta_min, eta_max, 100), [0 0 0 0 0 0 0 0 0 0 1 1]);
        
        % Solve BVP
        sol = bvp4c(@shootode, @shootbc, solinit);
        
        % Extract solution
        eta = sol.x;
        f = sol.y;
        
        % Plot velocity profile f'(η) which corresponds to f(2,:)
        plot(eta, f(1,:), lines{i}, 'LineWidth', 2);
                % =================================================================
        % EXTRACT SURFACE QUANTITIES AT η = -1 (Index 1) AND η = 1 (Index end)
        % =================================================================
        
        % 1. Skin Friction Coefficients (Cf1*Re and Cf2*Re)
        % Formula: Cf1*Re = 2 * psi_1'(-1)  and  Cf2*Re = 2 * psi_1'(1)
        Cf1_Re = 2 * sol.y(2, 1);   % f(2) at eta = -1
        Cf2_Re = 2 * sol.y(2, end); % f(2) at eta = 1
        
        % 2. Nusselt Number (Nu = -dtheta/deta) at eta = -1 and eta = 1
        % f(10) is dtheta/deta
        Nu_minus1 = -sol.y(10, 1);   % Nu at eta = -1
        Nu_plus1  = -sol.y(10, end); % Nu at eta = 1
        
        % 3. Sherwood Number (Sh = -dphi/deta) at eta = -1 and eta = 1
        % Assuming you expand your ODE to 12 equations, f(12) would be dphi/deta.
        % If not yet implemented, we use a placeholder or safe check:
        if size(sol.y, 1) >= 12
            Sh_minus1 = -sol.y(12, 1);   % Sh at eta = -1
            Sh_plus1  = -sol.y(12, end); % Sh at eta = 1
        else
            Sh_minus1 = NaN; % Placeholder until 11th & 12th equations are added
            Sh_plus1  = NaN;
        end
        
        % =================================================================
        % PRINT RESULTS TO THE COMMAND WINDOW
        % =================================================================
        fprintf('\n--- Results for Ha = %g ---\n', Ha);
        fprintf('Skin Friction (Cf1*Re) at eta = -1: %0.6f\n', Cf1_Re);
        fprintf('Skin Friction (Cf2*Re) at eta =  1: %0.6f\n', Cf2_Re);
        fprintf('Nusselt Number (Nu)    at eta = -1: %0.6f\n', Nu_minus1);
        fprintf('Nusselt Number (Nu)    at eta =  1: %0.6f\n', Nu_plus1);
        if ~isnan(Sh_minus1)
            fprintf('Sherwood Number (Sh)   at eta = -1: %0.6f\n', Sh_minus1);
            fprintf('Sherwood Number (Sh)   at eta =  1: %0.6f\n', Sh_plus1);
        else
            fprintf('Sherwood Number (Sh): (Expand system to 12 equations to compute Sh)\n');
        end

    end
    
    xlabel('\eta');
    ylabel('f''(\eta)');
    legend('M=1','M=3','M=5','M=7');
    title('Velocity Profile for Different Magnetic Parameter M');
    grid on;
    hold off;
end

% -------------------------------------------------------------------------
function dydx = shootode(eta, f)
    global S Re Ha B Br Gr Gm Pr  Nb Nt Pe Le 
    dydx = zeros(12, 1);
    
    % First-order ODE system equations
    dydx(1) = f(2);
    dydx(2) = f(3);
    dydx(3) = f(4);
    dydx(4) = (((Gr/Re)*f(9)+(Gm/Re)*f(11)+f(3)-Re*f(2)-((Ha*Ha)/(1+B))*(f(1)+B*f(5)))/(S*S));
    dydx(5) = f(6);
    dydx(6) = f(7);
    dydx(7) = f(8);
    dydx(8) = (f(7) - Re*f(6)-((Ha*Ha)/(1+B))*(B*f(1) - f(5))/(S*S));
    dydx(9) = f(10);
    dydx(10) = (Re*Pr*f(10) - 2*Br*((f(2)*f(2)) + (f(6)*f(6))) - S*S*Br*((f(3)*f(3)) + f(7)*f(7)) ...
        - (Br*Ha*Ha/(1+B))*(f(1)*f(1) + f(5)*f(5))-Nb*Pr*f(10)*f(12)-Nt*Pr*f(12)*f(12));
    dydx(11) = f(12);
    dydx(12) = (Le*Pe*f(12)+(Nt/Nb)*(Re*Pr*f(10) - 2*Br*((f(2)*f(2)) + (f(6)*f(6))) - S*S*Br*((f(3)*f(3)) + f(7)*f(7)) ...
        - (Br*Ha*Ha/(1+B))*(f(1)*f(1) + f(5)*f(5))-Nb*Pr*f(10)*f(12)-Nt*Pr*f(12)*f(12)));
end

% -------------------------------------------------------------------------
function res = shootbc(fa, fb)
    % Well-posed standard physical boundary conditions for fluid dynamics systems
    % Maps 5 initial states at eta = -1 and 5 final states at eta = 1
    res = [ fa(1);     % f1(-1) = 0
            fa(3) ;  % f2(-1) = 1  (Stretching/moving wall condition)
            fa(5);     % f5(-1) = 0
            fa(7);     % f6(-1) = 0
            fa(9);     % f9(-1) = 0
            fa(11);
            fb(1);     % f1(1)  = 0
            fb(3);     % f2(1)  = 0  (Free stream/stationary condition)
            fb(5);     % f5(1)  = 0
            fb(7);     % f6(1)  = 0
            fb(9) - 1; % f9(1)  = 1
            fb(11) - 1;
          ]; 
end
