function [e_fastsum] = hiFastsum(x_fastsum, y_fastsum, dim, eps_B, B3, weighting, boolFast)
%HIFASTSUM Calculates Fastums solution for 3D-Hilbert Transform
%   x_fastsum   start points for fastsum
%   y_fastsum   finish points for fastsum (actual positions for result are
%               x_fastsum)
%   dim         the dimension of the result (x=1, y=2)
%   eps_B       the outer Boundary
%   B3          the coefficient that shall be transformed
%   weighting   the geometrical wheigting of the coefficient
%   boolFast    1: fast compution, 2: direct computation
%RETURN:
%   e_fastsum   B1/B2, solution of the fastsum algorithm


    % coefficients
    alpha1 = B3.*weighting;
    alpha2 = B3.*x_fastsum(:,dim).*weighting;
    
    % Initialize Fastsum parameters
    d = 2;       % number of dimensions
    N = length(x_fastsum);       % number of source knots
    M = length(y_fastsum);       % number of target knots
    kernel = 'one_over_cube';
    c = 0;          % kernel parameter (not needed)
    p = 2;          % degree of smoothness of regularization  %3 recommended?
    flags = 1;      % flags (could be EXACT_NEARFIELD or NEARFIELD_BOXES)
    n = 64;         % bandwidth in frequency domain for NFFT (every interval of 1/n shall at most contain 2N/n points) 600 recommended?
    % n = round(sqrt(p*N))*1;        % bandwidth in frequency domain for NFFT
    eps_I = max(p/n,eps_B);    % inner boundary, nearfield radius
    % eps_I = eps_I / 10;
    m = p;          % window cut-off parameter for NFFT (for p<=5 m=p, for large p m=5 sufficient
    nn_oversampled = 2*n; % oversampled bandwidth in frequency domain for NFFT
    
    %% Perform fastsum
    time = 0;
    
    % First Sum
    disp(['Start fastsum 1 (time: ', num2str(time),'s)']);
    
    tic;
    plan=fastsum_init(d,kernel,c,flags,n,p,eps_I,eps_B);
    fastsum_set_x(plan,x_fastsum,nn_oversampled,m)
    fastsum_set_alpha(plan,alpha1)
    fastsum_set_y(plan,y_fastsum,nn_oversampled,m)

    if boolFast == 0% direct computation
        fastsum_trafo_direct(plan)   
        e1_fastsum = fastsum_get_f(plan);
    else % fast computation
        tic;
        fastsum_trafo(plan)         
        e1_fastsum = fastsum_get_f(plan);
    end
    fastsum_finalize(plan)
    time = time + toc;

    % Second Sum
    disp(['Start fastsum 2 (time: ', num2str(time), 's)']);
    tic;
    plan=fastsum_init(d,kernel,c,flags,n,p,eps_I,eps_B);
    fastsum_set_x(plan,x_fastsum,nn_oversampled,m)
    fastsum_set_y(plan,y_fastsum,nn_oversampled,m)
    fastsum_set_alpha(plan,alpha2)
    if boolFast == 0% direct computation
        fastsum_trafo_direct(plan)   % direct computation
        e2_fastsum = fastsum_get_f(plan);
    else
        fastsum_trafo(plan)         % fast computation
        e2_fastsum = fastsum_get_f(plan);
    end
    fastsum_finalize(plan)
    time = time + toc;

    %% Post Calculation
    e_fastsum = (y_fastsum(:,dim).*real(e1_fastsum) - real(e2_fastsum))/(-2*pi);
    disp(['End fastsums    (time: ', num2str(time), 's)']);

    % When target nods differ from start nodes, interpolate back to start
    % nodes:
    if ~isequal(x_fastsum,y_fastsum)
        tic;
        F = scatteredInterpolant(y_fastsum(:,1),y_fastsum(:,2),e_fastsum);
        F.Method = 'linear';
        e_fastsum = F(x_fastsum(:,1),x_fastsum(:,2));   
        time = time + toc;
        disp(['End interpolation    (time: ', num2str(time), 's)']);
    end
end

