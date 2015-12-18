function [U] = NeumannSolution(p, t, e)
%
% Purpose:  Solves  - u_tt - laplace u = f
%           in a domain  D,  
%           with boundary conditions
%           dn grad u = 0, where  n  is the outer
%           normal to the boundary.
%   
%        p  2 x n  vector of node coordinates, that is  p(1, nn)
%           and  p(2, nn)  are the  x1  and  x2  coordinates,
%           respectively, of node number  nn. 
%
%        t  4 x m  matrix of element node numbers  t(1, el) = nn1,
%           t(2, el) = nn2,  and  t(3, el) = nn3, 1 <= el <= m,
%           pointing to corresponding node coordinates  p(:, nn1),
%           p(:, nn2)  and  p(:, nn3),  and subdomain reference
%           tags  t(4, el),  useful if data are given by different
%           expressions in different parts of the domain.
%
%        e  7 x q  matrix where  e(:, bel)  holds the following
%           information about bdry element  bel:  e(1:2, bel)  are
%           the node numbers (pointers to  p) of the two (endpoint)
%           nodes of the bdry element,  e(3:4, bel)  are the relative
%           positions of these two nodes on their bdry segment, with
%           reference number  e(5, bel),  and  e(6:7, bel) are the
%           reference tags of the domains/triangle elements tags to
%           the left and right, respectively, of the bdry element.
%           The exterior region has reference number 0.
%
%        U  n x 1  matrix , the solution values
%           at node points  p . 
%           
%
    
    %% Initialize
    
    % Settings
    duration = 10;    % Duration in second
    f = @(x, time) 0;   % Load function
    forceframes = 5;    % Number of time steps to apply forces


    % Precompute triangle areas for better performance
    t = [t; zeros(4,size(t,2))];
    for el = 1:size(t,2)
        NodeCoords = p(:,t(1:3,el));
        t(5,el) = area(NodeCoords);
    end

    % Initialize
    % To get an automatically adjusting timestep, let the time constant k
    % be the square of the smallest distance along the boundary.
    k = min(sum((p(:,e(1,:))-p(:,e(2,:))).^2));  % Minimum of squared distances
    k = k;

    iters = floor(duration/k)  % Duration in frames

    n = size(p, 2);            % number of points
    x1 = p(1,:)';              % x-parameter
    x2 = p(2,:)';              % y-parameter


    % Result matrix
    U = zeros(n,iters);        % Initialize result matrix
    
    %  Compute matrices and vectors  %
    M = MassMatrix(p,t);     % Compute mass matrix
    S = StiffMatrix(p,t);    % Compute stiffness matrix
    F = LoadVector(p,t,f,0); % Compute initial load vector

    Minv = pinv(M);          % Calculate the pseudo inverse once (performance)

    M = sparse(M);           % Convert M to a sparse matrix
    S = sparse(S);           % Convert S to a sparse matrix 

    [u_0, v_0] = bounds(x1,x2);  % Set boundary conditions

    % Computer initial values (from instructions), a backwards step.
    a_0    = M\(F - S*u_0);
    U(:,1) = u_0 - k*v_0 + k^2/2*a_0;
    U(:,2) = u_0;
    
    
    %% Animate interactive solution
    figure(1); set(gcf, 'DoubleBuffer','on');

    ph = pdeplot(p, e, t, 'zdata', U(:,1), 'xydata', U(:,1),...
            'mesh','off','colorbar','off','colormap','cool');
    grid on;
    %axis equal;
    set(gca,'ZLim',[-4 4]);
    view(0,90);
    
    loadTimeOut = 0;  % Countdown to see if load is active
    updateGraph = 1;  % Sould we update graphics?
    
    for i = 2:iters;
        if mod(i,100)==0
            disp(iters-i);    % Display iterations left
        end
        
        loadTimeOut = loadTimeOut - 1;
        
        % Update graph (use mod to change how often, affects performance)
        if updateGraph == 1 && mod(i,100)==0
            % Replace data
            Z = U(:,i);
            %Z(p(1,:)<0.02) = 0;

            set(ph, 'zdata', [U(t(1,:),i)';U(t(2,:),i)';U(t(3,:),i)']);
            set(ph, 'cdata', [Z(t(1,:))';Z(t(2,:))';Z(t(3,:))']);
            drawnow;
            
            % Interactive stuff
            %
            %  Click/grab point with the rotate tool and click
            %    f - Very small force
            %    d - medium force
            %    s - very strong force
            %
            %    q - exit loop
            %    v - disable graphics
            %    
            %    a - used as a hack to avoid missing input
            %
            pt = ptinxy();  % Last point clicked (in xy-plane)

            kkey = get(gcf, 'CurrentCharacter');
            
            if ~strcmp(kkey, 'a')
                loadTimeOut = forceframes+1;
            end

            if strcmp(kkey, 'f')
                f=@(x,T) 0.05*(T<((forceframes+i)*k))/k* ...
                         (((x(1)-pt(1,1)).^2+(x(2)-pt(1,2)).^2)<0.1);
            elseif strcmp(kkey,'d')
                f=@(x,T) 0.5*(T<((forceframes+i)*k))/k* ...
                         (((x(1)-pt(1,1)).^2+(x(2)-pt(1,2)).^2)<0.1);
            elseif strcmp(kkey,'s')
                f=@(x,T) 5*(T<((forceframes+i)*k))/k* ...
                         (((x(1)-pt(1,1)).^2+(x(2)-pt(1,2)).^2)<0.1);
            elseif strcmp(kkey,'q')
                set(gcf,'CurrentCharacter', 'a');
                break;
            elseif strcmp(kkey,'v')
                updateGraph = 0; % Ignore graphics
            end
            set(gcf,'CurrentCharacter', 'a');
        end
        
        % If no force, ignore force calculations
        if loadTimeOut > 0
            time = k*i;
            F = LoadVector(p,t,f,time);
        end
        
        % Do calculations
        b = k^2*F + M*(2*U(:,i)-U(:,i-1)) - k^2*S*U(:,i);
        %U(:,i+1) = M\b;
        U(:,i+1) = Minv*b; % Multiplying with (pseudo)inverse is MUCH faster!

        % Subtract mean to keep the surface in place
        U(:,i+1) = U(:,i+1) - sum(U(:,i+1))/size(U(:,i+1),1);
    end
end


%%% Subroutines %%%

%% Matrix assembly
function M = MassMatrix(p,t)
% Purpose:  Calculate the mass matrix
    n = size(p, 2);                       % Number of points
    M = zeros(n, n);                      % Initiate matrix
    for el = 1 : size(t, 2)
        dM = ElemMassMatrix(p, t, el);    % Calc submatrix
        nn = t(1:3, el);
        M(nn, nn) = M(nn, nn) + dM;       % Add submatrix
    end
end
    
function S = StiffMatrix(p,t)
% Purpose:  Calculate the stiffness matrix
    n = size(p, 2);                       % Number of points
    S = zeros(n, n);                      % Initiate matrix
    for el = 1 : size(t, 2)
        dS = ElemStiffMatrix(p, t, el);   % Calc submatrix
        nn = t(1:3, el);
        S(nn, nn) = S(nn, nn) + dS;       % Add submatrix
    end
end
    
function F = LoadVector(p,t,f,time)
% Purpose:  Calculate the load vector
    n = size(p, 2);                             % Number of points
    F = zeros(n, 1);                            % Initiate vector
    for el = 1 : size(t, 2)
        dF = ElemLoadVector(p, t, f, time, el); % Calc subvector
        nn = t(1:3, el);
        F(nn) = F(nn) + dF;                     % Add submatrix
    end
end

%% Elemens for matrix assembly
function dM = ElemMassMatrix(p, t, el)
% Purpose:  Compute the contribution to mass matrix M from a triangle
    % Quadrature
    %   area/3*[1/2 1/4 1/4; 1/4 1/2 1/4; 1/4 1/4 1/2]
    
    %     area
    dM = t(5,el)*(ones(3)+eye(3))/12;
end

function dS = ElemStiffMatrix(p, t, el)
% Purpose:  Compute the contribution to stiffness matrix M from a triangle
    NodeCoords = p(:, t(1:3, el));
    Dphi = ElementPhiGradients(NodeCoords);
    dx = t(5,el);               % element area
    dS = (Dphi' * Dphi) * dx;
end
    
function dF = ElemLoadVector(p, t, f, time, el)
% Purpose:  Compute the contribution to the load vector F from a triangle
%           using quadrature
    NodeCoords = p(1:2,t(1:3,el));
    MP = .5*(NodeCoords + NodeCoords(:,[2 3 1])); % Midpoint coords
    weight = t(5,el)/6; % area/6
    
    f_abc = [f(MP(:,1), time); f(MP(:,2), time); f(MP(:,1), time)];
    
    % Some products expressed in matrix form
    dF = weight*[1 0 1;...
                 1 1 0;...
                 0 1 1] * f_abc;
end

%% Helper functions

function Dphi = ElementPhiGradients(Nodes)
% Purpose:  Returns the gradients of the three element basis functions phi
%           on a triangle with nodes Nodes.
  v = Nodes(:,3)-Nodes(:,1);
  w = Nodes(:,2)-Nodes(:,1);
  gr = v - dot(v,w)*w/norm(w)^2;
  Dphi3 = gr/norm(gr)^2;
  gr = w - dot(w,v)*v/norm(v)^2;
  Dphi2 = gr/norm(gr)^2;
  Dphi1 = - (Dphi2 + Dphi3);
  Dphi = [Dphi1 Dphi2 Dphi3];
end

function ar=area(NodeCoords)
% Purpose:  Returns the area of a triangular element with given
%           node coordinates
  x1 = NodeCoords(1,:);
  x2 = NodeCoords(2,:);
  ar = (x1(2)-x1(1))*(x2(3)-x2(1))-(x1(3)-x1(1))*(x2(2)-x2(1));
  ar = ar/2;
end

function f = loadFun(x, time)
% Purpose:  Defines the load function.
    x1 = x(1,:);
    x2 = x(2,:);
    
    rad = 0.1;
    
    f = zeros(length(x1),1);
    
    %if (sqrt(x1.^2 + x2.^2) < rad)
    %    f = f + 10*cos(10*time);
    %else
    %    f = f + 0;
    %end
end

function [u_0,v_0] = bounds(x1,x2)
% Purpose: Defines the boundary (initial) condition.
    u_0 = zeros(length(x1),1);
    %u_0 = x1<-3.1;%(sqrt(sum(X.^2))<0.5)' * 1;
    %u_0(floor(length(x1)/2)) = 1;
    %v_0 = x1<-3.9;%(sqrt(sum(X.^2))<0.5)' * 1;
    v_0 = zeros(length(x1),1);
end

function pt = ptinxy()
% Purpose:  Return the point on xy-plane the mouse is on
    pt = get(gca, 'CurrentPoint');

    p1 = pt(1,:);
    p2 = pt(2,:);
    dv = p2 - p1;
    pt = p1 - p1(3)/dv(3) * dv;
end
