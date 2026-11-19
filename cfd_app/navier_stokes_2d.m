function [U,V,P] = navier_stokes_2d( ...
    Lx,Ly,Mx,My,nu,dt,maxIter,rho,lidVelocity,pressureIter, ...
    xperiodic,yperiodic)
% NAVIER_STOKES_2D Incompressible flow on a staggered MAC grid.
%
% Optional periodic flags:
%   xperiodic = 1: periodic left/right boundaries
%   yperiodic = 1: periodic bottom/top boundaries
%
% Nonperiodic boundaries are stationary walls, except the top wall,
% which moves at lidVelocity when yperiodic = 0.
%
% U,V,P are Mx-by-My arrays interpolated to grid nodes.

    if nargin < 11
        xperiodic = 0;
    end
    if nargin < 12
        yperiodic = 0;
    end

    validateattributes(Lx,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(Ly,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(Mx,{'numeric'},{'scalar','integer','>=',3});
    validateattributes(My,{'numeric'},{'scalar','integer','>=',3});
    validateattributes(nu,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(dt,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(maxIter,{'numeric'}, ...
        {'scalar','integer','nonnegative'});
    validateattributes(rho,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(lidVelocity,{'numeric'}, ...
        {'scalar','real','finite'});
    validateattributes(pressureIter,{'numeric'}, ...
        {'scalar','integer','positive'});
    validateattributes(xperiodic,{'numeric','logical'},{'scalar','binary'});
    validateattributes(yperiodic,{'numeric','logical'},{'scalar','binary'});

    Nx = Mx-1;
    Ny = My-1;

    dx = Lx/Nx;
    dy = Ly/Ny;
    idx2 = 1/dx^2;
    idy2 = 1/dy^2;

    if nu*dt*(idx2+idy2) > 0.5
        error('navier_stokes_2d:UnstableTimeStep', ...
            'Diffusion stability limit exceeded. Reduce dt.');
    end

    Uface = zeros(Nx+1,Ny);
    Vface = zeros(Nx,Ny+1);
    p = zeros(Nx,Ny);

    % Include the periodic boundary face in momentum updates.
    uStart = 2;
    vStart = 2;

    if xperiodic
        uStart = 1;
    end
    if yperiodic
        vStart = 1;
    end

    for it = 1:maxIter
        Us = Uface;
        Vs = Vface;

        %% Horizontal momentum prediction
        for i = uStart:Nx
            for j = 1:Ny
                uc = Uface(i,j);

                % Horizontal neighbours and adjacent pressure-cell index.
                if xperiodic
                    left = mod(i-2,Nx)+1;
                    right = mod(i,Nx)+1;
                else
                    left = i-1;
                    right = i+1;
                end

                ul = Uface(left,j);
                ur = Uface(right,j);

                % Vertical neighbours or wall ghost values.
                if j > 1
                    ub = Uface(i,j-1);
                elseif yperiodic
                    ub = Uface(i,Ny);
                else
                    ub = -uc;
                end

                if j < Ny
                    ut = Uface(i,j+1);
                elseif yperiodic
                    ut = Uface(i,1);
                else
                    ut = 2*lidVelocity-uc;
                end

                lap = (ur-2*uc+ul)*idx2 ...
                    + (ut-2*uc+ub)*idy2;

                vc = 0.25*(Vface(left,j)+Vface(i,j) ...
                    +Vface(left,j+1)+Vface(i,j+1));

                if uc >= 0
                    dudx = (uc-ul)/dx;
                else
                    dudx = (ur-uc)/dx;
                end

                if vc >= 0
                    dudy = (uc-ub)/dy;
                else
                    dudy = (ut-uc)/dy;
                end

                Us(i,j) = uc+dt*(nu*lap-uc*dudx-vc*dudy);
            end
        end

        %% Vertical momentum prediction
        for i = 1:Nx
            for j = vStart:Ny
                vc = Vface(i,j);

                if i > 1
                    vl = Vface(i-1,j);
                elseif xperiodic
                    vl = Vface(Nx,j);
                else
                    vl = -vc;
                end

                if i < Nx
                    vr = Vface(i+1,j);
                elseif xperiodic
                    vr = Vface(1,j);
                else
                    vr = -vc;
                end

                if yperiodic
                    bottom = mod(j-2,Ny)+1;
                    top = mod(j,Ny)+1;
                else
                    bottom = j-1;
                    top = j+1;
                end

                vb = Vface(i,bottom);
                vt = Vface(i,top);

                lap = (vr-2*vc+vl)*idx2 ...
                    + (vt-2*vc+vb)*idy2;

                uc = 0.25*(Uface(i,bottom)+Uface(i+1,bottom) ...
                    +Uface(i,j)+Uface(i+1,j));

                if uc >= 0
                    dvdx = (vc-vl)/dx;
                else
                    dvdx = (vr-vc)/dx;
                end

                if vc >= 0
                    dvdy = (vc-vb)/dy;
                else
                    dvdy = (vt-vc)/dy;
                end

                Vs(i,j) = vc+dt*(nu*lap-uc*dvdx-vc*dvdy);
            end
        end

        % Duplicate periodic faces or enforce impermeable walls.
        if xperiodic
            Us(end,:) = Us(1,:);
        else
            Us(1,:) = 0;
            Us(end,:) = 0;
        end

        if yperiodic
            Vs(:,end) = Vs(:,1);
        else
            Vs(:,1) = 0;
            Vs(:,end) = 0;
        end

        cfl = dt*(max(abs(Us(:)))/dx+max(abs(Vs(:)))/dy);
        if ~isfinite(cfl) || cfl > 1
            error('navier_stokes_2d:AdvectionCFL', ...
                'Advection CFL exceeded 1. Reduce dt.');
        end

        %% Pressure projection
        divergence = (Us(2:end,:)-Us(1:end-1,:))/dx ...
                   + (Vs(:,2:end)-Vs(:,1:end-1))/dy;

        rhs = (rho/dt)*divergence;

        % Remove round-off incompatibility for the singular Poisson system.
        rhs = rhs-mean(rhs(:));

        omega = 1.6;

        for sweep = 1:pressureIter
            for i = 1:Nx
                for j = 1:Ny
                    neighbourSum = 0;
                    diagonal = 0;

                    if i > 1 || xperiodic
                        left = mod(i-2,Nx)+1;
                        neighbourSum = neighbourSum+p(left,j)*idx2;
                        diagonal = diagonal+idx2;
                    end

                    if i < Nx || xperiodic
                        right = mod(i,Nx)+1;
                        neighbourSum = neighbourSum+p(right,j)*idx2;
                        diagonal = diagonal+idx2;
                    end

                    if j > 1 || yperiodic
                        bottom = mod(j-2,Ny)+1;
                        neighbourSum = neighbourSum+p(i,bottom)*idy2;
                        diagonal = diagonal+idy2;
                    end

                    if j < Ny || yperiodic
                        top = mod(j,Ny)+1;
                        neighbourSum = neighbourSum+p(i,top)*idy2;
                        diagonal = diagonal+idy2;
                    end

                    pNew = (neighbourSum-rhs(i,j))/diagonal;
                    p(i,j) = (1-omega)*p(i,j)+omega*pNew;
                end
            end

            p = p-p(1,1);
        end

        Uface(2:Nx,:) = Us(2:Nx,:) ...
            -(dt/rho)*(p(2:Nx,:)-p(1:Nx-1,:))/dx;

        Vface(:,2:Ny) = Vs(:,2:Ny) ...
            -(dt/rho)*(p(:,2:Ny)-p(:,1:Ny-1))/dy;

        if xperiodic
            Uface(1,:) = Us(1,:)-(dt/rho)*(p(1,:)-p(end,:))/dx;
            Uface(end,:) = Uface(1,:);
        else
            Uface(1,:) = 0;
            Uface(end,:) = 0;
        end

        if yperiodic
            Vface(:,1) = Vs(:,1)-(dt/rho)*(p(:,1)-p(:,end))/dy;
            Vface(:,end) = Vface(:,1);
        else
            Vface(:,1) = 0;
            Vface(:,end) = 0;
        end
    end

    %% Interpolate velocities to nodes
    U = zeros(Mx,My);
    V = zeros(Mx,My);

    U(:,2:My-1) = 0.5*(Uface(:,1:Ny-1)+Uface(:,2:Ny));
    V(2:Mx-1,:) = 0.5*(Vface(1:Nx-1,:)+Vface(2:Nx,:));

    if yperiodic
        U(:,1) = 0.5*(Uface(:,end)+Uface(:,1));
        U(:,end) = U(:,1);
    else
        U(:,1) = 0;
        U(:,end) = lidVelocity;
    end

    if xperiodic
        V(1,:) = 0.5*(Vface(end,:)+Vface(1,:));
        V(end,:) = V(1,:);
    else
        V(1,:) = 0;
        V(end,:) = 0;

        % Stationary side walls take precedence at lid corners.
        U(1,:) = 0;
        U(end,:) = 0;
    end

    %% Interpolate pressure to nodes
    % Periodic boundaries wrap cell values. Wall boundaries use
    % constant extrapolation, consistent with zero normal gradient.
    leftCells = [1, 1:Nx];
    rightCells = [1:Nx, Nx];
    bottomCells = [1, 1:Ny];
    topCells = [1:Ny, Ny];

    if xperiodic
        leftCells(1) = Nx;
        rightCells(end) = 1;
    end

    if yperiodic
        bottomCells(1) = Ny;
        topCells(end) = 1;
    end

    P = 0.25*( ...
        p(leftCells,bottomCells)+p(rightCells,bottomCells) ...
        +p(leftCells,topCells)+p(rightCells,topCells));
end
