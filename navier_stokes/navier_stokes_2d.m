function [U,V,P] = navier_stokes_2d(Lx,Ly,Mx,My,nu,dt,maxIter,rho,lidVelocity,pressureIter)
% NAVIER_STOKES_2D  Incompressible 2-D lid-driven cavity flow.
%
% [U,V,P] = navier_stokes_2d(Lx,Ly,Mx,My,nu,dt,maxIter, ...
%                           rho,lidVelocity,pressureIter)
%
% Solves the constant-density, laminar Navier-Stokes equations in a
% rectangular, closed cavity. The top wall moves at lidVelocity (m/s),
% while the other walls are stationary. No-slip walls are assumed.
%
% METHOD
%   Finite differences on a staggered (MAC) grid:
%   - u on vertical cell faces, v on horizontal cell faces
%   - pressure at cell centres
%   - explicit momentum prediction with first-order upwind advection
%   - pressure Poisson solve with successive over-relaxation (SOR)
%   - projection to approximately enforce zero velocity divergence
%
% OUTPUTS
%   U,V,P: Mx-by-My matrices at grid nodes, for use with the app.
%   The underlying staggered velocities and cell-centre pressure are
%   interpolated to nodes for plotting. These returned node values are
%   for VISUALISATION; use face velocities for a discrete divergence test.
%
% NOTE
%   An educational laminar cavity solver, not a general-purpose CFD tool.
%   Verify timestep, pressure convergence and grid independence before
%   interpreting numerical results. Pressure accuracy depends on the
%   projection and wall discretisation.

    validateattributes(Lx,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(Ly,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(Mx,{'numeric'},{'scalar','real','finite','integer','>=',3});
    validateattributes(My,{'numeric'},{'scalar','real','finite','integer','>=',3});
    validateattributes(nu,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(dt,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(maxIter,{'numeric'},{'scalar','real','finite','integer','nonnegative'});
    validateattributes(rho,{'numeric'},{'scalar','real','finite','positive'});
    validateattributes(lidVelocity,{'numeric'},{'scalar','real','finite'});
    validateattributes(pressureIter,{'numeric'},{'scalar','real','finite','integer','positive'});

    Nx = Mx-1;
    Ny = My-1;
    dx = Lx/Nx;
    dy = Ly/Ny;
    idx2 = 1/dx^2;
    idy2 = 1/dy^2;

    diffNumber = nu*dt*(idx2+idy2);
    if diffNumber > 0.5
        error('navier_stokes_2d:UnstableTimeStep', ...
              'Diffusion stability limit exceeded. Reduce dt.');
    end

    Uface = zeros(Nx+1,Ny);
    Vface = zeros(Nx,Ny+1);
    p = zeros(Nx,Ny);

    for it = 1:maxIter
        Us = Uface;
        Vs = Vface;

        for i = 2:Nx
            for j = 1:Ny
                uc = Uface(i,j);
                if j == 1
                    ub = -uc;
                else
                    ub = Uface(i,j-1);
                end
                if j == Ny
                    ut = 2*lidVelocity-uc;
                else
                    ut = Uface(i,j+1);
                end

                lap = (Uface(i+1,j)-2*uc+Uface(i-1,j))*idx2 ...
                    + (ut-2*uc+ub)*idy2;

                vc = 0.25*(Vface(i-1,j)+Vface(i,j) ...
                          +Vface(i-1,j+1)+Vface(i,j+1));

                if uc >= 0
                    dudx = (uc-Uface(i-1,j))/dx;
                else
                    dudx = (Uface(i+1,j)-uc)/dx;
                end
                if vc >= 0
                    dudy = (uc-ub)/dy;
                else
                    dudy = (ut-uc)/dy;
                end

                Us(i,j) = uc + dt*(nu*lap - uc*dudx - vc*dudy);
            end
        end

        for i = 1:Nx
            for j = 2:Ny
                vc = Vface(i,j);
                if i == 1
                    vl = -vc;
                else
                    vl = Vface(i-1,j);
                end
                if i == Nx
                    vr = -vc;
                else
                    vr = Vface(i+1,j);
                end

                lap = (vr-2*vc+vl)*idx2 ...
                    + (Vface(i,j+1)-2*vc+Vface(i,j-1))*idy2;

                uc = 0.25*(Uface(i,j-1)+Uface(i+1,j-1) ...
                          +Uface(i,j)+Uface(i+1,j));

                if uc >= 0
                    dvdx = (vc-vl)/dx;
                else
                    dvdx = (vr-vc)/dx;
                end
                if vc >= 0
                    dvdy = (vc-Vface(i,j-1))/dy;
                else
                    dvdy = (Vface(i,j+1)-vc)/dy;
                end

                Vs(i,j) = vc + dt*(nu*lap - uc*dvdx - vc*dvdy);
            end
        end

        Us(1,:) = 0;
        Us(end,:) = 0;
        Vs(:,1) = 0;
        Vs(:,end) = 0;

        cfl = dt*(max(abs(Us(:)))/dx + max(abs(Vs(:)))/dy);
        if ~isfinite(cfl) || cfl > 1
            error('navier_stokes_2d:AdvectionCFL', ...
                  'Advection CFL exceeded 1 or velocity diverged. Reduce dt.');
        end

        divergence = (Us(2:end,:)-Us(1:end-1,:))/dx ...
                   + (Vs(:,2:end)-Vs(:,1:end-1))/dy;
        rhs = (rho/dt)*divergence;

        omega = 1.6;
        for sweep = 1:pressureIter
            for i = 1:Nx
                for j = 1:Ny
                    neighbourSum = 0;
                    diagonal = 0;
                    if i > 1
                        neighbourSum = neighbourSum + p(i-1,j)*idx2;
                        diagonal = diagonal + idx2;
                    end
                    if i < Nx
                        neighbourSum = neighbourSum + p(i+1,j)*idx2;
                        diagonal = diagonal + idx2;
                    end
                    if j > 1
                        neighbourSum = neighbourSum + p(i,j-1)*idy2;
                        diagonal = diagonal + idy2;
                    end
                    if j < Ny
                        neighbourSum = neighbourSum + p(i,j+1)*idy2;
                        diagonal = diagonal + idy2;
                    end
                    pNew = (neighbourSum-rhs(i,j))/diagonal;
                    p(i,j) = (1-omega)*p(i,j)+omega*pNew;
                end
            end
        end

        p = p-p(1,1);

        Uface(2:Nx,:) = Us(2:Nx,:) ...
            - (dt/rho)*(p(2:Nx,:)-p(1:Nx-1,:))/dx;
        Vface(:,2:Ny) = Vs(:,2:Ny) ...
            - (dt/rho)*(p(:,2:Ny)-p(:,1:Ny-1))/dy;
        Uface(1,:) = 0;
        Uface(end,:) = 0;
        Vface(:,1) = 0;
        Vface(:,end) = 0;
    end

    U = zeros(Mx,My);
    V = zeros(Mx,My);
    P = zeros(Mx,My);

    U(2:Mx-1,2:My-1) = 0.5*(Uface(2:Nx,1:Ny-1) ...
                                   +Uface(2:Nx,2:Ny));
    V(2:Mx-1,2:My-1) = 0.5*(Vface(1:Nx-1,2:Ny) ...
                                   +Vface(2:Nx,2:Ny));

    U(2:end-1,end) = lidVelocity;

    P(2:Mx-1,2:My-1) = 0.25*(p(1:Nx-1,1:Ny-1) ...
                                +p(2:Nx,1:Ny-1) ...
                                +p(1:Nx-1,2:Ny) ...
                                +p(2:Nx,2:Ny));
    P(1,:) = P(2,:);
    P(end,:) = P(end-1,:);
    P(:,1) = P(:,2);
    P(:,end) = P(:,end-1);
end
