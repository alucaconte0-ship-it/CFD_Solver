function u = solve_unsteady_diff(uinitial, Lx, Ly, Mx, My, nu, maxIter, dt, xperiodic, yperiodic, tBC, bBC, lBC, rBC, showplot, varargin)
%Check inputs are correct format
inputs = [Lx, Ly, Mx, My, nu, maxIter, dt, xperiodic, yperiodic, tBC, bBC, lBC, rBC, showplot];
classes = {'numeric'};
attributes = {'size',[1,14]};
validateattributes(inputs, classes, attributes)
validateattributes(Lx, classes, {'scalar','positive'})
validateattributes(Ly, classes, {'scalar','positive'})
validateattributes(dt, classes, {'scalar','positive'})
validateattributes(Mx, classes, {'scalar','integer','>=',3})
validateattributes(My, classes, {'scalar','integer','>=',3})
validateattributes(nu, classes, {'scalar','positive'})
validateattributes(maxIter, classes, {'scalar','integer','nonnegative'})
validateattributes(xperiodic, classes, {'scalar','binary'}) %must be 1 or 0
validateattributes(yperiodic, classes, {'scalar','binary'})
validateattributes(showplot, classes, {'scalar','binary'})
%define dx and dy
dx = Lx/(Mx - 1);
dy = Ly/(My - 1);
%Error check uinitial is of size Mx by My
if size(uinitial,1) ~= Mx || size(uinitial,2) ~= My
    error("solve_unsteady_diff:DomainSizeError", "uinitial should be an array of size Mx by My")
end
%Define u to be uinitial
u = uinitial;
% Print any optional arguments
nVarargs = length(varargin);
for k = 1:nVarargs
    fprintf(' %d\n', varargin{k})
end
%Loop over all timesteps
for it = 1:maxIter
    %Enforce Boundary Conditions
    if xperiodic == true  %xperiodic = 1 — roll the rectangle into a cylinder left-to-right.
        u(1,:) = u(end-1,:); 
        u(end,:) = u(2,:);
    else
        u(1,:)   = lBC;  % Left boundary
        u(end,:) = rBC;  % Right boundary
    end
    if yperiodic ==true %yperiodic = 1 — roll it into a cylinder top-to-bottom.
        u(:,1) = u(:,end-1);
        u(:,end) = u(:,2);
    else
        u(:,1)   = bBC;  % Bottom boundary
        u(:,end) = tBC;  % Top boundary
    end
    %Set u at next timestep (un) to include boundary values
    utp1 = u;
    %Loop over all points and update un from 5 point stencil of u values
    for i = 2:Mx-1
        for j = 2:My - 1
            utp1(i,j) = u(i,j) + nu*dt*( (u(i+1,j)-2*u(i,j)+u(i-1,j))/dx^2+(u(i,j+1)-2*u(i,j)+u(i,j-1))/dy^2);
        end
    end
    %Overwrite u array with new u at t+1
    u = utp1;
    %Plotting
    if (showplot)
        %Note transpose as u arrays defined with Mx x My so flip to plot
        %(matlab interprets columns as horizontal and rows as vertical)
        contourf(u'); %surf also possible
        colorbar;
        pause(0.01)
    end
end
