function fig = unsteady_diffusion_app
    % Create app window
    fig = uifigure('Name','2D Unsteady Diffusion Solver', ...
        'Position',[100 100 1100 700], ...
        'Color',[0.95 0.96 0.98]);

    % Title
    uilabel(fig,'Text','2D Unsteady Diffusion Solver', ...
        'Position',[30 650 500 35], ...
        'FontSize',22,'FontWeight','bold', ...
        'FontColor',[0.12 0.25 0.40]);

    % Panel for inputs; components use coordinates relative to this panel
    panel = uipanel(fig,'Title','Simulation Parameters', ...
        'Position',[15 50 440 560], ...
        'BackgroundColor',[1 1 1]);

    % Inputs inside the panel
    Lx = addInput('Lx',1,510);
    Ly = addInput('Ly',1,475);
    Mx = addInput('Mx',31,440);
    My = addInput('My',31,405);
    nu = addInput('Diffusivity',0.1,370);
    dt = addInput('Time step',0.001,335);
    maxIter = addInput('Iterations',100,300);

    tBC = addInput('Top BC',0,265);
    bBC = addInput('Bottom BC',0,230);
    lBC = addInput('Left BC',0,195);
    rBC = addInput('Right BC',0,160);

    xp = uicheckbox(panel,'Text','Periodic X', ...
        'Position',[290 510 120 25]);
    yp = uicheckbox(panel,'Text','Periodic Y', ...
        'Position',[290 475 120 25]);

    uilabel(panel,'Text','Initial condition', ...
        'Position',[275 425 145 25]);
    mode = uidropdown(panel, ...
        'Items',{'Gaussian','Single hot cell','Uniform'}, ...
        'Position',[275 395 145 25]);

    initvalue = addInput('Initial Value',100,125);

    % Run button
    uibutton(panel,'Text','Run Simulation','Position',[105 55 230 45], 'FontSize',14,'FontWeight','bold', 'BackgroundColor',[0.13 0.32 0.52], 'FontColor',[1 1 1], 'ButtonPushedFcn',@runSimulation);

    % Plot area
    ax = uiaxes(fig,'Position',[485 95 590 520]);
    title(ax,'Simulation Output');
    xlabel(ax,'x');
    ylabel(ax,'y');
    ax.FontSize = 12;

    % Helper function for numeric inputs
    function field = addInput(name,value,y)
        uilabel(panel,'Text',name, 'Position',[20 y 115 25], 'FontSize',12,'FontColor',[0.25 0.30 0.38]);

        field = uieditfield(panel,'numeric', 'Value',value, 'Position',[140 y 115 28], 'FontSize',12, 'BackgroundColor',[0.97 0.98 1]);
    end

    % Run simulation when button is pressed
    function runSimulation(~,~)
        lx = Lx.Value;
        ly = Ly.Value;
        mx = Mx.Value;
        my = My.Value;
        diffusion = nu.Value;
        timestep = dt.Value;
        iterations = maxIter.Value;
        A = initvalue.Value;

        % Create grid
        x = linspace(0,lx,mx);
        y = linspace(0,ly,my);
        [X,Y] = ndgrid(x,y);

        % Initial condition
        switch mode.Value
            case 'Gaussian'
                sigma = 0.12*min(lx,ly);
                uinitial = A*exp(-((X-lx/2).^2 + ...
                    (Y-ly/2).^2)/(2*sigma^2));
            case 'Single hot cell'
                uinitial = zeros(mx,my);
                uinitial(round((mx+1)/2),round((my+1)/2)) = A;
            case 'Uniform'
                uinitial = A*ones(mx,my);
        end

        % Call existing solver; 0 disables its separate plot
        u = solve_unsteady_diff(uinitial,lx,ly,mx,my, ...
            diffusion,iterations,timestep, ...
            double(xp.Value),double(yp.Value), ...
            tBC.Value,bBC.Value,lBC.Value,rBC.Value,0);

        % Plot result
        contourf(ax,x,y,u.',30,'LineColor','none');
        colormap(ax,'turbo');
        colorbar(ax);
        title(ax,'2D Diffusion');
        xlabel(ax,'x');
        ylabel(ax,'y');
    end
end
