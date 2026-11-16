function cfd_app_V4
% GUI for 2D diffusion and 2D lid-driven cavity flow.
% Requires solve_unsteady_diff.m and navier_stokes_2d.m.

    % Main window
    fig = uifigure('Name','2D Unsteady Diffusion and Navier Stokes', ...
        'Position',[100 100 1100 780]);

    main = uigridlayout(fig,[1 2]);
    main.ColumnWidth = {350,'1x'};
    main.Padding = [12 12 12 12];

    % Input panel
    panel = uipanel(main,'Title','Simulation parameters');
    panel.Layout.Column = 1;

    controls = uigridlayout(panel,[25 2]);
    controls.ColumnWidth = {175,'1x'};
    controls.RowHeight = [repmat({28},1,22),{38},{28},{'1x'}];
    controls.RowSpacing = 3;
    controls.Padding = [8 8 8 8];

    % Solver selection
    solverLabel = uilabel(controls,'Text','Solver');
    solverLabel.Layout.Row = 1;
    solverLabel.Layout.Column = 1;

    solverDrop = uidropdown(controls, ...
        'Items',{'2D Diffusion','2D Navier-Stokes'}, ...
        'Value','2D Diffusion', ...
        'ValueChangedFcn',@solverChanged);

    solverDrop.Layout.Row = 1;
    solverDrop.Layout.Column = 2;

    % Shared inputs
    LxField   = addNumber('Domain length Lx',1,2);
    LyField   = addNumber('Domain length Ly',1,3);
    MxField   = addNumber('Grid points Mx',31,4);
    MyField   = addNumber('Grid points My',31,5);
    nuField   = addNumber('Diffusivity / viscosity',0.1,6);
    dtField   = addNumber('Time step dt',0.001,7);
    iterField = addNumber('Time steps',100,8);

    % Diffusion inputs
    [topField,topLabel] = ...
        addNumber('Top boundary',0,9);

    [bottomField,bottomLabel] = ...
        addNumber('Bottom boundary',0,10);

    [leftField,leftLabel] = ...
        addNumber('Left boundary',0,11);

    [rightField,rightLabel] = ...
        addNumber('Right boundary',0,12);

    xCheck = uicheckbox(controls,'Text','Periodic in x');
    xCheck.Layout.Row = 13;
    xCheck.Layout.Column = [1 2];

    yCheck = uicheckbox(controls,'Text','Periodic in y');
    yCheck.Layout.Row = 14;
    yCheck.Layout.Column = [1 2];

    modeLabel = uilabel(controls,'Text','Initial condition');
    modeLabel.Layout.Row = 15;
    modeLabel.Layout.Column = 1;

    modeDrop = uidropdown(controls, ...
        'Items',{'Single hot cell','Uniform'}, ...
        'ValueChangedFcn',@initialConditionChanged);

    modeDrop.Layout.Row = 15;
    modeDrop.Layout.Column = 2;

    [amplitudeField,amplitudeLabel] = ...
        addNumber('Initial amplitude',100,16);

    % Hotspot grid-point inputs
    % Grid point (1,1) is bottom-left
    [hotspotXField,hotspotXLabel] = ...
        addNumber('Hotspot grid X',16,17);

    [hotspotYField,hotspotYLabel] = ...
        addNumber('Hotspot grid Y',16,18);

    hotspotXField.ValueChangedFcn = @hotspotGridChanged;
    hotspotYField.ValueChangedFcn = @hotspotGridChanged;

    % Navier-Stokes inputs
    [rhoField,rhoLabel] = addNumber('Fluid density rho',1,19);
    [lidField,lidLabel] = addNumber('Lid velocity',1,20);
    [pressureField,pressureLabel] = addNumber('Pressure iterations',100,21);

    % Buttons
    runButton = uibutton(controls,'push', ...
        'Text','Run simulation', ...
        'FontWeight','bold', ...
        'ButtonPushedFcn',@runSimulation);

    runButton.Layout.Row = 22;
    runButton.Layout.Column = 1;

    resetButton = uibutton(controls,'push', ...
        'Text','Reset', ...
        'ButtonPushedFcn',@resetSimulation);

    resetButton.Layout.Row = 22;
    resetButton.Layout.Column = 2;

    % Information
    info = uilabel(controls, ...
        'WordWrap','on');

    info.Layout.Row = 23;
    info.Layout.Column = [1 2];

    % Status
    status = uilabel(controls, ...
        'Text','Ready.', ...
        'WordWrap','on');

    status.Layout.Row = 24;
    status.Layout.Column = [1 2];


    % -------------------------------------------------------------
    % Output panel
    % -------------------------------------------------------------
    right = uipanel(main,'Title','Simulation output');
    right.Layout.Column = 2;

    ax = uiaxes(right, ...
        'Position',[45 55 600 600]);

    xlabel(ax,'x');
    ylabel(ax,'y');


    % -------------------------------------------------------------
    % Main output selector
    % -------------------------------------------------------------
    viewLabel = uilabel(right, ...
        'Text','View:', ...
        'Position',[45 675 35 22]);

    viewDrop = uidropdown(right, ...
        'Items',{'Plot','Matrix'}, ...
        'Value','Plot', ...
        'Position',[85 675 120 22], ...
        'ValueChangedFcn',@viewChanged);


    % -------------------------------------------------------------
    % Navier-Stokes matrix quantity selector
    % -------------------------------------------------------------
    quantityLabel = uilabel(right, ...
        'Text','Matrix quantity:', ...
        'Position',[225 675 95 22], ...
        'Visible','off');

    quantityDrop = uidropdown(right, ...
        'Items',{ ...
        'Velocity magnitude', ...
        'U velocity', ...
        'V velocity', ...
        'Pressure'}, ...
        'Value','Velocity magnitude', ...
        'Position',[325 675 150 22], ...
        'Visible','off', ...
        'ValueChangedFcn',@quantityChanged);


    % Matrix view occupies same area as graph
    matrixTable = uitable(right, ...
        'Position',[45 55 600 600], ...
        'Visible','off');


    % -------------------------------------------------------------
    % Store latest results
    % -------------------------------------------------------------
    lastMatrix = [];

    lastU = [];
    lastV = [];
    lastP = [];

    lastX = [];
    lastY = [];

    currentColorbar = gobjects(0);


    % -------------------------------------------------------------
    % Initial hotspot
    % -------------------------------------------------------------
    x0 = LxField.Value/2;
    y0 = LyField.Value/2;

    marker = gobjects(0);

    ax.ButtonDownFcn = @selectHotspot;

    LxField.ValueChangedFcn = @domainChanged;
    LyField.ValueChangedFcn = @domainChanged;

    MxField.ValueChangedFcn = @gridChanged;
    MyField.ValueChangedFcn = @gridChanged;

    solverChanged([],[]);


    % =============================================================
    % Create numeric input
    % =============================================================
    function [field,lbl] = addNumber(caption,value,row)

        lbl = uilabel(controls,'Text',caption);
        lbl.Layout.Row = row;
        lbl.Layout.Column = 1;

        field = uieditfield(controls, ...
            'numeric', ...
            'Value',value);

        field.Layout.Row = row;
        field.Layout.Column = 2;

    end


    % =============================================================
    % Switch between solvers
    % =============================================================
    function solverChanged(~,~)

        isDiffusion = ...
            strcmp(solverDrop.Value,'2D Diffusion');

        diffControls = { ...
            topField,topLabel, ...
            bottomField,bottomLabel, ...
            leftField,leftLabel, ...
            rightField,rightLabel, ...
            xCheck,yCheck, ...
            modeLabel,modeDrop, ...
            amplitudeField,amplitudeLabel, ...
            hotspotXField,hotspotXLabel, ...
            hotspotYField,hotspotYLabel};

        nsControls = { ...
            rhoField,rhoLabel, ...
            lidField,lidLabel, ...
            pressureField,pressureLabel};

        for k = 1:numel(diffControls)

            diffControls{k}.Visible = ...
                onoff(isDiffusion);

        end

        for k = 1:numel(nsControls)

            nsControls{k}.Visible = ...
                onoff(~isDiffusion);

        end

        updateHotspotControlVisibility();


        % Show matrix quantity selector only for Navier-Stokes
        quantityLabel.Visible = ...
            onoff(~isDiffusion);

        quantityDrop.Visible = ...
            onoff(~isDiffusion);


        if isDiffusion

            info.Text = ...
                ['Diffusion: click the graph or enter a hotspot grid point. ' ...
                'Grid point (1,1) is the bottom-left.'];

        else

            info.Text = ...
                ['Navier-Stokes: lid-driven cavity flow. ' ...
                'Matrix view can display speed, U, V or pressure.'];

        end

        status.Text = 'Ready.';

        resetPlot();

    end


    % =============================================================
    % Reset all inputs
    % =============================================================
    function resetSimulation(~,~)

        solverDrop.Value = '2D Diffusion';

        LxField.Value = 1;
        LyField.Value = 1;

        MxField.Value = 31;
        MyField.Value = 31;

        nuField.Value = 0.1;
        dtField.Value = 0.001;
        iterField.Value = 100;

        topField.Value = 0;
        bottomField.Value = 0;
        leftField.Value = 0;
        rightField.Value = 0;

        xCheck.Value = false;
        yCheck.Value = false;

        amplitudeField.Value = 100;

        hotspotXField.Value = 16;
        hotspotYField.Value = 16;

        rhoField.Value = 1;
        lidField.Value = 1;
        pressureField.Value = 100;

        quantityDrop.Value = ...
            'Velocity magnitude';

        x0 = 0.5;
        y0 = 0.5;

        solverChanged([],[]);

        status.Text = ...
            'Reset complete.';

    end


    % =============================================================
    % Reset hotspot when domain changes
    % =============================================================
    function domainChanged(~,~)

        if LxField.Value > 0 && ...
                LyField.Value > 0 && ...
                MxField.Value >= 3 && ...
                MyField.Value >= 3

            hotspotXField.Value = ...
                round((MxField.Value + 1)/2);

            hotspotYField.Value = ...
                round((MyField.Value + 1)/2);

            updateHotspotFromGrid();

            resetPlot();

            status.Text = ...
                'Hotspot reset to centre.';

        end

    end


    % =============================================================
    % Reset hotspot when grid changes
    % =============================================================
    function gridChanged(~,~)

        if isfinite(MxField.Value) && ...
                isfinite(MyField.Value) && ...
                MxField.Value >= 3 && ...
                MyField.Value >= 3

            hotspotXField.Value = ...
                round((MxField.Value + 1)/2);

            hotspotYField.Value = ...
                round((MyField.Value + 1)/2);

            updateHotspotFromGrid();

            resetPlot();

            status.Text = ...
                'Hotspot reset to centre grid point.';

        end

    end


    % =============================================================
    % Initial condition changed
    % =============================================================
    function initialConditionChanged(~,~)

        updateHotspotControlVisibility();

        resetPlot();

    end


    % =============================================================
    % Show hotspot coordinates when required
    % =============================================================
    function updateHotspotControlVisibility()

        showHotspot = ...
            strcmp(solverDrop.Value,'2D Diffusion') && ...
            strcmp(modeDrop.Value,'Single hot cell');

        hotspotXField.Visible = ...
            onoff(showHotspot);

        hotspotXLabel.Visible = ...
            onoff(showHotspot);

        hotspotYField.Visible = ...
            onoff(showHotspot);

        hotspotYLabel.Visible = ...
            onoff(showHotspot);

    end


    % =============================================================
    % Hotspot grid coordinate changed
    % =============================================================
    function hotspotGridChanged(~,~)

        if ~strcmp( ...
                solverDrop.Value,'2D Diffusion') || ...
                strcmp(modeDrop.Value,'Uniform')

            return

        end

        updateHotspotFromGrid();

        drawHotspotMarker();

        status.Text = sprintf( ...
            'Hotspot grid point: (%d, %d)', ...
            round(hotspotXField.Value), ...
            round(hotspotYField.Value));

    end


    % =============================================================
    % Convert grid point to physical coordinate
    % =============================================================
    function updateHotspotFromGrid()

        if MxField.Value < 3 || ...
                MyField.Value < 3 || ...
                LxField.Value <= 0 || ...
                LyField.Value <= 0

            return

        end

        i = round(hotspotXField.Value);
        j = round(hotspotYField.Value);

        i = min( ...
            MxField.Value, ...
            max(1,i));

        j = min( ...
            MyField.Value, ...
            max(1,j));

        hotspotXField.Value = i;
        hotspotYField.Value = j;

        dxLocal = ...
            LxField.Value/(MxField.Value-1);

        dyLocal = ...
            LyField.Value/(MyField.Value-1);

        % Grid point (1,1) corresponds to (0,0)
        x0 = (i-1)*dxLocal;
        y0 = (j-1)*dyLocal;

    end


    % =============================================================
    % Change Plot / Matrix view
    % =============================================================
    function viewChanged(~,~)

        updateOutputView();

    end


    % =============================================================
    % Change Navier-Stokes matrix quantity
    % =============================================================
    function quantityChanged(~,~)

        if strcmp(viewDrop.Value,'Matrix')

            updateOutputView();

        end

    end


    % =============================================================
    % Update output view
    % =============================================================
    function updateOutputView()

        if strcmp(viewDrop.Value,'Matrix')

            ax.Visible = 'off';

            matrixTable.Visible = 'on';

            % Hide colourbar in Matrix view
            if ~isempty(currentColorbar) && ...
                    isgraphics(currentColorbar)

                currentColorbar.Visible = 'off';

            end


            if ~isempty(lastMatrix)

                % ---------------------------------------------
                % Choose which matrix to display
                % ---------------------------------------------
                if strcmp( ...
                        solverDrop.Value, ...
                        '2D Navier-Stokes')

                    switch quantityDrop.Value

                        case 'Velocity magnitude'

                            matrixToShow = ...
                                lastMatrix;


                        case 'U velocity'

                            matrixToShow = ...
                                lastU;


                        case 'V velocity'

                            matrixToShow = ...
                                lastV;


                        case 'Pressure'

                            matrixToShow = ...
                                lastP;

                    end

                else

                    % Diffusion
                    matrixToShow = ...
                        lastMatrix;

                end


                % Flip vertically:
                % y=max at top
                % y=0 at bottom
                matrixTable.Data = ...
                    flipud(matrixToShow.');


                % x coordinates left -> right
                matrixTable.ColumnName = ...
                    arrayfun( ...
                    @(value) sprintf('%.3g',value), ...
                    lastX, ...
                    'UniformOutput',false);


                % y coordinates top -> bottom
                matrixTable.RowName = ...
                    arrayfun( ...
                    @(value) sprintf('%.3g',value), ...
                    flip(lastY), ...
                    'UniformOutput',false);

            else

                matrixTable.Data = [];

                matrixTable.ColumnName = ...
                    'numbered';

                matrixTable.RowName = ...
                    'numbered';

            end


        else

            % Plot view
            matrixTable.Visible = 'off';

            ax.Visible = 'on';

            if ~isempty(currentColorbar) && ...
                    isgraphics(currentColorbar)

                currentColorbar.Visible = 'on';

            end

        end

    end


    % =============================================================
    % Clear graph
    % =============================================================
    function resetPlot()

        cla(ax);

        colorbar(ax,'off');

        currentColorbar = ...
            gobjects(0);

        lastMatrix = [];

        lastU = [];
        lastV = [];
        lastP = [];

        lastX = [];
        lastY = [];

        matrixTable.Data = [];

        matrixTable.ColumnName = ...
            'numbered';

        matrixTable.RowName = ...
            'numbered';

        title(ax,'Simulation Output');

        xlabel(ax,'x');
        ylabel(ax,'y');

        if LxField.Value > 0 && ...
                LyField.Value > 0

            xlim(ax,[0 LxField.Value]);
            ylim(ax,[0 LyField.Value]);

        end

        drawHotspotMarker();

        updateOutputView();

    end


    % =============================================================
    % Select hotspot with mouse
    % =============================================================
    function selectHotspot(~,~)

        if ~strcmp( ...
                solverDrop.Value,'2D Diffusion') || ...
                strcmp(modeDrop.Value,'Uniform')

            return

        end

        point = ax.CurrentPoint;

        xc = point(1,1);
        yc = point(1,2);

        if ~isfinite(xc) || ...
                ~isfinite(yc) || ...
                xc < 0 || ...
                xc > LxField.Value || ...
                yc < 0 || ...
                yc > LyField.Value

            return

        end

        dxLocal = ...
            LxField.Value/(MxField.Value-1);

        dyLocal = ...
            LyField.Value/(MyField.Value-1);

        % Convert mouse click to nearest grid point
        i = min( ...
            MxField.Value, ...
            max(1,round(xc/dxLocal)+1));

        j = min( ...
            MyField.Value, ...
            max(1,round(yc/dyLocal)+1));

        hotspotXField.Value = i;
        hotspotYField.Value = j;

        x0 = (i-1)*dxLocal;
        y0 = (j-1)*dyLocal;

        drawHotspotMarker();

        status.Text = sprintf( ...
            'Hotspot grid point: (%d, %d)', ...
            i,j);

    end


    % =============================================================
    % Draw hotspot marker
    % =============================================================
    function drawHotspotMarker()

        if isgraphics(marker)

            delete(marker);

        end

        marker = ...
            gobjects(0);

        if ~strcmp( ...
                solverDrop.Value,'2D Diffusion') || ...
                strcmp(modeDrop.Value,'Uniform')

            return

        end

        hold(ax,'on');

        marker = plot( ...
            ax, ...
            x0,y0, ...
            'rx', ...
            'MarkerSize',12, ...
            'LineWidth',2);

        marker.HitTest = 'off';

        marker.PickableParts = ...
            'none';

        hold(ax,'off');

    end


    % =============================================================
    % Run simulation
    % =============================================================
    function runSimulation(~,~)

        runButton.Enable = 'off';
        resetButton.Enable = 'off';

        status.Text = ...
            'Calculating...';

        drawnow;


        try

            % -----------------------------------------------------
            % Shared inputs
            % -----------------------------------------------------
            Lx = LxField.Value;
            Ly = LyField.Value;

            Mx = MxField.Value;
            My = MyField.Value;

            nu = nuField.Value;
            dt = dtField.Value;

            maxIter = ...
                iterField.Value;


            % -----------------------------------------------------
            % Validate inputs
            % -----------------------------------------------------
            validateattributes( ...
                Lx, ...
                {'numeric'}, ...
                {'scalar','real','finite','positive'});

            validateattributes( ...
                Ly, ...
                {'numeric'}, ...
                {'scalar','real','finite','positive'});

            validateattributes( ...
                Mx, ...
                {'numeric'}, ...
                {'scalar','integer','>=',3});

            validateattributes( ...
                My, ...
                {'numeric'}, ...
                {'scalar','integer','>=',3});

            validateattributes( ...
                nu, ...
                {'numeric'}, ...
                {'scalar','real','finite','positive'});

            validateattributes( ...
                dt, ...
                {'numeric'}, ...
                {'scalar','real','finite','positive'});

            validateattributes( ...
                maxIter, ...
                {'numeric'}, ...
                {'scalar','integer','nonnegative'});


            x = linspace(0,Lx,Mx);
            y = linspace(0,Ly,My);

            dx = Lx/(Mx-1);
            dy = Ly/(My-1);


            % =====================================================
            % DIFFUSION
            % =====================================================
            if strcmp( ...
                    solverDrop.Value, ...
                    '2D Diffusion')

                % Stability check
                if nu*dt*( ...
                        1/dx^2 + ...
                        1/dy^2) > 0.5

                    error( ...
                        'Diffusion time step too large: reduce dt.');

                end


                amplitude = ...
                    amplitudeField.Value;


                % -------------------------------------------------
                % Initial condition
                % -------------------------------------------------
                switch modeDrop.Value

                    case 'Single hot cell'

                        uinitial = ...
                            zeros(Mx,My);

                        validateattributes( ...
                            hotspotXField.Value, ...
                            {'numeric'}, ...
                            {'scalar','integer','>=',1,'<=',Mx});

                        validateattributes( ...
                            hotspotYField.Value, ...
                            {'numeric'}, ...
                            {'scalar','integer','>=',1,'<=',My});

                        i = ...
                            hotspotXField.Value;

                        j = ...
                            hotspotYField.Value;

                        x0 = ...
                            (i-1)*dx;

                        y0 = ...
                            (j-1)*dy;

                        uinitial(i,j) = ...
                            amplitude;


                    case 'Uniform'

                        uinitial = ...
                            amplitude*ones(Mx,My);

                end


                % -------------------------------------------------
                % Run diffusion solver
                % -------------------------------------------------
                u = solve_unsteady_diff( ...
                    uinitial, ...
                    Lx,Ly, ...
                    Mx,My, ...
                    nu, ...
                    maxIter, ...
                    dt, ...
                    double(xCheck.Value), ...
                    double(yCheck.Value), ...
                    topField.Value, ...
                    bottomField.Value, ...
                    leftField.Value, ...
                    rightField.Value, ...
                    0);


                % -------------------------------------------------
                % Plot diffusion
                % -------------------------------------------------
                cla(ax);

                [~,h] = contourf( ...
                    ax, ...
                    x,y, ...
                    u.', ...
                    20, ...
                    'LineColor','none');

                h.HitTest = 'off';

                h.PickableParts = ...
                    'none';

                colormap(ax,'parula');

                title(ax, ...
                    sprintf( ...
                    'Diffusion at t = %.4g', ...
                    maxIter*dt));


                cb = colorbar(ax);

                cb.Label.String = ...
                    'Value';

                currentColorbar = cb;


                % Keep graph original size
                ax.Position = ...
                    [45 55 600 600];


                % Put colourbar beside graph
                cb.Position = ...
                    [0.925 0.13 0.025 0.74];


                % -------------------------------------------------
                % Store diffusion result
                % -------------------------------------------------
                lastMatrix = u;

                lastU = [];
                lastV = [];
                lastP = [];

                lastX = x;
                lastY = y;


                drawHotspotMarker();

                updateOutputView();

                status.Text = ...
                    'Diffusion simulation complete.';


            % =====================================================
            % NAVIER-STOKES
            % =====================================================
            else

                % -------------------------------------------------
                % Navier-Stokes inputs
                % -------------------------------------------------
                rho = ...
                    rhoField.Value;

                lidVelocity = ...
                    lidField.Value;

                pressureIter = ...
                    pressureField.Value;


                validateattributes( ...
                    rho, ...
                    {'numeric'}, ...
                    {'scalar','real','finite','positive'});

                validateattributes( ...
                    lidVelocity, ...
                    {'numeric'}, ...
                    {'scalar','real','finite'});

                validateattributes( ...
                    pressureIter, ...
                    {'numeric'}, ...
                    {'scalar','integer','positive'});


                % -------------------------------------------------
                % Run Navier-Stokes solver
                %
                % Keep U, V AND P
                % -------------------------------------------------
                [U,V,P] = navier_stokes_2d( ...
                    Lx,Ly, ...
                    Mx,My, ...
                    nu,dt, ...
                    maxIter, ...
                    rho, ...
                    lidVelocity, ...
                    pressureIter);


                if ~isequal(size(U),[Mx My]) || ...
                        ~isequal(size(V),[Mx My]) || ...
                        ~isequal(size(P),[Mx My])

                    error( ...
                        ['Navier-Stokes solver must return ' ...
                        'Mx-by-My arrays.']);

                end


                % Velocity magnitude
                speed = ...
                    sqrt(U.^2 + V.^2);


                % -------------------------------------------------
                % Plot velocity magnitude
                % -------------------------------------------------
                cla(ax);

                contourf( ...
                    ax, ...
                    x,y, ...
                    speed.', ...
                    20, ...
                    'LineColor','none');

                colormap(ax,'parula');


                cb = colorbar(ax);

                cb.Label.String = ...
                    'Velocity magnitude';

                currentColorbar = cb;


                % Keep graph same size
                ax.Position = ...
                    [45 55 600 600];


                % Put colourbar beside graph
                cb.Position = ...
                    [0.925 0.13 0.025 0.74];


                % -------------------------------------------------
                % Store ALL Navier-Stokes outputs
                % -------------------------------------------------
                lastMatrix = speed;

                lastU = U;
                lastV = V;
                lastP = P;

                lastX = x;
                lastY = y;


                % -------------------------------------------------
                % Velocity arrows
                % -------------------------------------------------
                hold(ax,'on');


                stride = ...
                    max( ...
                    1, ...
                    ceil(max(Mx,My)/15));

                ix = ...
                    1:stride:Mx;

                iy = ...
                    1:stride:My;


                quiver( ...
                    ax, ...
                    x(ix), ...
                    y(iy), ...
                    U(ix,iy).', ...
                    V(ix,iy).', ...
                    'k');


                hold(ax,'off');


                title(ax, ...
                    sprintf( ...
                    'Velocity magnitude at t = %.4g', ...
                    maxIter*dt));


                updateOutputView();

                status.Text = ...
                    'Navier-Stokes simulation complete.';

            end


            % -----------------------------------------------------
            % Axis labels
            % -----------------------------------------------------
            xlabel(ax,'x');
            ylabel(ax,'y');

            xlim(ax,[0 Lx]);
            ylim(ax,[0 Ly]);


        catch ME

            status.Text = ...
                'Simulation could not run.';

            uialert( ...
                fig, ...
                ME.message, ...
                'Simulation error');

        end


        runButton.Enable = 'on';
        resetButton.Enable = 'on';

    end


    % =============================================================
    % Convert true/false to on/off
    % =============================================================
    function s = onoff(value)

        if value

            s = 'on';

        else

            s = 'off';

        end

    end

end
