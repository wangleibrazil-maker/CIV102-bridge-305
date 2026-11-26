% CIV102 Appendix I, Sections 0-1: a moving train on a simply supported span.
% Run the script. SFDi(i,:) and BMDi(i,:) are the diagrams at train_shift(i).
% All distances are in mm, forces in N, and bending moments in N mm.

%% 0. Bridge, train and display settings
L = 1200;
dx = 1;                              % Bridge stations and train moves use the same 1 mm grid
x = 0:dx:L;
x_train = [52 228 392 568 732 908];  % Axle coordinates when train_shift = 0
P = 400;
P_train = (P / 6) * ones(size(x_train));   % Load Case 1, downward load magnitudes
% For the first pass of Load Case 2, an alternative axle load vector is
% P_train = [135 135 135 135 182 182] / 2; % Locomotive at the right end
% P_train = fliplr(P_train);                % Reverse the axle order for another run

show_shift = [];                      % Optional extra view: 0 for the reference, 120 for centered
animate_train = false;                % Set true to play every 1 mm step

% First frame: rightmost axle at A. Last frame: leftmost axle at B.
train_shift = -x_train(end):dx:(L - x_train(1));
n_train = numel(train_shift);
SFDi = zeros(n_train, numel(x));
BMDi = zeros(n_train, numel(x));
RA = zeros(n_train, 1);
RB = zeros(n_train, 1);

%% 1. Calculate one SFD and one BMD for each train position
for i = 1:n_train
    axle_x = x_train + train_shift(i);
    on_bridge = (axle_x >= 0) & (axle_x <= L);
    a = axle_x(on_bridge);            % Axle coordinates on this span
    F = P_train(on_bridge);           % Matching downward forces

    RB(i) = sum(F .* a) / L;          % Sum of moments about A = 0
    RA(i) = sum(F) - RB(i);           % Sum of vertical forces = 0

    % At each grid point, add the signed point forces (including the supports).
    point_force = zeros(size(x));
    point_force(1) = RA(i);
    point_force(end) = RB(i);
    for k = 1:numel(a)
        j = a(k) / dx + 1;            % Grid starts at x = 0; MATLAB indices start at 1
        point_force(j) = point_force(j) - F(k);
    end

    V = cumsum(point_force);          % Shear just to the right of each node; V(end) is outside
    Vseg = V(1:end-1);                % Constant shear in the intervals between nodes
    SFDi(i,:) = [Vseg, Vseg(end)];    % At B store the shear just inside the span
    BMDi(i,:) = [0, cumsum(Vseg * dx)]; % M(0) = 0; each interval adds V times dx
end


% At each x, compare all train positions (down the column).
SFD = max(abs(SFDi), [], 1);
% At an interior node the preceding segment gives the other side of a shear jump.
SFD(2:end) = max(SFD(2:end), SFD(1:end-1)); % Envelope of the shear magnitude on either side (N)
BMD = max(BMDi, [], 1);               % Positive sagging moment envelope (N mm)

%% 2. Global peaks across all train positions and bridge sections
Vmax = max(SFD);                     % Maximum shear magnitude on the grid
Mmax = max(BMD);                     % Maximum sagging moment on the grid


% A tolerance includes mirror-image cases that differ only by roundoff.
tolV = 1e-10 * max(1, Vmax);
tolM = 1e-10 * max(1, Mmax);

% Columns of each result: train row, then bridge column.
% Shear column j acts between x(j) and x(j+1). The last SFDi column repeats B-.
[trainRow, bridgeCol] = find(abs(abs(SFDi(:,1:end-1)) - Vmax) <= tolV);
shearPeaks = sortrows([trainRow, bridgeCol]);
[trainRow, bridgeCol] = find(abs(BMDi - Mmax) <= tolM);
momentPeaks = sortrows([trainRow, bridgeCol]);

fprintf('\n%g N train | %d train positions | %g mm step\n', sum(P_train), n_train, dx);
fprintf('Maximum |V| = %.2f N\n', Vmax);
for k = 1:size(shearPeaks,1)
    i = shearPeaks(k,1);
    j = shearPeaks(k,2);
    fprintf('  shift = %+g mm, leftmost axle = %g mm, bridge x = %g to %g mm, V = %+.2f N\n', ...
            train_shift(i), x_train(1) + train_shift(i), x(j), x(j+1), SFDi(i,j));
end
fprintf('Maximum M = %.2f N mm\n', Mmax);
for k = 1:size(momentPeaks,1)
    i = momentPeaks(k,1);
    j = momentPeaks(k,2);
    fprintf('  shift = %+g mm, leftmost axle = %g mm, bridge x = %g mm, M = %.2f N mm\n', ...
            train_shift(i), x_train(1) + train_shift(i), x(j), BMDi(i,j));
end

%% 3. Section properties (Design 0)
b = [80 1.27 1.27 5 5 100]
h = [1.27 73.73 73.73 1.27 1.27 1.27]
yc = [0.635 38.135 38.135 74.365 74.365 75.635]
ybar = sum(b.*h.*yc)/sum(b.*h)
I = sum(b.*h.^3/12 + b.*h.*(yc-ybar).^2)
Qglue = 100*1.27*(75.635 - ybar)
Qcent = 80*1.27*(ybar - 0.635) + 2*1.27*(ybar - 1.27)*(ybar - 1.27)/2


%% 4. Stresses, capacities and FOS
H = 76.27;
S_top = BMD .* (H - ybar) / I;
S_bot = BMD .* ybar / I;
T_cent = SFD .* Qcent / (I * 2 * 1.27);
T_glue = SFD .* Qglue / (I * 10);

t = 1.27
k0 = pi^2*4000/(12*(1-0.2^2))

S_buck1 = 4 * k0 * (t/78.73)^2
S_buck2 = 0.425 * k0 * (t/10)^2
S_buck3 = 6 * k0 * (t/(H - t - ybar))^2
T_buck = 5 * k0 * ((t/73.73)^2 + (t/400)^2)


%% LEI DO FOS PLEASE


%% 5. Draw both diagrams at the first shear-critical and moment-critical positions
% All tied positions are printed above. Each column of this figure is one train position.
figure('Name', 'Critical train positions', 'NumberTitle', 'off', 'Color', 'w', ...
       'Units', 'normalized', 'Position', [0.12 0.15 0.74 0.68]);
criticalRows = [shearPeaks(1,1), momentPeaks(1,1)];
caseNames = {'Peak shear', 'Peak moment'};
for col = 1:2
    i = criticalRows(col);
    axV = subplot(2,2,col);
    axM = subplot(2,2,col+2);
    curves = makeDiagramPair(axV, axM, x, Vmax, Mmax);
    showFrame(curves, SFDi(i,:), BMDi(i,:), x_train + train_shift(i), L);
    caption = sprintf('%s case  |  shift = %+g mm', caseNames{col}, train_shift(i));
    title(axV, caption);
    if col == 1
        j = shearPeaks(1,2);
        sectionMid = (x(j) + x(j+1)) / 2;
        plot(axV, sectionMid, SFDi(i,j), 'ko', 'MarkerFaceColor', 'y', ...
             'Clipping', 'off', 'HandleVisibility', 'off');
        title(axV, {caption, sprintf('V = %+.2f N on x = %g to %g mm', ...
                                     SFDi(i,j), x(j), x(j+1))});
    else
        j = momentPeaks(1,2);
        plot(axM, x(j), BMDi(i,j), 'ko', 'MarkerFaceColor', 'y');
        title(axM, sprintf('M = %.2f N mm at x = %g mm', BMDi(i,j), x(j)));
    end
end

%% 6. Optional: view a specified train position or play the crossing
if animate_train || ~isempty(show_shift)
    if animate_train
        frames = 1:n_train;
    else
        [~, frames] = min(abs(train_shift - show_shift));
    end
    fig = figure('Name', 'CIV102 moving train', 'Color', 'w');
    axV = subplot(2,1,1);
    axM = subplot(2,1,2);
    curves = makeDiagramPair(axV, axM, x, Vmax, Mmax);
    for i = frames
        if ~ishandle(fig)
            break
        end
        showFrame(curves, SFDi(i,:), BMDi(i,:), x_train + train_shift(i), L);
        title(axV, sprintf('Train shift = %+g mm    R_A = %.2f N    R_B = %.2f N', ...
                          train_shift(i), RA(i), RB(i)));
        drawnow
    end
end

function curves = makeDiagramPair(axV, axM, x, Vmax, Mmax)
    L = x(end);
    hold(axV, 'on');
    plot(axV, [0 L], [0 0], 'k:', 'HandleVisibility', 'off');
    % Extra endpoints draw both support jumps to zero outside the span.
    curves.shear = stairs(axV, [0 x L], zeros(1, numel(x)+2), 'b', 'LineWidth', 1.5);
    curves.axles = plot(axV, NaN, NaN, 'kv', 'MarkerFaceColor', 'k', 'LineStyle', 'none');
    legend(axV, [curves.shear curves.axles], {'SFD', 'Axles'}, 'Location', 'northeast');
    xlim(axV, [0 L]);
    ylim(axV, [-1.12 1.12] * max(1, Vmax));
    ylabel(axV, 'Shear V (N)');
    grid(axV, 'on');

    hold(axM, 'on');
    curves.moment = plot(axM, x, zeros(size(x)), 'r', 'LineWidth', 1.5);
    xlim(axM, [0 L]);
    ylim(axM, [-0.04 1.1] * max(1, Mmax));
    set(axM, 'YDir', 'reverse');      % Positive sagging is drawn downward
    xlabel(axM, 'Bridge coordinate x (mm)');
    ylabel(axM, 'Moment M (N mm)');
    title(axM, 'BMD (sagging positive, drawn down)');
    grid(axM, 'on');
end

function showFrame(curves, shear, moment, axle_x, L)
    a = axle_x(axle_x >= 0 & axle_x <= L);
    set(curves.axles, 'XData', a, 'YData', zeros(size(a)));
    set(curves.shear, 'YData', [0 shear 0]);
    set(curves.moment, 'YData', moment);
end
