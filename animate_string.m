function animate_string(t,V,p,filename,writer,mode_number,omega)
own_writer = nargin < 5;
if own_writer
    writer = VideoWriter(filename,'MPEG-4');
    writer.FrameRate = 30;
    writer.Quality = 100;
    open(writer);
    video_cleanup = onCleanup(@() close(writer));
end
fig = figure('Units','pixels','Position',[50 50 1920 1080], ...
    'Color','w','Visible','off','MenuBar','none','ToolBar','none');
fig_cleanup = onCleanup(@() close(fig));
ax = axes(fig,'Position',[0.085 0.24 0.84 0.49], ...
    'TickLabelInterpreter','latex','FontSize',22,'LineWidth',1.2);
hold(ax,'on');
x = linspace(0,p.L,p.n+2);
plot(ax,[0 p.L],[0 0],'--','Color',[0.78 0.81 0.84],'HandleVisibility','off');
line_plot = plot(ax,x,zeros(size(x)),'-','Color',[0.12 0.32 0.52],'LineWidth',3);
mass_plot = plot(ax,x(2:end-1),zeros(1,p.n),'o','MarkerSize',14, ...
    'MarkerFaceColor',[0.10 0.55 0.68],'MarkerEdgeColor',[0.07 0.24 0.34],'LineWidth',1.4);
fixed_plot = plot(ax,0,0,'s','MarkerSize',14,'MarkerFaceColor',[0.25 0.28 0.32], ...
    'MarkerEdgeColor',[0.25 0.28 0.32]);
driven_plot = plot(ax,p.L,0,'d','MarkerSize',15,'MarkerFaceColor',[0.91 0.46 0.14], ...
    'MarkerEdgeColor',[0.7 0.28 0.08]);
axis(ax,'equal');
xlim(ax,[-0.06 1.06]*p.L);
ylim(ax,[-0.14 0.14]*p.L);
set(ax,'XLimMode','manual','YLimMode','manual','DataAspectRatioMode','manual');
grid(ax,'on');
ax.GridAlpha = 0.12;
xlabel(ax,'Horizontal position, $x$ (m)','Interpreter','latex','FontSize',24);
ylabel(ax,'Displacement, $u$ (m)','Interpreter','latex','FontSize',24);
if nargin >= 6
    heading = sprintf('Vibrating string: $n=%d$, mode %d',p.n,mode_number);
    subheading = sprintf('$u_f(t)=A\\cos(\\omega t),\\quad \\omega=%.4f$ rad/s, $A=%.3g$ m',omega,p.Uf_func(0));
else
    heading = sprintf('Vibrating string: $n=%d$ masses',p.n);
    subheading = 'Fixed left end; prescribed right-end motion';
end
title(ax,heading,'Interpreter','latex','FontSize',32);
annotation(fig,'textbox',[0.06 0.84 0.88 0.065],'String',subheading, ...
    'Interpreter','latex','FontSize',24,'EdgeColor','none','HorizontalAlignment','center');
annotation(fig,'textbox',[0.06 0.77 0.88 0.06], ...
    'String',sprintf('$M=%.1f$ kg, $T_f=%.1f$ N, $L=%.1f$ m, $c=%.3f$ N s',p.M,p.Tf,p.L,p.c), ...
    'Interpreter','latex','FontSize',22,'EdgeColor','none','HorizontalAlignment','center');
legend(ax,[line_plot mass_plot fixed_plot driven_plot], ...
    {'String','Masses','Fixed endpoint','Driven endpoint'}, ...
    'Interpreter','latex','Orientation','horizontal','Location','southoutside','FontSize',20,'Box','off');
clock_label = annotation(fig,'textbox',[0.1 0.09 0.8 0.06],'String','', ...
    'Interpreter','latex','FontSize',25,'EdgeColor','none','HorizontalAlignment','center');
if nargin >= 6
    annotation(fig,'textbox',[0.08 0.02 0.84 0.05], ...
        'String','Cosine-driven response at a natural frequency; steady-state initial conditions. Equal axis scaling.', ...
        'Interpreter','latex','FontSize',17,'EdgeColor','none','HorizontalAlignment','center');
end
frame_times = (t(1):1/writer.FrameRate:t(end)-0.5/writer.FrameRate)';
states = interp1(t,V(:,1:p.n),frame_times);
[~,preview_index] = max(sum(states.^2,2));
for i = 1:numel(frame_times)
    uf = p.Uf_func(frame_times(i));
    line_plot.YData = [0 states(i,:) uf];
    mass_plot.YData = states(i,:);
    driven_plot.YData = uf;
    clock_label.String = sprintf('$t=%.2f$ s',frame_times(i));
    drawnow;
    frame = getframe(fig);
    if i == 1
        assert(size(frame.cdata,2) >= 1920 && size(frame.cdata,1) >= 1080,'Video frame is smaller than 1080p.');
    end
    assert(isequal(ax.XLim,[-0.06 1.06]*p.L) && isequal(ax.YLim,[-0.14 0.14]*p.L));
    assert(isequal(ax.DataAspectRatio,[1 1 1]));
    if i == preview_index && nargin >= 6
        imwrite(frame.cdata,sprintf('mode_n%d_%d.png',p.n,mode_number));
    end
    writeVideo(writer,frame);
end
end
