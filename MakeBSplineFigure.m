% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

t_knot = [0;1;3;4;5;8;10];

t = linspace(0,10,1001)';

maxK = 4;
maxD = 4;
iSpline = 2;
iSubplot = 1;
h = zeros(maxK*maxD,1);

FigureSize = [50 50 figure_width_2col 400*scaleFactor];

figure('Units', 'points', 'Position', FigureSize)
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

for K=1:maxK
    B = bspline(t,t_knot,K);
    for D=1:maxK
        index=(K-1)*maxD+D;
        if D > min(K,maxD)
            % We create an empty subplot and set its visibility to zero,
            % otherwise packfig (below) creates a whitebox.
            h(index) = subplot(maxK,maxD,index);
            set(h(index),'visible','off')
        else
            h(index) = subplot(maxK,maxD,iSubplot);
            
            plot(t,squeeze(B(:,iSpline+K-1,D)), 'LineWidth', 1.0*scaleFactor, 'Color', 'k');
            ylim([-2 2])
            hold on
            y1=get(gca,'ylim');
            for i=1:length(t_knot)
                plot([t_knot(i) t_knot(i)],y1, 'LineWidth', 0.5*scaleFactor, 'Color', [0.5 0.5 0.5]);
            end
            set( gca, 'FontSize', figure_axis_tick_size);
            
            if (K<maxK)
                set(gca, 'XTick', []);
            end
            
            if (D>1)
                set(gca, 'YTick', []);
            else
                ylabel(sprintf('$K=%d$',K),'Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
            end
            
            if (D==K)
                if D==2
                    title(sprintf('$1^{st}$-deriv'),'Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
                elseif D==3
                    title(sprintf('$2^{nd}$-deriv'),'Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
                elseif D==4
                    title(sprintf('$3^{rd}$-deriv'),'Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
                end
                iSubplot = K*maxD + 1;
            else
                iSubplot = iSubplot + 1;
            end
        end
    end
end

packfig(maxK,maxD)
hfig = tightfig;
set(hfig,'position',FigureSize)

print('-depsc2', 'figures/bsplines.eps')