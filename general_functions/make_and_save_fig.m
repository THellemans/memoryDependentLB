function [marker_list, line_list] = make_and_save_fig(x, y, legend_entries, legend_loc, xLab, yLab, xLim, yLim, filename, logplot, dist_points, personal_xvals)
%

[m,N]=size(y);
close all
h=figure;
marker_list=['d','x','o','+','s','*'];
line_list={'-','--','-.',':','-','--'};

if ~exist('dist_points', 'var')
    dist_points=ceil(N/20);
end
if ~exist('personal_xvals', 'var')
    personal_xvals=false;
end

for k=1:m
    if logplot
        if personal_xvals
            semilogy(x(k,:),y(k,:), 'k', 'LineWidth', 1.1,'Marker',marker_list(k),...
            'LineStyle', line_list{k},'MarkerSize',6, 'MarkerIndices',1:dist_points:N)
        else
            semilogy(x,y(k,:), 'k', 'LineWidth', 1.1,'Marker',marker_list(k),...
            'LineStyle', line_list{k},'MarkerSize',6, 'MarkerIndices',1:dist_points:N)
        end
    else
        if personal_xvals
            plot(x(k,:),y(k,:), 'k', 'LineWidth', 1.1,'Marker',marker_list(k),...
            'LineStyle', line_list{k},'MarkerSize',6, 'MarkerIndices',1:dist_points:N)
        else
            plot(x,y(k,:), 'k', 'LineWidth', 1.1,'Marker',marker_list(k),...
            'LineStyle', line_list{k},'MarkerSize',6, 'MarkerIndices',1:dist_points:N)
        end
    end
    hold on
end

legend(legend_entries,'Interpreter','latex','Location',legend_loc)
xlabel(xLab,'Interpreter','latex')
ylabel(yLab,'Interpreter','latex')
xlim(xLim)
ylim(yLim)
set(gca,'fontsize',14)
box on
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,char(filename),'-dpdf','-r0')

end