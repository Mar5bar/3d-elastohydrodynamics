function plot_ans(X,ts,skip,start,stop)
%% Basic plotter for output.
if nargin < 3
    skip = 1;
end
if nargin < 4
    start = 1;
end
if nargin < 5
    stop = length(ts);
end

f1=figure;
hold on
colormap('jet')
for i = start : skip : min(stop,length(ts))
    % Stop the drawing if the figure is closed.
    if ishghandle(f1) ~= true
        break
    end
    try delete(h1)
    end
    h1 = plot3(X(1,:,i),X(2,:,i),X(3,:,i),'black','LineWidth',2);
    axis equal
    grid on
    title(['T = ',num2str(ts(i))]);
    drawnow
end
end