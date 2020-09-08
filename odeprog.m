function status = odeprog(t,y,flag,varargin)
% Adapted from:
% Tim Franklin
% Virginia Tech
% Jesse Norris
% Wake Forrest
% May 2006

global odeprogglobvar
global timeouts

if nargin < 3 || isempty(flag) 
    if(etime(clock,odeprogglobvar(8:13))>0.5)
        tfin=odeprogglobvar(1);
        sstrt=odeprogglobvar(2:7);
        figure(95);
        new_t = t(end);
        perc=new_t/tfin;
        area([new_t tfin-new_t;new_t tfin-new_t]);
        title([num2str(perc*100) '%'],'Interpreter','none');
        odeprogglobvar(8:13)=clock;
    end
    % Stop the computation after some elapsed time, in seconds.
    % Should be v large for v stiff problems.
    if (etime(clock,odeprogglobvar(2:7)) > 20000 & ishandle(95))
        % disp('Exiting rotation due to time elapsed...')
        timeouts = timeouts + 1;
        clear init_time
        if (ishandle(95))
            close(95)
        end
    end
else
    switch(flag)
    case 'init'
        odeprogglobvar=zeros(1,13);
        odeprogglobvar(1)=t(end);
        odeprogglobvar(2:7)=clock;
        odeprogglobvar(8:13)=clock;
        tfin=odeprogglobvar(1);
        sstrt=odeprogglobvar(2:7);
        figure(95); 
        set(gcf,'Position',[4,40,100,500],'MenuBar','none','NumberTitle','off');
        axes('Position',[0.5,0.05,0.25,0.8]);
        axis([1,2,0,tfin]);
        set(gca,'XTickLabel',[],'NextPlot','replacechildren');
        ylabel('Simulation Progress - Time (s)');
        title('0%','Interpreter','none');
        area([0 tfin;0 tfin]);
        uicontrol('Style', 'pushbutton', 'String', 'Abort','Position', [7 460 90 30], 'Callback', 'close(gcf)')
        pause(0.1);

    case 'done'    
        if(ishandle(95))
            close(95);
        end
    end
end
status = 0;
drawnow;