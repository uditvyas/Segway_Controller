function drawcartpend(y,m,M,L)
x = y(1);
th = y(3);

% dimensions
L = L/3;  % pendulum length
R = 0.2*sqrt(M);
mr = .1*sqrt(m); % mass radius

% positions
y = R; % cart vertical position

px = x + L*sin(th);
py = y + L*cos(th);

plot([-10 20],[0 0],'k','LineWidth',2)
hold on
rectangle('Position',[x-R,y-R,2*R,2*R],'Curvature',[1 1],'FaceColor',[1 0.1 0.1])

plot([x px],[y py],'k','LineWidth',2)
plot(x, y, 'b*')
plot(px, py, 'b*')
rectangle('Position',[px-mr/2,py-mr/2,mr,mr],'Curvature',1,'FaceColor',[.1 0.1 1])

% set(gca,'YTick',[])
% set(gca,'XTick',[])
xlim([-1 16]);
ylim([-2 2.5]);
set(gcf,'Position',get(0, 'Screensize'))
% box off
drawnow
hold off