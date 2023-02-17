N = length(exactReachSet);

% find top X classes and set up bounds
im_lb = [];
im_ub = [];
for i=1:N
    [lb, ub] = exactReachSet(i).getRanges;
    
    % get top-5 classes (standard imagenet comparison for accuracy, i.e., 
    % one of the top 5 classes is correct)
    %[val_m,idx_m] = max(exactReachSet(i).getRanges);
    topX = 5;
    [topXvals,topXidx]=maxk(exactReachSet(i).getRanges(), topX);
    
    % get class names
    topXclasses = net.Layers(end).Classes(topXidx);
    topXclasses_str = cellstr(topXclasses); % string names in cell array
    
    lb = reshape(lb, [1000, 1]);
    ub = reshape(ub, [1000, 1]);
    im_lb = [im_lb lb];
    im_ub = [im_ub ub];
end

point1 = [1 1 topXidx(1)]; % top class
point2 = [1 1 topXidx(2)]; % 2nd class
S1 = [];
for i=1:N
    S1 = [S1 exactReachSet(i).project2D(point1, point2)];
end

point1 = [1 1 topXidx(1)]; % top class
point2 = [1 1 topXidx(3)]; % 3rd class

S2 = [];
for i=1:N
    S2 = [S2 exactReachSet(i).project2D(point1, point2)];
end

im_lb = min(im_lb, [], 2);
im_ub = max(im_ub, [], 2);


im_center = (im_lb + im_ub)/2;
err = (im_ub - im_ub)/2;

x = 1:1:net.Layers(end).OutputSize; % all classes
y = im_center;

figure;

% this plot can be a bit misleading, it isn't showing all the sets, but 
% rather only the centers, so comparing ranges on this plot is inaccurate,
% the plots below do this however as they plot the reach sets
subplot(2,2,[1 2]);
e = errorbar(x,y,err);
hold on;
e.LineStyle = 'none';
e.Color = 'red';
% show top-X classes in different colors and also with a dot
e = errorbar(x(topXidx(1:topX)),y(topXidx(1:topX)),err(topXidx(1:topX)),'b.');
e.LineStyle = 'none';
xlabel('Output Category ID', 'FontSize', 13);
ylabel('Range', 'FontSize', 13);

topXplots_str = {};
for i_str = 1 : topX
    topXplots_str{i_str} = [topXclasses_str{i_str}, ' (', num2str(topXidx(i_str)), ')'];
end
                   

dim = [0.2 0.5 0.3 0.3];
box_str = {['Top-', num2str(topX), ' classes'], topXplots_str{1:topX}};
annotation('textbox',dim,'String',box_str,'FitBoxToText','on');

% todo: put these on the same plot so it is easier to compare
subplot(2,2,3)
Star.plotBoxes_2D(S1, 1, 2, 'blue');
xlabel(topXplots_str{1}, 'FontSize', 13);
ylabel(topXplots_str{2}, 'FontSize', 13);

% % lower and upper bounds of each class
% xvs_l = exactReachSet(1).getRanges();
% xvs_u = exactReachSet(2).getRanges();
% xlb = xvs_l(topXidx(1));
% xub = xvs_u(topXidx(1));
% ylb = xvs_l(topXidx(2));
% yub = xvs_u(topXidx(2));
% 
% xlim([xlb * 0.95, xub * 1.05]);
% ylim([ylb * 0.95, yub * 1.05]); % todo: something wrong here, need
% projection or something possibly, or to increase arithmetic precision
% (vpa) since these values are so tiny

subplot(2,2,4)
Star.plotBoxes_2D(S2, 1, 2, 'blue');
xlabel(topXplots_str{1}, 'FontSize', 13);
ylabel(topXplots_str{3}, 'FontSize', 13);
%xlim([6, 6.5]);
%ylim([6, 6.5]);