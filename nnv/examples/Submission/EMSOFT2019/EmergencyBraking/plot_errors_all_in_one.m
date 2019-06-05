load errors.mat;

n1 = length(lb_poly);
n2 = length(lb_interval);
n3 = length(lb_star_exact);
n4 = length(lb_star_approx);

poly_vel_errs = [];
poly_dis_errs = [];
for i=1:n1
    err1 = max(abs(lb_poly(1, i)-lb_star_exact(1,i)), abs(ub_poly(1,i) - ub_star_exact(1,i)));
    err2 = max(abs(lb_poly(2, i)-lb_star_exact(2,i)), abs(ub_poly(2,i) - ub_star_exact(2,i)));
    
    %err1 = 100*err1/(ub_star_exact(1,i) - lb_star_exact(1,i));
    %err2 = 100*err2/(ub_star_exact(2,i) - lb_star_exact(2,i));
    
    poly_dis_errs = [poly_dis_errs err1];
    poly_vel_errs = [poly_vel_errs err2];
end

interval_vel_errs = [];
interval_dis_errs = [];
for i=1:n2
    err1 = max(abs(lb_interval(1, i)-lb_star_exact(1,i)), abs(ub_interval(1,i) - ub_star_exact(1,i)));
    err2 = max(abs(lb_interval(2, i)-lb_star_exact(2,i)), abs(ub_interval(2,i) - ub_star_exact(2,i)));
    
    %err1 = 100*err1/(ub_star_exact(1,i) - lb_star_exact(1,i));
    %err2 = 100*err2/(ub_star_exact(2,i) - lb_star_exact(2,i));
    
    interval_dis_errs = [interval_dis_errs err1];
    interval_vel_errs = [interval_vel_errs err2];
end

star_approx_vel_errs = [];
star_approx_dis_errs = [];
for i=1:n4
    err1 = max(abs(lb_star_approx(1, i)-lb_star_exact(1,i)), abs(ub_star_approx(1,i) - ub_star_exact(1,i)));
    err2 = max(abs(lb_star_approx(2, i)-lb_star_exact(2,i)), abs(ub_star_approx(2,i) - ub_star_exact(2,i)));
    
    %err1 = 100*err1/(ub_star_exact(1,i) - lb_star_exact(1,i));
    %err2 = 100*err2/(ub_star_exact(2,i) - lb_star_exact(2,i));
    
    star_approx_dis_errs = [star_approx_dis_errs err1];
    star_approx_vel_errs = [star_approx_vel_errs err2];
end


figure; % distance errors
subplot(2,2,1)
steps = 1:1:n1;
plot(steps, poly_dis_errs, '--*');
hold on;
plot(steps, poly_vel_errs, '--x');
legend('Distance error (m)', 'Velocity error (m/s)', 'fontsize', 8, 'Location', 'northwest');
xlabel('Step Number');
ylabel('Overapprox Error');
title('Polyhedron');
set(gca,'FontSize',14)
subplot(2,2,2)
steps = 1:1:n2;
plot(steps, interval_dis_errs, '--*', 'color', 'blue');
hold on;
plot(steps, interval_vel_errs, '--x', 'color', 'red');
legend('Distance error (m)', 'Velocity error (m/s)', 'fontsize', 8, 'Location', 'northwest');
xlabel('Step Number');
ylabel('Overapprox Error');
title('Interval');
set(gca,'FontSize',14)

subplot(2,2,[3, 4])
steps = 1:1:n4;
plot(steps, star_approx_dis_errs, '--*', 'color', 'black');
hold on;
plot(steps, star_approx_vel_errs, '--x', 'color', 'green');
legend('Distance error (m)', 'Velocity error (m/s)', 'fontsize', 8, 'Location', 'northwest');
xlabel('Step Number');
ylabel('Overapprox Error');
title('Approximate-Star');
set(gca,'FontSize',14)
