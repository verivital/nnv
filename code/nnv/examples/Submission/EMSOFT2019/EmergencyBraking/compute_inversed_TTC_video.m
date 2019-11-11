load reachSet.mat;
n = length(S);
start = tic;
B = [];
for i=1:n
    if isa(S(i), 'Star')
        B = [B S(i).getBox]; % get box bounds the star set
    else
        B = [B S(i)];
    end
end

% Estimate minimum value of the function g = v^2 + 2ad
g_min =  zeros(n,1);
g_max = zeros(n,1);

for i=1:n
    
    fun_min = @(x)x(2)^2 + 2*x(1)*x(3);
    fun_max = @(x)-x(2)^2 - 2*x(1)*x(3);
    lb = B(i).lb;
    ub = B(i).ub;
    x0 = lb;   
    [~,fval] = fmincon(fun_min,x0,[],[],[],[],lb,ub);
    g_min(i) = fval;
    [~,fval] = fmincon(fun_max,x0,[],[],[],[],lb,ub);
    g_max(i) = -fval;
end

% compute TTC^-1
inv_TTC_min = zeros(n,1);
inv_TTC_max = zeros(n,1);

for i=1:n
    
    if g_max(i) <= 0
        
        inv_TTC_min(i) = 0;
        inv_TTC_max(i) = 0;
        
    elseif g_min(i) >= 0
              
        fun_min = @(x)(-x(3)/(x(2) - sqrt(x(2)^2 + 2*x(1)*x(3))));
        fun_max = @(x)(x(3)/(x(2) - sqrt(x(2)^2 + 2*x(1)*x(3))));
        lb = B(i).lb;
        ub = B(i).ub;
        x0 = lb; 
        if lb(3) == 0 && ub(3) == 0
            inv_TTC_min(i) = lb(2)/ub(1);
            inv_TTC_max(i) = ub(2)/lb(1);
        else
            [~,fval] = fmincon(fun_min,x0,[],[],[],[],lb,ub);
            inv_TTC_min(i) = fval;
            [~,fval] = fmincon(fun_max,x0,[],[],[],[],lb,ub);
            inv_TTC_max(i) = -fval;
        end

    elseif g_min(i) < 0 && g_max(i) > 0
        
        nonlcon = @feasiblecon;
        fun_max = @(x)(x(3)/(x(2) - sqrt(x(2)^2 + 2*x(1)*x(3))));
        lb = B(i).lb;
        ub = B(i).ub;
        x0 = lb;
        j = 0;
        while (j < 1000)
            for k=1:3
                    x0(k) = (ub(k) - lb(k)).*rand(1, 1) + lb(k);
            end
            x1 = fun_max(x0);
            if isreal(x1)
                break;
            end
        end
        [~,fval] = fmincon(fun_max,x0,[],[],[],[],lb,ub,nonlcon);
        inv_TTC_max(i) = -fval;
        inv_TTC_min(i) = 0;       
        
    end
    
end


% plot TTC^-1 reachable set
inv_TTC_d = [];
inv_TTC_acc = [];
inv_TTC_v = [];
for i=1:n
    lb = [inv_TTC_min(i); B(i).lb(1)];
    ub = [inv_TTC_max(i); B(i).ub(1)];
    inv_TTC_d = [inv_TTC_d Star(lb, ub)];
    lb = [inv_TTC_min(i); B(i).lb(3)];
    ub = [inv_TTC_max(i); B(i).ub(3)];
    inv_TTC_acc = [inv_TTC_acc Star(lb, ub)];
    lb = [inv_TTC_min(i); B(i).lb(2)];
    ub = [inv_TTC_max(i); B(i).ub(2)];
    inv_TTC_v = [inv_TTC_v Star(lb, ub)];
end

get_inv_TTC_time = toc(start);

N = length(inv_TTC_v);

inv_tau = zeros(N, 1);
for i=1:N
    tau = 2*(B(i).ub(2)/25);
    inv_tau(i) = 1/tau;
end

video = VideoWriter('ReachSet');
video.FrameRate = 2; 

open(video);

figure; 
for i=2:N
    subplot(2,1,1);
    times = 1:1:i;
    inv_TTC_v1 = inv_TTC_v(1:i);
    Star.plotRanges_2D(inv_TTC_v1, 2, times, 'b');
    xlabel('Time steps');
    ylabel('Velocity');
    xlim([1 i]);
    set(gca,'FontSize',16);
    
    subplot(2,1,2);
    plot(times, inv_tau(1:i), 'red');
    hold on;
    Star.plotRanges_2D(inv_TTC_v1, 1, times, 'b');
    xlabel('Time steps');
    ylabel('$$TTC^{-1}$$', 'interpreter', 'latex');
    xlim([1 i]);
    title('$$TTC^{-1}$$ over time', 'interpreter', 'latex');
    set(gca,'FontSize',16);

    % write video 
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(video, frame);
    
end

close(video);



function [c, ceq] = feasiblecon(x)
    c = -x(2)^2 - 2*x(1)*x(3);
    ceq = [];
end
