function [Np , Nt , Ns]  = CP_specification(delta, confidence_LB, dv , train_mode, type)


if ispc
    % Windows systems
    [~, sys] = memory;
    ramGB = sys.PhysicalMemory.Total / 1024^3;
elseif isunix
    % Linux/macOS systems: use Java
    os = java.lang.management.ManagementFactory.getOperatingSystemMXBean();
    method = os.getClass().getMethod('getTotalPhysicalMemorySize', []);
    totalRamBytes = method.invoke(os, []);
    ramGB = double(totalRamBytes) / (1024^3);
else
    error('Unsupported operating system');
end

if strcmp(train_mode, 'gpu')


    gpm = gpuDevice().TotalMemory / 1024^3;
    Np = computeNumVectors(gpm, dv, type);


elseif strcmp(train_mode, 'cpu')


    if ispc
        % Windows systems
        Np = computeNumVectors(ramGB, dv, type);
    elseif isunix
        % Linux/macOS systems: use Java
        Np = computeNumVectors(ramGB, dv, type);
    else
        error('Unsupported operating system');
    end
    


else
    error('only gpu or cpu.')
end

Ns  = findMinimumM_binary(delta, confidence_LB);

Yolo_found = 8112;

Ntdv = min( [8000*Yolo_found,  Yolo_found*floor( 0.75 * ramGB*(1024^3) / (dv*8)  ),  Yolo_found*max(floor(Ns/10), 300)] );

Nt = floor(Ntdv/dv);

Nt = min(Nt , floor(Ns/10));

Npdv = max(100*Yolo_found, Yolo_found*min(Np , floor(Nt/3)));
Np = min(Np , floor(Npdv / dv));


if Np > Nt

    Np = Nt;

else
    
    r = mod(Nt, Np);
    if r < Nt/4
        Nt = Nt - r;
    else
        r1 = findMinimumR(Nt, Np);
        Np = Np - r1;
    end

end



end




function m = findMinimumM_binary(delta, confidence_LB)

    % Validate inputs
    if delta <= 0 || delta >= 1
        error('delta must be in (0,1)');
    end
    if confidence_LB <= 0 || confidence_LB >= 1
        error('confidence_LB must be in (0,1)');
    end

    % Define the function
    confFun = @(m) 1 - betacdf(delta, m - 1, 2);

    % Lower and upper bounds for m
    m_min = 2;
    m_max = 1e6;

    % Check feasibility at m_max
    if confFun(m_max) < confidence_LB
        error('No feasible m found up to m = 1e6. Increase m_max.');
    end

    % Binary search
    while m_min < m_max
        m_mid = floor((m_min + m_max)/2);
        conf_mid = confFun(m_mid);
        if conf_mid < confidence_LB
            % Need larger m
            m_min = m_mid + 1;
        else
            % Try smaller m
            m_max = m_mid;
        end
    end

    m = m_min;
end





function Np = computeNumVectors(gpm, dv, type)

    % Reference constants from System A
    refGpuMem_GB = 48;                 % GB
    refVectorDim = [720, 960, 12];     % vector dimension
    refVectorCount = 150;              % number of vectors

    % Compute elements
    refElems = prod(refVectorDim);
    currElems = prod(dv);

    % Effective memory per element in GB
    memPerElem_GB = (refGpuMem_GB / refVectorCount) / refElems;

    % Estimate number of vectors
    nVectors = floor(gpm / (currElems * memPerElem_GB));

    if strcmp(type , 'single')
        Np = nVectors;
    elseif strcmp(type , 'double')
        Np = floor(0.5*nVectors);
    else
        error('We only support single and double arrays.')
    end

end


function r = findMinimumR(Nt, Np)
% Find smallest r such that mod(Nt, Np - r) = 0
%
% Nt - large integer
% Np - smaller integer
%
% Returns:
% r - smallest nonnegative integer r < Np satisfying the condition

    if ~isscalar(Nt) || ~isscalar(Np) || Nt <= 0 || Np <= 0
        error('Nt and Np must be positive scalars.');
    end

    % Get all divisors of Nt
    divisors = unique(factorAll(Nt));

    % Keep only divisors <= Np
    divisors = divisors(divisors <= Np);

    % Compute r = Np - d
    rCandidates = Np - divisors;

    % Keep only r >=0
    rCandidates = rCandidates(rCandidates >=0);

    if isempty(rCandidates)
        error('No feasible r found.');
    end

    % Return minimal r
    r = min(rCandidates);
end
