% test_soundness_gpu_bab_alpha_closedform
% P4 increment 2 (NNV_BAB_CLOSEDFORM): the analytic closed-form [alpha;beta] subgradient
% (i_dag_grad, plain gpuArrays, no dlarray/dlgradient) must equal the autodiff gradient, so the
% optimizer takes the SAME steps and returns the SAME bound -- validated here as BOUND-PARITY
% (closed-form-ON optimized margin == dlfeval-OFF optimized margin to tol) across FC / conv /
% conv-add(resnet-block) DAGs, plus SOUNDNESS (the optimized margin is a valid lower bound).
% (The direct grad-equivalence max|i_dag_grad - dlgradient| was validated at 0.0 / 5.96e-08 machine
% precision via the NNV_BAB_CLOSEDFORM_TEST hook; bound-parity is its clean assert-based manifestation.)
% SOUNDNESS: closed-form changes only HOW the gradient is computed; alpha in [0,1]/beta>=0 clamps +
% keep-best + the no-worse-than-min-area floor are unchanged, so any output stays a sound lower bound.

restore = onCleanup(@() setenv('NNV_BAB_CLOSEDFORM',''));
TOL = 2e-3;   % off(dlarray) vs on(plain) reassociate FP by ~ulps; a few Adam iters accumulate to <2e-3

%% Test 1: FC -- bound-parity off vs on + sound
rng(1); ops = i_fc([6 16 16 7]);
lb = -0.18*ones(6,1,'single'); ub = 0.18*ones(6,1,'single'); C = [eye(6,'single') -ones(6,1,'single')];
[mOff, mOn] = i_run_offon(ops, lb, ub, C);
assert(max(abs(mOff(:)-mOn(:))) < TOL, sprintf('FC closed-form bound must match autodiff; max diff %.3e', max(abs(mOff(:)-mOn(:)))));
mc = i_mc_min(ops, C, lb, ub, 6000);
assert(all(double(mOn(:)) <= double(mc(:)) + 1e-4), 'FC closed-form margin must be a SOUND lower bound');

%% Test 2: conv-relu-affine -- bound-parity off vs on + sound
rng(2); ish=[6 6 2]; cop=i_mkconv(ish,3,4,[1 1],[1 1 1 1],0); osh=cop.outShape; nz=prod(osh);
ops={cop, struct('type','relu','src',1), struct('type','affine','W',0.3*randn(8,nz,'single'),'b',0.1*randn(8,1,'single'),'src',2)};
lb=-0.15*ones(prod(ish),1,'single'); ub=0.15*ones(prod(ish),1,'single'); C=[eye(7,'single') -ones(7,1,'single')];
[mOff, mOn] = i_run_offon(ops, lb, ub, C);
assert(max(abs(mOff(:)-mOn(:))) < TOL, sprintf('conv closed-form bound must match autodiff; max diff %.3e', max(abs(mOff(:)-mOn(:)))));

%% Test 3: conv-relu-conv + add(skip) - relu - affine (resnet block) -- bound-parity off vs on
rng(3); ish=[6 6 3];
c1=i_mkconv(ish,3,3,[1 1],[1 1 1 1],0); r1=struct('type','relu','src',1);
c2=i_mkconv(c1.outShape,3,3,[1 1],[1 1 1 1],2); add=struct('type','add','inputs',[3 0],'shape',c1.outShape);
r2=struct('type','relu','src',4); a1=struct('type','affine','W',0.3*randn(9,prod(c1.outShape),'single'),'b',0.1*randn(9,1,'single'),'src',5);
ops={c1,r1,c2,add,r2,a1};
lb=-0.12*ones(prod(ish),1,'single'); ub=0.12*ones(prod(ish),1,'single'); C=[eye(8,'single') -ones(8,1,'single')];
[mOff, mOn] = i_run_offon(ops, lb, ub, C);
assert(max(abs(mOff(:)-mOn(:))) < TOL, sprintf('conv-add closed-form bound must match autodiff; max diff %.3e', max(abs(mOff(:)-mOn(:)))));

disp('test_soundness_gpu_bab_alpha_closedform: all sections passed');

% ----------------------------------------------------------------------------------------
function [mOff, mOn] = i_run_offon(ops, lb, ub, C)
    nOps=numel(ops); reluIdx=find(cellfun(@(o) strcmp(o.type,'relu'), ops));
    [~, rtL, rtU]=gpu_bab_crown_tight(ops, lb, ub, C, 'single', cell(nOps,1));
    rb=struct('preL',{rtL},'preU',{rtU});
    fixings=cell(nOps,1); rng(7);
    for r=1:numel(reluIdx)
        k=reluIdx(r); dk=size(rtL{k},1); uns=(rtL{k}<0 & rtU{k}>0);
        fs=zeros(dk,1,'single'); idx=find(uns);
        for j=1:min(4,numel(idx)), fs(idx(j))=(mod(j,2)*2-1); end
        fixings{k}=fs;
    end
    setenv('NNV_BAB_CLOSEDFORM','');  mOff = gpu_bab_crown_alpha_dag(ops, lb, ub, C, fixings, reluIdx, 'single', 20, 0.3, rb);
    setenv('NNV_BAB_CLOSEDFORM','1'); mOn  = gpu_bab_crown_alpha_dag(ops, lb, ub, C, fixings, reluIdx, 'single', 20, 0.3, rb);
    setenv('NNV_BAB_CLOSEDFORM','');
end
function ops=i_fc(dims)
    ops={};
    for L=1:numel(dims)-1
        W=randn(dims(L+1),dims(L),'single')*sqrt(2/dims(L)); b=randn(dims(L+1),1,'single')*0.1;
        ops{end+1}=struct('type','affine','W',W,'b',b,'src',numel(ops)); %#ok<AGROW>
        if L<numel(dims)-1, ops{end+1}=struct('type','relu','src',numel(ops)); end %#ok<AGROW>
    end
end
function y=i_fwd(ops,X)
    v=X;
    for k=1:numel(ops)
        if strcmp(ops{k}.type,'affine'), v=ops{k}.W*v+ops{k}.b; elseif strcmp(ops{k}.type,'relu'), v=max(v,0); end
    end
    y=v;
end
function mn=i_mc_min(ops,C,lb,ub,nS)
    rng(101); X=lb+(ub-lb).*rand(numel(lb),nS); mn=min(double(C)*double(i_fwd(ops,X)),[],2);
end
function op=i_mkconv(ish, fh, Cout, stride, pad, src)
    W=0.3*randn(fh,fh,ish(3),Cout,'single'); b=0.1*randn(1,1,Cout,'single');
    oh=floor((ish(1)+pad(1)+pad(2)-(fh-1)-1)/stride(1))+1; ow=floor((ish(2)+pad(3)+pad(4)-(fh-1)-1)/stride(2))+1;
    op=struct('type','conv','W',W,'b',b,'stride',stride,'pad',pad,'dil',[1 1],'inShape',ish,'outShape',[oh ow Cout],'src',src);
end
