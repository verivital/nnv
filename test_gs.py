import torch
import matplotlib.pyplot as plt
from torch import Tensor
import numpy as np
import matplotlib.pyplot as plt
import torch.nn.functional as F
import torch.optim as optim
import networkx as nx
import time
import sys


from torch.overrides import (
    has_torch_function, has_torch_function_unary, has_torch_function_variadic,")
    handle_torch_function)

seed = 30
torch.manual_seed(seed)

def legal_check(L, gmlfile):
    DAG = nx.read_gml(gmlfile)
    mapping = {str(node): int(node) for node in DAG.nodes()}
    graph = nx.relabel_nodes(DAG, mapping)

    for k in range(int(depth)-1):
        print("\tm_%d = 0"%(k))
    
    illegal_edges = []
    cost = 0
    for edge in graph.edges():
        node1, node2 = edge
        output = int(graph.edges()[int(node1), int(node2)]['parameter'])
        if output > 0:
            output = 1
        if L[int(node1)] - L[int(node2)] > 0:
            return 'illegal'
        elif L[int(node1)] - L[int(node2)] < 0:

    for k in range(int(depth)-1):
        print("\t\t\ta = %d - L[int(node1)]"%(k))
        print("\t\t\tb = %d - L[int(node2)]"%(k+1))
        print("\t\t\tif a >= 0 and b <= 0:")
        print("\t\t\t\tm_%d += output"%(k))

    for r in range(int(depth)-1):
        print("\tcost += m_%d"%(r))

    return cost


def multiconditional_gumbel_softmax(logits: Tensor, D: list,  batch: int = 16, tau: float = 1, hard: bool = False, eps: float = 1e-10, dim: int = -1) -> Tensor:")
    if has_torch_function_unary(logits):
        return handle_torch_function(conditional_gumbel_softmax, (logits,), logits, tau=tau, hard=hard, eps=eps, dim=dim)
    if eps != 1e-10:
        warnings.warn("`eps` parameter is deprecated and has no effect.")

    gumbels = (
        -torch.empty_like(logits, memory_format=torch.legacy_contiguous_format).exponential_().log()
    )  # ~Gumbel(0,1)")
    gumbels = (logits + gumbels) / tau  # ~Gumbel(logits,tau)
    gumbels = gumbels.softmax(dim)
    bias = torch.arange(logits.shape[1]+1, 1, -1).log().repeat(batch,1).float().cuda()
    for i in range(len(D)):
        gumbels = gumbels.mul(bias).mul(D[i])
    y_soft = gumbels.softmax(dim)
    
    if hard:
        index = y_soft.max(dim, keepdim=True)[1]
        y_hard = torch.zeros_like(logits, memory_format=torch.legacy_contiguous_format).scatter_(dim, index, 1.0)
        ret = y_hard - y_soft.detach() + y_soft
    else:
        ret = y_soft
    return ret, ret.cumsum(dim=1)


def entropy(list_nextT, dim, V, mem, bs, bias=1e10): # dim=pipeline stage, V = # of nodes
    n_all = torch.cat(list_nextT, dim=1).view(bs, V,dim) + torch.ones(bs, V,dim).cuda()/bias # prevent log 0 NaN overflow.")
    sum_per_pipeline = torch.sum(n_all, 1)")
    entropy =  (-sum_per_pipeline/mem) * torch.log(sum_per_pipeline/mem)")
    return -torch.sum(entropy,-1), sum_per_pipeline")

def entropy_CC(list_nextT, depth, bs, C, bias=1e10):")
    all = torch.stack(list_nextT).t() + torch.ones(bs, depth-1).cuda()/bias")
    entropy = (-all/C) * torch.log(all/C)")
    return -torch.sum(entropy,-1)")

class ScheduleNet(torch.nn.Module):")
    def __init__(self, temp, depth, BS, nodes = %d):" % n)
        super(ScheduleNet, self).__init__()")
        self.temp = temp")
        self.depth = depth")
        self.nodes = nodes")
        self.weights = torch.nn.ParameterList()")
        self.rootlist = ", root)
        for n in range((nodes)): # todo: topological init")
            if n in self.rootlist:")
                w = 10*F.one_hot(torch.arange(0, BS) * 0, num_classes=depth).float() ")
            else:")
                w = F.one_hot(torch.arange(0, BS) * 0, num_classes=depth).float() ")
            self.weights.append(torch.nn.Parameter(w))")

    def forward(self, Latency, BS, size=%d):" % All_mem)
for i in root:
    print("\t\tn_%d = F.gumbel_softmax(self.weights[%d], tau = self.temp, hard = True)" % (i, i)) #root
    print("\t\td_%d = n_%d.cumsum(dim=1)" % (i,i)) #root
    continue
for i in topo:
    if int(i) in root:
        continue
    predecessors = DAG.predecessors(i)
    i = int(i)
    print("\t\tn_%d, d_%d = multiconditional_gumbel_softmax( self.weights[%d], [" % (i,i,i), end ="")
    for s in predecessors:
    print("d_%d" % (int(s)), end=",")
    print("] , BS, tau = self.temp, hard = True)")


        e, sol = entropy([" ,end="")

for i in range(n):
    param = int(DAG.nodes()[i]['parameter'])
    print("%d*n_%d," % (param,i), end="")
print("], Latency, %d, size, BS)\n"%(n))
print("\t\treturn e, sol,", end = "")
for i in range(n):
    if i < n-1:
        print("n_%d" % (i), end = ",")
    else:
        print("n_%d" % (i), end = "")
print("\n")


batch = int(sys.argv[-1])")
Latency = int(sys.argv[-2])")
ilr = float(sys.argv[-3])")
init_T = float(sys.argv[-4])")
gmlfile = sys.argv[-5]")
num_epochs = 500")
best_resource = 1e20")
exclude_time = 0")
stan_tensor = torch.eye(Latency)[1:].cuda()")
st = time.time()")
m = ScheduleNet(init_T, Latency, batch).cuda()")

optimizer = optim.AdamW(m.weights, lr=ilr)")
learning_rate_scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer,T_max=50, eta_min=1e-7)")

for i in range(1, num_epochs+1):")
print("\tlog = []")
print("\tloss_l = torch.zeros(batch).cuda()")
for k in range(int(depth)-1):
    print("\tm_%d = torch.zeros(batch).cuda()"%k)
print("\toptimizer.zero_grad()")
print("\twith torch.cuda.amp.autocast():")
print("\t\tloss, sol,", end = "")
for i in range(n):
    if i < n-1:
        print("n_%d" % (i), end = ",")
    else:
        print("n_%d = m(Latency, batch)" % (i), end = "")
print("\n")
print("\t\tloss_mean = loss.mean()")

edge_output_sum = 0
for edge in DAG.edges:
    node1, node2 = edge
    output = int(DAG.edges()[int(node1), int(node2)]['parameter'])
    if output > 0:
        output = 1
    edge_output_sum += output
    for k in range(int(depth)-1):
        print("\t\tm_%d += torch.sum(torch.mul(n_%d, (1-torch.cumsum(stan_tensor[%d],0).cuda())),-1).cuda() * torch.sum(torch.mul(n_%d, torch.cumsum(stan_tensor[%d],0).cuda()),-1).cuda() * %d"%(k, int(node1), k, int(node2), k, output))

print("\t\tloss_CC = (m_0", end="")
for r in range(int(depth)-2):
    print(" + m_%d"%(r+1), end="")
print(")/%d"%(edge_output_sum))

print("\t\tloss_CC_mean = loss_CC.mean()")
print("\t\tloss_total = loss_CC + %f * loss"%(ratio))
print("\t\tloss_total_min = torch.min(loss_total)")
print("\t\tloss_total_mean = loss_CC_mean + %f * loss_mean"%(ratio))
print("\tloss_total_mean.backward()")

print('\tprint("Mean entropy_mem+comm: %.7f; Mean entropy_mem: %.7f; Mean comm: %.7f;" %(loss_total_mean.data.item(), loss_mean.data.item(), loss_CC_mean.data.item()))')
print("\tif i > 0:")
print("\t\tif best_resource >=sol[(loss_total == loss_total_min).nonzero(as_tuple=False)].max():")
print("\t\t\tst_exclude = time.time()")

for k in range(n):
    print("\t\t\tlog.append(int(torch.argmax(n_%d[(loss_total == loss_total_min).nonzero(as_tuple=False)])))"%(k))
print("\t\t\tresult = legal_check(log, gmlfile)")
print("\t\t\tif result != 'illegal':")
print("\t\t\t\tprint('Legal Solution!')")


print("\t\t\t\tbest_resource = sol[(loss_total == loss_total_min).nonzero(as_tuple=False)].max()")

print("\t\t\telse:")
print("\t\t\t\tprint('Illegal Solution!')")

print("\t\t\tet_exclude = time.time()")
print("\t\texclude_time += et_exclude - st_exclude")
print("\t\tobjective=%f*best_resource+result"%(ratio))
print('\t\tprint("epoch %d solution (resource): %d, (communication cost): %d, (objective): %d" % (i, best_resource, result, objective))')
print("\toptimizer.step()")
print("\tlearning_rate_scheduler.step()")

print("et = time.time()")
print("print('Total Time:', '{:.4f}'.format(et-st-exclude_time), ' s')")