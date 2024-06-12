res_dir = "results_approx_acasxu/";

ss = dir(res_dir + "*.txt");

nP = length(ss); % number of properties

vRes(nP) = string;
vTime(nP) = string;

for i=1:nP
    fileName = res_dir + ss(i).name;
    iN = split(fileName, "_");
    iN = split(iN{end}, '.');
    iN = str2double(iN{1});
    fid = fopen(fileName, "r");
    vRes(iN) = fgetl(fid);
    vTime(iN) = fgetl(fid);
end

resTable = table(vRes', vTime');

unsat = sum(count(vRes, "unsat"));
sat = sum(count(vRes, "sat", "IgnoreCase",true));
unknown = sum(count(vRes, "unknown"));
sat = sat - unsat;
