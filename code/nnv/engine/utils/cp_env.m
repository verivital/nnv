function python_path= cp_env()
%python environment path for cp methods
nnv_path = nnvroot();
python_path = nnv_path + "/.venv/bin/python";
end

