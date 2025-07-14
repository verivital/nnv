classdef ProbReach_ImageStar


    properties
        model
        LB
        de
        indices
        original_dim
        output_dim
        mode
        params
    end

    methods
        function obj = ProbReach_ImageStar(model,LB, UB,indices,SizeOut,mode,params)
            obj.model = model;
            obj.LB = LB;
            obj.de = UB-LB;
            obj.indices = indices;

            SizeIn = size(LB);            
            lenIn = length(SizeIn);
            if lenIn == 1 
                obj.original_dim = [SizeIn , 1, 1];
            elseif lenIn == 2 
                obj.original_dim = [SizeIn , 1];
            elseif lenIn == 3
                obj.original_dim = SizeIn;
            else
                obj.original_dim = SizeIn(1:3);
            end

            lenOut = length(SizeOut);
            if lenOut == 1 
                obj.output_dim = [SizeOut , 1, 1];
            elseif lenOut == 2 
                obj.output_dim = [SizeOut , 1];
            elseif lenOut == 3
                obj.output_dim = SizeOut;
            else
                obj.output_dim = SizeOut(1:3);
            end

            obj.mode = mode;

            thisFile = mfilename('fullpath');
            [dir0, ~, ~] = fileparts(thisFile);
            params.files_dir = fullfile(dir0, 'Temp_files_mid_run');
            obj.params = params;
        end

        function Matrix = mat_generator_no_third(obj, values )
            
            height = obj.original_dim(1);
            width = obj.original_dim(2);
            n_channel = obj.original_dim(3);
            Matrix = zeros(height, width, n_channel);

            N_perturbed = size(obj.indices,1);

            t = 0;
            for c = 1: n_channel
                for i = 1:N_perturbed
                    index = obj.indices(i,:);
                    t = t+1;
                    Matrix(index(1) , index(2) , c) = values(t);
                end
            end


        end

        function out = forward(obj, x)

            model_source = class(obj.model);

            switch model_source

                case 'SeriesNetwork'
                    out = obj.model.predict(x);

                case 'DAGNetwork'
                    out = obj.model.predict(x);

                case 'dlnetwork'
                    dlX = dlarray(x, obj.params.dlarrayType);
                    out = obj.model.predict(dlX);

                case 'NN'
                    out = obj.model.evaluate(x);

                otherwise
                    error("Unknown model source: " + model_source + ". We only cover NN, SeriesNetwork, dlnetwork and DAGNetwork.");
            end
        end

        function needs = Provide(obj)
            
            out_height = obj.output_dim(1);
            height = obj.original_dim(1);
            out_width = obj.output_dim(2);
            width = obj.original_dim(2);
            n_class = obj.output_dim(3);
            n_channel = obj.original_dim(3);
            N = obj.params.Nt;
            N_dir = obj.params.N_dir;
            trn_batch = obj.params.trn_batch;

            N_perturbed = size(obj.indices , 1);

            rng(0)
            Y = zeros(out_height,out_width,n_class,N);
            X = zeros(n_channel*N_perturbed , N);
            
            %%%%%%%%%%%%%%
            parfor i=1:N
                disp(i)
                Rand = rand(n_channel*N_perturbed,1);
                Rand_matrix = obj.mat_generator_no_third(Rand);
                d_at = zeros(height,width,n_channel);
                for c=1:n_channel
                    d_at(:,:,c) = obj.de(:,:,c) .* Rand_matrix(:,:,c) ;
                end
                Inp = single(obj.LB + d_at);
                X(:,i) = single(Rand);
                Y(:,:,:,i) = obj.forward(Inp);
            end
            %%%%%%%%%%%%%
            n1 = numel(Y(:,:,:,1));
            if n1 < 8000  %%% SVD Algorithm in MATLAB has shown scalability upto this extent.

                Y = reshape(Y, [n1 , N]);
                [Directions , ~] = eig(Y*Y'/N);
                cd(obj.params.files_dir)

            else



                %%%%  dtype = 'single';
                SizePerElement = 4;
                NumElement = out_height*out_width*n_class*N_dir;
                TotalSize = (NumElement * SizePerElement) / 1024^3;

                if TotalSize > 1.9
                    leng = floor(N_dir * 1.9 / TotalSize);
                    r = mod(N_dir , leng);
                    N_dir = N_dir - r;
                else
                    leng = N_dir;
                end

                ind = 0;
                last_i = N_dir/leng ;

                for i=1:last_i
                    eval(['Y' num2str(i) ' = Y(:,:,:,ind+1:ind+leng);' ]);
                    ind = ind+leng;
                end

                cd(obj.params.files_dir)
                Text = ' save("Direction_data.mat" ';
                for i=1:last_i
                    Text = [ Text  ' ,"Y' num2str(i) '" ']; %#ok<*AGROW>
                end
                Text = [ Text ');'];
                eval(Text);

                Text = 'clear Y1'; %%% You always need to have something here even if it is nonsense, otherwise it clears everything
                for i=1:last_i
                    Text = [ Text  ' Y' num2str(i) '  '];
                end
                eval(Text);

                % mat_file_path = fullfile(obj.params.files_dir, 'Direction_data.mat');
                mat_file_path = obj.params.files_dir;
                cd ..
                command = sprintf([ 'python Direction_trainer.py --mat_file_path "%s" '...
                    '--num_files %d --N_dir %d --batch_size %d --height %d --width %d '...
                    '--n_class %d '], ...
                    mat_file_path,  last_i, N_dir , trn_batch, out_height, out_width, n_class);

                status = system(command);

                cd(obj.params.files_dir)
                delete('Direction_data.mat')

                % load('directions.mat')
                % pyenv
                pyenv('Version', obj.params.py_dir)
                npz = py.numpy.load('directions.npz');
                Directions_py = py.numpy.array(npz{'Directions'});

                directions_list = cell(Directions_py.tolist());
                rows = cellfun(@(row) double(cellfun(@double, cell(row))), directions_list, 'UniformOutput', false);
                Directions = single(vertcat(rows{:}));


                clear Directions_py  npz
                delete('directions.npz')


                Y = reshape(Y, [n1 , N]);


            end

            C = 20* (  0.001*(mean(Y,2))   +  (0.05-0.001) * 0.5*( min(Y , [], 2) + max(Y , [] , 2) ));

            dYV  = Directions' * (Y - C );
            dims = obj.params.dims;
            epochs = obj.params.epochs;
            lr = obj.params.lr;

            % Text = ' save("Reduced_dimension.mat", "dYV" ';
            % for i=1:last_i
            %     Text = [ Text  ' ,"X' num2str(i) '" ']; %#ok<*AGROW>
            % end
            % Text = [ Text ' , "dims", "epochs", "lr");'];
            % eval(Text);

            save("Reduced_dimension.mat", "dYV", "X", "dims", "epochs", "lr");
            
            cd ..

            if strcmp(obj.mode, 'relu')
                system('python Trainer_ReLU.py')
            elseif strcmp(obj.mode, 'Linear')
                system('python Trainer_Linear.py')
            end


            cd(obj.params.files_dir)
            delete('Reduced_dimension.mat')

            if strcmp(obj.mode, 'relu')
                load("trained_relu_weights_2h_norm.mat")

                Layers = cell(1,3);
                W = double(W1);
                b = double(b1)';
                L = LayerS(W, b, 'poslin');
                Layers{1} = L;

                W = double(W2);
                b = double(b2)';
                L = LayerS(W, b, 'poslin');
                Layers{2} = L;

                W = double(W3);
                b = double(b3)';
                L = LayerS(W, b, 'purelin');
                Layers{3} = L;

                small_net = NN(Layers);

                delete('trained_relu_weights_2h_norm.mat')
                %%%%%%%%%%%%%%%%%%%%%%%%
            elseif strcmp(obj.mode, 'Linear')
                load("trained_Linear_weights_norm.mat")
                                
                W = double(W);
                b = double(b)';
                L = LayerS(W, b, 'purelin');
                Layers = {L};

                small_net = NN(Layers);

                delete('trained_Linear_weights_norm.mat')
                %%%%%%%%%%%%%%%%%%%%%%%%
            end



            dim2 = out_height*out_width*n_class;

            res_trn = abs( Y - ( Directions * evaluate(small_net , X)  + C ) );


            threshold_normal = obj.params.threshold_normal;
            res_max = max( res_trn ,[] ,2 );
            indexha = find(res_max < threshold_normal);
            res_max(indexha,1) = threshold_normal;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


            clear X  Y  res_trn

            Ns = obj.params.Ns;
            Nsp = N_dir;


            if Ns > Nsp
                thelen = Nsp;
            else
                thelen = Ns;
            end


            if Ns > thelen
                chunck_size = thelen;
                Num_chuncks = floor(Ns / chunck_size);
                remainder = mod(Ns , chunck_size);
            else
                chunck_size = Ns;
                Num_chuncks = 1;
                remainder = 0;
            end
            chunck_sizes = chunck_size * ones(1, Num_chuncks);
            test_data_run = zeros(1,Num_chuncks);
            res_test_time = zeros(1,Num_chuncks);
            if remainder ~= 0
                chunck_sizes = [chunck_sizes , remainder]; %#ok<NASGU>
                test_data_run = zeros(1, Num_chuncks+1);
                res_test_time = zeros(1,Num_chuncks+1);
            end



            Rs = zeros(1, Ns);
            ind = 0;

            for nc=1:length(chunck_sizes)

                rng(nc)

                len = chunck_sizes(nc);
                Y_test_nc = zeros(out_height,out_width,n_class,len);
                X_test_nc = zeros(n_channel*N_perturbed , len);

                tic
                %%%%%%%%%%%%%%
               parfor i=1:len
                    disp(['Part ' num2str(nc) ' Part ' num2str(i) '.'])
                    Rand = rand(n_channel*N_perturbed,1);
                    Rand_matrix = obj.mat_generator_no_third(Rand);
                    d_at = zeros(height,width,n_channel);
                    for c=1:n_channel
                        d_at(:,:,c) = obj.de(:,:,c) .* Rand_matrix(:,:,c) ;
                    end
                    % Inp = single(obj.LB + d_at);
                    Inp = obj.LB + d_at;
                    X_test_nc(:,i) = Rand;
                    Y_test_nc(:,:,:,i) = obj.forward(Inp);
                end
                %%%%%%%%%%%%%

                test_data_run(nc) = toc;



                Y_test = reshape(Y_test_nc, [dim2 , len] );

                clear Y_test_nc

                tic
                res_tst = abs(Y_test - ( Directions * evaluate( small_net , X_test_nc) + C ));
                Rs(ind+1:ind+len) = max( res_tst ./ res_max);
                res_test_time(nc) = toc;

                clear Y_test X_test_nc res_tst

                ind = ind + len;

            end



            Rs_sorted = sort(Rs);
            R_star = Rs_sorted(:,Ns-1);
            Conf = R_star*res_max;





            needs = struct;
            needs.Conf = Conf;
            needs.small_net = small_net;
            needs.Directions = Directions;
            needs.C = C;

        end


    

        function OS = ProbReach(obj)

            needs = obj.Provide();

            Conf = needs.Conf;
            N_perturbed = size(obj.indices,1);
            out_height = obj.output_dim(1);
            out_width = obj.output_dim(2);
            n_class = obj.output_dim(3);
            n_channel = obj.original_dim(3);

            dim2 = out_height*out_width*n_class;
            
            H = Star();
            H.C = sparse(1,dim2);
            H.d = 0;
            H.predicate_lb = -Conf;
            H.predicate_ub =  Conf;
            H.dim = dim2;
            H.nVar= dim2;

            I = Star( zeros(n_channel*N_perturbed,1) , ones(n_channel*N_perturbed,1) );

            if strcmp(obj.mode , 'relu')
                needs.small_net.reachMethod = 'approx-star';
            end
            Principal_reach = reach(needs.small_net , I);
            Surrogate_reach = affineMap(Principal_reach , needs.Directions , needs.C);
            % Out = MinkowskiSum(Surrogate_reach , H);

            NumElement = numel(Surrogate_reach.V) + dim2^2;
            SizePerElement = 8; % double
            TotalSize = (NumElement * SizePerElement) / 1024^3;

            if ispc
                % Windows systems
                [~, sys] = memory;
                ramGB = sys.PhysicalMemory.Total / 1024^3;
            elseif isunix
                % Linux/macOS systems: use Java
                runtime = java.lang.Runtime.getRuntime();
                ramBytes = runtime.totalMemory();
                ramGB = ramBytes / 1024^3;
            else
                error('Unsupported operating system');
            end

            if TotalSize < 0.5*ramGB
                H.V = [zeros(dim2,1) eye(dim2)];
                Out = Sum(Surrogate_reach, H);


                OS = ImageStar();
                OS.numChannel = n_class;
                OS.height = out_height;
                OS.width = out_width;
                OS.V = reshape(double(Out.V) , [out_height, out_width, n_class, Out.nVar+1]);
                OS.C = double(Out.C);
                OS.d = double(Out.d);
                OS.numPred = Out.nVar;
                OS.pred_lb = double(Out.predicate_lb);
                OS.pred_ub = double(Out.predicate_ub);
                if ~isempty(Out.state_lb)
                    OS.im_lb = reshape(double(Out.state_lb) , [out_height, out_width, n_class]);
                    OS.im_ub = reshape(double(Out.state_ub) , [out_height, out_width, n_class]);
                end
            else
                disp(' The Image Star is large for your memory and should be presented in sparse format.')
                disp('Unfortunately matlab does not support sparse representation for (N>2)D arrays.')
                disp('Thus we provide the vectorized format of ImageStar() that is a Star() via sparse 2D arrays. ')
                p = input('Do you want to continue? Yes <-- 1 / No <-- 0     ');
                
                if p==1
                    H.V = [sparse(dim2,1) speye(dim2)];
                    OS = Star();
                    P_Center = sparse(double(Surrogate_reach.V(:,1)));
                    P_lb = double([Surrogate_reach.predicate_lb ; H.predicate_lb]);
                    P_ub = double([Surrogate_reach.predicate_ub ; H.predicate_ub]);

                    P_V = [double(Surrogate_reach.V(:,2:end))   double(H.V(:,2:end))];
                    OS.V = [P_Center, P_V];
                    OS.C = blkdiag(sparse(Surrogate_reach.C) , sparse(H.C));
                    OS.d = [Surrogate_reach.d; H.d];
                    OS.predicate_lb = P_lb;
                    OS.predicate_ub = P_ub;
                    OS.nVar = Surrogate_reach.nVar + H.nVar;
                    OS.dim = H.dim;


                else
                    error('The Output ImageStar is too large and can not be presented as ImageStar().')
                end
                
            end
            
        dir = obj.params.currentDir;
        cd(dir)

     end
    end
end