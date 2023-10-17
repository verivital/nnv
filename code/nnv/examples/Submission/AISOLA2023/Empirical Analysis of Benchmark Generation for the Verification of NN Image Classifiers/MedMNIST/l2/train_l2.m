%% Train all models using l2 regularization

disp("Training all models with L2 regularization...");

disp("... with glorot initialization...")
cd glorot;
medNist_training;
cd ..

disp("... with he initialization...")
cd he;
medNist_training;
cd ..

disp("... with narrow-normal initialization...")
cd narrow-normal;
medNist_training;
cd ..

disp("Finished training all L2 models.")