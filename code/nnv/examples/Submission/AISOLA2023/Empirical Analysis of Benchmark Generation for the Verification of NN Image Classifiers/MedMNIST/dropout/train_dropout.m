%% Train all models using dropout for regularization

disp("Training all models with dropout layer...");

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

disp("Finished training all dropout models.")