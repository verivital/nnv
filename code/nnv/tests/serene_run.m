function [all_results, failed]=serene_run(test_list)
  if isempty(test_list)
    test_list=["test_attacks", "test_vgg_attacks", "test_video_to_image", "test_c_code_gen","test_glpk_feasibility","test_io","test_cnn_cnn","test_DAGNN","test_fnn_ffnns","test_funcs_LogSig","test_funcs_poslin","test_funcs_relu","test_funcs_SatLin","test_funcs_SatLins","test_funcs_TanSig","test_layers_averagePooling","test_layers_batchNorm","test_layers_Conv2D","test_layers_fullyConnected","test_layers_layer","test_layers_layerS","test_layers_maxpool","test_layers_maxunpool","test_layers_relu","test_layers_TransposedConv2dLayer","test_SEGNET"];
  end

  all_results=[];
  failed=[];%not used yet
  for f=1:length(test_list)
    results=runtests(test_list(f));
    all_results=[all_results, results];
  end




end
