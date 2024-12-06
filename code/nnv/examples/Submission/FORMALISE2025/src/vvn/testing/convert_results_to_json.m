results_files = ["ZoomOut/C3D_small_robustness_results_FPV", "ZoomIn/C3D_small_robustness_results_FPV", "ScalabilityZoomOut/C3D_small_robustness_results_FPV", "ScalabilityZoomIn/C3D_small_robustness_results_FPV"];

for fp = 1:length(results_files)

    % Load the mat file
    filename = results_files(fp);
    filepath = "results/" + filename + ".mat";
    results = load(filepath);
    
    % Convert to JSON
    json_file = jsonencode(results, PrettyPrint=true);

    % Save as JSON file
    new_filename = "results/" + filename + ".json";
    fid = fopen(new_filename, "w");
    if fid == -1
        error("Cannot create JSON file.");
    end

    fwrite(fid, json_file, "char");
    fclose(fid);
end
    
