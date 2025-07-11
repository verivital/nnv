% Name of the source file
sourceFile = 'Yolo_results_timed_out_1.txt';

% Loop to create copies (start from 2 to avoid overwriting itself)
for i = 2:72
    % Generate target filename
    targetFile = ['Yolo_results_timed_out_' num2str(i) '.txt'];
    
    % Copy the file
    copyfile(sourceFile, targetFile);
end

disp('All copies created.');
