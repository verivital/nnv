% Test InstanceNormalization2DLayer
% To run: results = runtests('test_InstanceNormalization2DLayer')

%% Test 1: Layer Not Yet Implemented
% InstanceNormalization2DLayer is not yet implemented in NNV
% This test is a placeholder for future implementation

% Check if class exists (exist returns 8 for class, 2 for .m file)
% Note: exist(..., 'file') returns 7 for folders, so we must be specific
classExists = exist('InstanceNormalization2DLayer', 'class') == 8 || ...
              (exist('InstanceNormalization2DLayer', 'file') == 2);

if classExists
    % Layer exists - run actual tests
    error('InstanceNormalization2DLayer now exists! Please implement proper tests.');
else
    % Layer doesn't exist yet - skip gracefully
    warning('test_InstanceNormalization2DLayer:NotImplemented', ...
        'Skipping test: InstanceNormalization2DLayer not yet implemented in NNV. This is a placeholder for future work.');
    assert(true, 'Test skipped - layer not yet implemented');
end
