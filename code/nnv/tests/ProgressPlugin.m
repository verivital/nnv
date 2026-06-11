classdef ProgressPlugin < matlab.unittest.plugins.TestRunnerPlugin
    % PROGRESSPLUGIN  Append a START/DONE line per test to a file, flushed each write.
    %
    % If MATLAB crashes mid-test (e.g. out-of-memory / segfault -> exit 255 on a CI
    % runner), the last "START <name>" with no matching "DONE <name>" identifies the
    % offending test -- which the crashed job's stdout would not reveal. Used by
    % run_shard so CI shard crashes are self-diagnosing via the uploaded artifact.

    properties
        ProgressFile (1,1) string
    end

    methods
        function plugin = ProgressPlugin(progressFile)
            plugin.ProgressFile = progressFile;
        end
    end

    methods (Access = protected)
        function runTest(plugin, pluginData)
            plugin.write("START", pluginData.Name);
            runTest@matlab.unittest.plugins.TestRunnerPlugin(plugin, pluginData);
            plugin.write("DONE ", pluginData.Name);
        end
    end

    methods (Access = private)
        function write(plugin, tag, name)
            fid = fopen(plugin.ProgressFile, 'a');
            if fid > 0
                fprintf(fid, '%s %s\n', tag, name);
                fclose(fid);   % close each write to force a flush to disk before the test runs
            end
        end
    end
end
