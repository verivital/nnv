for ($i = 1 ; $i -le 399; $i++) {
    Write-Host "Running iteration $i"

    $matlabCmd = "try, Run_iteration($i); disp('MATLAB finished iteration $i.'); exit; catch ex, disp(getReport(ex)); exit; end;"

    

    # Start stopwatch
    $sw = [System.Diagnostics.Stopwatch]::StartNew()

    # Start MATLAB
    $matlabPath = "C:\Program Files\MATLAB\R2024b\bin\matlab.exe"
    $process = Start-Process -FilePath $matlabPath `
        -ArgumentList "-nosplash", "-nodesktop", "-wait", "-r", "`"$matlabCmd`"" `
        -PassThru

    while (-not $process.HasExited) {
        Start-Sleep -Seconds 2
        if ($sw.Elapsed.TotalSeconds -ge 100) {
            Write-Host "Iteration $i timed out. Killing MATLAB..."
            Stop-Process -Id $process.Id -Force
            Set-Content -Path ("tinyimage_results_timed_out_" + $i + ".txt") -Value "timed out"
            break
        }
    }

    if ($process.HasExited -and $sw.Elapsed.TotalSeconds -lt 100) {
        Write-Host "Iteration $i completed in $($sw.Elapsed.TotalSeconds) seconds."
    }
}

