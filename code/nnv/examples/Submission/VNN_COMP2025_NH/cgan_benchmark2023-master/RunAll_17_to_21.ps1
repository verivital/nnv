for ($i = 17 ; $i -le 21; $i++) {
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
        if ($sw.Elapsed.TotalSeconds -ge 1200) {
            Write-Host "Iteration $i timed out. Killing MATLAB..."
            Stop-Process -Id $process.Id -Force
            Set-Content -Path ("cgan_results_timed_out_" + $i + ".txt") -Value "timed out"
            break
        }
    }

    if ($process.HasExited -and $sw.Elapsed.TotalSeconds -lt 1200) {
        Write-Host "Iteration $i completed in $($sw.Elapsed.TotalSeconds) seconds."
    }
}

