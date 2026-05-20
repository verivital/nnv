# Locate MATLAB license server details (port@host) on the local Windows machine.
# Searches env vars, PATH, registry, and common install / preferences locations.

$found = @()
$searched = @()
$matlabRoots = @()

# --- 1. Environment variables ---
foreach ($var in 'MLM_LICENSE_FILE','LM_LICENSE_FILE') {
    foreach ($scope in 'Process','User','Machine') {
        $val = [Environment]::GetEnvironmentVariable($var, $scope)
        if ($val) {
            $found += [PSCustomObject]@{ Source = "env:$var ($scope)"; Value = $val }
        }
    }
}

# --- 2. matlab.exe on PATH ---
$cmd = Get-Command matlab.exe -ErrorAction SilentlyContinue
if ($cmd) {
    $bin = $cmd.Source
    $root = Split-Path -Parent (Split-Path -Parent $bin)
    $matlabRoots += $root
    Write-Host "matlab.exe on PATH:  $bin  (MATLABROOT=$root)"
}

# --- 3. Registry: HKLM and HKCU MathWorks keys ---
foreach ($hive in 'HKLM:\SOFTWARE\MathWorks\MATLAB','HKLM:\SOFTWARE\WOW6432Node\MathWorks\MATLAB','HKCU:\SOFTWARE\MathWorks\MATLAB') {
    if (Test-Path $hive) {
        Get-ChildItem $hive -ErrorAction SilentlyContinue | ForEach-Object {
            $props = Get-ItemProperty $_.PSPath -ErrorAction SilentlyContinue
            foreach ($name in 'MATLABROOT','InstallationFolder','LASTINSTALLDIR') {
                if ($props.$name -and (Test-Path $props.$name)) {
                    $matlabRoots += $props.$name
                    Write-Host "Registry $($_.PSChildName).$name -> $($props.$name)"
                }
            }
        }
    }
}

# --- 4. Common install / preferences roots (glob for any R20xxx) ---
$globRoots = @(
    "$env:ProgramFiles\MATLAB",
    "${env:ProgramFiles(x86)}\MATLAB",
    "C:\MATLAB",
    "$env:ProgramData\MathWorks\MATLAB",
    "$env:APPDATA\MathWorks\MATLAB",
    "$env:LOCALAPPDATA\MathWorks\MATLAB"
)
foreach ($r in $globRoots) {
    if (Test-Path $r) {
        Get-ChildItem $r -Directory -ErrorAction SilentlyContinue |
            Where-Object { $_.Name -match '^R20\d{2}[ab]$' } |
            ForEach-Object { $matlabRoots += $_.FullName }
    }
}

$matlabRoots = $matlabRoots | Sort-Object -Unique

# --- 5. Hunt under each root: *.lic (FlexLM SERVER lines) and
#       license_info.xml (modern MATLAB activation manifest). ---
$licModes = @()

foreach ($root in $matlabRoots) {
    $licDirs = @("$root\licenses", "$root")
    foreach ($d in $licDirs) {
        if (Test-Path $d) {
            $searched += $d

            # FlexLM-style .lic files (network license manager).
            Get-ChildItem -Path $d -Filter '*.lic' -File -ErrorAction SilentlyContinue | ForEach-Object {
                Select-String -Path $_.FullName -Pattern '^\s*SERVER\s+(\S+)\s+\S+\s+(\d+)' -AllMatches |
                    ForEach-Object {
                        foreach ($m in $_.Matches) {
                            $found += [PSCustomObject]@{
                                Source = $_.Path
                                Value  = "$($m.Groups[2].Value)@$($m.Groups[1].Value)"
                            }
                        }
                    }
            }

            # Modern MATLAB license manifests (license_info.xml).
            Get-ChildItem -Path $d -Filter 'license_info.xml' -File -ErrorAction SilentlyContinue | ForEach-Object {
                $file = $_.FullName
                try {
                    [xml]$xml = Get-Content $file -Raw
                    foreach ($entry in @($xml.root.ActivationEntry)) {
                        if (-not $entry) { continue }
                        $mode = "$($entry.licmode)".Trim()
                        if (-not $mode) { $mode = '(unspecified)' }
                        $licModes += [PSCustomObject]@{ Source = $file; Mode = $mode }

                        # Network-license entries may include host/port fields under
                        # several historical names; harvest whichever are present.
                        $host_ = $null; $port = $null
                        foreach ($n in 'lichost','serverhost','host','server') {
                            if ($entry.$n) { $host_ = "$($entry.$n)".Trim(); break }
                        }
                        foreach ($n in 'licport','serverport','port') {
                            if ($entry.$n) { $port = "$($entry.$n)".Trim(); break }
                        }
                        if ($entry.server -and "$($entry.server)" -match '^\s*(\d+)@(\S+)\s*$') {
                            $port = $matches[1]; $host_ = $matches[2]
                        }
                        if ($host_ -and $port) {
                            $found += [PSCustomObject]@{ Source = $file; Value = "$port@$host_" }
                        }
                    }
                } catch {
                    Write-Host "  (warning: could not parse ${file}: $_)"
                }
            }
        }
    }
}

# --- Report ---
Write-Host ""
Write-Host "=== Searched MATLAB roots ==="
if ($matlabRoots) { $matlabRoots | ForEach-Object { Write-Host "  $_" } } else { Write-Host "  (none found)" }

Write-Host ""
Write-Host "=== Searched license directories ==="
if ($searched) { $searched | Sort-Object -Unique | ForEach-Object { Write-Host "  $_" } } else { Write-Host "  (none)" }

if ($licModes.Count -gt 0) {
    Write-Host ""
    Write-Host "=== License modes declared in license_info.xml ==="
    $licModes | Format-Table -AutoSize
}

Write-Host ""
if ($found.Count -eq 0) {
    Write-Host "No MATLAB network license server (port@host) found on disk."
    if ($licModes.Count -gt 0) {
        $modes = ($licModes.Mode | Sort-Object -Unique) -join ', '
        Write-Host "Declared license mode(s): $modes"
        if ($modes -match 'onlinelicensing|standalone') {
            Write-Host ""
            Write-Host "Your MATLAB is configured for $modes -- a sign-in / individual"
            Write-Host "activation that does NOT use a port@host server. The Docker build's"
            Write-Host "--build-arg LICENSE_SERVER=<port>@<host> flag will not work with this"
            Write-Host "license type. To build the image you need either:"
            Write-Host "  (a) a network license (port@host) from your institution's MATLAB admin, or"
            Write-Host "  (b) a different deployment path (e.g. MATLAB Online, Code Ocean)."
        }
    } elseif ($matlabRoots.Count -gt 0) {
        Write-Host "MATLAB is installed but no license_info.xml or *.lic file was found."
        Write-Host "MATLAB may not have been activated yet. Launch it once and complete"
        Write-Host "activation, then re-run this script."
    } else {
        Write-Host "MATLAB does not appear to be installed on this machine."
    }
    Write-Host ""
    Write-Host "To query MATLAB itself for its license number, run:"
    Write-Host '  matlab -batch "disp(license); disp(getenv(''MLM_LICENSE_FILE''))"'
    exit 1
}

Write-Host "=== Found license server details ==="
$found | Format-Table -AutoSize
Write-Host "Use a port@host value above with:"
Write-Host "  docker build --build-arg LICENSE_SERVER=<value> ."
