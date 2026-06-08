#!/usr/bin/env bash
# Locate MATLAB license server details (port@host) on the local machine.
# Works on Linux, macOS, WSL, and Git Bash on Windows. Searches env vars,
# PATH, and common install / preferences locations.

set -u
found=0
declare -a roots=()
declare -a searched=()

echo "=== Environment variables ==="
env_hit=0
for var in MLM_LICENSE_FILE LM_LICENSE_FILE; do
    val="${!var:-}"
    if [ -n "$val" ]; then
        echo "  $var=$val"
        env_hit=1
        found=1
    fi
done
[ $env_hit -eq 0 ] && echo "  (none set)"

# --- matlab on PATH ---
if command -v matlab >/dev/null 2>&1; then
    bin="$(command -v matlab)"
    # Follow symlinks (e.g. /usr/local/bin/matlab -> .../R2024b/bin/matlab).
    real="$(readlink -f "$bin" 2>/dev/null || echo "$bin")"
    root="$(dirname "$(dirname "$real")")"
    echo
    echo "matlab on PATH:  $bin  ->  MATLABROOT=$root"
    roots+=("$root")
fi

# --- Common install + preferences locations ---
case "$(uname -s)" in
    Darwin)
        for d in /Applications/MATLAB_R*.app \
                 "$HOME/Library/Application Support/MathWorks/MATLAB"/*; do
            [ -d "$d" ] && roots+=("$d")
        done
        ;;
    Linux|*BSD*)
        for d in /usr/local/MATLAB/* /opt/matlab/* /opt/MATLAB/* \
                 "$HOME/.matlab"/*; do
            [ -d "$d" ] && roots+=("$d")
        done
        ;;
    MINGW*|MSYS*|CYGWIN*)
        # Git Bash / MSYS on Windows.
        for d in "/c/Program Files/MATLAB"/* \
                 "/c/Program Files (x86)/MATLAB"/* \
                 "/c/MATLAB"/* \
                 "$APPDATA/MathWorks/MATLAB"/* \
                 "$LOCALAPPDATA/MathWorks/MATLAB"/* \
                 "$PROGRAMDATA/MathWorks/MATLAB"/*; do
            [ -d "$d" ] && roots+=("$d")
        done
        ;;
esac

# Deduplicate.
if [ ${#roots[@]} -gt 0 ]; then
    mapfile -t roots < <(printf '%s\n' "${roots[@]}" | awk '!seen[$0]++')
fi

# --- Scan each root for FlexLM *.lic and modern license_info.xml ---
declare -a lic_modes=()

xml_extract() {
    # $1=file, $2=tag -- return inner text of first <tag>...</tag> in file.
    sed -n "s|.*<$2[^>]*>\([^<]*\)</$2>.*|\1|p" "$1" | head -n1
}

for root in "${roots[@]:-}"; do
    [ -n "$root" ] || continue
    for d in "$root/licenses" "$root"; do
        [ -d "$d" ] || continue
        searched+=("$d")

        # FlexLM-style .lic files (network license manager).
        for f in "$d"/*.lic; do
            [ -f "$f" ] || continue
            while IFS= read -r host hostid port; do
                if [ -n "${port:-}" ]; then
                    echo "  $f  ->  ${port}@${host}"
                    found=1
                fi
            done < <(awk 'tolower($1)=="server" {print $2, $3, $4}' "$f")
        done

        # Modern MATLAB license_info.xml manifest.
        for f in "$d"/license_info.xml; do
            [ -f "$f" ] || continue
            mode=$(xml_extract "$f" licmode)
            [ -z "$mode" ] && mode='(unspecified)'
            lic_modes+=("$f|$mode")

            # Network entries may carry host/port under historical names.
            xhost=$(xml_extract "$f" lichost)
            [ -z "$xhost" ] && xhost=$(xml_extract "$f" serverhost)
            [ -z "$xhost" ] && xhost=$(xml_extract "$f" host)
            xport=$(xml_extract "$f" licport)
            [ -z "$xport" ] && xport=$(xml_extract "$f" serverport)
            [ -z "$xport" ] && xport=$(xml_extract "$f" port)

            # Combined "<server>port@host</server>" form.
            srv=$(xml_extract "$f" server)
            if [ -n "$srv" ] && [[ "$srv" =~ ^[[:space:]]*([0-9]+)@([^[:space:]]+)[[:space:]]*$ ]]; then
                xport="${BASH_REMATCH[1]}"
                xhost="${BASH_REMATCH[2]}"
            fi

            if [ -n "$xhost" ] && [ -n "$xport" ]; then
                echo "  $f  ->  ${xport}@${xhost}"
                found=1
            fi
        done
    done
done

echo
echo "=== Searched MATLAB roots ==="
if [ ${#roots[@]} -eq 0 ]; then
    echo "  (none found)"
else
    printf '  %s\n' "${roots[@]}"
fi

echo
echo "=== Searched license directories ==="
if [ ${#searched[@]} -eq 0 ]; then
    echo "  (none)"
else
    printf '  %s\n' "${searched[@]}" | awk '!seen[$0]++'
fi

if [ ${#lic_modes[@]} -gt 0 ]; then
    echo
    echo "=== License modes declared in license_info.xml ==="
    printf '  %s\n' "${lic_modes[@]}" | awk -F'|' '{printf "  %-70s  %s\n", $1, $2}'
fi

echo
if [ $found -eq 0 ]; then
    echo "No MATLAB network license server (port@host) found on disk."
    if [ ${#lic_modes[@]} -gt 0 ]; then
        modes=$(printf '%s\n' "${lic_modes[@]}" | awk -F'|' '{print $2}' | sort -u | paste -sd, -)
        echo "Declared license mode(s): $modes"
        case "$modes" in
            *onlinelicensing*|*standalone*)
                echo
                echo "Your MATLAB is configured for $modes -- a sign-in / individual"
                echo "activation that does NOT use a port@host server. The Docker build's"
                echo "--build-arg LICENSE_SERVER=<port>@<host> flag will not work with this"
                echo "license type. To build the image you need either:"
                echo "  (a) a network license (port@host) from your institution's MATLAB admin, or"
                echo "  (b) a different deployment path (e.g. MATLAB Online, Code Ocean)."
                ;;
        esac
    elif [ ${#roots[@]} -gt 0 ]; then
        echo "MATLAB is installed but no license_info.xml or *.lic file was found."
        echo "MATLAB may not have been activated yet. Launch it once and complete"
        echo "activation, then re-run this script."
    else
        echo "MATLAB does not appear to be installed on this machine."
    fi
    echo
    echo "To query MATLAB itself for its license number, run:"
    echo "  matlab -batch \"disp(license); disp(getenv('MLM_LICENSE_FILE'))\""
    exit 1
fi

echo "Use a port@host value above with:"
echo "  docker build --build-arg LICENSE_SERVER=<value> ."
