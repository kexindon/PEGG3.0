#!/usr/bin/env bash
#
# PEGG 3.0 -- environment setup, installation, and the web visualizer.
#
#   ./pegg.sh install [ENV] [PYVER]   fresh conda env + PEGG        (default: pegg 3.9)
#   ./pegg.sh update  [ENV]           upgrade PEGG in an existing env
#   ./pegg.sh check   [ENV]           verify the install works
#   ./pegg.sh viz     [ENV]           run the web visualizer
#   ./pegg.sh notebook [ENV]          register a Jupyter kernel for VS Code
#
# PEGG needs Python 3.9 or 3.10. Newer Pythons cannot install the scikit-learn
# version the on-target scoring models were pickled with.

set -euo pipefail

ENV_NAME="${2:-pegg}"
PY_VERSION="${3:-3.9}"

# The visualizer lives in its own repository, expected beside this one.
VIZ_REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/prime-editing-visualization"

say()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
die()  { printf '\n\033[31mError: %s\033[0m\n' "$*" >&2; exit 1; }

# conda's shell function is not exported to non-interactive shells, so source
# its hook rather than relying on `conda activate` being already available.
load_conda() {
    local base
    base="$(conda info --base 2>/dev/null)" || die "conda not found. Install Miniconda first: https://docs.conda.io/en/latest/miniconda.html"
    # shellcheck disable=SC1091
    source "$base/etc/profile.d/conda.sh"
}

require_env() {
    conda env list | awk '{print $1}' | grep -qx "$ENV_NAME" \
        || die "No conda environment named '$ENV_NAME'. Create one:  ./pegg.sh install $ENV_NAME"
    conda activate "$ENV_NAME"
}

cmd_install() {
    case "$PY_VERSION" in
        3.9|3.10) ;;
        *) die "Python $PY_VERSION is not supported. Use 3.9 or 3.10 -- the scoring models are pickled with scikit-learn 1.1.1, which does not build on newer Pythons." ;;
    esac

    load_conda
    say "Creating conda environment '$ENV_NAME' (Python $PY_VERSION)"
    conda create -n "$ENV_NAME" python="$PY_VERSION" -y

    conda activate "$ENV_NAME"
    say "Installing PEGG from PyPI"
    pip install pegg

    cmd_check_body
    say "Done. Activate it with:  conda activate $ENV_NAME"
}

cmd_update() {
    load_conda
    require_env
    say "Upgrading PEGG in '$ENV_NAME'"
    # --upgrade-strategy only-if-needed: leave working dependencies alone, so an
    # upgrade cannot quietly pull numpy 2 (which breaks cyvcf2's compiled parts).
    pip install --upgrade --upgrade-strategy only-if-needed pegg
    cmd_check_body
}

# Import from a temporary directory: running inside a PEGG checkout would import
# ./pegg/ instead of the installed package, and the check would prove nothing.
cmd_check_body() {
    say "Verifying the installation"
    ( cd "$(mktemp -d)" && python - <<'PY'
import os, sys
import pegg
from pegg import prime, base, library, bystander
import numpy, sklearn, pandas

print('python       :', sys.version.split()[0])
print('pegg         :', os.path.dirname(pegg.__file__))
print('numpy        :', numpy.__version__)
print('scikit-learn :', sklearn.__version__)
print('pandas       :', pandas.__version__)

if 'site-packages' not in pegg.__file__:
    sys.exit('FAIL: imported a source checkout, not the installed package.')

# The .pkl/.pickle scorers and the canonical-transcript tables are package data;
# if packaging drops them, imports still succeed and only design fails.
data = [f for f in os.listdir(os.path.dirname(pegg.__file__))
        if not f.endswith('.py') and f != '__pycache__']
print('data files   :', len(data))
if len(data) < 11:
    sys.exit('FAIL: expected 11 data files, found %d.' % len(data))

df = pandas.DataFrame({'SEQ': ['ATGGCTAGCACCGGTG(C/T)CATCGGATCGGGCTAGCTAGGCTAAGCTTAGGCAT']})
out = prime.run(df, 'PrimeDesign', chrom_dict=None, silent_bystander=False)
print('pegRNAs      :', len(out))
scores = [c for c in out.columns if 'Score' in c]
if len(scores) != 3:
    sys.exit('FAIL: expected 3 score columns, got %r' % scores)
print('scorers      :', ', '.join(scores))

# ORF_start is required here: PrimeDesign input carries no genomic coordinates,
# so the reading frame cannot be looked up and has to be declared.
byst = prime.run(df, 'PrimeDesign', chrom_dict=None,
                 silent_bystander=True, ORF_start=0, seed=1)
print('bystander    : %d pegRNAs, %d carry a silent bystander'
      % (len(byst), int(byst['has_silent_bystander'].sum())))
print('\nOK')
PY
    ) || die "Verification failed -- see the output above."
}

cmd_check() { load_conda; require_env; cmd_check_body; }

cmd_notebook() {
    load_conda
    require_env
    say "Registering a Jupyter kernel"
    pip install --quiet ipykernel
    python -m ipykernel install --user \
        --name "$ENV_NAME" \
        --display-name "Python $(python -c 'import sys;print("%d.%d"%sys.version_info[:2])') (PEGG)"
    cat <<'EOF'

In VS Code: open the notebook, click the kernel picker (top right), then
  Select Another Kernel... -> Jupyter Kernel... -> "Python 3.x (PEGG)"

If it is not listed, run Developer: Reload Window from the command palette --
VS Code caches the environment list and will not notice a new one on its own.
EOF
}

cmd_viz() {
    load_conda
    require_env
    [ -d "$VIZ_REPO" ] || die "Visualizer not found at $VIZ_REPO
Clone it beside this repository:
  git clone https://github.com/kexindon/prime-editing-visualization.git \"$VIZ_REPO\""

    say "Installing the visualizer's requirements"
    pip install -r "$VIZ_REPO/requirements.txt"

    local port="${PORT:-5050}"
    say "Starting the visualizer on http://127.0.0.1:$port"

    # The visualizer prefers a sibling PEGG3.0 checkout over the installed
    # package, which is right while developing PEGG but wrong here: this script
    # sets up released versions. Point PEGG_PATH at a directory with no pegg/
    # inside it so that preference finds nothing and the import falls through to
    # the environment. Respect PEGG_PATH if the caller set one deliberately.
    local pegg_path="${PEGG_PATH:-$(mktemp -d)}"

    cat <<EOF
The app imports the PEGG installed in '$ENV_NAME'; /api/health reports which
copy is actually in use. Port $port rather than Flask's usual 5000, which
macOS occupies with its AirPlay Receiver.

  PORT=8080 ./pegg.sh viz $ENV_NAME              use a different port
  HOST=0.0.0.0 ./pegg.sh viz $ENV_NAME           reachable from other machines
  PEGG_PATH=/path/to/PEGG3.0 ./pegg.sh viz       run against a source checkout

Press Ctrl-C to stop.
EOF
    ( cd "$VIZ_REPO" && PEGG_PATH="$pegg_path" python run.py )
}

usage() {
    sed -n '3,13p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//'
    cat <<'EOF'

Examples:
  ./pegg.sh install                 create env 'pegg' on Python 3.9
  ./pegg.sh install pegg310 3.10    create env 'pegg310' on Python 3.10
  ./pegg.sh update pegg             upgrade PEGG in 'pegg'
  ./pegg.sh viz pegg                run the web visualizer
EOF
}

case "${1:-}" in
    install)  cmd_install  ;;
    update)   cmd_update   ;;
    check)    cmd_check    ;;
    viz)      cmd_viz      ;;
    notebook) cmd_notebook ;;
    ''|-h|--help|help) usage ;;
    *) die "Unknown command '${1}'. Run ./pegg.sh --help" ;;
esac
