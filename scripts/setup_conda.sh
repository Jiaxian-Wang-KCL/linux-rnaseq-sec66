#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."
[[ "$(uname -s)" == Linux && "$(uname -m)" == aarch64 ]] || {
  echo 'This installer is pinned to Linux aarch64.' >&2; exit 1;
}
installer=config/Miniforge3-Linux-aarch64.sh
if [[ ! -x "$HOME/miniforge3/bin/conda" ]]; then
  curl -fL --retry 3 -o "$installer" \
    https://github.com/conda-forge/miniforge/releases/download/26.7.2-0/Miniforge3-26.7.2-0-Linux-aarch64.sh
  echo "89b786c8d2c8b0fda7553914c1314ae4ddaa094503802f279377b19ac4463cb2  $installer" | sha256sum -c -
  bash "$installer" -b -p "$HOME/miniforge3"
fi
source "$HOME/miniforge3/etc/profile.d/conda.sh"
if [[ -f config/conda-linux-aarch64.lock.txt ]]; then
  conda create -y -n rnaseq --file config/conda-linux-aarch64.lock.txt
else
  CONDA_CHANNEL_PRIORITY=strict conda env create -f environment.yml
fi
