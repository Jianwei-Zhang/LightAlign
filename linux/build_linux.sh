#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
command -v g++ >/dev/null
mkdir -p "$here/bin"
g++ -O3 -D_GNU_SOURCE -std=c++17 \
    "$here/src/help.cpp" "$here/src/step4_5.cpp" \
    "$here/src/step3.cpp" "$here/src/step2.cpp" \
    "$here/src/step1.cpp" "$here/src/functions.cpp" \
    "$here/src/MurmurHash3.cpp" "$here/src/preprocess_fasta.cpp" \
    -o "$here/bin/LightAlign" -pthread -lstdc++fs
"$here/bin/LightAlign" -t 12 -h
