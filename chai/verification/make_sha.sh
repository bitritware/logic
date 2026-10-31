#!/bin/sh
# Hash every file of this directory except SHA256SUMS itself.
cd "$(dirname "$0")"
ls | grep -v '^SHA256SUMS$' | LC_ALL=C sort | xargs sha256sum > SHA256SUMS
echo "wrote SHA256SUMS ($(wc -l < SHA256SUMS) files)"
