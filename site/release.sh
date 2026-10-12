#!/bin/bash
# Archive the posters of this version as a GitHub release of the tag v<version>
# (site/_variables.yml), with the A0 PDFs and full PNGs of site/assets/posters,
# so they stay available once the website shows a later version.
# Usage (from the repo root, with the tag pushed and the posters of this
# version in site/assets/posters): site/release.sh
set -euo pipefail
cd "$(dirname "$0")"
version=$(sed -n 's/^version: *"\(.*\)"/\1/p' _variables.yml)
tag="v$version"
git show "$tag:site/_variables.yml" | grep -q "^version: *\"$version\"" ||
    { echo "tag $tag missing or not at version $version"; exit 1; }
gh release create "$tag" assets/posters/*.pdf assets/posters/*.png --verify-tag \
    --title "Cryosphere maps $tag" --notes "A0 posters (PDF) and 100 dpi PNGs of version $version."
