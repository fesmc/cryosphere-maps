#!/bin/bash
# Render the site and publish it to the gh-pages branch of origin as a single
# commit without history (the branch is replaced each time), so the map tiles
# and poster PDFs do not accumulate in the repository.
# Usage (from the repo root, after scripts/web.jl): site/publish.sh
set -euo pipefail
cd "$(dirname "$0")"
quarto render
touch _site/.nojekyll
gitdir=$(git rev-parse --absolute-git-dir)
export GIT_INDEX_FILE="$(mktemp -d)/index"          # a separate index: the working tree stays untouched
tree=$(cd _site && git --git-dir="$gitdir" --work-tree=. add -A && git --git-dir="$gitdir" write-tree)
commit=$(git commit-tree "$tree" -m "Publish site ($(git rev-parse --short HEAD))")
git push --force origin "${commit}:refs/heads/gh-pages"
echo "published $commit to gh-pages"
