#!/bin/bash
version=$(curl -s https://api.github.com/repos/marianski-lab/artdep/tags | \
jq -r 'first(.[].name | select(test("^v?[0-9]")))')
echo $version
