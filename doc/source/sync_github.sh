#!/bin/bash

make latexpdf
make html
cp -r _build/html/* docs/
git add -A
git commit -m 'update'
git push
