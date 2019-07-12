#!/bin/bash
git config --global user.email xiaohu@iastate.edu
git config --global user.name Rinoahu


git remote rm origin

git add -A .
git commit -m 'modify requirement'
git remote add origin https://github.com/Rinoahu/camelina_pangenome_analysis

git pull origin master
git push origin master

git checkout master
