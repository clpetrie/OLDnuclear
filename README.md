# nuclear
Here are some notes on how to do some stuff.

git clone https://github.com/clpetrie/nuclear
clone the archive

git fetch
to get information about the origin archive

git pull
to integrate the origin archive into yours

git push
to put your changes on the origin archive

git push origin --tags
pushes your tags to the origin

git add file
adds or stages a file for commit

git commit
commits staged file to local copy of archive

git tag -n
lists current tags

git tag -a v1.01
tags current version with version v1.01
and prompts for a longer description

git branch newbranchname
create a new branch

git checkout newbranchname
switch to branch newbranchname (it is simplest if everything is
committed in current branch)

git checkout master
back to master branch

git checkout master
git merge newbranchname
merges newbranchname into master

git branch -d newbranchname
deletes branch newbranchname (useful after merge to get rid duplicate
branch)

git config --global color.diff never
permanently gets rid of diff colors

git diff master..bugfix
view all differences between the branches master and bugfix

git diff --name-status master..bugfix
show file names that are different between the branches master and bugfix
