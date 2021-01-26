#!/usr/bin/env bash
# Creates a new branch, pulls changes from ddh.com, pushes the branch
# and creates a github pull request.
# This script requires the gh command can be installed.
# See https://cli.github.com/ to install gh.
# You should run `gh auth login` once to authenticate with github before
# running this script.

# turn on command exit status checking
set -e

# branch to hold changes
DESTBRANCH=ddh-com-$(date +"%Y-%m-%d_%H%M")

# upstream branch to pull changes from
UPBRANCH=org-reorg

unstagefiles () {
  PATTERN=$1
  echo "Unstaging files for pattern $PATTERN"
  for ADDED_FILE in $(git diff --staged --name-only --diff-filter=A | grep "$PATTERN")
  do
    echo "Removing $ADDED_FILE"
    git rm -f $ADDED_FILE
  done
  for MODIFIED_FILE in $(git diff --staged --name-only --diff-filter=M | grep "$PATTERN")
  do
    echo "Removing changes for $MODIFIED_FILE"
    git checkout HEAD -- $MODIFIED_FILE
  done
}

echo "Switching to master branch and pulling"
# make sure we are on the ddh.org master branch and up to date
git checkout master
git pull

# turn off command exit status checking for git ls-remote
set +e

# add a remote named upstream pointing to ddh.com
git ls-remote --exit-code upstream 2>/dev/null >/dev/null
RETVAL=$?
if [ $RETVAL -eq 0 ]
then
   echo "Upstream remote already exists."
else
   echo "Adding upstream remote to ddh.com"
   git remote add upstream git@github.com:matthewhirschey/ddh.com.git
fi

# turn on command exit status checking
set -e

echo "Fetching $UPBRANCH remote changes"
git fetch upstream $UPBRANCH

echo "Creating and switching to new $DESTBRANCH branch"
git checkout -b $DESTBRANCH

echo "Merging all changes from upstream/$UPBRANCH"
git merge --squash --no-commit -X theirs upstream/$UPBRANCH

echo "Removing private files from staged changes"
unstagefiles '_private.R$'

echo "Removing openshift files from staged changes"
unstagefiles "openshift/*"

echo "Removing tests/data files from staged changes"
unstagefiles "tests/data/*"

echo "Removing README changes from staged changes"
unstagefiles "README.md$"

echo "Committing changes to $DESTBRANCH"
MSGFILE=$(mktemp)
MSGDATE=$(date +"%Y-%m-%d")
echo "pulled changes from matthewhirschey/ddh.com on $MSGDATE" > $MSGFILE
echo "" >> $MSGFILE
echo "This commit was created using pullchanges.sh" >> $MSGFILE
echo "" >> $MSGFILE
echo "latest commit from ddh.com $UPBRANCH" >> $MSGFILE
git log -1 --stat upstream/$UPBRANCH >> $MSGFILE
echo "" >> $MSGFILE
git commit --file=$MSGFILE
rm $MSGFILE

echo "Pushing branch to origin(github)"
git push -u origin $DESTBRANCH

echo "Creating a pull request on github"
gh pr create --base master -w

echo "Done"

