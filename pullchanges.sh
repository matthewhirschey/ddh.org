#!/usr/bin/env bash
# Creates a new branch and pulls changes from ddh.com

# branch to hold changes
DESTBRANCH=ddh-com-$(date +"%Y-%m-%d_%H%M")

# upstream branch to pull changes from
UPBRANCH=org-reorg

echo "Switching to master branch and pulling"
# make sure we are on the ddh.org master branch and up to date
git checkout master
git pull

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

echo "Fetching $UPBRANCH remote changes"
git fetch upstream $UPBRANCH

echo "Creating and switching to new $DESTBRANCH branch"
git checkout -b $DESTBRANCH

echo "Merging all changes from upstream/$UPBRANCH"
git merge --squash --no-commit -X theirs upstream/$UPBRANCH

echo "Removing private files from staged changes"
PRIVATE_FILES=$(git diff --staged --name-only | grep '_private.R$')
for PRIVATE_FILE in $PRIVATE_FILES
do
   git rm -f $PRIVATE_FILE
done

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

echo "Done"

