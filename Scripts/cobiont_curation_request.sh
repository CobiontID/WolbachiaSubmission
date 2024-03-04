#!/usr/bin/env bash

# cobiont_curation_request.sh

grit_email="grit-jira@sanger.ac.uk"
sender_email=$2
yaml=$1

name=`basename $yaml .yaml`
echo "" \
	| mailx -s "$name is ready for curation" \
		-A $yaml $grit_email $sender_email

echo "Please review the following JIRA request."
echo $1
echo $2