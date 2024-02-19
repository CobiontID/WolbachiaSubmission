#!/usr/bin/env bash

# request_taxids.sh
# taxon_request_file format -> "proposed_name\tname_type\thost\tproject_id\tdescription"
# Name type options:
#	- "Environmental Name"
#	- "Sythetic Name"
#	- "Novel Species"
#	- "Unidentified Species"
#	- "Published Species"

ena_email="ena-dtol@ebi.ac.uk"
sender_email=$2
taxon_request_file=$1

echo "Please review the attached taxonomy ID requests." \
	| mailx -s "Taxon ID requests for Wolbachia and other endosymbionts from dToL samples" \
		-A $taxon_request_file $ena_email $sender_email

#echo "Please review the attached taxonomy ID requests."
#echo $1
#echo $2
