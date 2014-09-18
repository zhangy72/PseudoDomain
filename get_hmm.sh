#!/bin/bash -login
# this bash is used to get any pfam hmm file as $2 and $3 indicate
if [ $# -ne 2 ];then
  echo "get_hmm.sh -a pfam_accession"
  exit 0
fi
usage() {
	echo The options are: [-a] arg [-n] arg
	exit 0
}

while getopts ':a:n:' OPT
do
	case $OPT in
	a)  acc=$OPTARG
		flag=0
		#echo $acc
		;;
	n)  pfam_name=$OPTARG 
		flag=1
		#echo $pfam_name
		;;
	?) usage
		;;
	esac
done
if [ $flag -eq 0 ]; then #the user inputs accession number
	awk -v acc="$acc" '
								{
									if ($1 == "HMMER3/b") {
										head_line = $0
									} else if ($1 == "NAME") {
										pfam_name_line = $0
									} else if ($1 == "ACC" && $2 ~ acc) {
										flag = 1
										print head_line 
										print pfam_name_line
										print $0
									} else if (flag == 1) {
										print
										if ($1 == "//") {
											flag = 0
											exit 0
										}
									} 
								}
							' ~/Pfam-A-26.hmm
else
	echo "CAT $3"
	cat $3 | awk -v pfam_name="$pfam_name" '
								{
									if ($1 == "HMMER3/b") {
										head_line = $0
									} else if ($1 == "NAME" && $2 ~ pfam_name) {
										flag = 1
										print head_line
										print $0
									} else if (flag == 1) {
										print
										if ($1 == "//") {
											flag = 0
											exit 0
										}
									} 
								}
							' ~/Pfam-A-26.hmm
fi
										
						
	
