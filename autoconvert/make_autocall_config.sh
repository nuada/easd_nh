#!/bin/bash

output=/data/easd_nh/autocall_config.txt

cat autocall_config.txt > ${output}

for array in $(tail -n +12 /data/easd_nh/sample_sheet.csv | cut -d , -f 5); do
	echo "D:\\Illumina\\iScan Control Software\\DataOutput\\${array}" >> ${output}
done

sudo unix2dos ${output}
