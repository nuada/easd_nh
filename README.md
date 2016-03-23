# EASD New Horizon grant.

Scripts and utilities for GWAS data analysis.

## Docker
Create & tag image:
```
docker build -t nuada/easd_nh .
docker tag nuada/easd_nh:latest nuada/easd_nh:$(date +%F)
```

Create container:
```
mkdir tmp
chmod 1777 tmp
docker run -d --name easd_nh \
	--volume=/mnt/mfs/data/easd_nh:/data \
	--volume=/mnt/mfs/rawdata/HWI-H245/Beeline/easd_nh:/rawdata \
	--volume=/mnt/mfs/resources:/resources \
	--volume=$(pwd)/tmp:/tmp
	nuada/easd_nh
```
