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
docker run -it --name easd_nh \
	--volume=/data:/data \
	--volume=/rawdata:/rawdata \
	--volume=/resources:/resources \
	--volume=$(pwd)/tmp:/tmp \
	nuada/easd_nh
```
