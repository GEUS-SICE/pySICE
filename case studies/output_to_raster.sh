#!/usr/bin/env bash 
#PATH=~/local/snap/bin:$PATH #./S3_proc.sh -i ./dat_S3A -o ./out_S3A
start=`date +%s`
start0=`date +%s`
RED='\033[0;31m'
ORANGE='\033[0;33m'
GREEN='\033[0;32m'
NC='\033[0m' # No Color
MSG_OK() { printf "${GREEN}${1}${NC}\n"; }
MSG_WARN() { printf "${ORANGE}WARNING: ${1}${NC}\n"; }
MSG_ERR() { printf "${RED}ERROR: ${1}${NC}\n"; }

timing() { if [[ $TIMING == 1 ]]; then MSG_OK "$(date)"; fi; }

while [[ $# -gt 0 ]]
do
    key="$1"
    
    case $key in
	-h|--help)
	    echo "./csv2tiff.sh -i inpath "
	    echo "  -i: Path to the csv file"
	    exit 1;;
	-i)
	    f="$2"
	    shift # past argument
	    shift # past value
	    ;;

    esac
done

# output to raster

grid_width=$(head -n 1 latitude.csv | grep -o -E '[0-9]+')
grid_height=$(( $(wc -l latitude.csv| awk '{print $1} ') / $grid_width ))
if [[ ! $grid_height =~ ^-?[0-9]+$ ]]; then 
	MSG_WARN "Grid width: ${grid_width}"
	MSG_WARN "Grid height: ${grid_height}"
	MSG_ERR "Width or height of csv files not integer"
	exit 1
fi

echo $grid_width
echo $grid_height

# if [ ! -f ${f}.vrt ]; then
cat <<-EOF > ${f}.vrt
<OGRVRTDataSource>
<OGRVRTLayer name="${f}_tmp">
<SrcDataSource>${f}_tmp.csv</SrcDataSource>
<GeometryType>wkbPoint</GeometryType>
<GeometryField encoding="PointFromColumns" x="alon" y="alat"/>
</OGRVRTLayer>
</OGRVRTDataSource>
EOF
echo VRT GENERATED
# fi

# if [ ! -f ${f}_3413.vrt ]; then
cat <<-EOF > ${f}_3413.vrt
<OGRVRTDataSource>
<OGRVRTLayer name="${f}_3413">
<SrcDataSource>${f}_3413.csv</SrcDataSource>
<GeometryType>wkbPoint</GeometryType>
<GeometryField encoding="PointFromColumns" x="X" y="Y"/>
</OGRVRTLayer>
</OGRVRTDataSource>
EOF
# fi

# adding header		
{ printf 'alat,alon,var\n'; cat ${f}.csv; } > ${f}_tmp.csv

# reprojecting
ogr2ogr -f CSV -t_srs EPSG:3413 -s_srs EPSG:4326 -lco GEOMETRY=AS_XY -nln ${f}_3413 ${f}_3413.csv ${f}.vrt
echo REPROJECTED

#converting to tiff			
gdal_grid -q -a count:radius1=10000:radius2=10000  -outsize ${grid_width} ${grid_height}  -a_srs EPSG:3413 -l ${f}_3413 ${f}_3413.vrt  concavehull_${f}_tmp.tif		
gdal_calc.py -A concavehull_${f}_tmp.tif --calc="A>1" --outfile scene_extent.tif

gdal_grid -q -l ${f}_3413 -zfield "var" -outsize ${grid_width} ${grid_height} -a_srs EPSG:3413 -a nearest:radius1:5000:radius2:5000:nodata=-999  ${f}_3413.vrt var_tmp.tif
gdal_calc.py --quiet -A var_tmp.tif -B scene_extent.tif --outfile=${f}.tif --calc="A*(B==1) -999*(B==0)" --NoDataValue=-999
	

# rm *_tmp.tif
# rm *_tmp.csv
# rm *.vrt