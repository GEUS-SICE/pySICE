# output to raster

parallel --bar --verbose \
	"cat olci_toa.dat | sed 's/\t/,/g' | cut -d, -f2,3,{1}|grep -v NaN > olci_toa_{2}.csv" \
	::: $(seq 4 31) \
	:::+ sza vza saa vaa height toa1 toa2 toa3 toa4 toa5 toa6 toa7 toa8 toa9 toa10 \
	toa11 toa12 toa13 toa14 toa15 toa16 toa17 toa18 toa19 toa20 toa21 ozone water

grass -e -c ../../mask.tif ~/tmp/G_CSVll2GeoTIFFxy
grass ~/tmp/G_CSVll2GeoTIFFxy/PERMANENT/ --exec ../../CSVll2GeoTIFFxy.sh $(ls olci_toa_*.csv)