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

while [[ $# -gt 0 ]]; do
    key="$1"
    
    case $key in
	-h|--help)
	    echo "./sice_f.sh -i inpath "
	    echo "  -i: Path to folder containing csv inputs to sice.f"
	    exit 1;;
	-i)
	    INPATH="$2"
	    shift # past argument
	    shift # past value
	    ;;
    esac
done

if [ -z $INPATH ] ;then
    echo "-i option not set"
    echo " "
    $0 -h
    exit 1
fi

DEST=${INPATH}

timing
MSG_OK ${DEST}
# input for sice
# ns,alat,alon,sza,vza,saa,vaa,height,(toa(iks),iks=1,21), ozone, water
# The number of lines in this file MUST be equal to the number of lines in the file 'nlines.dat'
wc -l < ${DEST}/olci_toa_newformat.dat > ${DEST}/nlines.dat

# moving files to processor folder
cp ${DEST}/olci_toa_newformat.dat ./fortran/olci_toa_newformat.dat
cp ${DEST}/nlines.dat ./fortran/nlines.dat

end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
start=`date +%s`
# ===========  Running fortran SICE ======================
MSG_OK "Running sice.exe"
#gfortran ./fortran/sice.f -o ./fortran/sice.exe

cd ./fortran
./sice.exe
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
start=`date +%s`

# =========== copying back to home folder =========================

cd ..
cp ./fortran/bba_alex_reduced.dat		${DEST}/bba_alex_reduced.dat
cp ./fortran/bba.dat					${DEST}/bba.dat
cp ./fortran/boar.dat					${DEST}/boar.dat
cp ./fortran/planar_albedo.dat		${DEST}/planar_albedo.dat
cp ./fortran/spherical_albedo.dat		${DEST}/spherical_albedo.dat
cp ./fortran/impurity.dat				${DEST}/impurity.dat
cp ./fortran/nlines.dat				${DEST}/nlines.dat
cp ./fortran/notsnow.dat				${DEST}/notsnow.dat
cp ./fortran/size.dat					${DEST}/size.dat
cp ./fortran/retrieved_O3.dat 		${DEST}/retrieved_O3.dat

rm ./fortran/{bba_alex_reduced,bba,boar}.dat 
rm ./fortran/planar_albedo.dat
rm ./fortran/spherical_albedo.dat 
rm ./fortran/impurity.dat
rm ./fortran/nlines.dat 
rm ./fortran/olci_toa_newformat.dat 
rm ./fortran/notsnow.dat
rm ./fortran/size.dat
rm ./fortran/retrieved_O3.dat

timing
MSG_OK "Finished: ${folder}"
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
echo Total execution time was `expr $end - $start0` seconds.

	