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

# input for sice
# ns,alat,alon,sza,vza,saa,vaa,height,(toa(iks),iks=1,21), ozone, water
# The number of lines in this file MUST be equal to the number of lines in the file 'nlines.dat'
wc -l < ${DEST}/olci_toa_newformat.dat > ${DEST}/nlines.dat

# moving files to processor folder
cp ${DEST}/olci_toa_newformat.dat ./Fortran/olci_toa_newformat.dat
cp ${DEST}/nlines.dat ./Fortran/nlines.dat

end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
start=`date +%s`
# ===========  Running FORTRAN SICE ======================
MSG_OK "Running sice.exe"
# gfortran ./Fortran/sice.f -o ./Fortran/sice.exe

cd ./Fortran
./sice.exe
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
start=`date +%s`

# =========== copying back to home folder =========================

cd ..
cp ./Fortran/bba_alex_reduced.dat		${DEST}/bba_alex_reduced.dat
cp ./Fortran/bba.dat					${DEST}/bba.dat
cp ./Fortran/boar.dat					${DEST}/boar.dat
cp ./Fortran/planar_albedo.dat		${DEST}/planar_albedo.dat
cp ./Fortran/spherical_albedo.dat		${DEST}/spherical_albedo.dat
cp ./Fortran/impurity.dat				${DEST}/impurity.dat
cp ./Fortran/nlines.dat				${DEST}/nlines.dat
cp ./Fortran/notsnow.dat				${DEST}/notsnow.dat
cp ./Fortran/size.dat					${DEST}/size.dat
cp ./Fortran/retrieved_O3.dat 		${DEST}/retrieved_O3.dat

rm ./Fortran/{bba_alex_reduced,bba,boar}.dat 
rm ./Fortran/planar_albedo.dat
rm ./Fortran/spherical_albedo.dat 
rm ./Fortran/impurity.dat
rm ./Fortran/nlines.dat 
rm ./Fortran/olci_toa_newformat.dat 
rm ./Fortran/notsnow.dat
rm ./Fortran/size.dat
rm ./Fortran/retrieved_O3.dat


timing
MSG_OK "Finished: ${folder}"
end=`date +%s`
echo Execution time was `expr $end - $start` seconds.
echo Total execution time was `expr $end - $start0` seconds.

	