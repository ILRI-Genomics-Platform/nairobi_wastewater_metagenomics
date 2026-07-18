#!/usr/bin/env bash
#SBATCH -p batch
#SBATCH -J centrifugedb
#SBATCH --qos vip
#SBATCH -n 16

set -uE
trap ' echo Error $? occured on $LINENO && exit 1' ERR

module load centrifuge/1.0.4

INDIR=$PWD
#DATEID="`date +"%Y%m%d"`_${src_dir}"
DATEID="20240222_"
OUTDIR=$INDIR/${DATEID}centrifugedb/
mkdir -p ${OUTDIR}/{taxonomy,library,centrifugeDb}

centrifugev1041="/home/gkibet/bioinformatics/github/metagenomics/scripts/centrifuge-1.0.4.1/"
centrifugeDownload=${centrifugev1041}/centrifuge-download
centrifugeBuild=${centrifugev1041}/centrifuge-build
#Create a custom centrifuge database that has:

# step 1
# Download NCBI taxonomy
echo -e "\tProceeding with building a custom database.\n\tStep 1: Downloading taxonomy information from NCBI..."
${centrifugeDownload} -P 30 -o ${OUTDIR}/taxonomy taxonomy
if [ $? -eq 0 ]
then
	echo -e "\tStep 1 successfully completed..."
elif [ $? -ne 0 ]
then
	echo -e "\tDownloading Taxonomy data from NCBI NOT successful..."
fi

# Step 2
# Download genomes: archaea,bacteria,fungi,viral,protozoa
echo -e "\nStep 2: Downloading genomes from NCBI..."
echo -e "\tAll genomes from NCBI..."
${centrifugeDownload} -P 16 \
	-o ${OUTDIR}/library \
	-m \
	-d "archaea,bacteria,fungi,viral,protozoa" \
	refseq >> ${OUTDIR}/seqid2taxid.map

echo -e "\tHuman genomes from NCBI..."
${centrifugeDownload} -P 16 -o ${OUTDIR}/library -d "vertebrate_mammalian" -a "Chromosome" -t 9606 -c 'reference genome' refseq >> ${OUTDIR}/seqid2taxid.map
if [ $? -eq 0 ]
then
	echo -e "\tStep 2 successfully completed..."
elif [ $? -ne 0 ]
then
	echo -e "\tDownloading genomes from NCBI NOT successful..."
fi


# Step 3
# Build the database
echo -e "\tStep 3: Building custom kraken2 database from NCBI data..."
#cat ${OUTDIR}/library/*/*.fna > ${OUTDIR}/input-sequences.fna
#https://askubuntu.com/questions/1028197/file-list-command-line-hidden-and-subfolders
echo -e "Compiling a list of all Genome files..."
find ${OUTDIR}/library/ \
	-maxdepth 2 \
	-type f \
	-name *.fna \
	! -name ${OUTDIR}/genomeFlist.txt \
	-printf  '%f,%h,%s,%TD %Tr\n' > ${OUTDIR}/genomeFlist.txt
sed -i '1s/^/genomeFileName,subDirPath,fileSizeKB,lastModDateTime\n/' ${OUTDIR}/genomeFlist.txt
#copy FASTA files/sequences to one file:
#cat library/*/*.fna > input-sequences.fna
{
        read
        while IFS="" read -r p || [ -n "$p" ];
        do
                genomeFile=$(echo $p | cut -f1 -d,);
                dirPath=$(echo $p | cut -f2 -d,);
		src_dir=${dirPath##*/};
                fileSize=$(echo $p | cut -f3 -d,);
                dateTime=$(echo $p | cut -f4 -d,);
                echo -e "Processing ${genomeFile} in ${src_dir}";
	cat ${dirPath}/${genomeFile} >> ${OUTDIR}/input-sequences.fna
	done
}<${OUTDIR}/genomeFlist.txt
#Build a centrifuge database
${OUTDIR}/library/*/*.fna
echo -e "Building the centrifuge database..."
cd ${OUTDIR}
centrifuge-build -p 16 \
	--bmax 1342177280 \
	--conversion-table ${OUTDIR}/seqid2taxid.map \
	--taxonomy-tree ${OUTDIR}/taxonomy/nodes.dmp \
	--name-table ${OUTDIR}/taxonomy/names.dmp \
	${OUTDIR}/input-sequences.fna ${OUTDIR}/centrifugeDb/centrifugeDb
if [ $? -eq 0 ]
then
        echo -e "\tBuilding the centrifuge database successfully completed..." 
elif [ $? -ne 0 ]
then
        echo -e "\tBuilding the centrifuge database from nt database NOT successful..." 
fi
cd ${INDIR}
# Step 4
# Remove intermediate files to save on space
#echo -e "\tStep 4: Removing intermediate files to save on space..."
# rm ${OUTDIR}/library/*/*.fna
# if [ $? -eq 0 ]
# then
# 	echo -e "\tStep 4 successfully completed..."
# elif [ $? -ne 0 ]
# then
# 	echo -e "\tRemoving intermediate files NOT successful..."
# fi
