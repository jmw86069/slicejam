#!/bin/bash
#
# version 0.0.41.900
#
# merge, slop, split, count peaks
#
#
# PEAKS="peakfile1.bed peakfile2.bed peakfile3.bed"
# BAMS="bamfile1.bam bamfile2.bam bamfile3.bam"
#
# CHROMSIZES="hg19.chromsizes.txt"
# BEDTOOLS="/path/to/bin/bedtools"
# BEDOPS="/path/to/bedops-v2.4.35/bin"
# DATE="17jun2021"

SLICEJAM_VERSION="0.0.41.900"

## BED peak files
#PEAKS=($(ls *VEH*_25*sing*bed));
PEAKFILES=(${PEAKS})

## BAM alignment files
#BAMS=($(ls bams/*VEH*_25*sing*bam));
BAMFILES=(${BAMS});

####################################################################
## Define ANSI colors for output
export    RESET="\e[0m"    # reset
export     BOLD="\e[1m"    # hicolor/bold
export   FBLACK="\e[30m"   # foreground black
export     FRED="\e[31m"   # foreground red
export   FGREEN="\e[32m"   # foreground green
export  FYELLOW="\e[33m"   # foreground yellow
export    FBLUE="\e[34m" # foreground blue
export FMAGENTA="\e[35m"   # foreground magenta
export    FCYAN="\e[36m"   # foreground cyan

####################################################################
## Print documentation when PEAKS is empty
if [ "." == ".${PEAKFILES}" ]; then
   echo -e "${BOLD}mergeSplitCountBed.sh${RESET}: (version: ${SLICEJAM_VERSION})";
   echo -e "=====================";
   echo -e "This program:";
   echo -e " - takes bed peak files (${BOLD}BEDS${RESET})";
   echo -e " - merges peaks within ${BOLD}PEAKGAP${RESET} distance together";
   echo -e " - extends the edges of the merged peaks by ${BOLD}PEAKSLOP${RESET} bases";
   echo -e " - splits merged regions larger than ${BOLD}PEAKMAX${RESET} width into slices";
   echo -e "     of ${BOLD}PEAKMAX${RESET} width."
   echo -e " - counts peaks per merged peak for each file in ${BOLD}BEDS${RESET}";
   echo -e " - if ${BOLD}BAMS${RESET} is defined, runs featureCounts for coverage per";
   echo -e "     merged peak for each file in ${BOLD}BAMS${RESET}";
   echo -e "";
   echo -e "Usage:";
   echo -e 'PEAKS=$(ls *bed) BAMS=$(ls *bam) DEBUG=0 DRYRUN=1 ./mergeSplitCountBed.sh';
   echo -e 'PEAKS=$(ls *VEH*_25*sing*bed) BAMS=$(ls bams/*VEH*_25*sing*bam) DRYRUN=1 ./mergeSplitCountBed.sh';
   echo -e "";
   echo -e "Parameters:";
   echo -e "   PEAKS       (required) list of BED format files";
   echo -e "   BAMS        (optional) list of BAM files used for coverage matrix";
   echo -e "   PROJECT     (optional) prefix project name, default ${BOLD}PROJECT='temp_'${RESET}";
   echo -e "               It is recommended to add an underscore '_' at the end, but not required.";
   echo -e "   DRYRUN      1 performs a dry-run, prints commands without running ${BOLD}(default)${RESET}";
   echo -e "               0 runs the full script";
   echo -e "   DEBUG       0 does no debug function ${BOLD}(default)${RESET}";
   echo -e "               1 keeps intermediate files without deleting them";
   echo -e "               2 stops before processing";
   echo -e "   PEAKGAP     the distance below which peaks are merged, default ${BOLD}PEAKGAP=100${RESET}";
   echo -e "   PEAKSLOP    the width to expand merged peak edges, default ${BOLD}PEAKSLOP=0${RESET}";
   echo -e "   PEAKMAX     peaks are sliced to this maximum width, default ${BOLD}PEAKMAX=1000${RESET}";
   echo -e "   CHROMSIZES  path to tab-delimited chromosome sizes, default";
   echo -e "               ${BOLD}/ddn/gs1/shared/dirib/reference_genomes/hg19/hg19.chromSizes${RESET}";
   echo -e "   FCTHREADS   number of threads for featureCounts, largely ignored, ${BOLD}FCTHREADS=32${RESET}";
   echo -e "   BEDTOOLS    path to the binary executable file 'bedtools', optional";
   echo -e "   BEDOPS      path to folder that contains bedops binaries, optional";
   echo -e "   DO_FC       1 runs featureCounts if BAMS is defined, default ${BOLD}DO_FC=1${RESET}";
   echo -e "   DO_AUC      1 runs bedtools multicov if BAMS is defined, default ${BOLD}DO_AUC=0${RESET}";
   echo "";
   echo -e "Files are named using file stems defined with this format:";
   echo -e "   PROJECTmergeBed_[num_BED_files]_[run_number]_[date]"
   echo -e "for example:";
   echo -e "   PROJECTmergeBed_4files_num1_30may2018";
   echo -e "If the script is re-run with 4 BED files, it will create a new file stem:"
   echo -e "   PROJECTmergeBed_4files_num2_30may2018";
   echo "";
   echo -e "Files can be defined inside quotes like this:";
   echo -e 'PEAKS="UL3_VEH_25_singfrag.Lee_summits.bed p2w5_VEH_25_singleFragment.Lee_summits.bed p4w4_VEHrep1_25_singleFragment.Lee_summits.bed p4w4_VEHrep2_25_2_singleFragment.Lee_summits.bed" ./mergeSplitCountBed.sh';
   exit;
fi;


## User-defined parameters
PEAKGAP=${PEAKGAP:-"150"};
PEAKSLOP=${PEAKSLOP:-"75"};
PEAKMAX=${PEAKMAX:-"300"};
DRYRUN=${DRYRUN:-"1"};
DEBUG=${DEBUG:-"0"};
DO_FC=${DO_FC:-"1"};
DO_AUC=${DO_AUC:-"0"};

PROJECT=${PROJECT:-"temp_"}
ERROR="0";
## Genome-specific chromsizes
## path to tab-delimited chromosome names and lengths
CHROMSIZES=${CHROMSIZES:-"/ddn/gs1/shared/dirib/reference_genomes/hg19/hg19.chromSizes"};
if [ ! -f "${CHROMSIZES}" ]; then
   echo -e "${FRED}CHROMSIZES file is not available at path:${BOLD}${CHROMSIZES}${RESET}";
   ERROR="1";
fi;

## featureCounts threads, in this case use 32 threads
FCTHREADS=${FCTHREADS:-"32"};

## path to bedtool executable
BEDTOOLS_DEFAULT="/ddn/gs1/biotools/bedtools/bin/bedtools";
BEDTOOLS_FOUND=$(dirname $(which bedtools));
BEDTOOLS_FOUND=${BEDTOOLS_FOUND:-"${BEDTOOLS_DEFAULT}"};
BEDTOOLS=${BEDTOOLS:-"${BEDTOOLS_FOUND}"};
if [ ! -f "${BEDTOOLS}" ]; then
   echo -e "${FRED}bedtools is not available at path:${BOLD}${BEDTOOLS}${RESET}";
   ERROR="1";
fi;

## path to bedops executable
BEDOPS_DEFAULT="/ddn/gs1/biotools/bedops/bin";
BEDOPS_FOUND=$(dirname $(which bedops));
BEDOPS_FOUND=${BEDOPS_FOUND:-"${BEDOPS_DEFAULT}"};
echo ${BEDOPS_FOUND}
BEDOPS=${BEDOPS:-"${BEDOPS_FOUND}"};
if [[ ( ! "." == ".${BEDOPS}" ) ]]; then
   BEDOPS="${BEDOPS}/"
fi;

## current date used in naming output files
DATE1=$(date | perl -p -e 'chomp;@v=split(/[ ]+/);$date=lc("$v[2]$v[1]$v[5]\n");print $date;$_=""')
DATE=${DATE:-${DATE1}}

####################################################################
## Define a couple useful functions
function join_by { local IFS="$1"; shift; echo "$*"; }
strdither() {
   STROUT="";
   NUM=0;
   for STR in ${@}; do
      NUM=$((1-NUM));
      if [ ${NUM} == 1 ]; then
         STROUT="${STROUT} ${FYELLOW}${BOLD}${STR}${RESET}${FYELLOW}";
      else
         STROUT="${STROUT} ${STR}";
      fi;
   done;
   echo "${STROUT}";
}

####################################################################
## Print user-defined parameters
echo -e "${FCYAN}${BOLD}mergeSplitCountBed parameters:${RESET}";
echo -e "${FCYAN}   VERSION:${BOLD}${SLICEJAM_VERSION}${RESET}";
echo -e "${FCYAN}   PEAKGAP:${BOLD}${PEAKGAP}${RESET}";
echo -e "${FCYAN}  PEAKSLOP:${BOLD}${PEAKSLOP}${RESET}";
echo -e "${FCYAN}   PEAKMAX:${BOLD}${PEAKMAX}${RESET}";
echo -e "${FCYAN}CHROMSIZES:${BOLD}${CHROMSIZES}${RESET}";
echo -e "${FCYAN}     DO_FC:${BOLD}${DO_FC}${RESET}";
echo -e "${FCYAN}    DO_AUC:${BOLD}${DO_AUC}${RESET}";

####################################################################
## list the BED peak files being used
INCT=${#PEAKFILES[@]};
if [[ ${INCT} -gt 1 ]]; then
   INCOUNT="${INCT}files";
else
   INCOUNT="${INCT}file";
fi;
echo -e "${BOLD}Peak files (${INCOUNT}):${RESET}";
for IN in ${PEAKFILES[*]}; do
   echo -e "   ${FCYAN}${IN}${RESET}";
done;

####################################################################
## list the BED peak files being used
BAMCT=${#BAMFILES[@]};
if [[ ${BAMCT} -gt 0 ]]; then
   echo -e "${BOLD}BAM files:${RESET}";
   for BAM in ${BAMFILES[*]}; do
      echo -e "   ${FCYAN}${BAM}${RESET}";
   done;
fi;

####################################################################
## Exit if ERROR="1"
if [[ ( "1" == "${ERROR}" ) ]]; then
   echo -e "${FRED}Exiting due to errors.${RESET}";
   exit;
fi;

####################################################################
## Make sure BAM indices are present
if [[ ${BAMCT} -gt 0 ]]; then
   echo -e "${FMAGENTA}==================================================${RESET}";
   echo -e "${BOLD}Creating BAM index files if needed${RESET}";
   for BAM in ${BAMFILES[*]}; do
      BAI=${BAM/.bam/.bam.bai};
      if [ ! -f "${BAI}" ]; then
         echo -e "   ${FCYAN}samtools index ${BOLD}${BAM}${RESET}";
         if [ ! "1" == "${DRYRUN}" ]; then
            samtools index ${BAM};
         else
            echo -e "      ${FRED}DRYRUN${RESET}";
         fi;
      fi;
   done;
fi;

####################################################################
## First define a file name to use for the new files
TEMPNUM=1;
TEMPBED="${PROJECT}mergeBed_${INCOUNT}_gap${PEAKGAP}_slop${PEAKSLOP}_max${PEAKMAX}_num${TEMPNUM}_${DATE}.bed";
#TEMPBED="${PROJECT}mergeBed_${INCOUNT}_num${TEMPNUM}_${DATE}.bed";
echo -e "${FMAGENTA}==================================================${RESET}";
echo -e "${BOLD}Creating temporary BED file using ${FYELLOW}num${TEMPNUM} ${FCYAN}${TEMPBED}${RESET}";
while [ -f ${TEMPBED} ]
do
   TEMPNUM=$((TEMPNUM+1));
   echo -e "      ${FYELLOW}File already exists, trying ${BOLD}num${TEMPNUM}${RESET}";
   TEMPBED="${PROJECT}mergeBed_${INCOUNT}_gap${PEAKGAP}_slop${PEAKSLOP}_max${PEAKMAX}_num${TEMPNUM}_${DATE}.bed";
   #TEMPBED="${PROJECT}mergeBed_${INCOUNT}_num${TEMPNUM}_${DATE}.bed";
done;

####################################################################
## if DEBUG=2 then exit here
if [ "2" == "${DEBUG}" ]; then
   exit;
fi;

####################################################################
## First combine all BED files into one temp file
## Note this command uses "sort -k1,1V" where the "V" does version sort,
## ordering chromosomes chr1, chr2, ..., chr10, chr11, etc.
## Not all sort offer "V" in which case it should be omitted.
echo -e "${FMAGENTA}==================================================${RESET}";
echo -e "${BOLD}Combining BED files${RESET}:";
#echo -e "${FYELLOW}cat $(strdither ${PEAKFILES[*]}) | sort -k1,1V -k2,2n > ${BOLD}${TEMPBED}${RESET}";
echo -e "${FYELLOW}cat $(strdither ${PEAKFILES[*]}) | sort -k1,1 -k2,2n > ${BOLD}${TEMPBED}${RESET}";
if [ ! "1" == "${DRYRUN}" ]; then
   #cat ${PEAKFILES[*]} | grep -v ^track | sort -k1,1V -k2,2n > ${TEMPBED};
   cat ${PEAKFILES[*]} | grep -v ^track | sort -k1,1 -k2,2n > ${TEMPBED};
else
   echo -e "   ${FRED}DRYRUN${RESET}";
fi;


####################################################################
## Merge peaks
MERGEBED=${TEMPBED/.bed/.mergeBed};
MERGECMD="${BEDTOOLS} merge -d ${PEAKGAP} -i ${TEMPBED} | bedtools slop -b ${PEAKSLOP} -g ${CHROMSIZES} -i -";
echo -e "${FMAGENTA}==================================================${RESET}";
echo -e "${BOLD}Merging peaks${RESET}:";
echo -e "${FYELLOW}${MERGECMD} > ${MERGEBED}${RESET}";
if [ ! "1" == "${DRYRUN}" ]; then
   ${BEDTOOLS} merge -d ${PEAKGAP} -i ${TEMPBED} | bedtools slop -b ${PEAKSLOP} -g ${CHROMSIZES} -i - > ${MERGEBED};
else
   echo -e "   ${FRED}DRYRUN${RESET}";
fi;


####################################################################
## Split peaks to a maximum size
SPLITBED=${TEMPBED/.bed/.splitBed};
SPLITCMD="${BEDOPS}bedops --chop ${PEAKMAX} ${MERGEBED}";
echo -e "${FMAGENTA}==================================================${RESET}";
echo -e "${BOLD}Splitting peaks${RESET}:";
echo -e "${FYELLOW}${SPLITCMD} > ${SPLITBED}${RESET}";
if [ ! "1" == "${DRYRUN}" ]; then
   ${SPLITCMD} > ${SPLITBED};
else
   echo -e "   ${FRED}DRYRUN${RESET}";
fi;


####################################################################
## Multi-intersect using bedops
MULTIBED=${TEMPBED/.bed/.multiBed};
export TAB=`echo -e '\011'`;
echo -e "${FMAGENTA}==================================================${RESET}";
echo -e "${BOLD}Counting peak overlaps${RESET}:";
for NUM in $(seq 1 ${INCT}); do
   FILENUM=$((NUM-1));
   MULTINUMBED=${TEMPBED/.bed/.bedops_${NUM}};
   #MULTICMD="grep -v ^track ${PEAKFILES[${FILENUM}]} | ${BEDOPS}bedmap --echo --count ${MERGEBED} -";
   MULTICMD="grep -v ^track ${PEAKFILES[${FILENUM}]} | ${BEDOPS}sort-bed - | ${BEDOPS}bedmap --echo --count <(${BEDOPS}sort-bed ${SPLITBED}) -";
   echo -e "   ${BOLD}Counting overlaps${RESET}:";
   echo -e "   ${FYELLOW}${MULTICMD} > ${MULTINUMBED}${RESET}"
   if [ ! "1" == "${DRYRUN}" ]; then
      grep -v ^track ${PEAKFILES[${FILENUM}]} | ${BEDOPS}sort-bed - | ${BEDOPS}bedmap --echo --count --delim "${TAB}" <(${BEDOPS}sort-bed ${SPLITBED}) - > ${MULTINUMBED};
   else
      echo -e "      ${FRED}DRYRUN${RESET}";
   fi;
done;
MULTINUMBEDALL=${TEMPBED/.bed/.bedops_all};
MULTINUMBEDTMP=${TEMPBED/.bed/.bedops_all_tmp};
MULTINUMBED1=${TEMPBED/.bed/.bedops_1};
echo -e "${FMAGENTA}==================================================${RESET}";
echo -e "${BOLD}Pasting each set of counts to one file${RESET}";
if [ ! "1" == "${DRYRUN}" ]; then
   if [ ! "1" == "${DEBUG}" ]; then
      mv ${MULTINUMBED1} ${MULTIBED};
   else
      cp ${MULTINUMBED1} ${MULTIBED};
   fi;
   HEADER=$(join_by "${TAB}" chrom start end ${PEAKFILES[*]});
fi;
for NUM in $(seq 2 ${INCT}); do
   MULTINUMBED=${TEMPBED/.bed/.bedops_${NUM}};
   echo -e "${FYELLOW}cat ${MULTINUMBED} | cut -f4 | paste ${MULTIBED} - > ${MULTINUMBEDTMP}${RESET}";
   echo -e "${FYELLOW}mv ${MULTINUMBEDTMP} ${MULTIBED}${RESET}";
   if [ ! "1" == "${DRYRUN}" ]; then
      cat ${MULTINUMBED} | cut -f4 | paste ${MULTIBED} - > ${MULTINUMBEDTMP};
      if [ ! "1" == "${DEBUG}" ]; then
         echo -e "   ${FCYAN}rm ${MULTINUMBED}${RESET}";
         rm ${MULTINUMBED};
      fi;
      mv ${MULTINUMBEDTMP} ${MULTIBED};
   else
      echo -e "   ${FRED}DRYRUN${RESET}";
   fi;
done;
if [ ! "1" == "${DRYRUN}" ]; then
   sed -i 1i"${HEADER}" ${MULTIBED};
fi;


####################################################################
## area under the curve for each merged-split peak
## featureCounts took about 9 minutes
## for 4 BAM files and 165,000 features
if [[ ( "1" == "${DO_FC}" ) ]]; then
if [[ ${BAMCT} -gt 0 ]]; then
   ## First create SAF format
   SPLITSAF=${TEMPBED/.bed/.splitSaf};
   SPLITSAFCMD="cat ${SPLITBED} | bed2saf.pl > ${SPLITSAF}";
   echo -e "${FMAGENTA}==================================================${RESET}";
   echo -e "${BOLD}Running bed2saf.pl${RESET}:";
   echo -e "${FYELLOW}${SPLITSAFCMD}${RESET}";
   if [ ! "1" == "${DRYRUN}" ]; then
      cat ${SPLITBED} | \
      perl -p -e '
         while(<>) {
            chomp;
            @v=split(/\t/);
            if ($v[0] =~ /^#/ || $v[0] =~ /^\s*$/) {
               next;
            }
            if (@v > 4) {
               $strand = $v[4];
            } else {
               $strand = "*";
            }
            if (@v > 3) {
               print "$v[3]\t$v[0]\t$v[1]\t$v[2]\t$strand\n";
            } else {
               $peakNum = $peakNum + 1;
               print "peak_$peakNum\t$v[0]\t$v[1]\t$v[2]\t$strand\n";
            }
         }' > ${SPLITSAF};
   else
      echo -e "   ${FRED}DRYRUN${RESET}";
   fi;

   ## Run featureCounts
   FCBED=${TEMPBED/.bed/.fc};
   FCBEDSTDOUT=${TEMPBED/.bed/.fc.stdout};
   FCBEDSTDERR=${TEMPBED/.bed/.fc.stderr};
   FCCMD="featureCounts -T ${FCTHREADS} -F SAF -M -s 0 --donotsort -f -O -a ${SPLITBED} -o ${FCBED}";
   echo -e "${FMAGENTA}==================================================${RESET}";
   echo -e "${BOLD}Running featureCounts${RESET}:";
   echo -e "${FMAGENTA}${FCCMD} ${FYELLOW}$(strdither ${BAMFILES[*]})${RESET}";
   if [ ! "1" == "${DRYRUN}" ]; then
      featureCounts -T ${FCTHREADS} -F SAF -M -s 0 --donotsort -f -O -a ${SPLITSAF} -o ${FCBED} \
      ${BAMFILES[*]} \
      1>${FCBEDSTDOUT} 2>${FCBEDSTDERR};
   else
      echo -e "   ${FRED}DRYRUN${RESET}";
   fi;
fi;
fi;


####################################################################
## area under the curve for each merged-split peak
## multicov took about 45 minutes
if [[ ( "1" == "${DO_AUC}" ) ]]; then
if [[ ${BAMCT} -gt 0 ]]; then
   MULTICOVBED=${TEMPBED/.bed/.multiCovBed};
   #MULTICOVCMD="bedtools multicov -bed ${MERGEBED} -bams ${BAMFILES[*]}";
   MULTICOVCMD="bedtools multicov -bed ${SPLITBED} -bams ${BAMFILES[*]}";
   echo -e "${FMAGENTA}==================================================${RESET}";
   echo -e "${BOLD}Running bedtools multicov${RESET}:";
   echo -e "${FYELLOW}${MULTICOVCMD} > ${BOLD}${MULTICOVBED}${RESET}";
   if [ ! "1" == "${DRYRUN}" ]; then
      bedtools multicov -bed ${SPLITBED} -bams ${BAMFILES[*]} > ${MULTICOVBED};
   else
      echo -e "   ${FRED}DRYRUN${RESET}";
   fi;
fi;
fi;
