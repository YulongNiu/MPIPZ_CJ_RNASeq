
## originally by Yulong Niu
## yulong.niu@hotmail.com

date

RAW_PATH=/netscratch/dep_psl/grp_rgo/yniu/CJFe/raw_data_1stadd
CLEAN_PATH=/netscratch/dep_psl/grp_rgo/yniu/CJFe/clean_data_1stadd

FASTP_PATH=/home/yniu/Biotools

CORENUM=16

cd ${RAW_PATH}

## sample names
fq=($(ls | grep fq.gz))
fqnames=($(echo "${fq[@]%_*}" | tr ' ' '\n' | sort -u | tr '\n' ' '))

for i in ${fqnames[@]}
do
    echo "Trimming ${i}_R1.fq.gz ${i}_R2.fq.gz."

    ${FASTP_PATH}/fastp -w ${CORENUM} \
                 -z 6 \
                 -p -c \
                 -h ${i}.html \
                 -i ${i}_R1.fq.gz -I ${i}_R2.fq.gz \
                 -o ${CLEAN_PATH}/${i}_R1.fq.gz -O ${CLEAN_PATH}/${i}_R2.fq.gz
done

date
