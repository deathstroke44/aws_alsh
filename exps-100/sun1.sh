cd ../build

src=astro1m
id=${src}1
dataset=/data/kabir/similarity-search/dataset/${src}/base.fvecs
query=/data/kabir/similarity-search/dataset/${src}/query.fvecs
groundtruth=/data/kabir/similarity-search/dataset/${src}/groundtruth.ivecs
N=79106
DD=512
QN=200
DATA_HASH=${id}_data_hash.dh
QUERY_HASH=${id}_query_hash.dh
OUTPUT=${id}_output.out

./alsh -A fraction_recall_s2alsh -n ${N} -q ${QN} -d ${DD} -D ${dataset} -Q ${query} -W ../data/Mnist784/Mnist784_normal.w -G ${groundtruth} -O ${OUTPUT} -U 3.140000 --data_hash_filename ${DATA_HASH} --query_hash_filename ${QUERY_HASH}