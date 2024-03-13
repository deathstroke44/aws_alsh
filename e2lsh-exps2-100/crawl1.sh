cd ../build

src=crawl
id=${src}1
dataset=/data/kabir/similarity-search/dataset/${src}/base.fvecs
query=/data/kabir/similarity-search/dataset/${src}/query.fvecs
groundtruth=/data/kabir/similarity-search/dataset/${src}/groundtruth.ivecs
N=1989995
DD=300
QN=10000
DATA_HASH=${id}_data_hash_e2lsh.dh
QUERY_HASH=${id}_query_hash_e2lsh.dh
OUTPUT=${id}_output_e2lsh.out

./alsh -A fraction_recall_e2kl -n ${N} -q ${QN} -d ${DD} -D ${dataset} -Q ${query} -W ../data/Mnist784/Mnist784_normal.w -G ${groundtruth} -O ${OUTPUT} -U 3.140000 --data_hash_filename ${DATA_HASH} --query_hash_filename ${QUERY_HASH}