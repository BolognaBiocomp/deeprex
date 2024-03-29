#!/bin/bash

export DEEPREX_ROOT=$(realpath ..)
export PATH=${DEEPREX_ROOT}:${PATH}

bts_list=blind_test_set/blind_test_set.lst

tmpworkdir=$(mktemp --tmpdir=. -d)

output_file=bts_deeprex.tsv
echo -e "PROTEIN\tPOSITION\tRES\tPRED_EXP\tRI\tHYDRO\tCONS\tTRUE_EXP\tRSA" > ${output_file}
echo "Benchmarking DeepREx on blind test set (it may take hours)"
for pdbid in $(cat ${bts_list}); do
  fas_tmp_file=$(mktemp --tmpdir=${tmpworkdir} --suffix=.fa)
  a3m_tmp_file=$(mktemp --tmpdir=${tmpworkdir} --suffix=.a3m)
  hhm_tmp_file=$(mktemp --tmpdir=${tmpworkdir} --suffix=.hhm)
  out_tmp_file=$(mktemp --tmpdir=${tmpworkdir} --suffix=.tsv)
  rsa_tmp_file=$(mktemp --tmpdir=${tmpworkdir} --suffix=.rsa)
  gunzip -c blind_test_set/fasta/${pdbid}.fasta.gz > ${fas_tmp_file}
  gunzip -c blind_test_set/aln/${pdbid}.a3m.gz > ${a3m_tmp_file}
  gunzip -c blind_test_set/aln/${pdbid}.hhm.gz > ${hhm_tmp_file}
  gunzip -c blind_test_set/rsa/${pdbid}.rsa.gz > ${rsa_tmp_file}
  deeprex.py aln -f ${fas_tmp_file} -a ${a3m_tmp_file} -m ${hhm_tmp_file} -o ${out_tmp_file} &> /dev/null
  paste ${out_tmp_file} ${rsa_tmp_file} | cut -f 1-7,12,13 >> ${output_file}
  rm ${fas_tmp_file} ${a3m_tmp_file} ${hhm_tmp_file} ${out_tmp_file} ${rsa_tmp_file}
done
echo "DeepREx performance on blind test set"
python eval_performance.py ${output_file}
rm -rf ${tmpworkdir}
