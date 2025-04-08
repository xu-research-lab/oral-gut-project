phylophlan_metagenomic \
    -i test \
    -o test_output \
--nproc 22 \
--database_folder  /home/huangxc/db/SGB.Jul20 \
-n 1 \
    -d SGB.Jul20 \
    --verbose 2>&1 | tee logs/phylophlan_metagenomic.log