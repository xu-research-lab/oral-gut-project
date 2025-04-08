#./SparCC.py feces_nonPMA.tsv -c feces_nonPMA_cor_spacc_out.txt >feces_nonPMA_sparcc.log
#./MakeBootstraps.py feces_nonPMA.tsv -n 100 t_bootstrap_#.txt -p feces_nonPMA_pvals/ >>feces_nonPMA_sparcc.log
for n in {0..99};do ./SparCC.py feces_nonPMA_pvals/feces_nonPMA.tsv.permuted_${n}.txt -i 20 -c feces_nonPMA_bootstrap_cor_${n}.txt>>feces_nonPMA_sparcc.log;done
./PseudoPvalues.py feces_nonPMA_cor_spacc_out.txt feces_nonPMA_pvals/feces_nonPMA_bootstrap_cor_#.txt -o feces_nonPMA_pvals/pvals.two_sided.txt -t two_sided >>feces_nonPMA_sparcc.log
#./SparCC.py feces_PMA.tsv -c feces_PMA_cor_spacc_out.txt >feces_PMA_sparcc.log
#./MakeBootstraps.py feces_PMA.tsv -n 100 t_bootstrap_#.txt -p feces_PMA_pvals/ >>feces_PMA_sparcc.log
for n in {0..99};do ./SparCC.py feces_PMA_pvals/feces_nonPMA.tsv.permuted_${n}.txt -i 20 -c feces_PMA_bootstrap_cor_${n}.txt>>feces_PMA_sparcc.log;done
./PseudoPvalues.py feces_PMA_cor_spacc_out.txt feces_PMA_pvals/feces_PMA_bootstrap_cor_#.txt -o feces_PMA_pvals/pvals.two_sided.txt -t two_sided >>feces_PMA_sparcc.log
echo "first"
#./SparCC.py saliva_nonPMA.tsv -c saliva_nonPMA_cor_spacc_out.txt >saliva_nonPMA_sparcc.log
#./MakeBootstraps.py saliva_nonPMA.tsv -n 100 t_bootstrap_#.txt -p saliva_nonPMA_pvals/ >>saliva_nonPMA_sparcc.log
for n in {0..99};do ./SparCC.py saliva_nonPMA_pvals/feces_nonPMA.tsv.permuted_${n}.txt -i 20 -c saliva_nonPMA_bootstrap_cor_${n}.txt>>saliva_nonPMA_sparcc.log;done
./PseudoPvalues.py saliva_nonPMA_cor_spacc_out.txt saliva_nonPMA_pvals/saliva_nonPMA_bootstrap_cor_#.txt -o saliva_nonPMA_pvals/pvals.two_sided.txt -t two_sided >>saliva_nonPMA_sparcc.log

echo "second"
