for d in trees/*; do
    if [[ -d $d ]]; then
        pushd $d
        refAcc=$(grep ^refseq_acc config.toml | cut -d \' -f 2)
        if [[ ! -e $refAcc.gbff ]]; then
            time python ../../scripts/get_gbff.py $refAcc angie@soe.ucsc.edu
        fi
        if [[ ! -e $refAcc.fa ]]; then
            time python ../../scripts/gbff_to_fasta.py $refAcc.gbff $refAcc.fa
        fi
        popd
    fi
done
