for i in `seq 1 256`
do 
    echo "Rscript --vanilla ScriptSensi_M1_jan2021.R $i" > "res/run$i.sh"
    qsub -cwd ./res/run$i.sh
    echo $i
done
