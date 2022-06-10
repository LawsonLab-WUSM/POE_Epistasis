#!/bin/bash
#SBATCH --mem=8G
#SBATCH -n 1
#SBATCH -N 1

#[Number of ASE Genes] [Number of Expressed Genes] [ASE file] [General File]

Wa=10
Wg=1000

for Ai in `seq 1 $Wa $1`
	do
		for Gi in `seq 1 $Wg $2`
			do
				Af=$(($Ai+$Wa-1))
				Gf=$(($Gi+$Wg-1))
				

				if (($1 < $Af))
					 then
						Af=$1
				fi

				if (($2 < $Gf))
					 then
						Gf=$2
				fi

				echo "ASE: $Ai $Af Gen: $Gi $Gf"
				sbatch run_GeneCorrelation.sh $3 $4 $Ai $Af $Gi $Gf
			done
done

