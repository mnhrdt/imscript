for i in {00..99}; do
	for j in r g b; do
		SEED=$[i+1] ./bin/genk 10 8 9 2|fftshift - /tmp/k_$i.pgm
		qnm convoluciona /tmp/lena_${j}.pgm /tmp/k_$i.pgm - | qeasy 0 2200 - /tmp/lspeck_${i}_${j}.png
	done
	plambda /tmp/lspeck_${i}_{r,g,b}.png "x y z join3" |downsa v 2|qeasy 0 255 - /tmp/dlspeck_$i.png
done
