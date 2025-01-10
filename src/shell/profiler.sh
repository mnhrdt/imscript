for i in `seq 11 -1 1`; do
for j in {0..9}; do
	make clean
	/usr/bin/time -f "%e" make full -j $i >/dev/null 2>/tmp/oe_${i}_${j}.txt
done
done

for i in {1..11}; do
	cat /tmp/oe_${i}_*.txt|xargs -L 1 echo $i
done | gnuplot -e 'set term png;plot [0:12] [0:150] "-"' >ref_`hostname`_.png
