help2man --help >/dev/null && for i in `ls ../src/*.c|grep -v help_stuff`; do
	X=`basename $i .c`
	grep -q string_oneliner $i && ../bin/$X --manraw |\
		sed 's/by help2man 1.*/by help2man/' > man/$X.1
done
true
