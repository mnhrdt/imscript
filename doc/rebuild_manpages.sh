help2man --help >/dev/null && for i in `ls ../src/*.c|grep -v help_stuff`; do
	X=`basename $i .c`
	grep -q string_oneliner $i && ../bin/$X --manraw > man/$X.1
done
true
