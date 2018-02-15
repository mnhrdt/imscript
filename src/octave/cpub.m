function cpub(x)
	f = [ tempname()  ".tiff" ];
	c = [ "cpu "  f  " 2>/dev/null &" ];
	r = [ "rm "  f ];
	iio_write(f, x);
	system(c);
	system(r);
endfunction
