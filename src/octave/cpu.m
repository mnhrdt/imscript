function cpu(x)
	f = [ tempname()  ".tiff" ];
	c = [ "cpu "   f   " 2>/dev/null ; rm "   f];
	iio_write(f, x);
	system(c);
endfunction
