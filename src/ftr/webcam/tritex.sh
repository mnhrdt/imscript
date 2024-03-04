#plambda zero:1512x1512 "randg randg randg rgb"|fft|plambda ':R 1 ^ /'|ifft|GETPIXEL=periodic plambda 'x[0] x[0](100,0) - x[1] x[1](50,87) - x[2] x[2](-50,87) - rgb rgb2xyz'|qauto -p 0.1 - /tmp/tritex.png
plambda zero:1512x1512 randg|fft|plambda ':R 1 ^ /'|ifft|GETPIXEL=periodic plambda 'x 3 * x(100,0) x(50,87) x(-50,87) - - -'|qauto -p 0.1 - /tmp/tritex2.png
