// derivative of a planar homography

void homaff(float A[6], float H[9], float x0[2])
{
	//        / a b p \        / A B P |
	//    H = | c d q |    A = \ C D Q |
	//        \ r s t /
	//
	float a = H[0]; float b = H[1]; float p = H[2];
	float c = H[3]; float d = H[4]; float q = H[5];
	float r = H[6]; float s = H[7]; float t = H[8];
	float x = x0[0]; float y = x0[1];
	float n = r*x + s*x + t;
	float A = ( (a*s - b*r)*y + (a*t - p*r) ) / (n*n);
	float B = ( (b*r - a*s)*x + (b*t - p*s) ) / (n*n);
	float C = ( (c*s - d*r)*y + (c*t - q*r) ) / (n*n);
	float D = ( (d*r - c*s)*x + (d*t - q*s) ) / (n*n);
	float P = (a*x + b*y + p)/n - (A*x + B*x);
	float Q = (c*x + d*y + q)/n - (C*x + D*x);
	A[0] = A; A[1] = B; A[2] = P;
	A[3] = C; A[4] = D; A[5] = Q;
}
