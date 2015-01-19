#!/usr/bin/octave -qf

format loose

if nargin != 3
	fprintf(stderr, "usage:\n\tbasify p q r > P.txt\n");
	exit(1);
endif

% read input vector
p = str2double(argv())';

% define canonical basis
e1 = [1 0 0];
e2 = [0 1 0];
e3 = [0 0 1];

% normalize input vector
p /= norm(p);

% prevent gimbal lock around the vertical direction (ad-hoc)
if norm(p - e3) < 1e-5
	d = e1;
else
	d = e3;
endif

% build an orthonormal basis {p, q, r}
q = cross(p, d);
q /= norm(q);
r = cross(p, q);


% project the canonical basis
pe1 = e1 - dot(p, e1)*p;
pe2 = e2 - dot(p, e2)*p;
pe3 = e3 - dot(p, e3)*p;

% build matrices with the two bases
mpei = [pe1 ; pe2 ; pe3];
mpqr = [p ; q ; r];

% compute the coefficients of the projection with this pair of basis
P = mpqr * mpei';

% remove the first row
P = P(2:3,:);

printf("%g %g %g 0\n%g %g %g 0\n", P(1,1),P(1,2),P(1,3), P(2,1),P(2,2),P(2,3));
