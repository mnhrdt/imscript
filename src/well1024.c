/* ************************************************************************* */
/* Copyright:  Francois Panneton and Pierre L'Ecuyer, University of Montreal */
/*             Makoto Matsumoto, Hiroshima University                        */
/* Notice:     This code can be used freely for personal, academic,          */
/*             or non-commercial purposes. For commercial purposes,          */
/*             please contact P. L'Ecuyer at: lecuyer@iro.UMontreal.ca       */
/* ************************************************************************* */

// Note: some variables moved and renamed by Enric Meinhardt-Llopis
// (this modification is probably illegal, please do not use or distribute it)
//

#define MAT0POS(t,v) (v^(v>>t))
#define MAT0NEG(t,v) (v^(v<<(-(t))))
#define Identity(v) (v)

#define V0            well1024_STATE[ well1024_state_i                   ]
#define VM1           well1024_STATE[(well1024_state_i+3) & 0x0000001fU]
#define VM2           well1024_STATE[(well1024_state_i+24) & 0x0000001fU]
#define VM3           well1024_STATE[(well1024_state_i+10) & 0x0000001fU]
#define VRm1          well1024_STATE[(well1024_state_i+31) & 0x0000001fU]
#define newV0         well1024_STATE[(well1024_state_i+31) & 0x0000001fU]
#define newV1         well1024_STATE[ well1024_state_i                   ]

static unsigned int well1024_state_i = 0;
static unsigned int well1024_STATE[32] = {0};
static unsigned int well1024_z0;
static unsigned int well1024_z1;
static unsigned int well1024_z2;

void well1024_init (unsigned int *init)
{
	well1024_state_i = 0;
	for (int j = 0; j < 32; j++)
		well1024_STATE[j] = init[j];
}

void well1024_seed(unsigned int seed)
{
	// silly seed function
	unsigned int t[32];
	for (int i = 0; i < 32; i++)
		t[i] = 0xdeadbeef ^ (0xc0ffee + i*i*(seed-1));
	well1024_init(t);
}

double well1024 (void)
{
	well1024_z0 = VRm1;
	well1024_z1 = Identity(V0)      ^ MAT0POS(8, VM1);
	well1024_z2 = MAT0NEG(-19, VM2) ^ MAT0NEG(-14,VM3);
	newV1       = well1024_z1       ^ well1024_z2;
	newV0       = MAT0NEG(-11,well1024_z0)
	            ^ MAT0NEG(-7, well1024_z1)
	            ^ MAT0NEG(-13,well1024_z2);
	well1024_state_i = (well1024_state_i + 31) & 0x0000001fU;
	return 2.32830643653869628906e-10 * well1024_STATE[well1024_state_i];
}
#undef MAT0POS
#undef MAT0NEG
#undef Identity
#undef V0
#undef VM1
#undef VM2
#undef VM3
#undef VRm1
#undef newV0
#undef newV1
