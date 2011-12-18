
//static const char *wisdomfilename = "/home/eml/.wisdom";
static char *wisdomfilename = NULL;
static char wisdomfilenamestatic[1000];

void evoke_wisdom(void)
{
	// XXX: race condition
	if (!wisdomfilename) {
		char *homedir = getenv("HOME");
		if (!homedir)
			return; // no wisdom for you
		wisdomfilename = wisdomfilenamestatic;
		size_t n = strlen(homedir);
		assert(n + 10 < 1000);
		memcpy(wisdomfilename, homedir, n);
		memcpy(wisdomfilename + n, "/.wisdomf", 9);
	}

	FILE *f = fopen(wisdomfilename, "r");
	if (f) {
		fftwf_import_wisdom_from_file(f);
		fclose(f);
	}
	else {
		f = fopen(wisdomfilename, "w");
		if (!f) error("could not create wisdom file");
		fclose(f);
		fprintf(stderr, "created wisdom file \"%s\"\n", wisdomfilename);
	}
}

void bequeath_wisdom(void)
{
	assert(wisdomfilename);
	FILE *f = fopen(wisdomfilename, "w");
	if (f)
	{
		fftwf_export_wisdom_to_file(f);
		fclose(f);
	}
	else
		error("could not bequeath wisdom");
}
