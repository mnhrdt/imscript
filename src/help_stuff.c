

static int print_help(char *v, int verbosity)
{
	if (verbosity == 0) printf("%s\n", help_string_oneliner);
	if (verbosity == 1) puts(help_string_usage);
	if (verbosity == 2) puts(help_string_long);
	return 0;
}

static int print_version(void)
{
	puts(help_string_version);
	return 0;
}

static int do_man(void)
{
#ifdef __OpenBSD__
#define MANPIPE "|mandoc -a"
#else
#define MANPIPE "|man -l -"
#endif
	char buf[0x200];
	snprintf(buf, 0x200, "help2man -N -S imscript -n \"%s\" %s" MANPIPE,
			help_string_oneliner, help_string_name);
	return system(buf);
}


