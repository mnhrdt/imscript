#include <string.h> // strcmp
#include <stdio.h>  // puts, snprintf
#include <stdlib.h> // system, exit

static int do_man(int raw)
{
#ifdef __OpenBSD__
#define MANPIPE "|mandoc -a"
#else
#define MANPIPE "|man -l -"
#endif
	char buf[0x200];
	snprintf(buf, 0x200, "help2man -N -S imscript -n \"%s\" %s %s",
			help_string_oneliner, help_string_name,
			raw ? "" : MANPIPE
			);
	return system(buf);
}

static void if_help_is_requested_print_it_and_exit_the_program(char *s)
{
	if (!strcmp(s, "-?"))        exit(0*puts(help_string_oneliner));
	if (!strcmp(s, "-h"))        exit(0*puts(help_string_usage));
	if (!strcmp(s, "--help"))    exit(0*puts(help_string_long));
	if (!strcmp(s, "--version")) exit(0*puts(help_string_version));
	if (!strcmp(s, "--man"))     exit(do_man(0));
	if (!strcmp(s, "--manraw"))  exit(do_man(1));
	if (!strcmp(s, "--help-oneliner")) puts(help_string_oneliner);
}

