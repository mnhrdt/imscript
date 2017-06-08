#include <stdio.h>     /* for printf */
#include <stdlib.h>    /* for exit */
#include <unistd.h>    /* for getopt */
int main_bad (int argc, char **argv) {
    int c;
    int digit_optind = 0;
    int aopt = 0, bopt = 0;
    char *copt = 0, *dopt = 0;
    while ( (c = getopt(argc, argv, "abc:d:012")) != -1) {
        int this_option_optind = optind ? optind : 1;
        switch (c) {
        case '0':
        case '1':
        case '2':
            if (digit_optind != 0 && digit_optind != this_option_optind)
              printf ("digits occur in two different argv-elements.\n");
            digit_optind = this_option_optind;
            printf ("option %c\n", c);
            break;
        case 'a':
            printf ("option a\n");
            aopt = 1;
            break;
        case 'b':
            printf ("option b\n");
            bopt = 1;
            break;
        case 'c':
            printf ("option c with value '%s'\n", optarg);
            copt = optarg;
            break;
        case 'd':
            printf ("option d with value '%s'\n", optarg);
            dopt = optarg;
            break;
        case '?':
            break;
        default:
            printf ("?? getopt returned character code 0%o ??\n", c);
        }
    }
    if (optind < argc) {
        printf ("non-option ARGV-elements: ");
        while (optind < argc)
            printf ("%s ", argv[optind++]);
        printf ("\n");
    }
    exit (0);
}

int main_good(int argc, char **argv)
{
	int aopt = (int)pick_option(&argc, &argv, "a", NULL);
	int bopt = (int)pick_option(&argc, &argv, "b", NULL);
	char *copt = pick_optionp(&argc, &argv, "c", "");
	char *dopt = pick_optionp(&argc, &argv, "d", "");
}
