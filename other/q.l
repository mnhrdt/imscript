%option noyywrap

%{
#include <stdio.h>
%}

number [1-9][0-9]*
name [a-zA-Z][a-zA-Z0-9]*


%%

{number}     { printf("a number \"%s\"\n", yytext); }
{name}       { printf("an identifier \"%s\"\n", yytext); }
{name}[{]{number}[}]       { printf("a numbered identifier \"%s\"\n", yytext); }
[ \t\n]+     /* */

%%

int main()
{
	yyin = stdin;
	yylex();
	return 0;
}
