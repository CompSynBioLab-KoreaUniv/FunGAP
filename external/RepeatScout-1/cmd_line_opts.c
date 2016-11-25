#include "cmd_line_opts.h"

#include <stdlib.h>
#include <string.h>

int co_get_int(int argc, char** argv, const char* text, int* res) {
	int opt_len = strlen(text);
	int found = 0;

	int ii;
	for(ii = 0; ii < argc; ++ii) {
		if( argv[ii][0] == '-' ) {
			if(!strncmp(text, argv[ii], opt_len)) {
				*res = atoi(argv[ii+1]);
				found = 1;
				break;
			}

			else {
				++ii;
			}
		}
	}

	return found;
}

int co_get_bool(int argc, char** argv, const char* text, int* res) {
	int opt_len = strlen(text);
	int found = 0;

	int ii;
	for(ii = 0; ii < argc; ++ii) {
		if( argv[ii][0] == '-' ) {
			if(!strncmp(text, argv[ii], opt_len)) {
				*res = 1;
				found = 1;
				break;
			}

			else {
				++ii;
			}
		}
	}

	return found;
}

int co_get_float(int argc, char** argv, const char* text, float* res) {
	int opt_len = strlen(text);
	int found = 0;

	int ii;
	for(ii = 0; ii < argc; ++ii) {
		if( argv[ii][0] == '-' ) {
			if(!strncmp(text, argv[ii], opt_len)) {
				*res = atof(argv[ii+1]);
				found = 1;
				break;
			}

			else {
				++ii;
			}
		}
	}

	return found;
}

int co_get_string(int argc, char** argv, const char* text, char* * res) {
	int opt_len = strlen(text);
	int found = 0;

	int ii;
	for(ii = 0; ii < argc; ++ii) {
		if( argv[ii][0] == '-' ) {
			if(!strncmp(text, argv[ii], opt_len)) {
				*res = argv[ii+1];
				found = 1;
				break;
			}

			else {
				++ii;
			}
		}
	}

	return found;
}
