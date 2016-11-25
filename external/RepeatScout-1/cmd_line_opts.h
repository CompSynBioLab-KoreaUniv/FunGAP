#ifndef __CMD_LINE_OPTS_H__
#define __CMD_LiNE_OPTS_H__

int co_get_int(int argc, char** argv, const char* text, int* );
int co_get_bool(int argc, char** argv, const char* text, int* );
int co_get_float(int argc, char** argv, const char* text, float* );
int co_get_string(int argc, char** argv, const char* text, char* * );

#endif
