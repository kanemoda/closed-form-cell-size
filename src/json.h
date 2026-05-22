#ifndef JSON_H
#define JSON_H

#include <stdint.h>

/*
 * Minimal flat-object JSON parser. Sufficient for our config files:
 * top-level object with string keys and number/string/bool values.
 * Supports `//` line comments as a convenience extension.
 *
 * Not supported (deliberately): nested objects/arrays, scientific
 * notation longer than strtod handles, unicode escapes.
 */

#define JSON_MAX_KEY     64
#define JSON_MAX_STR    256
#define JSON_MAX_ENTRIES 64

typedef enum { JSON_TYPE_NUM, JSON_TYPE_STR, JSON_TYPE_BOOL } json_type_t;

typedef struct {
    char        key[JSON_MAX_KEY];
    json_type_t type;
    double      num;
    int         b;
    char        str[JSON_MAX_STR];
} json_entry_t;

typedef struct {
    json_entry_t entries[JSON_MAX_ENTRIES];
    int          n_entries;
} json_object_t;

int json_parse_file (const char *path, json_object_t *obj);

int json_get_double (const json_object_t *obj, const char *key, double   *out);
int json_get_int    (const json_object_t *obj, const char *key, int      *out);
int json_get_u64    (const json_object_t *obj, const char *key, uint64_t *out);
int json_get_string (const json_object_t *obj, const char *key, char *buf, int buflen);
int json_get_bool   (const json_object_t *obj, const char *key, int      *out);

#endif
