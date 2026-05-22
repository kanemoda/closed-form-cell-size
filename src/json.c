#include "json.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    const char *p;
    const char *end;
} parser_t;

static void skip_ws(parser_t *ps) {
    while (ps->p < ps->end) {
        char c = *ps->p;
        if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
            ps->p++;
        } else if (c == '/' && ps->p + 1 < ps->end && ps->p[1] == '/') {
            while (ps->p < ps->end && *ps->p != '\n') ps->p++;
        } else {
            break;
        }
    }
}

static int parse_string(parser_t *ps, char *buf, int buflen) {
    skip_ws(ps);
    if (ps->p >= ps->end || *ps->p != '"') return -1;
    ps->p++;
    int n = 0;
    while (ps->p < ps->end && *ps->p != '"') {
        if (n >= buflen - 1) return -1;
        char c = *ps->p++;
        if (c == '\\' && ps->p < ps->end) {
            char e = *ps->p++;
            switch (e) {
                case 'n':  c = '\n'; break;
                case 't':  c = '\t'; break;
                case 'r':  c = '\r'; break;
                case '"':  c = '"';  break;
                case '\\': c = '\\'; break;
                default:   c = e;    break;
            }
        }
        buf[n++] = c;
    }
    if (ps->p >= ps->end) return -1;
    ps->p++; /* closing quote */
    buf[n] = '\0';
    return 0;
}

static int parse_number(parser_t *ps, double *out) {
    skip_ws(ps);
    char *endp;
    *out = strtod(ps->p, &endp);
    if (endp == ps->p) return -1;
    ps->p = endp;
    return 0;
}

static int parse_value(parser_t *ps, json_entry_t *e) {
    skip_ws(ps);
    if (ps->p >= ps->end) return -1;
    char c = *ps->p;
    if (c == '"') {
        e->type = JSON_TYPE_STR;
        return parse_string(ps, e->str, JSON_MAX_STR);
    }
    if (c == 't') {
        if (ps->p + 4 <= ps->end && strncmp(ps->p, "true", 4) == 0) {
            e->type = JSON_TYPE_BOOL; e->b = 1; ps->p += 4; return 0;
        }
        return -1;
    }
    if (c == 'f') {
        if (ps->p + 5 <= ps->end && strncmp(ps->p, "false", 5) == 0) {
            e->type = JSON_TYPE_BOOL; e->b = 0; ps->p += 5; return 0;
        }
        return -1;
    }
    e->type = JSON_TYPE_NUM;
    return parse_number(ps, &e->num);
}

int json_parse_file(const char *path, json_object_t *obj) {
    FILE *f = fopen(path, "rb");
    if (!f) return -1;
    if (fseek(f, 0, SEEK_END) != 0) { fclose(f); return -1; }
    long sz = ftell(f);
    if (sz < 0) { fclose(f); return -1; }
    rewind(f);
    char *buf = (char *)malloc((size_t)sz + 1);
    if (!buf) { fclose(f); return -1; }
    size_t got = fread(buf, 1, (size_t)sz, f);
    fclose(f);
    if (got != (size_t)sz) { free(buf); return -1; }
    buf[sz] = '\0';

    parser_t ps = { buf, buf + sz };
    obj->n_entries = 0;
    skip_ws(&ps);
    if (ps.p >= ps.end || *ps.p != '{') { free(buf); return -1; }
    ps.p++;

    while (1) {
        skip_ws(&ps);
        if (ps.p >= ps.end) { free(buf); return -1; }
        if (*ps.p == '}') { ps.p++; break; }

        if (obj->n_entries >= JSON_MAX_ENTRIES) { free(buf); return -1; }
        json_entry_t *e = &obj->entries[obj->n_entries];
        if (parse_string(&ps, e->key, JSON_MAX_KEY) != 0) { free(buf); return -1; }
        skip_ws(&ps);
        if (ps.p >= ps.end || *ps.p != ':') { free(buf); return -1; }
        ps.p++;
        if (parse_value(&ps, e) != 0) { free(buf); return -1; }
        obj->n_entries++;

        skip_ws(&ps);
        if (ps.p < ps.end && *ps.p == ',') { ps.p++; continue; }
        if (ps.p < ps.end && *ps.p == '}') { ps.p++; break; }
        free(buf); return -1;
    }

    free(buf);
    return 0;
}

static const json_entry_t *find(const json_object_t *obj, const char *key) {
    for (int i = 0; i < obj->n_entries; i++) {
        if (strcmp(obj->entries[i].key, key) == 0) return &obj->entries[i];
    }
    return NULL;
}

int json_get_double(const json_object_t *obj, const char *key, double *out) {
    const json_entry_t *e = find(obj, key);
    if (!e || e->type != JSON_TYPE_NUM) return -1;
    *out = e->num;
    return 0;
}

int json_get_int(const json_object_t *obj, const char *key, int *out) {
    const json_entry_t *e = find(obj, key);
    if (!e || e->type != JSON_TYPE_NUM) return -1;
    *out = (int)e->num;
    return 0;
}

int json_get_u64(const json_object_t *obj, const char *key, uint64_t *out) {
    const json_entry_t *e = find(obj, key);
    if (!e || e->type != JSON_TYPE_NUM) return -1;
    /* Doubles exactly represent integers up to 2^53; seeds above that are
     * not faithfully round-tripped through this minimal parser. */
    *out = (uint64_t)e->num;
    return 0;
}

int json_get_string(const json_object_t *obj, const char *key, char *buf, int buflen) {
    const json_entry_t *e = find(obj, key);
    if (!e || e->type != JSON_TYPE_STR) return -1;
    int n = (int)strlen(e->str);
    if (n >= buflen) n = buflen - 1;
    memcpy(buf, e->str, n);
    buf[n] = '\0';
    return 0;
}

int json_get_bool(const json_object_t *obj, const char *key, int *out) {
    const json_entry_t *e = find(obj, key);
    if (!e || e->type != JSON_TYPE_BOOL) return -1;
    *out = e->b;
    return 0;
}
