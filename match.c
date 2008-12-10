#include <stdlib.h>
#include <libgen.h>
#include <fnmatch.h>

int match_ (const char *str, const char *pattern) {
  /*  printf("str %s\npattern %s\n",str,pattern); */
  /*
#ifdef 1
  */

  return !fnmatch(pattern,str,0);
  /*#else 
    return gmatch (str,pattern);
#endif */
}

