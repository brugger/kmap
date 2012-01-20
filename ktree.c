
#include "kmap.h" 



struct node * create_node() {
  struct node * new = malloc(sizeof(struct node));
  new->bases  = '\0';
  new->np[0]  = 0;
  new->np[1]  = 0;
  new->np[2]  = 0;
  new->np[3]  = 0;
  new->count  = 0;

  return new;
}
  
