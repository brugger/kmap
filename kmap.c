

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <strings.h>
#include <malloc.h>

#include "kmap.h"

struct node {
  char *bases;
  struct node *np[4];
  
  int count;
};

struct node *root;

struct fasta {
  char * name, * seq;
  int length;
};

int kmer_length = 36;

struct fasta * FastaIn();

int main(int argc, char **argv ) {
  root = malloc(sizeof(struct node));

  printf("Counting kmers, for fun an profit\n");


  struct fasta *fs;

  while (fs = FastaIn() ) {
    if (! fs ) 
      break;

    printf("%s with %d bases\n", fs->name, fs->length);
    split_and_build(fs->seq, fs->length);
    break;
  }

  struct node *nn = root->np[1];
  printf("%s\n", nn->bases);

  return 0;
}


split_and_build(char* seq, int length ) {

  int i;


  printf("%s\n", seq);

  for(i=0;i <= length - kmer_length; i++) {
    char *kmer = malloc(sizeof(char)*(kmer_length + 1));
    strncpy(kmer, &seq[i], kmer_length);
    kmer[kmer_length] = 0;

    
    //    printf("[%s]\n", kmer);
    if (! ATGC_only(kmer, kmer_length)) {
      free(kmer);
      continue;
    }


    add2tree(root, kmer, kmer_length);
  }
}


inline int base2pos(const char b) {
  switch (b) {
    case 'A':
    case 'a': return 0;

    case 'C':
    case 'c': return 1;

    case 'G':
    case 'g':  return 2;

    case 'T':
    case 't':  return 3;

    }
}


int ATGC_only( char* bases, int length ) {
  int i;
  for(i=0;i<length;i++) {
    if ( bases[i] == 'N' || bases[i] == 'n')
      return 0;
  }

  return 1;
}

void print_node(struct node *node) {
  
  printf("node seq: %s\t[%p,%p,%p,%p]\t%d\n", node->bases, node->np[0], node->np[1], node->np[2], node->np[3], node->count);

}


int add2tree(struct node *node, char *string, int length ) {

  print_node(node);

  printf("adding : |%s|\n", string);

  // Only for the initial unassigned root node.
  if (node->bases == 0 && 
      node->np[0] == 0 &&
      node->np[1] == 0 &&
      node->np[2] == 0 &&
      node->np[3] == 0 ) {
    printf("Init the root node\n");
    int basepos = base2pos(string[0]);
    node->np[ basepos ] =  malloc(sizeof(struct node));
    struct node * next_node = node->np[ basepos ];
    next_node->bases = string;
    next_node->count++;
    return 1;
  }

  int basepos = base2pos(string[0]);
  // at the root, and the next level down node does not exist, create it and return.
  if (node->bases == 0 &&  ! node->np[ basepos ] ) {
    node->np[ basepos ] = malloc(sizeof(struct node));
    struct node * next_node = node->np[ basepos ];
    next_node->bases = malloc(sizeof(char)*length+1);
    next_node->bases = string;
    next_node->count++;
    return 1;
  }
  
  if (node->bases == 0) {
    printf("Stepping down one level\n");
    add2tree(node->np[ basepos ], string, length);
  }
  else {
    int shared_length = common_string(node->bases, string, length);
    printf("SHARED :: %d bases\n", shared_length);

    int post_length = length - shared_length;    

    // no nodes further down, split the strings into the shared bit, and create two sub-nodes with the rest...
    if (node->np[0] == 0 &&
	node->np[1] == 0 &&
	node->np[2] == 0 &&
	node->np[3] == 0 ) {

      char *shared_bases = malloc(sizeof(char)*shared_length+1);
      strncpy(shared_bases, string, shared_length);
      shared_bases[shared_length]= '\0';

      char *post_bases1 = malloc(sizeof(char)*post_length);
      strcpy(post_bases1, &string[shared_length]);
      char *post_bases2 = malloc(sizeof(char)*post_length);
      post_bases2 = strcpy(post_bases2, &node->bases[shared_length]);
      
      node->count++;
      
      free(string);
      free(node->bases);

      node->bases = shared_bases;
      
      int basepos1 = base2pos(post_bases1[0]);
      int basepos2 = base2pos(post_bases2[0]);
      
      node->np[ basepos1 ] = malloc(sizeof(struct node));
      struct node * next_node = node->np[ basepos1 ];
      next_node->bases = post_bases1;
      next_node->count++;
      
      node->np[ basepos2 ] = malloc(sizeof(struct node));
      next_node = node->np[ basepos2 ];
      next_node->bases = post_bases2;
      next_node->count++;
      //      printf("'%s'\n'%s'\n shares %d bases [%s] [%s] - [%s] from the 5' end\n", node->bases, string, shared_length, shared_bases, post_bases1, post_bases2);
      printf("shares %d bases [%s] [%s] - [%s] from the 5' end\n", shared_length, shared_bases, post_bases1, post_bases2);
      print_node(node);

      return;
    }
    else { 
      
      printf("node sequence %s\n", node->bases);

      char *post_bases = malloc(sizeof(char)*post_length);
      strcpy(post_bases, &string[shared_length]);
      free(string);
      int post_base_pos = base2pos(post_bases[0]);
      printf("trimmed string %s to %s, index %d\n", string, post_bases, post_base_pos);
      node->np[ post_base_pos ] = malloc(sizeof(struct node));
      add2tree(node->np[ post_base_pos ], post_bases, post_length);
      print_node(node);
    }

  }
  return;
}




int common_string(char *s1, char *s2, int length) {

  //  printf("%s\n%s\n", s1, s2);

  int i;
  
  for(i=0;i<length;i++) {
    //printf("CS: %c-%c\n",s1[i],s2[i]);
    if (!s1[i] || !s2[i] || s1[i] != s2[i])
      return i;
  }



  return ;
}

struct fasta * FastaIn() {
/******************************************************************************
 * Reads a sequence from stdin into a fasta structure.
 *****************************************************************************/
  long p, ps;
  size_t i;
  char c;
  char *s;
  struct fasta *fs;

  fs = malloc(sizeof(struct fasta));

  // First the header/name (> Name ......)
  i = 0;

  p = getline(&fs->name, &i, stdin);
  fs->name[strlen(fs->name) - 1] = 0; 
  if (p == -1) {
    perror("FastaIn ");
    exit(4);
  }
  else if (p == EOF) {
    free(fs->name);
    free(fs);
    return NULL;
  }
  else if (fs->name[0] != '>') {
    fprintf(stderr, "Input error : The stdin should be fasta formated %s\n", 
            fs->name);
    exit(1);
  }

  i = 0;
  s = 0;
  fs->seq = NULL;//malloc(1 * sizeof(char));
  ps = 0;
  while ((c = getchar()) != '>' || c == EOF) {
    //while not EOF or new name
    ungetc(c, stdin);
    p = getline(&s, &i, stdin);

    if (p == EOF)
      break;

    if (p < 0) {
      perror("getline");
      break;
    }

    ps += p;

    fs->seq = (char * ) realloc(fs->seq, (ps+p)  * sizeof(char));
    strcpy(&fs->seq[ps - p], s);

    fs->seq[--ps] = 0;

    //free(s);
    p = 0;
    //s = 0;
  }
  ungetc(c, stdin);  
  free(s);
  fs->length = strlen(fs->seq);

  return(fs);
}
