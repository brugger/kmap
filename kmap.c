

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

  char kmer[kmer_length + 1];
  kmer[kmer_length] = 0;

  printf("%s\n", seq);

  for(i=0;i <= length - kmer_length; i++) {
    strncpy(kmer, &seq[i], kmer_length);
    
    printf("[%s]\n", kmer);
    if (! ATGC_only(kmer, kmer_length))
	continue;


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




int add2tree(struct node *node, char *string, int length ) {

  printf("|%s|\n", string);

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
    next_node->bases = malloc(sizeof(char)*length+1);
    strcpy(next_node->bases,string);
    next_node->count++;
    return 1;
  }

  int basepos = base2pos(string[0]);
  if ( ! node->np[ basepos ] ) {
    node->np[ basepos ] = malloc(sizeof(struct node));
    struct node * next_node = node->np[ basepos ];
    next_node->bases = malloc(sizeof(char)*length+1);
    strcpy(next_node->bases,string);
    next_node->count++;
    return 1;
  }

  int shared = common_string(node->bases, string);
  printf("%s\n%s\n shares %d bases from the 5' end\n", node->bases, string, shared);

  add2tree(node->np[ basepos ], string, length);

  return;
}

int common_string(char *s1, char *s2, int length) {

  int i;
  
  for(i=0;i<length;i++) {
    if (s1[i] != s2[i])
      return i-1;
  }

  return length;
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
