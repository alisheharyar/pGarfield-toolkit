#include <iostream>
#include <string.h>
#include <cstdio>

#include "wcpplib/stream/findmark.h"
#include "wcpplib/stream/prstream.h"

/*
Copyright (c) 2000 I. B. Smirnov

Permission to use, copy, modify, distribute and sell this file for any purpose
is hereby granted without fee, provided that the above copyright notice,
this permission notice, and notices about any modifications of the original
text appear in all copies and in supporting documentation.
The file is provided "as is" without express or implied warranty.
*/

int findmark(std::istream &file, const char *s) {
  int ic;
  int l = strlen(s);  // length does not include end symbol
  int n;
  char *fs = new char[l + 1];
  for (n = 0; n < l; n++) {
    if ((ic = file.get()) == EOF) {
      delete[] fs;
      return 0;
    } else {
      fs[n] = ic;
    }
  }
  fs[l] = '\0';
  while (strcmp(fs, s) != 0) {
    for (n = 1; n < l; n++)
      fs[n - 1] = fs[n];
    if ((ic = file.get()) == EOF) {
      delete[] fs;
      return 0;
    }
    fs[l - 1] = ic;
    fs[l] = '\0';
  }
  delete[] fs;
  return 1;
}

int findmark_other_repeat(std::istream &file, std::ostream &outfile,
                          const char *s) {
  int ic;
  int l = strlen(s);  // length does not include end symbol
  int n;
  char *fs = new char[l + 1];
  for (n = 0; n < l; n++) {
    if ((ic = file.get()) == EOF) {
      fs[n] = ic;
      outfile.write(fs, n + 1);
      delete[] fs;
      return 0;
    } else {
      fs[n] = ic;
    }
  }
  fs[l] = '\0';
  int s_not_new_line = 0;
  while (strcmp(fs, s) != 0) {
    if (fs[0] != '\n') {
      if ((fs[0] == '\t' || fs[0] == ' ') && s_not_new_line == 0) {
        ;
      } else {
        s_not_new_line = 1;
        outfile.put(fs[0]);
      }
    } else {
      if (s_not_new_line == 1) outfile.put(fs[0]);  // to finish valid line
      s_not_new_line = 0;  // to avoid switching to new line with empty old one
    }
    for (n = 1; n < l; n++)
      fs[n - 1] = fs[n];
    if ((ic = file.get()) == EOF) {
      delete[] fs;
      return 0;
    }
    fs[l - 1] = ic;
    fs[l] = '\0';
  }
  delete[] fs;
  return 1;
}

int find1ofnmark(std::istream &file, int q, char *s[]) {
  //mcout<<"find1ofnmark is started\n";
  //char c;
  int ic;
  int *l = new int[q];
  int *pos_fs = new int[q];
  int s_init_pos_fs = 0;
  int i;
  int l_max = -1;
  for (i = 0; i < q; i++) {
    l[i] = strlen(s[i]);  // length does not include end symbol
    if (l[i] > l_max) l_max = l[i];
  }
  //l_max++;
  //Iprintn(mcout, q);
  //Iprintn(mcout, l_max);
  char *fs = new char[l_max + 1];
  //int qfs=0;                   // number of symbols in fs
  for (i = 0; i < q; i++) {
    pos_fs[i] = l_max;
  }
  fs[l_max] = '\0';
  //Iprintn(mcout, file.good());
  //Iprintn(mcout, file.eof());
  //Iprintn(mcout, file.fail());
  //Iprintn(mcout, file.bad());
  //mcout<<"State:"<< file.rdstate() <<'\n';
  while ((ic = file.get()) != EOF) {
    //Iprintn(mcout, ic);
    for (i = 1; i < l_max; i++) {
      fs[i - 1] = fs[i];
    }
    fs[l_max - 1] = ic;        // new symbol
    if (s_init_pos_fs == 0) {  // shift positions
      int ss = 1;
      for (i = 0; i < q; i++) {
        if (l_max - pos_fs[i] < l[i]) {
          pos_fs[i]--;
          ss = 0;
          //mcout<<"s[i]="<<s[i]<<'\n';
          //mcout<<"i="<<i<<" l_max="<<l_max<<" pos_fs[i]="<<pos_fs[i]
          //     <<" l[i]="<<l[i]<<" "<< &(fs[ pos_fs[i] ])<<'\n';
        }
      }
      if (ss == 1) s_init_pos_fs = 1;
    }
    for (i = 0; i < q; i++) {
      if (strcmp(&(fs[pos_fs[i]]), s[i]) == 0) {
        delete[] l;
        delete[] fs;
        delete[] pos_fs;
        return i;
      }
    }
  }
  //Iprintn(mcout, ic);
  delete[] l;
  delete[] fs;
  delete[] pos_fs;
  return -1;
}

int find1ofnmark(std::istream &file, int q,  // number of strings
                 const String str[])         // addresses
    {
  char **s = new char *[q];
  int i;
  for (i = 0; i < q; i++) {
#ifdef USE_STLSTRING
    s[i] = new char[strlen(str[i].c_str()) + 1];
    strcpy(s[i], str[i].c_str());
//Iprint(mcout, s[i]);
#else
    s[i] = new char[strlen(str[i]) + 1];
    strcpy(s[i], str[i]);
#endif
  }
  int iret = find1ofnmark(file, q, s);
  for (i = 0; i < q; i++) {
    delete[] s[i];
  }
  delete[] s;
  return iret;
}
