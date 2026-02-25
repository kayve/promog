/*

This project aims to simplify the picture of proteomic studies without losing fine details. These studies are defined in the medical literature by data from myriad quantitative techniques that are difficult to distil holistically. The original thrust was determining which exact proteomic gene products were confined within plasma membranes.  According to the work of Singer and Nicolson from Science 175; 720-731; 1972, these proteins could be thought of as being constrained in space along folded sheets confined to two dimensional diffusion only, as opposed to having complete freedom to diffuse in three dimensions.  This idea was coined the Fluid Mosaic Model (FMM) of the Structure of Cell Membranes.  The raw proteomic data chosen for this study was obtained from the UniProt Knowledgebase provided publicly at http://www.uniprot.org/uniprotkb. As this work evolved, a four compartment model proposed by Satoh et al from Multiple Sclerosis; 15: 531-541; doi:10.1177/1352458508101943; 2009 was used.  This four compartment model was 1) nuclear, 2) cytosolic, 3) membrane, and 4) extracellular proteins.


        Copyright (C)  2026     Kayven Riese
                                kayvey@gmail.com
                                (415) 902-5513
                                3591 Quail Lakes Drive Unit 84
                                Stockton, CA   95207

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, see
<https://www.gnu.org/licenses/>.

*/

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <time.h>
#include <sys/mman.h>


int main (int argc, char *argv[]) {

   const int TRUE = 1;
   const int FALSE = 0;
   const int STDOUT = 1;
   const int STDIN = 0;
   const int BAD_INFILE = -1;
   const int OUT_OF_MEMORY = -2;
   const int TIME_ERR = -3;
   const int BAD_ARGC  = -4;
   const int MAXLINE = getpagesize();
   const int PAGESIZE = getpagesize();
   //const long long BLOCKSIZE = 1<<30;

   int fd, done, block_count, bytes_read,status,line_len;
   int biggest_line, file_arg = 1, bs_arg = 2,  opt;
   char  err_msg[MAXLINE], line[MAXLINE], *block, alloc_type='v';
   long long i = 0, inf_size, numlines = 0, j = 0, last_nl = 0, k = 0;
   long long line_begin = 0, last_line_begin = 0;
   char mem_method[20];

   FILE * fp;
   struct stat fstat;
   struct timespec t_begin, t_end, t_res;

   if (argc < 3) {
      sprintf(err_msg, "USAGE: demog_gets_sf [-mvap] <datafile> <log base 2 of BLOCKSIZE> O:=} Not");
      perror(err_msg);
      return  BAD_ARGC;
   }

   while ((opt =  getopt(argc, argv, "mvap"))  != EOF)  {
      switch (opt)  {
        case 'm':
          alloc_type = 'm';
          file_arg++;
          bs_arg++;
          break;
        case 'a':
          alloc_type = 'a';
          file_arg++;
          bs_arg++;
          break;
        case 'v':
          alloc_type = 'v';
          file_arg++;
          bs_arg++;
          break;
        case 'p':
          alloc_type = 'p';
          file_arg++;
          bs_arg++;
          break;
        case '?':
          sprintf(err_msg,"invalid option to %s:",argv[0]);
          perror(err_msg);
      } //--- switch (opt)  ---//
   } //--- while ((opt =  getopt(argc, argv, "mvap"))  != EOF)  ---//
  
   const long long BLOCKSIZE = 1 <<atoi(argv[bs_arg]);   
 
   switch (alloc_type) { 
      case 'm':
         block = malloc(BLOCKSIZE);
         strcpy(mem_method,"malloc");
         break;
      case 'v':
         block = valloc(BLOCKSIZE);
         strcpy(mem_method,"valloc");
         break;
      case 'a':
         block = alloca(BLOCKSIZE);
         strcpy(mem_method,"alloca");
         break;
      case 'p':
         posix_memalign((void **)&block, BLOCKSIZE, PAGESIZE);
         strcpy(mem_method,"posix_memalign");
         break;
   } //--- switch (alloc_type) ---///
 
   if (block == NULL) {
      sprintf(err_msg, "CAN'T MALLOC %lld BYTES  \ncause:",BLOCKSIZE);
      perror(err_msg);
      return OUT_OF_MEMORY;
   } //----- if ((block = malloc(BLOCKSIZE)) < 0) -----// 

   if ((fp = fopen( argv[file_arg], "r")) < 0) {
      sprintf(err_msg, "CAN'T OPEN FILE: %s  \ncause:",argv[file_arg]);
      perror(err_msg);
      fclose(fp);
      return BAD_INFILE;
   } //----- if ((fd = open( argv[1], O_RDONLY )) < 0)  -----//

   done = FALSE;

   if (clock_gettime(CLOCK_REALTIME, &t_begin)) {
      sprintf(err_msg,"failed to get start time\n\0");
      perror(err_msg);
      return TIME_ERR;
   }
 

   /*if (done) printf ("DONE IS TRUE\n");
   else printf ("DONE IS FALSE\n");/**/


   while (!done)  {
      status = (int)fgets(block,MAXLINE,fp);
      line_len = strlen(block); 
      if (line_len > biggest_line)
        biggest_line = line_len;
      done = !status;
   } //----- while (!done)  -----//


   if (clock_gettime(CLOCK_REALTIME, &t_end)) {
      sprintf(err_msg,"failed to get end time\n");
      perror(err_msg);
      return TIME_ERR;
   }

   //printf("there are %d lines in the file\n",numlines);
   printf("--------------------------------------------\n");
   printf("processing %s.\n", argv[file_arg]);
   printf("the memory allocation method is %s\n",mem_method);
   printf("the biggest line has %d characters\n",biggest_line);
   printf ("BLOCKSIZE is %lld\n",BLOCKSIZE);
   printf("it took ");
   print_interval(&t_begin,&t_end);
   printf(" to run.\n");
   printf("--------------------------------------------\n");
   close(fd);
} //----- int main (int argc, char *argv[])  -----//

