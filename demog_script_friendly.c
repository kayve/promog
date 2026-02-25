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
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>


#define    BILLION   1000000000 
#define STDOUT 1
#define STDIN 0 


int main (int argc, char *argv[]) {

  const int MAXLINE  = getpagesize();
  const int PAGESIZE = getpagesize();
  const int GOOD_EXIT  = 0;
  const int  FALSE   = 0;
  const int  TRUE =  1;
  const int  BAD_ARGC  = -1;
  const int BAD_DATAFILE  = -2 ;
  const int BAD_FSTAT = -3 ;
  const int TIME_ERR = -4 ;
  //const int  BLOCKSIZE   = 1<<25; 
//  const long long BLOCKSIZE   = 1<<31; 

/*********************************************************************
 *********************************************************************/
 
  int fd, fd_plist, a,b,c, i, len, tot, j,k,l,m,n,got_EOL, line_WRAP;
  int tot_proteins, tot_human_proteins, line_num, wrap_SHIFT;
  int wrap_ANNEAL, done, block_done, n_prot_lines, corrupt_infile; 
  int max_prot_lines,max_prot_chars, this_prot_chars, bytes_read;
  int file_arg = 1, bs_arg = 2;
  float r;
  int  max_line, zero_line;
  long long this_seek, line_begin, seek_result, last_lbegin,status;
  long long char_count;
  char *line, *this_line, this_char, *block, *wrap_frag, err_msg[MAXLINE],opt;
  char alloc_type = 'v', mem_method[20];
  struct stat statbuf;
  struct timespec t_begin, t_end, t_res;

  if (argc < 3) {
    sprintf(err_msg,"USAGE: demog_script_friendly [-mvap] <datafile> <log base 2 of BLOCKSIZE> O:=} Not");
    perror(err_msg);
    return BAD_ARGC;
  }
 
  while ((opt = getopt(argc,argv,"mvap")) !=EOF) {
    switch (opt) {
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
    }//--- switch (opt) ---//
  }//--- while ((opt= getopt(argc,argv,"m")) !=EOF) ---// 

  const long long BLOCKSIZE   = 1 << atoi(argv[bs_arg]); 

  if ((fd = open( argv[file_arg], O_RDONLY )) < 0) {
    sprintf(err_msg,"CAN'T OPEN FILE: %s \ncause", argv[file_arg]);
    perror(err_msg);
    close(fd);
    return BAD_DATAFILE; 
  }

  if ((fd_plist = open( "prots_loc", O_WRONLY )) < 0) {
    sprintf(err_msg,"CAN'T OPEN FILE: %s \ncause", "prots_loc");
    perror(err_msg);
    close(fd_plist);
    return BAD_DATAFILE; 
  }

  if (fstat(fd, &statbuf) < 0 ) {
    sprintf(err_msg,"BAD FILE STATS: %s \ncause", argv[1]);
    perror(err_msg);
    close(fd);

    return BAD_FSTAT;
  }

  i=0;
  a=0;
  b=0;
  c=0;
  tot_proteins = 0;
  line_num = 0;
  char_count = 0;
  tot_human_proteins = 0;
  max_prot_lines = 0;
  char_count = 0;
  zero_line = FALSE;
  corrupt_infile = FALSE;
  max_prot_chars = 0;
  this_prot_chars = 0;
  n_prot_lines = 0;

  /********************************************************
   *
   *   Data parameterization: 
   *
   *    BLOCKSIZE chunks of the files are removed as this_seek
   *     is incremented by that amount in the outer  while
   *     loop with termination on read failure.  The next
   *     level of while loop iterates on a variable line_begin,
   *     which must be less than BLOCKSIZE. This index steps
   *     through the currently examined block by leaps determined
   *     by an iterator b that checks for newline. 
   *
   ********************************************************/
  this_seek = 0;
  max_line = 0;
  line_WRAP = FALSE;
  done = FALSE; 
  block_done = FALSE; 
  wrap_ANNEAL = FALSE;
  if ((seek_result = lseek(fd, (off_t) this_seek, SEEK_SET))<0)
    done = TRUE;  

  switch (alloc_type) { 
     case 'm':
       block = malloc(BLOCKSIZE);
       strcpy(mem_method, "malloc");
       break;
     case 'p':
       posix_memalign((void **)&block,BLOCKSIZE,PAGESIZE);
       strcpy(mem_method, "posix_memalign");
       break;
     case 'a':
       block = alloca(BLOCKSIZE);
       strcpy(mem_method, "alloca");
       break;
     case 'v':
       block = valloc(BLOCKSIZE);
       strcpy(mem_method, "valloc");
       break;
     default:
       block = valloc(BLOCKSIZE);
  } //--- switch (alloc_type) ---//

  line_begin = 0;
  last_lbegin = 0;
  n_prot_lines = 0;
  //system("freecolor -o");
  
  if (clock_gettime(CLOCK_REALTIME, &t_begin)) {
    sprintf(err_msg, "failed to get start time\n\0");
    perror(err_msg);
    return TIME_ERR;
  }

/*********************************************************************
 *   this_seek--ing BLOCKSIZE steps
 *********************************************************************/
  if ((bytes_read = read(fd, block, BLOCKSIZE)) <= 0)  
    done = TRUE;
  while (!done)  {
   /*****************************************************************
    *   hopping line by line through one block
    *   block_done indicates data in block has been exhausted. 
    *
    *****************************************************************/
    block_done =  FALSE;
    while (!block_done) {
       
      b=0;
      /**************************************************************
       *  finding next line_begin with b 
       **************************************************************/
      got_EOL =  FALSE;
      while (!got_EOL) { 
           if (line_begin == bytes_read) {
             this_seek += BLOCKSIZE;
             if ((status= (lseek(fd, (off_t) this_seek, SEEK_SET)) < 0)) {
               done = TRUE;  
               got_EOL = TRUE;
               block_done = TRUE;
               break;
             }
             if ((bytes_read = read(fd, block, BLOCKSIZE)) <= 0)  {
               got_EOL = TRUE;
               block_done = TRUE;
               done = TRUE;
               break;
             }
             else
               line_begin = 0;
             }
           b++; 
           if (line_begin + b >= bytes_read) {
           /**************************************************************
            *  this line has wrapped accross a block.
            *  shift the file pointer so that the beginning of the line
            *  is the beginning of the new block.
            **************************************************************/
             if (bytes_read  < BLOCKSIZE) {
               /**************************************************************
                *  If bytes_read is less than BLOCKSIZE, this should imply
                *  that we are reading the last block of the infile.  Therefore,
                *  all loops should end.  In the event the file does not terminate
                *  correctly, the corrupt_infile flag is set. 
                **************************************************************/
               if (block[bytes_read-1] != '\n' )
                 corrupt_infile = TRUE; 
               got_EOL = TRUE;
               block_done = TRUE;
               done = TRUE;
               break;
             }
             wrap_SHIFT = b;
             this_seek += BLOCKSIZE - wrap_SHIFT;
             if ((status= (lseek(fd, (off_t) this_seek, SEEK_SET)) < 0))
               done = TRUE;  
             line_begin = 0;
             last_lbegin = 0;
           if ((bytes_read = read(fd, block, BLOCKSIZE)) <= 0)  {
         /**************************************************************
          *
          *  there is no more file to read--end all loops.
          *
          **************************************************************/
             //printf("status is %lld\n",status);
             got_EOL = TRUE;
             block_done = TRUE;
             done = TRUE;
           }
         } // --- if (line_begin + b == bytes_read) ---// 

         if ( block[line_begin + b ] == '\n')  {
           n_prot_lines++;
           char_count += b;
           this_prot_chars += b;
           got_EOL = TRUE;
         }
         if ( b == MAXLINE) {
           b -= 2;
           got_EOL = TRUE;
         } //----- if ( b == MAXLINE) -----// 

      }//----- while (!got_EOL) -----// 

      /**************************************************************
       **************************************************************/
      if (line_begin == bytes_read)
        break;
      if ( b > max_line) {
          max_line = b;  
      } // --- if ( b > max_line) ---// 
      if (line_begin != bytes_read)
         line_num++;
      if (b != 0)
        printas(STDOUT,block,line_begin,line_begin+b-1);
      a = 0;
      if ((block[line_begin] == '/')&&(block[line_begin+1] == '/')) {
        tot_proteins++;
        if ( n_prot_lines > max_prot_lines)
          max_prot_lines = n_prot_lines;
        if (this_prot_chars > max_prot_chars)
          max_prot_chars = this_prot_chars;
        n_prot_lines = 0;
        this_prot_chars= 0;
      } // -- if ((block[line_begin] == '/')&&(block[line_begin+1] == '/')) --// 

      if ((block[line_begin] == 'O')&&(block[line_begin+1] == 'S')) {
    /*****************************************************************
     *
     *  Organism Species (OS) line
     *
     *  http://www.expasy.org/sprot/userman.html#OS_line
     *
     *****************************************************************/
        if ((block[line_begin+5] == 'H')&&(block[line_begin+7] == 'm') &&(block[line_begin+10] == 's')) {
    /*****************************************************************
     *
     *    Human protein (Homo Sapiens)
     *
     *****************************************************************/
          tot_human_proteins++;
     
        } // if ((block[line_begin+5] == 'H')&&(block[line_begin+7] == 'm') &&(block[line_begin+10] == 's')) {
      } // -- if ((block[line_begin] == 'O')&&(block[line_begin+1] == 'S')) --// 
        last_lbegin = line_begin;
        line_begin += b + 1;
        if (line_begin >= bytes_read) {
          block_done = TRUE;
        }
      }//----- while (!block_done)-----// 
    /*****************************************************************
     *****************************************************************/

    if (line_begin != bytes_read)
      this_seek += BLOCKSIZE;
    block_done = FALSE;
    lseek(fd, (off_t) this_seek, SEEK_SET);  
    
  }//----- while (!done)  was ((bytes_read = read(fd, block, BLOCKSIZE)) > 0) -----// 

  if (clock_gettime(CLOCK_REALTIME, &t_end)) {
    sprintf(err_msg,"failed to get end time\n\0");
    perror(err_msg);
    return TIME_ERR;
  }
    /*****************************************************************
     *****************************************************************/
  printf("----------------------------------------\n"); 
  printf("processing %s \n", argv[file_arg] ); 
  printf("the memory allocation method is %s\n",mem_method); 
  printf("The longest line has %d characters\n",max_line);
  printf("There are a total of %d lines\n",line_num);
  printf("subtotal human proteins: %d\n",tot_human_proteins);
  printf("There are %d total proteins \n",tot_proteins);
  printf("The protein with the most lines has %d lines\n",max_prot_lines);/**/
  printf("BLOCKSIZE IS %d\n",BLOCKSIZE);
  printf("it took");
  print_interval(&t_begin,&t_end);
  printf(" to run.\n");
  printf("----------------------------------------\n"); 
    /*****************************************************************
     *
     *****************************************************************/
  close(fd);
  return GOOD_EXIT;
} //---- int main (int argc, char *argv[]) -----//


int printas (int fd, char * s, long long begin, long long end) {
/*****************************************************************
 *
 *   PRINTAS--Print Array Segment
 *  
 *  Prints an array segment to file fd avoiding segmentation fault
 *    
 *    String s[...begin to end. ...]
 *
 *****************************************************************/
  long long i;
  int j,k;
  if (begin <= end) {
    i = begin; 
    j = 0;
  }
}


