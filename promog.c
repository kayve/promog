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
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <regex.h>
#include <sys/mman.h>
#include <sys/stat.h>


#define    BILLION   1000000000 
#define STDOUT 1
#define STDIN 0 

static char * line;

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
  const int REGEX_ERR = -5 ;
  const int  BLOCKSIZE   = 1<<14; 
//  const long long BLOCKSIZE   = 1<<31; 

/*********************************************************************
 *  Determinants of "membrane" count
 *********************************************************************/
  const int MEMB_INDEX = 1;
  const int CELL_MEMB_INDEX = 0;
  const int CELL_SURF_INDEX = 21;
  const int GO_PMEMB_INDEX = 6;
  const int GO_INT_MEMB_INDEX = 4;

/*********************************************************************
 *  Determinants of "cytoplasmic" count
 *********************************************************************/
  const int CYTOPLASM_INDEX = 2;
  const int CYTOSOL_INDEX = 3;
  const int SOLUBLE_INDEX = 41;
  const int GO_CYTOSOL_INDEX = 7;
  const int GO_CYTOPLASM_INDEX = 2;

/*********************************************************************
 *  Determinants of "extracellular" count
 *********************************************************************/
  const int EXTRACELLULAR_INDEX = 4;
  const int SECRETED_INDEX = 5;
  const int GO_EXTRACELLULAR_INDEX = 3;
  const int GO_ECM_INDEX = 5;

/*********************************************************************
 *  Determinants of "nucleus" count
 *********************************************************************/
  const int NUCLEUS_INDEX = 6;
  const int TELOMERE_INDEX = 30;
  const int GO_NUCLEUS_INDEX = 0;
  const int GO_DNA_BIND_INDEX = 8;

/*********************************************************************
 *  Regular expression data 
 *********************************************************************/
  const int REGEX_COUNT = 45;
  const char* REGEX_RAW_ARRAY[] = { "[Cc]ell [Mm]embrane", "[mM]embrane", "[cC]ytoplasm", 
               "[cC]ytosol", "[Ee]xtracellular", "[Ss]ecreted", "[Nn]ucleus", 
               "[Mm]itochondrion", "[Ee]ndoplasmic reticulum lumen", "[Cc]ell junction", 
               "[Pp]eriplasm", "[Vv]acuole", "[Pp]lastid", "[Cc]apsid", 
               "[Ee]ndoplasmic reticulum", "[Ee]ndosome", "[Ll]ysosome", "[Vv]irion", 
               "[Cc]entromere", "[Pp]eroxisome", "[Gg]olgi", "[Cc]ell [Ss]urface", 
               "[Gg]lyoxysome", "[Gg]lyocosome", "[Zz]ona pellucida", "[Kk]inetochore", 
               "[Ss]pore", "[Bb]acterial", "[Ff]imbrium", "[Mm]elanosome", "[Tt]elomere", 
               "[Pp]odosome", "[Cc]ilium", "[Tt]richocyst", "[Hh]ydrogenosome", 
               "[Ss]arcoplasmic [Rr]eticulum", "[Aa]xon", "[Mm]icrosome", "[Aa]ngiotensin",
               "[Cc]hlorosome", "[tT]hylakoid", "[Ss]oluble", "[bB]ud", "[Ff]lagellum",
               "[Vv]iral"}; 

  const char* NAMES_ARRAY[] = { "Cell Membrane", "Membrane", "Cytoplasm", "Cytosol", 
               "Extracellular", "Secreted", "Nucleus", "Mitochondrion", 
               "Endoplasmic reticulum lumen", "Cell junction", "Periplasm", "Vacuole", 
               "Plastid", "Capsid", "Endoplasmic reticulum", "Endosome", "Lysosome", 
               "Virion", "Centromere", "Peroxisome", "Golgi", "Cell Surface", 
               "Glyoxysome", "Glyocosome", "Zona pellucida", "Kinetochore", "Spore", 
               "Bacterial", "Fimbrium", "Melanosome", "Telomere", "Podosome", "Cilium",
               "Trichocyst", "Hydrogenosome", "Sarcoplasmic Reticulum", "Axon", "Microsome",
               "Angiotensin", "Chlorosome", "Thylakoid", "Soluble", "Bud", "Flagellum",
               "Viral"}; 

  const int GO_COUNT = 9;
  const char * GO_RAW_REGEX_ARRAY[] = {"GO:0005634","GO:0007165","GO:0005737","GO:0005576","GO:0016021","GO:0031012","GO:0005886","GO:0005829","GO:0003677"};

  const char * GO_NAMES_ARRAY[] = {"Nucleus","Signal Transduction","Cytoplasm","Extracellular","Integral to Membrane","Extracellular Matrix","Plasma Membrane","Cytosol","DNA binding"};

  const int GO_MINOR_COUNT = 5;
  const char * GO_RAW_REGEX_MINOR_ARRAY[] = {"GO:0009103","GO:0030573","GO:0055114","GO:0033644"};
  const char * GO_NAMES_MINOR_ARRAY[] = {"lipopolysaccharide biosynthetic process","Bile Aid Catabolic Process","Oxidation Reduction","Host Cell Membrane"
};
 
/*********************************************************************
 *  Essential counters & roll flaps 
 *********************************************************************/
  int a,b,c, i, len, tot, j,k,l,m,n,got_EOL, line_WRAP;
  int tot_proteins, tot_human_proteins, line_num, wrap_SHIFT;
  int wrap_ANNEAL, done, block_done, n_prot_lines, corrupt_infile; 
  int max_prot_lines,max_prot_chars, this_prot_chars, bytes_read;
  int  max_line, zero_line;

/*********************************************************************
 *   File Descriptors 
 *********************************************************************/
  int fd, fd_plist, fd_ERROR, fd_GO_REMAINDER; 
  int fd_FT_TOPO_DOM, fd_REMAINDER;

/*********************************************************************
 *  Boolean protein attribute flags
 *********************************************************************/
  int is_FT_TRANSMEM, is_FT_INTRAMEM, is_FT_TD_extracellular, is_FT_TD_cytoplasmic; 
  int is_me_DUPE, is_REMAINDER, has_FT_SIGNAL, has_FT_SIG_TRANSMEM, has_FT_DNA_BIND;
  int is_FT_LIPID, is_mc_DUPE, is_mn_DUPE, is_ce_DUPE, is_cn_DUPE, is_ne_DUPE;
  int is_FLAGGED, is_SCL_ARRAY[REGEX_COUNT], has_SCL, has_DR_GO, has_GO_ARRAY[GO_COUNT]; 
  int is_GO_REMAINDER, has_GO_MINOR_ARRAY[GO_MINOR_COUNT], is_brain, is_muscle;

/*********************************************************************
 *   Human protein tabulators
 *********************************************************************/
  int hum_transmem, hum_extracellular, hum_cytoplasmic, hum_SIGNAL, hum_SIG_TRANSMEM;
  int hum_DNA_BIND, hum_mem, hum_intramem, hum_itmem, hum_lipid_bind, hum_membrane; 
  int in_SCL, hum_REMAINDER, hum_SCL_ARRAY[REGEX_COUNT], hum_SCL_NULL, hum_DR_GO;
  int hum_GO_ARRAY[GO_COUNT], hum_GO_MINOR_ARRAY[GO_MINOR_COUNT], hum_nuclear; 

/*********************************************************************
 *   Total protein tabulators
 *********************************************************************/
  int tot_transmem, tot_extracellular, tot_cytoplasmic, tot_SIGNAL, tot_SIG_TRANSMEM;
  int tot_DNA_BIND, tot_mem, tot_intramem, tot_itmem, tot_lipid_bind, tot_REMAINDER;
  int tot_SCL_ARRAY[REGEX_COUNT], tot_membrane, tot_SCL_NULL, tot_DR_GO;
  int tot_GO_ARRAY[GO_COUNT], tot_GO_MINOR_ARRAY[GO_MINOR_COUNT], tot_nuclear; 
  int tot_muscle, tot_brain, brain_cytoplasmic, brain_nuclear, brain_membrane;
  int brain_extracellular, muscle_cytoplasmic, muscle_nuclear, muscle_membrane;
  int muscle_extracellular;

/*********************************************************************
 *   Miscellaneous, timing, memory 
 *********************************************************************/
  int file_arg = 1, bs_arg = 2, this_is_human = 0;
  float r;
  long long this_seek, line_begin, seek_result, last_lbegin,status;
  long long char_count;
  char *this_line, this_char, *block, *wrap_frag, err_msg[MAXLINE],opt;
  char alloc_type = 'v', mem_method[20];
  struct stat statbuf;
  struct timespec t_begin, t_end, t_res;

/*********************************************************************
 *  REGular EXpressions compiled variables   
 *********************************************************************/
  regex_t rgx_array[REGEX_COUNT], rgx_GO_array[GO_COUNT], 
      rgx_GO_minor_array[GO_MINOR_COUNT], rgx_brain, rgx_muscle;
  int regex_status;
/*********************************************************************
 *  END VARIABLES SECTION 
 *  *************  **************  **************   ************** 
 *  END VARIABLES SECTION 
 *  *************  **************  **************   ************** 
 *  END VARIABLES SECTION 
 *  *************  **************  **************   ************** 
 *  END VARIABLES SECTION 
 *********************************************************************/

  line = malloc((MAXLINE+2)*sizeof(char));

 #if 0
 for (i=0;i<REGEX_COUNT;i++)
   printf("index %d: RAW REGEX: %s NAME: %s\n",i,REGEX_RAW_ARRAY[i],NAMES_ARRAY[i]);
 #endif

/*********************************************************************
 *  REGular EXpression COMPilations
 *********************************************************************/
  for (i=0;i<REGEX_COUNT;i++) {
    regex_status = regcomp(&rgx_array[i],REGEX_RAW_ARRAY[i] , REG_EXTENDED|REG_NOSUB);
    if (regex_status) {
        fprintf(stderr, "Could not compile regex for %s\n",REGEX_RAW_ARRAY[i]);
        exit(REGEX_ERR);
    }
  }

  for (i=0;i<GO_COUNT;i++) {
    regex_status = regcomp(&rgx_GO_array[i],GO_RAW_REGEX_ARRAY[i] , REG_EXTENDED|REG_NOSUB);
    if (regex_status) {
        fprintf(stderr, "Could not compile regex for %s\n",GO_RAW_REGEX_ARRAY[i]);
        exit(REGEX_ERR);
    }
  }

  for (i=0;i<GO_MINOR_COUNT;i++) {
    regex_status = regcomp(&rgx_GO_minor_array[i],GO_RAW_REGEX_MINOR_ARRAY[i] 
         , REG_EXTENDED|REG_NOSUB);
    if (regex_status) {
        fprintf(stderr, "Could not compile regex for %s\n",GO_RAW_REGEX_MINOR_ARRAY[i]);
        exit(REGEX_ERR);
    }
  }

  regex_status = regcomp(&rgx_brain,"TISSUE=Brain", REG_EXTENDED|REG_NOSUB);
    if (regex_status) {
        fprintf(stderr, "Could not compile regex for %s\n","TISSUE=Brain");
        exit(REGEX_ERR);
    }

  regex_status = regcomp(&rgx_muscle,"TISSUE=Muscle", REG_EXTENDED|REG_NOSUB);
    if (regex_status) {
        fprintf(stderr, "Could not compile regex for %s\n","TISSUE=Muscle");
        exit(REGEX_ERR);
    }

  if (argc < 2) {
    sprintf(err_msg,"USAGE: demog_script_friendly [-mvap] <datafile> [log base 2 of BLOCKSIZE] O:=} Not");
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

#if 0 
  const long long BLOCKSIZE;
  if (argc >= 3)
    BLOCKSIZE   = 1 << atoi(argv[bs_arg]); 
  else
    BLOCKSIZE   = 1 << 14; 
#endif

  if ((fd = open( argv[file_arg], O_RDONLY )) < 0) {
    sprintf(err_msg,"CAN'T OPEN FILE: %s \ncause", argv[file_arg]);
    perror(err_msg);
    close(fd);
    return BAD_DATAFILE; 
  }

#if 0
  if ((fd_plist = open( "human_subc_loc.txt", O_WRONLY )) < 0) {
    sprintf(err_msg,"CAN'T OPEN FILE: %s \ncause", "human_subc_loc.txt");
    perror(err_msg);
    close(fd_plist);
    return BAD_DATAFILE; 
  }

  if ((fd_FT_TOPO_DOM = open( "human_FT_TOPO_DOM.txt", O_WRONLY )) < 0) {
    sprintf(err_msg,"CAN'T OPEN FILE: %s \ncause", "human_FT_TOPO_DOM.txt");
    perror(err_msg);
    close(fd_FT_TOPO_DOM);
    return BAD_DATAFILE; 
  }

  if ((fd_GO_REMAINDER = open( "GO_REMAINDER.txt", O_WRONLY )) < 0) {
    sprintf(err_msg,"CAN'T OPEN FILE: %s \ncause", "GO_REMAINDER.txt");
    perror(err_msg);
    close(fd_FT_TOPO_DOM);
    return BAD_DATAFILE; 
  }

  if ((fd_REMAINDER = open( "human_REMAINDER.txt", O_WRONLY )) < 0) {
    sprintf(err_msg,"CAN'T OPEN FILE: %s \ncause", "human_REMAINDER.txt");
    perror(err_msg);
    close(fd_REMAINDER);
    return BAD_DATAFILE; 
  }

  if ((fd_ERROR = open( "ERROR.txt", O_WRONLY )) < 0) {
    sprintf(err_msg,"CAN'T OPEN FILE: %s \ncause", "ERROR.txt");
    perror(err_msg);
    close(fd_ERROR);
    return BAD_DATAFILE; 
  }

#endif

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
  hum_transmem = 0;
  hum_mem = 0;
  hum_DR_GO = 0;
  hum_itmem = 0;
  hum_intramem = 0;
  hum_SIGNAL= 0;
  tot_transmem = 0;
  tot_mem = 0;
  tot_brain = 0;
  tot_muscle = 0;
  tot_lipid_bind = 0;
  hum_lipid_bind = 0;
  tot_itmem = 0;
  tot_intramem = 0;
  tot_SIGNAL= 0;
  hum_DNA_BIND = 0;
  hum_REMAINDER = 0;
  tot_REMAINDER = 0;
  hum_nuclear = 0;
  tot_nuclear = 0;
  brain_nuclear = 0;
  muscle_nuclear = 0;
  hum_extracellular = 0;
  brain_extracellular = 0;
  muscle_extracellular = 0;
  hum_cytoplasmic = 0;
  brain_cytoplasmic = 0;
  muscle_cytoplasmic = 0;
  hum_membrane = 0;
  brain_membrane = 0;
  muscle_membrane = 0;
  hum_SCL_NULL = 0;
  tot_SCL_NULL = 0;
  tot_DR_GO = 0;
  
  for(i=0;i<REGEX_COUNT;i++) {
     hum_SCL_ARRAY[i] = 0;
     tot_SCL_ARRAY[i] = 0;
     is_SCL_ARRAY[i] = FALSE;
  }

  for(i=0;i<GO_COUNT;i++) {
     hum_GO_ARRAY[i] = 0;
     tot_GO_ARRAY[i] = 0;
     has_GO_ARRAY[i] = FALSE;
  }

  for(i=0;i<GO_MINOR_COUNT;i++) {
     hum_GO_MINOR_ARRAY[i] = 0;
     tot_GO_MINOR_ARRAY[i] = 0;
     has_GO_MINOR_ARRAY[i] = FALSE;
  }

  tot_extracellular = 0;
  tot_membrane = 0;
  tot_cytoplasmic = 0;
  max_prot_lines = 0;
  char_count = 0;
  zero_line = FALSE;
  corrupt_infile = FALSE;
  max_prot_chars = 0;
  this_prot_chars = 0;
  n_prot_lines = 0;

  this_is_human = FALSE;
  is_FT_TRANSMEM = FALSE;
  is_brain = FALSE;
  is_muscle = FALSE;
  is_REMAINDER = TRUE;
  is_GO_REMAINDER = TRUE;
  is_FT_INTRAMEM = FALSE;
  is_FT_LIPID = FALSE;
  has_FT_SIGNAL = FALSE;
  has_DR_GO = FALSE;
  has_SCL = FALSE;
  has_FT_DNA_BIND = FALSE;
  has_FT_SIG_TRANSMEM = FALSE;
  is_FT_TD_extracellular = FALSE;
  is_FT_TD_cytoplasmic = FALSE;
  in_SCL = FALSE;
  is_FLAGGED = FALSE;

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
         //printf ("e!gE.");
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
         }// if (line_begin + b == bytes_read) {

         if ( block[line_begin + b ] == '\n')  {
           n_prot_lines++;
           char_count += b;
           this_prot_chars += b;
           got_EOL = TRUE;
         }
         if ( b == MAXLINE) {
           b -= 2;
           got_EOL = TRUE;
         }//----- if ( b == MAXLINE) -----// 
      }//----- while (!got_EOL) -----// 

      if (line_begin == bytes_read)
        break;
      if ( b > max_line)
          max_line = b;  
      if (line_begin != bytes_read)
         line_num++;

      #if 0
      if (b != 0)
        printas(STDOUT,block,line_begin,line_begin+b-1);
      #endif

      a = 0;

      if ((block[line_begin] == '/')&&(block[line_begin+1] == '/')) {
    /*****************************************************************
     *  END OF RECORD
     *****************************************************************/
        tot_proteins++;
        if ( n_prot_lines > max_prot_lines)
          max_prot_lines = n_prot_lines;
        if (this_prot_chars > max_prot_chars)
          max_prot_chars = this_prot_chars;
        if (this_is_human) {
       /*****************************************************************
        * HUMAN DATA 
        *****************************************************************/
          if (!has_SCL) {
            hum_SCL_NULL++;
            is_REMAINDER = FALSE;
          }
          if (is_REMAINDER) 
            hum_REMAINDER++;
          if (is_FT_TRANSMEM) 
            hum_transmem++;
          if (is_FT_INTRAMEM)
            hum_intramem++;
          if (is_FT_LIPID)
            hum_lipid_bind++;
          if ((is_FT_INTRAMEM)&&(is_FT_TRANSMEM))
            hum_itmem++;
          if (is_FT_TD_extracellular)
            hum_extracellular++;
          if (is_FT_TD_cytoplasmic)
            hum_cytoplasmic++;
          if (has_FT_SIGNAL)
            hum_SIGNAL++;
          if (has_FT_DNA_BIND)
            hum_DNA_BIND++;
          if (has_DR_GO)
            hum_DR_GO++;
          if ((has_FT_SIGNAL) && (is_FT_TRANSMEM))
            hum_SIG_TRANSMEM++;
          for(i=0;i<REGEX_COUNT;i++) 
             if(is_SCL_ARRAY[i])
                hum_SCL_ARRAY[i]++;

          for(i=0;i<GO_COUNT;i++) 
             if(has_GO_ARRAY[i])
                hum_GO_ARRAY[i]++;

          for(i=0;i<GO_MINOR_COUNT;i++) 
             if(has_GO_MINOR_ARRAY[i])
                hum_GO_MINOR_ARRAY[i]++;

          if ((is_FT_TRANSMEM) || (is_FT_INTRAMEM) || (is_FT_LIPID) || 
              (is_SCL_ARRAY[CELL_SURF_INDEX]) || (is_SCL_ARRAY[CELL_MEMB_INDEX])
              || (is_SCL_ARRAY[MEMB_INDEX]) )
            hum_membrane++;

          if ((is_SCL_ARRAY[CYTOPLASM_INDEX]) || (is_SCL_ARRAY[CYTOSOL_INDEX]) ||
	      (is_SCL_ARRAY[SOLUBLE_INDEX]) || (has_GO_ARRAY[GO_CYTOSOL_INDEX]) ||
              (has_GO_ARRAY[GO_CYTOPLASM_INDEX]))
            hum_cytoplasmic++;

          if ((has_FT_SIGNAL) || (is_SCL_ARRAY[EXTRACELLULAR_INDEX]) ||
              (is_SCL_ARRAY[SECRETED_INDEX]) || (has_GO_ARRAY[GO_EXTRACELLULAR_INDEX]) ||
              (has_GO_ARRAY[GO_ECM_INDEX]))
            hum_extracellular++;

          if ((is_SCL_ARRAY[NUCLEUS_INDEX]) || (is_SCL_ARRAY[TELOMERE_INDEX]) ||
              (has_GO_ARRAY[GO_NUCLEUS_INDEX]) || (has_GO_ARRAY[GO_DNA_BIND_INDEX]))
            hum_nuclear++;
        }//----  HUMAN DATA -----//

       /*****************************************************************
       * TOTAL DATA 
       *****************************************************************/
        if (is_FT_TRANSMEM)
          tot_transmem++;
        if (!has_SCL) {
          tot_SCL_NULL++;
          is_REMAINDER = FALSE;
        }
        if (is_REMAINDER) 
          tot_REMAINDER++;
        if (is_muscle) 
          tot_muscle++;
        if (is_brain) 
          tot_brain++;
        if (is_FT_LIPID)
          tot_lipid_bind++;
        if (has_DR_GO)
          tot_DR_GO++;
        if (is_FT_INTRAMEM)
          tot_intramem++;
        if ((is_FT_INTRAMEM)&&(is_FT_TRANSMEM))
          tot_itmem++;
        if (is_FT_TD_extracellular)
          tot_extracellular++;
        if (is_FT_TD_cytoplasmic)
          tot_cytoplasmic++;
        if (has_FT_SIGNAL)
          tot_SIGNAL++;
        if (has_FT_DNA_BIND)
          tot_DNA_BIND++;
        if ((has_FT_SIGNAL) && (is_FT_TRANSMEM))
          tot_SIG_TRANSMEM++;
        for(i=0;i<REGEX_COUNT;i++) 
           if(is_SCL_ARRAY[i])
              tot_SCL_ARRAY[i]++;

        for(i=0;i<GO_COUNT;i++) 
           if(has_GO_ARRAY[i])
              tot_GO_ARRAY[i]++;

        for(i=0;i<GO_MINOR_COUNT;i++) 
           if(has_GO_MINOR_ARRAY[i])
              tot_GO_MINOR_ARRAY[i]++;

        if ((is_FT_TRANSMEM) || (is_FT_INTRAMEM) || (is_FT_LIPID) || 
            (is_SCL_ARRAY[CELL_SURF_INDEX]) || (is_SCL_ARRAY[CELL_MEMB_INDEX])
            || (is_SCL_ARRAY[MEMB_INDEX]) )
        {
          tot_membrane++;
          if (is_brain)
            brain_membrane++;
          if (is_muscle)
            muscle_membrane++;
        }

        if ((is_SCL_ARRAY[CYTOPLASM_INDEX]) || (is_SCL_ARRAY[CYTOSOL_INDEX]) ||
            (is_SCL_ARRAY[SOLUBLE_INDEX]) || (has_GO_ARRAY[GO_CYTOSOL_INDEX]) ||
            (has_GO_ARRAY[GO_CYTOPLASM_INDEX]))
        {
          tot_cytoplasmic++;
          if (is_brain)
            brain_cytoplasmic++;
          if (is_muscle)
            muscle_cytoplasmic++;
        }

        if ((has_FT_SIGNAL) || (is_SCL_ARRAY[EXTRACELLULAR_INDEX]) ||
            (is_SCL_ARRAY[SECRETED_INDEX]) || (has_GO_ARRAY[GO_EXTRACELLULAR_INDEX]) ||
            (has_GO_ARRAY[GO_ECM_INDEX]))
        {
          tot_extracellular++;
          if (is_brain)
            brain_extracellular++;
          if (is_muscle)
            muscle_extracellular++;
        }

        if ((is_SCL_ARRAY[NUCLEUS_INDEX]) || (is_SCL_ARRAY[TELOMERE_INDEX]) ||
            (has_GO_ARRAY[GO_NUCLEUS_INDEX]) || (has_GO_ARRAY[GO_DNA_BIND_INDEX]))
        {
          tot_nuclear++;
          if (is_brain)
            brain_nuclear++;
          if (is_muscle)
            muscle_nuclear++;
        }

    /*****************************************************************
     *   RESET this protein data
     *****************************************************************/
        n_prot_lines = 0;
        this_prot_chars= 0;
        this_is_human = FALSE;
        is_FT_TRANSMEM = FALSE;
        is_REMAINDER = TRUE;
        is_GO_REMAINDER = TRUE;
        is_FT_INTRAMEM = FALSE;
        is_FT_LIPID = FALSE;
        is_brain = FALSE;
        is_muscle = FALSE;
        has_FT_SIGNAL = FALSE;
        has_SCL = FALSE;
        has_DR_GO = FALSE;
        has_FT_DNA_BIND = FALSE;
        has_FT_SIG_TRANSMEM = FALSE;
        is_FT_TD_extracellular = FALSE;
        is_FT_TD_cytoplasmic = FALSE;
        in_SCL = FALSE;
        is_FLAGGED = FALSE;

        for(i=0;i<REGEX_COUNT;i++) 
           is_SCL_ARRAY[i] = FALSE;

        for(i=0;i<GO_COUNT;i++) 
           has_GO_ARRAY[i] = FALSE;

        for(i=0;i<GO_MINOR_COUNT;i++) 
           has_GO_MINOR_ARRAY[i] = FALSE;

      }// ---if ((block[line_begin] == '/')&&(block[line_begin+1] == '/'))---// 

    /*****************************************************************
     *  Done with End Record processing 
     *****************************************************************/

#if 1
      if ((block[line_begin] == 'R')&&(block[line_begin+1] == 'C')) {
    /*****************************************************************
     *
     *  Reference Comment (RC) line
     *
     *  http://ca.expasy.org/sprot/userman.html#RC_line
     *
     *****************************************************************/
        for(i=0;i<b+1;i++)
          line[i] = block[line_begin+i];
        line[i] = '\0';
        if (!(regexec(&rgx_muscle, line, (size_t)0,NULL,0))) {
          is_muscle = TRUE;
        }
        if (!(regexec(&rgx_brain, line, (size_t)0,NULL,0))) {
          is_brain = TRUE;
        }
        for(i=0;i<b+2;i++)
          line[i] = '\0';

      }//---if ((block[line_begin] == 'R')&&(block[line_begin+1] == 'C'))---// 
#endif

      if ((block[line_begin] == 'O')&&(block[line_begin+1] == 'S')) {
    /*****************************************************************
     *
     *  Organism Species (OS) line
     *
     *  http://ca.expasy.org/sprot/userman.html#OS_line
     *
     *****************************************************************/
        if ((block[line_begin+5] == 'H')&&(block[line_begin+7] == 'm') 
             &&(block[line_begin+10] == 's')) {
    /*****************************************************************
     *    Human protein (Homo Sapiens)
     *****************************************************************/
          tot_human_proteins++;
          this_is_human = TRUE;
        } 
     } 

     if ((block[line_begin] == 'F')&&(block[line_begin+1] == 'T') ) {
        /*****************************************************************
        *
        *  Feature Table (FT) line
        *
        *  http://www.expasy.org/sprot/userman.html#FT_line
        *
        *****************************************************************/
          if ((block[line_begin+5] == 'T') && (block[line_begin+6] == 'R') &&
              (block[line_begin+7] == 'A') && (block[line_begin+8] == 'N') &&
              (block[line_begin+9] == 'S') && (block[line_begin+10] == 'M') &&
              (block[line_begin+11] == 'E') && (block[line_begin+12] == 'M')) {
              is_FT_TRANSMEM = TRUE;
              is_REMAINDER = FALSE;
          } //--- if FT  TRANSMEM ----//

          if ((block[line_begin+5] == 'L') && (block[line_begin+6] == 'I') &&
              (block[line_begin+7] == 'P') && (block[line_begin+8] == 'I') &&
              (block[line_begin+9] == 'D')) {
              is_FT_LIPID= TRUE;
              is_REMAINDER = FALSE;
          } //--- if FT LIPID ----//

          if ((block[line_begin+5] == 'I') && (block[line_begin+6] == 'N') &&
              (block[line_begin+7] == 'T') && (block[line_begin+8] == 'R') &&
              (block[line_begin+9] == 'A') && (block[line_begin+10] == 'M') &&
              (block[line_begin+11] == 'E') && (block[line_begin+12] == 'M')) {
              is_FT_INTRAMEM = TRUE;
              is_REMAINDER = FALSE;
          } //--- if FT  INTRAMEM ----//

          if ((block[line_begin+5] == 'S') && (block[line_begin+6] == 'I') &&
              (block[line_begin+7] == 'G') && (block[line_begin+8] == 'N') &&
              (block[line_begin+9] == 'A') && (block[line_begin+10] == 'L')) {
              has_FT_SIGNAL = TRUE;
              is_REMAINDER = FALSE;
          } //--- if FT  SIGNAL ----//
          if ((block[line_begin+5] == 'D') && (block[line_begin+6] == 'N') &&
              (block[line_begin+7] == 'A') && (block[line_begin+8] == '_') &&
              (block[line_begin+9] == 'B') && (block[line_begin+10] == 'I') &&
              (block[line_begin+11] == 'N') && (block[line_begin+12] == 'D')) {
              has_FT_DNA_BIND = TRUE;
              is_REMAINDER = FALSE;
          } //--- if FT  DNA_BIND----//

        } //---if ((block[line_begin] == 'F')&&(block[line_begin+1] == 'T')) {


        if ((block[line_begin] == 'C')&&(block[line_begin+1] == 'C')) {
        /*****************************************************************
        *
        *  Comment Block (CC) line
        *
        *  http://www.expasy.org/sprot/userman.html#CC_line
        *
        *****************************************************************/
           if ((block[line_begin+9] == 'S') && (block[line_begin+10] == 'U')&&
               (block[line_begin+11] == 'B') && (block[line_begin+12] == 'C')&&
               (block[line_begin+13] == 'E') && (block[line_begin+14] == 'L')&&
               (block[line_begin+15] == 'L') && (block[line_begin+16] == 'U')&&
               (block[line_begin+17] == 'L') && (block[line_begin+18] == 'A')&&
               (block[line_begin+19] == 'R') && (block[line_begin+21] == 'L')&&
               (block[line_begin+22] == 'O') && (block[line_begin+23] == 'C')&&
               (block[line_begin+24] == 'A') && (block[line_begin+25] == 'T')&&
               (block[line_begin+26] == 'I') && (block[line_begin+27] == 'O')&&
               (block[line_begin+5] == '-')&&(block[line_begin+6] == '!') &&
               (block[line_begin+7] == '-') && (block[line_begin+28] == 'N')) {
            /*****************************************************************
             *
             *  CC   -!-  SUBCELLULAR LOCATION
             *
              *****************************************************************/
                has_SCL = TRUE;
                in_SCL = TRUE;
             } //---  CC   -!-  SUBCELLULAR LOCATION ----//
             else {
               if (((block[line_begin+5] == '-')&&(block[line_begin+6] == '!') &&
                 (block[line_begin+7] == '-'))||
                   ((block[line_begin+5] == '-')&&(block[line_begin+6] == '-') &&
                  (block[line_begin+7] == '-')))
                    in_SCL = FALSE;
             }
             if (in_SCL)  {
                 #if 1
                 for(i=0;i<b+1;i++)
                    line[i] = block[line_begin+i];
                 line[i] = '\0';
                 #endif
                 
                 for(i=0;i<REGEX_COUNT;i++) {
                    if (!(regexec(&rgx_array[i], line, (size_t)0,NULL,0))) {
                       is_SCL_ARRAY[i] = TRUE;
                       is_REMAINDER = FALSE;
                    }
                 }
             } 

        } //---if ((block[line_begin] == 'C')&&(block[line_begin+1] == 'C')) {

     if ((block[line_begin] == 'D')&&(block[line_begin+1] == 'R') ) {
        /*****************************************************************
         *
         *  DR -!-  Database cross-Reference 
         *
         *****************************************************************/

         if ((block[line_begin+5] == 'G')&&(block[line_begin+6] == 'O') ) {
            /*****************************************************************
             *
             *  GO  -!-  Gene Ontology reference 
             *
             *****************************************************************/
             for(i=0;i<b+1;i++)
               line[i] = block[line_begin+i];
             line[i] = '\0';

             for(i=0;i<GO_COUNT;i++) {
                if (!(regexec(&rgx_GO_array[i], line, (size_t)0,NULL,0))) {
                   has_GO_ARRAY[i] = TRUE;
                   is_REMAINDER = FALSE;
                   is_GO_REMAINDER = FALSE;
                }
             }

             for(i=0;i<GO_MINOR_COUNT;i++) {
                if (!(regexec(&rgx_GO_minor_array[i], line, (size_t)0,NULL,0))) {
                   has_GO_MINOR_ARRAY[i] = TRUE;
                   is_REMAINDER = FALSE;
                   is_GO_REMAINDER = FALSE;
                }
             }
             
             if (is_GO_REMAINDER)
                write(fd_GO_REMAINDER,line,b+1);
         } //--- ((block[line_begin+5] == 'G')&&(block[line_begin+6] == 'O') ) ---//
 
     } //--- if ((block[line_begin] == 'D')&&(block[line_begin+1] == 'R') ) ---// 

        if (is_REMAINDER) 
            write(fd_REMAINDER,line,b+1);

        last_lbegin = line_begin;
        line_begin += b + 1;
        if (line_begin >= bytes_read) {
          block_done = TRUE;
        }

      }//----- while (!block_done)-----// 

    /*****************************************************************
     * wrap up flip
     *****************************************************************/
    if (line_begin != bytes_read)
      this_seek += BLOCKSIZE;
    block_done = FALSE;
    lseek(fd, (off_t) this_seek, SEEK_SET);  
    
  }//----- while (!done)  was ((bytes_read = read(fd, block, BLOCKSIZE)) > 0) -----// 

  hum_membrane -= hum_SIGNAL;
  tot_membrane -= tot_SIGNAL;

  if (clock_gettime(CLOCK_REALTIME, &t_end)) {
    sprintf(err_msg,"failed to get end time\n\0");
    perror(err_msg);
    return TIME_ERR;
  }
    /*****************************************************************
     *
     *  OUTPUT RESULTS
     *
     *****************************************************************/
  printf("----------------------------------------\n"); 
  printf("processing %s \n", argv[file_arg] ); 
  printf("the memory allocation method is %s\n",mem_method); 
  printf("The longest line has %d characters\n",max_line);
  printf("There are a total of %d lines\n",line_num);
  printf("--------HUMAN PROTEINS--------------------\n"); 
  printf("human proteins: %d\n",tot_human_proteins);
  printf("human FT TRANSMEM proteins: %d\n",hum_transmem);
  printf("human FT INTRAMEM proteins: %d\n",hum_intramem);
  printf("human proteins with covalent lipid binding (FT LIPID): %d\n",hum_lipid_bind);
  printf("human proteins both intra- & trans- membrane: %d\n",hum_itmem);
  printf("human FT SIGNAL signal peptide containing proteins: %d\n",hum_SIGNAL);
  printf("human proteins with both signal sequence and transmembrane: %d\n",hum_SIG_TRANSMEM);
  printf("human proteins with DNA_BIND: %d\n",hum_DNA_BIND);

  for(i=0;i<REGEX_COUNT;i++) 
     printf("%d: human proteins with CC SUBCELLULAR LOCATION \"%s\": %d\n",
         i, NAMES_ARRAY[i], hum_SCL_ARRAY[i]);

  for(i=0;i<GO_COUNT;i++) 
     printf("%d: human proteins with Gene Ontology \"%s\": %d\n",
         i, GO_NAMES_ARRAY[i], hum_GO_ARRAY[i]);

  printf("human proteins with no CC SUBCELLULAR LOCATION annotation: %d\n",hum_SCL_NULL);
  printf("human total membrane proteins: %d\n",hum_membrane);
  printf("human cytoplasmic proteins: %d\n",hum_cytoplasmic);
  printf("human extracellular proteins: %d\n",hum_extracellular);
  printf("human nuclear proteins: %d\n",hum_nuclear);
  printf("REMAINDER human proteins: %d\n",hum_REMAINDER);
  printf("----------------------------------------\n"); 

  printf("There are %d total proteins \n",tot_proteins);
  printf("total FT TRANSMEM proteins: %d\n",tot_transmem);
  printf("total FT INTRAMEM proteins: %d\n",tot_intramem);
  printf("total proteins with covalent lipid binding: %d\n",tot_lipid_bind);
  printf("total proteins both intra- & trans- membrane: %d\n",tot_itmem);
  printf("total FT SIGNAL signal peptide containing proteins: %d\n",tot_SIGNAL);
  printf("total proteins with both signal sequence and transmembrane: %d\n",tot_SIG_TRANSMEM);
  printf("total proteins with DNA_BIND: %d\n",tot_DNA_BIND);

  for(i=0;i<REGEX_COUNT;i++) 
     printf("%d: total proteins with CC SUBCELLULAR LOCATION \"%s\": %d\n",
         i, NAMES_ARRAY[i], tot_SCL_ARRAY[i]);

  for(i=0;i<GO_COUNT;i++) 
     printf("%d: total proteins with Gene Ontology \"%s\": %d\n",
         i, GO_NAMES_ARRAY[i], tot_GO_ARRAY[i]);

  printf("total proteins with no CC SUBCELLULAR LOCATION annotation: %d\n",tot_SCL_NULL);
  printf("total total membrane proteins: %d\n",tot_membrane);
  printf("total cytoplasmic proteins: %d\n",tot_cytoplasmic);
  printf("total extracellular proteins: %d\n",tot_extracellular);
  printf("total nuclear proteins: %d\n",tot_nuclear);
  printf("REMAINDER total proteins: %d\n",tot_REMAINDER);
  printf("----------------------------------------\n"); 
  printf("total brain proteins: %d\n",tot_brain);
  printf("brain nuclear proteins: %d\n",brain_nuclear);
  printf("brain cytoplasmic proteins: %d\n",brain_cytoplasmic);
  printf("brain membrane proteins: %d\n",brain_membrane);
  printf("brain extracellular proteins: %d\n",brain_extracellular);
  printf("----------------------------------------\n"); 
  printf("total muscle proteins: %d\n",tot_muscle);
  printf("muscle nuclear proteins: %d\n",muscle_nuclear);
  printf("muscle cytoplasmic proteins: %d\n",muscle_cytoplasmic);
  printf("muscle membrane proteins: %d\n",muscle_membrane);
  printf("muscle extracellular proteins: %d\n",muscle_extracellular);
  printf("----------------------------------------\n"); 
  printf("The protein with the most lines has %d lines\n",max_prot_lines);/**/
  printf("BLOCKSIZE IS %d\n",BLOCKSIZE);
  printf("it took");
  print_interval(&t_begin,&t_end);
  printf(" to run.\n");
  printf("----------------------------------------\n"); 
    /*****************************************************************
     *
     *****************************************************************/
  cellgram("human",argv[file_arg],(double)hum_nuclear,(double)hum_cytoplasmic,
      (double)hum_membrane,(double)hum_extracellular); 
  cellgram("total",argv[file_arg],(double)tot_nuclear,(double)tot_cytoplasmic,
      (double)tot_membrane,(double)tot_extracellular);
  cellgram("brain",argv[file_arg],(double)brain_nuclear,(double)brain_cytoplasmic,
      (double)brain_membrane,(double)brain_extracellular);
  cellgram("muscle",argv[file_arg],(double)muscle_nuclear,(double)muscle_cytoplasmic,
      (double)muscle_membrane,(double)muscle_extracellular);
  close(fd);
  return GOOD_EXIT;
}// int main (int argc, char *argv[]) -----//


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
  long long i = begin;

  if (begin <= end) {
     while (i <= end ) { 
       write(fd,&s[i],1);
       i++;
     }
  }
  write(fd,"\n",1);
  //flush(fd);
}


