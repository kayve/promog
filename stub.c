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
#include <stdlib.h>
#include <cairo.h>

int
main (int argc, char *argv[]) {
  
   char *title = argv[1];
   char *infile = argv[2];
   double nuclear = atof(argv[3]);
   double cytosolic = atof(argv[4]);
   double membrane = atof(argv[5]);
   double extracellular = atof(argv[6]);

   
 
   printf(" nuclear is %6.2f\n", (float)nuclear); 
   printf(" cytosolic is %6.2f\n", (float)cytosolic); 
   printf(" membrane is %6.2f\n", (float)membrane); 
   printf(" extracelluar is %6.2f\n", (float)extracellular); 

   cellgram(title,infile,nuclear,cytosolic,membrane,extracellular);  

} //----- main (int argc, char *argv[])-----// 
