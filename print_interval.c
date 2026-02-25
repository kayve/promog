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
#include <time.h>
#include <sys/time.h>

#define BILLION 1000000000
#define MILLION 1000000
#define THOUSAND 1000

int print_interval(struct timespec* start, struct timespec* end) {
  long long interval,ninterval;
  int uinterval, sinterval,minterval;
 
  interval = ((int)(end->tv_sec)-(int)(start->tv_sec));
  ninterval = ((int)(end->tv_nsec)-(int)(start->tv_nsec));
  ninterval += BILLION*interval;
  if (ninterval < THOUSAND) {
    printf(" %lld nsec",ninterval);  
  } // ----(ninterval < THOUSAND) -----//
  else {
    uinterval = ninterval/THOUSAND;
    if (uinterval < THOUSAND) {
       ninterval -= THOUSAND*uinterval;
       if (ninterval > 99)
         printf(" %d.%lld usec",uinterval,ninterval);
       else
         if (ninterval > 9)
           printf(" %d.%d%lld usec",uinterval,0,ninterval);
         else
           printf(" %d.%d%d%lld usec",uinterval,0,0,ninterval);
    } //-----(uinterval < THOUSAND) -----// 
    else  {
       minterval = uinterval/THOUSAND; 
       if (minterval < THOUSAND) {
         uinterval -= THOUSAND*minterval;
         if (uinterval > 99)
           printf(" %d.%d msec",minterval,uinterval);
         else
           if (uinterval > 9)
             printf(" %d.%d%d msec",minterval,0,uinterval);
           else
             printf(" %d.%d%d%d msec",minterval,0,0,uinterval);
       } //------ (minterval < THOUSAND) -----// 
       else {
         sinterval = minterval/THOUSAND; 
         minterval -= THOUSAND*sinterval;
         if (minterval > 99)
           printf(" %d.%d sec",sinterval,minterval);
         else
           if (minterval > 9)
             printf(" %d.%d%d sec",sinterval,0,minterval);
           else
             printf(" %d.%d%d%d sec",sinterval,0,0,minterval);
       } //------ ELSE (minterval >= THOUSAND) -----// 
    } //-----ELSE(uinterval >= THOUSAND) -----// 
  } // ----(ninterval < THOUSAND) -----//
} //----print_interval(struct timespec* start, struct timespec* end)----// 

