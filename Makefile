DEBUG_FLAG = -g -O0
CAIRO_FLAG = `pkg-config --cflags --libs cairo`


promog : promog.o cellgram.o print_interval.o
	gcc -o promog -lrt promog.o cellgram.o print_interval.o ${CAIRO_FLAG} -lm 

promog.o :  
	gcc -c promog.c ${DEBUG_FLAG} ${CAIRO_FLAG} -lm 

cellgram.o :  
	gcc -c cellgram.c ${DEBUG_FLAG} ${CAIRO_FLAG} -lm 

stub : stub.o cellgram.o
	gcc -o stub -lmrt stub.o cellgram.o ${CAIRO_FLAG} -lm 

stub.o :  
	gcc -c stub.c ${DEBUG_FLAG} ${CAIRO_FLAG} -lm 

dirty_demog : dirty_demog.o print_interval.o
	gcc -o dirty_demog -lmrt print_interval.o dirty_demog.o ${DEBUG_FLAG} -lm 

dirty_demog.o :
	gcc -c dirty_demog.c ${DEBUG_FLAG} -lm 

demog_script_friendly: demog_script_friendly.o print_interval.o
	gcc print_interval.o -o demog_script_friendly -lmrt demog_script_friendly.o ${DEBUG_FLAG} -lm 

demog_gets_sf:	demog_gets_sf.o print_interval.o
	gcc print_interval.o -o demog_gets_sf -lmrt demog_gets_sf.o ${DEBUG_FLAG} -lm 

print_interval.o :
	gcc -c print_interval.c ${DEBUG_FLAG} -lm 

demog_script_friendly.o :
	gcc -c demog_script_friendly.c ${DEBUG_FLAG} -lm 

demog_gets_sf.o :
	gcc -c demog_gets_sf.c ${DEBUG_FLAG} -lm 

clean :
	rm ./*.o
