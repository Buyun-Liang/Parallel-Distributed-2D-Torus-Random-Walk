CC	= mpicc 
CFLAGS  = -g -Wall 

main.ex: main.o aux.o  
	$(CC) -o main.ex $(CFLAGS) -lm main.o aux.o

clean:
	'rm' OUT/*  *.o *.ex *~ #*

.c.o:	
	$(CC) -c $(CFLAGS)  $< 
