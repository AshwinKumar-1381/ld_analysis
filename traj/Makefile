# Makefile for the ld_analysis dir

CC = g++

CFLAGS = -w -g

msd: msd.cpp analysis.cpp
	$(CC) $(CFLAGS) msd.cpp analysis.cpp -o ./bin/msd.o

corr: corr.cpp analysis.cpp
	$(CC) $(CFLAGS) corr.cpp analysis.cpp -o ./bin/corr.o

pdf_tavg: pdf_tavg.cpp analysis.cpp
	$(CC) $(CFLAGS) pdf_tavg.cpp analysis.cpp -o ./bin/pdf.o

order: laneOrder.cpp analysis.cpp
	$(CC) $(CFLAGS) laneOrder.cpp analysis.cpp -o ./bin/order.o

rdf1D: computeRDF_1D.cpp analysis.cpp
	$(CC) $(CFLAGS) computeRDF_1D.cpp analysis.cpp -o ./bin/rdf1D.o

rdf2D: computeRDF_2Db.cpp analysis.cpp
	$(CC) $(CFLAGS) computeRDF_2Db.cpp analysis.cpp -o ./bin/rdf2D.o

clean:
	rm -f ./bin/*.o
