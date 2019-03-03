presek: presek.o
	gcc -Wall -o $@ $<
presek.o: 1.c
	gcc -Wall -c -o $@ $<

.PHONY: clean beauty

clean:
	rm -f *~ presek *.o *.out resenje.txt poligoni.txt
beauty:
	-indent *.c
