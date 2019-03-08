presek: presek.o
	gcc -Wall -o $@ $<
presek.o: 1.c
	gcc -Wall -c -o $@ $<

napravi_poligone: napravi_poligone.c
	gcc -Wall $< -o  $@ -lm


.PHONY: clean beauty

clean:
	rm -f *~ presek napravi_poligone *.o resenje.txt poligoni.txt 
beauty:
	-indent *.c
