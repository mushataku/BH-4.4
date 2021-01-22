gcc_options = -std=c++17 -Wall --pedantic-error

a.out : header.h calc.o
	g++ $(gcc_options) -o $@ $^

calc.o : calc.cpp header.h header.h.gch
	g++ -c $<

header.h.gch : header.h
	g++ $(gcc_options) -x c++-header -o $@ $<

run : a.out
	./a.out

clean :
	rm -f header.h.gch
	rm -f calc.o

# animation :
# 	python3 ../analytic_codes/animation.py
# 	# python3 ../analytic_codes/animation_all.py
# 	cp ../figs/animation.mp4 ../output/steady
# 	# cp ../figs/animation_all.mp4 ../output/steady

rm_data :
	rm -f ./output/position.csv

.PHONY : run clean rm_dara
# .PHONY : run clean animation rm_dara