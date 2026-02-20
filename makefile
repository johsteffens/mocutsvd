mocutsvd_test: test.c mocutsvd.c mocutsvd.h
	gcc -o mocutsvd_test test.c mocutsvd.c -fopenmp -march=native -O3 -lm

# use test_vg for valgrind testing
# valgrind does not work with -march=native or -fopenmp
mocutsvd_test_vg: test.c mocutsvd.c mocutsvd.h
	gcc -o mocutsvd_test_vg test.c mocutsvd.c -O3 -lm

# test in debug mode
mocutsvd_test_dbg: test.c mocutsvd.c mocutsvd.h
	gcc -o mocutsvd_test_dbg test.c mocutsvd.c -g -lm

mocutsvd_example: example.c mocutsvd.c mocutsvd.h
	gcc -o mocutsvd_example example.c mocutsvd.c -fopenmp -march=native -O3 -lm

clean:
	rm -f mocutsvd_example
	rm -f mocutsvd_test
	rm -f mocutsvd_test_vg	
	rm -f mocutsvd_test_dbg
