#Created by Rui 5/17/20

CC = c++
CPPFLAGS =-g -Wall -std=c++14 -O3
INCLUDES = -I/usr/local/include -L/usr/local/lib -lboost_unit_test_framework -static -lpthread
LIB = -I/usr/local/include -L/usr/local/lib -lpthread

TEST_CASES := algorithm_base_test glycan_test io_test lsh_test sim_test lsh_clustering_test  
TEST_CASES_2 := protein_test search_test glycan_builder_test search_engine_test svm_test


search:
	$(CC) $(CPPFLAGS) -o searching \
	apps/search/searching.cpp model/glycan/nglycan_complex.cpp $(LIB)


#  test
glycan_test:
	$(CC) $(CPPFLAGS) -o test/glycan_test \
	 model/glycan/glycan_test.cpp model/glycan/nglycan_complex.cpp $(INCLUDES)

sim_test:
	$(CC) $(CPPFLAGS) -o test/sim_test \
	util/calc/spectrum_sim_test.cpp util/calc/calc.cpp $(INCLUDES)

algorithm_base_test:
	$(CC) $(CPPFLAGS) -o test/algorithm_base_test \
	algorithm/base/base_test.cpp $(INCLUDES)

io_test:
	$(CC) $(CPPFLAGS) -o test/io_test \
	util/io/io_test.cpp  $(INCLUDES)

protein_test:
	$(CC) $(CPPFLAGS) -o test/protein_test \
	engine/protein/protein_test.cpp  $(INCLUDES)

search_test:
	$(CC) $(CPPFLAGS) -o test/search_test \
	algorithm/search/search_test.cpp $(INCLUDES)

glycan_builder_test:
	$(CC) $(CPPFLAGS) -o test/glycan_builder_test \
	engine/glycan/builder_test.cpp model/glycan/nglycan_complex.cpp $(INCLUDES)

search_engine_test:
	$(CC) $(CPPFLAGS) -o test/search_engine_test \
	engine/search/search_engine_test.cpp model/glycan/nglycan_complex.cpp $(INCLUDES)

# test
test: ${TEST_CASES} ${TEST_CASES_2}

# clean up
clean:
	rm -f core test/* *.o clustering searching