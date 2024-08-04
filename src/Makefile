C=g++
CFLAGS= -c -Wall -Werror -Wextra -g -std=c++17
TEST_CFLAGS = -lgtest -pthread -lgmock -lstdc++

SOURCES=$(wildcard *.cpp)
OBJECTS=$(SOURCES:.cpp=.o)
TESTS_CFILES = $(wildcard tests/*.cpp)
TESTS_OBJ = $(TESTS_CFILES:.cpp=.o)

LIB=s21_matrix_oop.a

all: $(LIB)

$(LIB): $(OBJECTS)
	ar -rc $(LIB) $(OBJECTS)
	ranlib $(LIB)

test : $(TESTS_OBJ) $(LIB)
	$(CC) $^ -o test $(TEST_CFLAGS)
	./test

%.o: %.cpp
	$(C) $(CFLAGS) -c $^ -o $@

style:
	clang-format -style=Google -i *.cpp *.h

clean:
	rm -f $(LIB) $(OBJECTS) $(TESTS_OBJ) test

gcov_report: clean
	# $(CC) $(GCF) ${SOURCES} tests/s21_matrix_tests.cpp $(TEST_CFLAGS) -o gcov_report
	g++ -o s21_matrix_tests.o $(CFLAGS) $(GCF) tests/s21_matrix_test.cpp
	g++ -o s21_matrix_oop.o $(CFLAGS) $(GCF) ${SOURCES}
	g++ -o gcov_report s21_matrix_oop.o s21_matrix_tests.o $(TEST_CFLAGS) $(GCF)
	./gcov_report
	lcov  -t "test" -o gcov_report.info -c -d .
	lcov --remove gcov_report.info '/usr/*' -o filtered_coverage.info
	genhtml -o report filtered_coverage.info