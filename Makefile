OBJS = main.o

CC = gcc
CCFLAGS = -g -O0 -Wall
CXX = g++
CXXFLAGS = -g -O0
LD = gcc
LDFLAGS = -g -O0
LIBS = -L/usr/local/lib -lgsl -lgslcblas -lm
INCLUDE = -I/usr/local/include

TARGET = TOV_solver

%.o : %.c
	$(CC) $(CCFLAGS) $(INCLUDE) -c $<

all: $(OBJS)
	$(LD) $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET)

clean:
	rm -rf $(OBJS) $(TARGET)
