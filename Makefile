OBJS = main.o

CC = gcc
CCFLAGS = -O3
CXX = g++
CXXFLAGS = -O3
LD = gcc
LDFLAGS = -O3
LIBS = -L/home/friesen/lib -lgsl -lgslcblas -lm
INCLUDE = -I/home/friesen/include

TARGET = TOV_solver

%.o : %.c
	$(CC) $(CCFLAGS) $(INCLUDE) -c $<

all: $(OBJS)
	$(LD) $(LDFLAGS) $(LIBS) $(OBJS) -o $(TARGET)

clean:
	rm -rf $(OBJS) $(TARGET)
