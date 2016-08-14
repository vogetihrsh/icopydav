# MAKE FILE FOR PAIRED END METHOD 
#

CC=$(CXX)
MCXX=mpic++
CFLAGS= -g

CAL=calAvg
SEGMENT=segment

all: $(CAL) $(SEGMENT)

$(SEGMENT): source/cpp/segment.cpp
	    $(MCXX) -o $(SEGMENT) source/cpp/segment.cpp
$(CAL): source/cpp/calculate.cpp 
	$(CC) $(CFLAGS) -o $(CAL) source/cpp/calculate.cpp 

clean: 
	$(RM) $(CAL) $(SEGMENT)
