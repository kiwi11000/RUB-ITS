OUTPUT = MatrixTest
EXTENSION = "cc"
#EXTENSION = "cpp"
MATRIX_DIR = /home/$(USER)/Matrix/

LIBRARIES = 
		

FILES = $(shell find . -maxdepth 1 -name "*.$(EXTENSION)")
HEADER = $(shell find . -maxdepth 1 -name "*.h")
INCLUDES = -I$(MATRIX_DIR)

#definitions for standard make options
# GDB_OPTION = -ggdb
GDB_OPTION = -O3
CPPFLAGS = -Wall
OUTPUT_OPTION = $(GDB_OPTION) $(INCLUDES)

ifeq ($(EXTENSION),"cc")
	OBJECTS = $(FILES:%.cc=%.o)
endif 
ifeq ($(EXTENSION),"cpp")
	OBJECTS = $(FILES:%.cpp=%.o)
endif

#build the output
all: TEST_IF_C_FILES_ARE_THERE $(OUTPUT) 

${OUTPUT}: $(OBJECTS) $(HEADER) 
	g++ $(GDB_OPTION) $(LIBRARIES) $(OBJECTS) -o ${OUTPUT}

TEST_IF_C_FILES_ARE_THERE:
	@test "$(FILES)" != ""

#rebuild the objects, if one header file has changed
$(OBJECTS): $(HEADER)

clean:
	rm -f $(OUTPUT) $(OBJECTS)

new: clean $(OUTPUT)
