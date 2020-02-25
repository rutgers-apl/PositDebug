CC := clang -Xclang -load -Xclang $(PD_HOME)/llvm-pass/build/PSan/libPSanitizer.so 
GCC := clang 
SRC_DIR := .
OBJ_DIR := .
SRC_FILES := $(wildcard $(SRC_DIR)/*.c)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.o,$(SRC_FILES))
OBJ_FILES1 := $(patsubst $(SRC_DIR)/%.c,$(OBJ_DIR)/%.pd.o,$(SRC_FILES))
LDFLAGS := -I . -I/$(PD_HOME)/include -L$(PD_HOME)/runtime/obj -lpsanitizer -I/$(SOFTPOSIT_HOME)/source/include  $(SOFTPOSIT_HOME)/build/Linux-x86_64-GCC/softposit.a 
CPPFLAGS := -O1  -g -std=c11  -lmpfr -lm 
CXXFLAGS :=  -I . -I/$(SOFTPOSIT_HOME)/source/include  $(SOFTPOSIT_HOME)/build/Linux-x86_64-GCC/softposit.a 

all: $(OBJ_FILES) $(OBJ_FILES1)

%.pd.o: %.c
	$(CC) $(CPPFLAGS) -o $@ $< $(LDFLAGS)
	chmod u+x $@

%.o: %.c
	$(GCC) $(CPPFLAGS) -o $@ $^ $(LDFLAGS)
	chmod u+x $@

clean:
	rm -rf *.o *.ll *.bc *.log *.out *.o *.dSYM
