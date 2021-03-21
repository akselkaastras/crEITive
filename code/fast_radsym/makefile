-include ../arch.make

SRC = src
INCLUDE = include
OBJ = obj
BIN = fast_radsym

CXXFLAGS += -I$(INCLUDE) -fopenmp

SRCS =$(wildcard $(SRC)/*.cpp)
INCLUDES =$(wildcard $(INCLUDE)/*.h)
OBJS = $(addprefix $(OBJ)/,$(patsubst %.cpp,%.o,$(notdir $(SRCS))))

.PHONY: all clean mrproper

all: $(BIN)

$(OBJS): | $(OBJ)

$(OBJ):
	mkdir $(OBJ)

$(OBJ)/radsym.o: $(filter-out $(INCLUDE)/radsym.h,$(INCLUDES))

$(OBJ)/dnmap_eigenvalues.o: $(INCLUDE)/conductivity_on_radius.h $(INCLUDE)/dnmap_eigenvalues_pcc.h

$(OBJ)/%.o: $(SRC)/%.cpp $(INCLUDE)/radsym.h
	$(CXX) -o $@ -c $< $(CXXFLAGS)

$(BIN): $(OBJS)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS)

clean:
	@rm -f $(OBJ)/*

mrproper: clean
	@rm -f $(BIN)
