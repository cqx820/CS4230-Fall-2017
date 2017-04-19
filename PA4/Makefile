INCLUDE=-I/usr/local/cuda/8.0/cuda/include 

SOURCE=spmv.cu
EXECUTABLE=spmv

$(EXECUTABLE): $(SOURCE)
	nvcc $(INCLUDE) $< -o $@ 

clean:
