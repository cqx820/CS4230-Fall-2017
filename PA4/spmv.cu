#include <stdio.h>
#include <math.h>
//#include <cutil.h>
#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?)(x):(y))

extern int cudaMemcpy();
extern int cudaFree();
extern void _syncthreads();
extern int cudaMemcpyToSymbol();
extern void MV_GPU_wrapper(int*, float*, float*, float*, int*, int, int, int);
extern int cudaMalloc();
extern __global__ void mv_GPU(int, int, int, int*, float*, float*, float*, int*); 
//extern __shared__ float*;

int block_size = 0;
int grid_num = 0;
int threads_per_block = 0;
int max_non_zero_per_row = 0;

int compare(float* a, float* b, int size, double threshold)
{
	int i;
	for(i = 0; i < size; i++)
	{
		if(abs(a[i] - b[i]) > threshold) return 0;
	}
	return 1;
}

void normalMV(int nr, int* ptr, float* data, float* t, float* b, int* indices){
	int i, j;
	for(i = 0; i < nr; i++){
		for(j = ptr[i]; j < ptr[i + 1]; j++){
			t[i] = t[i] + data[j] * b[indices[j]];
		}
	}
}


extern void MV_GPU_wrapper(int* ptr, float* data, float* t, float* b, int* indices, int nr, int nc, int n){

	float* devO1Ptr;
	float* devI1Ptr;
	float* devI2Ptr;
	int* devIdPtr;
	int* devptr_ptr;
	
	cudaMalloc((void**)&devO1Ptr, 4 * nr);
	cudaMalloc((void**)&devI1Ptr, 4 * n);
	cudaMalloc((void**)&devI2Ptr, 4 * nc);
	cudaMemcpy(devO1Ptr, t, 4 * nr, cudaMemcpyHostToDevice);
	cudaMemcpy(devI1Ptr, data, 4 * n, cudaMemcpyHostToDevice);
	cudaMemcpy(devI2Ptr, b, 4 * nc, cudaMemcpyHostToDevice);
	cudaMalloc((void**)&devptr_ptr, 4 * (nr + 1));
	cudaMalloc((void**)&devIdPtr, 4 * n);
	cudaMemcpy(devptr_ptr, ptr, 4 * (nr + 1), cudaMemcpyHostToDevice);
	cudaMemcpy(devIdPtr, indices, 4 * n, cudaMemcpyHostToDevice); 
	
	dim3 dimGrid(grid_num, 1);
	dim3 dimBlock(block_size, 1);
	dim3 threadsPerBlock(threads_per_block, 1);

//	printf("we have %d grids, %d blocks and %d threads per block\n", grid_num, block_size, threads_per_block);

	mv_GPU<<<dimGrid, dimBlock, threads_per_block>>>(nr, max_non_zero_per_row, block_size, devptr_ptr, devI1Ptr, devO1Ptr, devI2Ptr, devIdPtr);
	
	cudaMemcpy(t, devO1Ptr, nr * 4, cudaMemcpyDeviceToHost);

	cudaFree(devO1Ptr);
	cudaFree(devI1Ptr);
	cudaFree(devI2Ptr);
	cudaFree(devIdPtr);
	cudaFree(devptr_ptr);

//data number of non zero
}


extern __global__ void mv_GPU(int nr, int mx, int blockSize, int* ptr, float* data, float* t, float* b, int* indices)
{
	int bx;
	int tx;
	float suif_tmp0;
//	__shared__ float _P1[];

//	__shared__ float* AS = _P1;
//	__shared__ float* BS = AS + (sizeof(float) * blockSize);

//	int blksz = blockSize;	
	int k, j;
	bx = blockIdx.x;
	tx = threadIdx.x;
	int ptr_cur;
	int ptr_next;
	if(tx <= -(blockSize * bx) + (nr - 1)){
//		suif_tmp0 = 0.0;
		suif_tmp0 = ((float* )(float(*)[])t)[tx + blockSize * bx];
		ptr_cur = ((int*)(int(*)[])ptr)[tx + blockSize * bx];
		ptr_next = ((int*)(int(*)[])ptr)[blockSize * bx + tx + 1];	
	}
	
//	for(k = 0; k < grid_num - 1; k++){
	//	if(tx <= -(block_size * k) + (nr - 1)){
	//		((float*)(float(*)[blksz])BS)[blksz * k + tx - blksz * k] = ((float*)(float(*)[])data)[blksz * k + tx];
	//	}
	
	//	__syncthreads();
		


	for(j = 0; j < mx; j++){
		
		if(tx <= -(blockSize * bx) + (nr - 1)){
			if(ptr_next > (ptr_cur + j)){
				//suif_tmp0 = suif_tmp0 + ((float*)(float(*)[blksz]BS)[(ptr_cur + j) - (blksz * k)] * b[indices[ptr_cur + j]];
				suif_tmp0 = suif_tmp0 + data[ptr_cur + j] * b[indices[ptr_cur + j]];
			}

		}
//	__syncthreads();

	}
	
	__syncthreads();
//}
	if(tx <= -(blockSize * bx) + (nr - 1)){
		((float*)(float(*)[])t)[tx + blockSize * bx] = suif_tmp0;	
	}		
}


int main(int argc, char** argv){

	FILE* fp;
	char line[1024];
	int* ptr, *indices;
	float *data, *b, *t_h, *t_d;
	int i, j;
	int n, nc, nr;
	

	if(argc < 2) abort();
	
	if((fp = fopen(argv[1], "r")) == NULL) abort();
	
	fgets(line, 128, fp);
	while(line[0] == '%'){
		fgets(line, 128, fp);
	}

	sscanf(line, "%d %d %d\n", &nr, &nc, &n);	
	ptr = (int*)malloc((nr + 1) * sizeof(int));
	indices = (int*)malloc(n * sizeof(int));
	data = (float*)malloc(n * sizeof(int));
	b = (float*)malloc(nc * sizeof(int));
	t_h = (float*)malloc(nr * sizeof(int));
	t_d = (float*)malloc(nr * sizeof(int));

	int lastr = 0;
	for(i = 0; i < n; i++)
	{
		int r;
		fscanf(fp, "%d %d %f\n", &r, &(indices[i]), &(data[i]));
		indices[i]--;
		if(r != lastr){
			ptr[r - 1] = i;
			lastr = r;
		}
		
	}
	ptr[nr] = n;
	int temp = 0;
	for(i = 0; i < nr; i++){
		temp = ptr[i + 1] - ptr[i];
		max_non_zero_per_row = max(temp, max_non_zero_per_row);
	}
	
	
	for(i = 0; i < nr; i++){
		t_h[i] = 0.0;
		t_d[i] = 0.0;
	}
	for(i = 0; i < nc; i++)
		b[i] = (float) rand() / 1111111111;

	fclose(fp);		

//	block_size = (nr + 31) / 32;
	block_size = sqrt(nr) + (sqrt(nr) / 2);
	grid_num = block_size / 2;	
	
	threads_per_block = 32;
		
	
	cudaEvent_t start_event, end_event;
	float elapsed_time_seq, elapsed_time_gpu;
	cudaEventCreate(&start_event);
	cudaEventCreate(&end_event);
	cudaEventRecord(start_event, 0);
	normalMV(nr, ptr, data, t_h, b, indices);
	cudaEventRecord(end_event, 0);
	cudaEventSynchronize(end_event);
	cudaEventElapsedTime(&elapsed_time_seq, start_event, end_event);

	cudaEventCreate(&start_event);
	cudaEventCreate(&end_event);
	cudaEventRecord(start_event, 0);
	MV_GPU_wrapper(ptr, data, t_d, b, indices, nr, nc, n);
//	cudaThreadSynchronize();
	cudaEventRecord(end_event, 0);
	cudaEventSynchronize(end_event);
	cudaEventElapsedTime(&elapsed_time_gpu, start_event,  end_event);

	int res = compare(t_h, t_d, nr, 0.01);
	 
	if(res == 1)
		printf("VALID!\n Sequential Time: %.2f mesc\n  Parallel Time: %.2f mesc\n  Speedup = %.2f\n", elapsed_time_seq, elapsed_time_gpu, elapsed_time_seq / elapsed_time_gpu);
	else
		printf("INVALID...\n");
	
	return 0;
	
}
