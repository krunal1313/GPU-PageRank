
//PageRank value calculate function

__global__ void PRAdd(float *PR, const float* Grap, const float * sumOfOutDegree)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numberOfVertex )
    {
        float sum = 0.0;
        int k = 0;
        for (int j = i; j < numberOfVertex *numberOfVertex  ; j += numberOfVertex )
        {
            if (*(Grap + j) && sumOfOutDegree[k])
            {
                sum += PR[k] / sumOfOutDegree[k];

            }
            k++;
            //printf("%f\n", sum);
        } 
        PR[i] = Alpha  + (1 - Alpha)*sum;
    }
}

//Calculate Sum of out_degree of each vertex.
__global__ void claculateSumOfOutDegree(float * sumOfOutDegree, const float* Grap)
{
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if (i < numberOfVertex )
    {
        sumOfOutDegree[i]  = 0;
        for (int j = 0; j < numberOfVertex ; ++j)
        {
            sumOfOutDegree[i] += *(Grap +i*numberOfVertex  +j);
        }
    }

}

__global__ void calcfreq(int *dva,int *dfreq,long long siz)
{
	
	int i = threadIdx.x + blockIdx.x * blockDim.x;
		int stride = blockDim.x * gridDim.x;
  while (i < siz)
  {
          atomicAdd( &(dfreq[dva[i]]), 1 );
              i += stride;
  }
	
}

__global__ void mul(float *dans,int *dfreq,float *d_PR)
		{
		int i=blockIdx.x * blockDim.x + threadIdx.x;

	if(i<numberOfVertex)
	{
		dans[i]=dfreq[i] * d_PR[i];
	}
}


__global__ void thresh(float *d_graph,float *dans,int thr)
{
	int i= blockIdx.x * blockDim.x + threadIdx.x;
	if(i<numberOfVertex)
	{
		if(dans[i]<thr)
		{
			for(int j=0;j<numberOfVertex;j++)
			d_graph[(i)*numberOfVertex+j]=0;
			
			for(int j=0;j<numberOfVertex;j++)
				d_graph[j*numberOfVertex+(i)]=0;
		}
	
	}
}
//END condition: when the PR value stable

bool END(float a[], float b[])
{
    float sum = 0;
    for (int i = 0; i < numberOfVertex ; ++i)
    {
        sum += abs(a[i] - b[i]);
    }
   // printf("The Deviation Between Two Iteration: %f \n",sum);
    if (sum < END_WEIGHT)
    {
        return true;
    }
    
    return false;
}
