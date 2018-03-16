#include "graph.h"
#include "cudahelper.h"

int main()
{
    clock_t beg,end;
    beg=clock();
    int source = 0;
    int dest = 0;
    cudaError_t err = cudaSuccess;
    
    size_t size = numberOfVertex  * sizeof(float);

    float *sumOfOutDegree=new float[numberOfVertex ];

    //Allocate the device memory
    float *d_Sum_Of_Degree = NULL;
    cudaMalloc((void **)&d_Sum_Of_Degree, size);
    if(d_Sum_Of_Degree == NULL)
    {
        cout << "Failed"<<endl;
    }

    float *d_PR = NULL;
    cudaMalloc((void**)&d_PR,size);
    if (d_PR == NULL)
    {
        cout << "Failed" << endl;
    }

    float *d_Graph = NULL;
    
    cudaMalloc((void **)&d_Graph, size * numberOfVertex );
    if (d_Graph == NULL)
    {
        cout <<"Failed" << endl;
    }
	
    //Read Graph file.
    fstream fp("f1.txt",ios::in);
    if(!fp.is_open())
    {
        printf("Failed to open file.\n");
    }

    //host memory allocate

    float *Grap=new float[numberOfVertex*numberOfVertex ];
	
    float *PR=new float[numberOfVertex ];
    float *PR_Temp=new float[numberOfVertex];


    //init
    for (int i = 0; i < numberOfVertex ; ++i)
    {
        PR[i] = InitPageRankValue;
        PR_Temp[i] = InitPageRankValue;
    }

    for (int i = 0; i < numberOfVertex ; ++i)
    {
        for (int j = 0; j < numberOfVertex ; ++j)
        {
            Grap[i*numberOfVertex+j] = 0;
        }
    }

    int edge = 0;
	
    //read from Graph.txt
    while (!fp.eof())
    {
        fp >> source >> dest;
        Grap[(source-1)*numberOfVertex+(dest-1)] = 1;
	Grap[(dest-1)*numberOfVertex+(source-1)] = 1;
	edge++;
    }
    printf("Graph build Done!\n");

    //copy
    err = cudaMemcpy(d_Graph, Grap, numberOfVertex *size, cudaMemcpyHostToDevice);

    if (err != cudaSuccess)
    {
        fprintf(stderr, "Failed to copy vector B from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    //invoke PageRank.
    //CPU Routine
    for (int i = 0; i < numberOfVertex; ++i)
    {
        PR[i] = InitPageRankValue;
    }

    int iter = 0;
    for (int m = 0; m < Max_Iteration_Number; ++m)
    {

        iter ++;

        //CUDA event timing
        
        cudaEvent_t start, stop;
        cudaEventCreate(&start);
        cudaEventCreate(&stop);
        cudaEventRecord(start, 0);

        
        //calculate sum of out degree

	claculateSumOfOutDegree<<<256,256>>>(d_Sum_Of_Degree, d_Graph);
        err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            fprintf(stderr, "Failed to launch vectorAdd kernel (error code %s)!\n", cudaGetErrorString(err));
            exit(EXIT_FAILURE);
        }

        //copy
        err = cudaMemcpy(d_PR, PR, size, cudaMemcpyHostToDevice);

        if (err != cudaSuccess)
        {
            fprintf(stderr, "Failed to copy vector B from host to device (error code %s)!\n", cudaGetErrorString(err));
            exit(EXIT_FAILURE);
        }

	// add the sum
	PRAdd<<<256,256>>>(d_PR, d_Graph, d_Sum_Of_Degree);
        err = cudaGetLastError();

        if (err != cudaSuccess)
        {
            fprintf(stderr, "Failed to launch vectorAdd kernel (error code %s)!\n", cudaGetErrorString(err));
            exit(EXIT_FAILURE);
        }
        cudaEventRecord(stop, 0);
        cudaEventSynchronize(stop);
        float elapsedTime;
        cudaEventElapsedTime(&elapsedTime, start, stop);

        cudaMemcpy(PR_Temp, d_PR, size, cudaMemcpyDeviceToHost);

       if (END(PR_Temp, PR))
        {
            break;
        }
        else
	{
            for (int i = 0; i < numberOfVertex ; ++i)
            {
                PR[i] = PR_Temp[i];
            }
        }	
    }

    cudaMemcpy(d_PR, PR, size, cudaMemcpyHostToDevice);    

    fp.close();
    string ip,tok;
    fp.open("genes");
    map<string,int>gen;
    
    while(getline(fp,ip))
    {
	istringstream iss(ip);
	iss>>source>>tok;
	gen.insert(pair<string,int>(tok,source-1));
    }
	
   fp.close();
   fp.open("da.tsv");
   map<string,int>::iterator it;
	

   vector<int> v;
   long long t=0;
   while(getline(fp,ip))
   {
	t++;
	istringstream iss(ip);
	iss>>tok;
	it=gen.find(tok);
	if(it!=gen.end())
	{
		v.push_back(gen.find(tok)->second);
	}	
   }

   printf("\n size:%lld",t);
   int *va = new int[v.size()];
   copy(v.begin(),v.end(),va);
	
   int *dva=0;
   cudaMalloc((void **)&dva,sizeof(int)*v.size());
   err = cudaMemcpy(dva, va, sizeof(int)*v.size(), cudaMemcpyHostToDevice);

   if (err != cudaSuccess)
   {
   	fprintf(stderr, "Failed to copy vector B from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
    }

    int *dfreq=0;
    float *dans=0;
    cudaMalloc((void **)&dfreq,sizeof(int)*numberOfVertex);
    cudaMalloc((void **)&dans,sizeof(float)*numberOfVertex);
		
    calcfreq<<<256,256>>>(dva,dfreq,v.size());
		
    mul<<<256,256>>>(dans,dfreq,d_PR);
		
   thresh<<<256,256>>>(d_Graph,dans,100);

		

   err = cudaMemcpy(Grap, d_Graph, numberOfVertex *size,cudaMemcpyDeviceToHost);

   if (err != cudaSuccess)
   {
   	fprintf(stderr, "Failed to copy vector B from host to device (error code %s)!\n", cudaGetErrorString(err));
        exit(EXIT_FAILURE);
   }

   freopen("subnets2.txt","w",stdout);
   Graph g(numberOfVertex);
   
   for(int i=0;i<numberOfVertex;i++)
   {
	for(int j=0;j<numberOfVertex;j++)
	{
		if(Grap[i*numberOfVertex+j]==1)
		g.addEdge(i,j);
	}			
   }	
	
  g.printSCCs();
				

  end=clock();
  printf("\n Time taken is : %.9f",(double)(end-beg)/CLOCKS_PER_SEC);
}
