#include "graph.h"

using namespace std;
clock_t b1,e1,b2,e2,b3,e3,b4,e4;

bool END(float a[], float b[])
{
    float sum = 0;
    for (int i = 0; i < numberOfVertex ; ++i)
    {
        sum += abs(a[i] - b[i]);
    }
   // cout << sum <<endl;
    if (sum < END_WEIGHT)
    {
        return true;
    }

    return false;
}

void PageRank(float *Grap, float PR[])
{
    
    float *PR_Temp=new float[numberOfVertex ];

    //begin = clock();
    int iter = 0;  //迭代次数
    for (int m = 0; m < Max_Iteration_Number; ++m)
    {
        iter++;
		float *sumOfOutDegree=new float[numberOfVertex ];
	//	#pragma omp parallel for
		for (int i = 0; i < numberOfVertex ; ++i)
        {
            sumOfOutDegree[i] = 0.0;
        }

        //Calculate the sum of degree of each vertex
  	//	 #pragma omp parallel for
		for (int i = 0; i < numberOfVertex ; ++i)
        {
            float sum = 0;
	//		
			for (int j = 0; j < numberOfVertex ; ++j)
            {
                sum += *(Grap +i*numberOfVertex  +j);
            }
            sumOfOutDegree[i] = sum;
        }

        //Calculate the PR value of every vertex.
     //   #pragma omp parallel for
		for (int i = 0; i < numberOfVertex ; ++i)
        {
            float sum = 0;
            int k = 0;
			for (int j = i; j < numberOfVertex *numberOfVertex  ; j += numberOfVertex )
            {
                if (*(Grap + j) == 1)
                {
                    if(sumOfOutDegree[k] != 0)
                        sum += PR[k] / sumOfOutDegree[k];
                }
                k++;
                //printf("%f\n", sum);
            }
            PR_Temp[i] = Alpha  + (1 - Alpha)*sum;
        }

        if (END(PR_Temp, PR))
        {
            break;
        }
        else{
		//	#pragma omp parallel for
			for (int i = 0; i < numberOfVertex ; ++i)
            {
                PR[i] = PR_Temp[i];
            }
        }

    }
    //end = clock();
  //  printf("Calculate %d iteration of PageRank value cost us:%d ms.\n",iter, end - begin);

}



int main(int argc,char *argv[])
{
    char ch;
    int source = 0;
    int dist = 0;
    //int vertex = 5000;


	  ifstream fp;
	  fp.open("f1.txt");
	ofstream op;
	  //of.open("Pagerank.txt");
	
    float *Grap=new float[numberOfVertex*numberOfVertex];
	//for(int i=0;i<numberOfVertex;i++)
	//	Graph[i]=new float[numberOfVertex];

    float *PR=new float[numberOfVertex];
//	#pragma omp parallel for
    for (int i = 0; i < numberOfVertex; ++i)
    {
        PR[i] = InitPageRankValue;
    }
//	#pragma omp parallel for
    for (int i = 0; i < numberOfVertex; ++i)
    {
        for (int j = 0; j < numberOfVertex; ++j)
        {
            Grap[i*numberOfVertex+j] = 0;
        }
    }

	string ip;
	string tok,tokk;
	
    while (getline(fp,ip)){

		istringstream iss(ip);
        iss >>source>>dist;
       // std::cout << source<< ' '<< dist<< std::endl;
		
        Grap[(source-1)*numberOfVertex+(dist-1)]=1;
    	Grap[(dist-1)*numberOfVertex+(source-1)]=1;
	}
    printf("Graph build complete.\n");

    //invoke PageRank.
    b1=clock();
	PageRank(Grap, PR);
	e1=clock();

    //output to file PageRankValue.txt
	//#pragma omp parallel
    /*for (int i = 0; i < numberOfVertex; ++i)
    {
        of<< i <<" "<< PR[i] << endl;
    }
    //printf("%f\n", Graph[0][0]);
	*/
	fp.close();
	fp.open("genes");
	map<string,int>gen;
	while(getline(fp,ip))
	{
		istringstream iss(ip);
		iss>>source>>tok;
		gen.insert(pair<string,int>(tok,source-1));
	}
	
	fp.close();
	
	for(int i=0;i<atoi(argv[1]);i++)
	{
	fp.open(argv[2+i]);
	map<string,int>::iterator it;
	
	float *freq=new float[numberOfVertex];
	vector<int>v;
//	ofstream tp;
//	tp.open("vf.txt");
	
	while(getline(fp,ip))
	{
		istringstream iss(ip);
		iss>>tok;
		it=gen.find(tok);
		if(it!=gen.end())
		{
			//freq[gen.find(tok)->second]++;
			v.push_back(gen.find(tok)->second);
		//	tp<<gen.find(tok)->second<<endl;
		}	
	}
	//tp.close();
	fp.close();
	cout<<"size:"<<v.size();
	int *va=new int[v.size()];
	copy(v.begin(),v.end(),va);
	
//	#pragma omp parallel for
	for(long i=0;i<numberOfVertex;i++)
	{
		freq[i]=0;
	}
	
	b2=clock();
//	#pragma omp parallel for
	for(long long i=0;i<v.size();i++)
	{
		#pragma omp atomic
		freq[va[i]]++;
	}
	e2=clock();
	
	float *ans=new float[numberOfVertex];
	
	b3=clock();
//	#pragma omp parallel for
	for(long i=0;i<numberOfVertex;i++)
	{
		ans[i]=PR[i]*freq[i];
		//cout<<ans[i]<<endl;
	}
	e3=clock();
	
	b4=clock();
//	#pragma omp parallel for
	for(long i=0;i<numberOfVertex;i++)
	{
		if(ans[i]<100)
		{
					
			for(int j=0;j<numberOfVertex;j++)
			Grap[(i)*numberOfVertex+j]=0;
			
			for(int j=0;j<numberOfVertex;j++)
			Grap[j*numberOfVertex+(i)]=0;
			
			
		}
	
	}

	free(va);
	free(ans);
	free(freq);

}
	e4=clock();
		/*int cnt;
		op.open("result1.txt");
		vector<int>vi;
		for(int i=0;i<numberOfVertex;i++)
		{
			cnt=0;
			vi.clear();
			
			for(int j=0;j<numberOfVertex;j++)
			{
				if(Graph[i*numberOfVertex+j]==1)
				{
					cnt++;
					vi.push_back(j);
				}
				
			}
			if(cnt>=3)
			{
				op<<i<<" : ";
				for(int k=0;k<vi.size();k++)
				op<<vi[k]<<" ";
				
				op<<endl;
			}
		 } 
	*/
	freopen("subnets.txt","w",stdout);
		Graph g(numberOfVertex);
//		#pragma omp parallel for
		for(int j=0;j<numberOfVertex;j++)
			{
				for(int i=0;i<numberOfVertex;i++)
				{
					if(Grap[i*numberOfVertex+j]==1)
						g.addEdge(i,j);
					
				}			
			}	
	g.printSCCs();
	cout<<"Total time :"<<(double)(e1-b1+e2-b2+e3-b3+e4-b4)/CLOCKS_PER_SEC;
    //getchar();
}
