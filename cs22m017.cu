#include <iostream>
#include <stdio.h>
#include <cuda.h>

#define max_N 100000
#define max_P 30
#define BLOCKSIZE 1024

using namespace std;


__device__ int remind;//device variable 

/* I am solving this problem in three stage
1.
In the first stage I am calculating the prefix sum over the facility array so that if I know the centre number and facility number'
we can get its maximum capacity in O(1) time. This is also used to get the number of request for any facility of any centre.
in my code.I am doing it via parallel prefix sum algorithm.

2.In the second stage I am sorting the array storing the request info so that i know the number of request for any facilty 
of any centre and I know the staring and last index of all the request corresponding to any facility of any centre.


3.in the third stage I am finally iterating over all the facility than serving its corresponding request if possible by launching kernel.
*/

//*******************************************

// Write down the kernels here




/*it is used to set the value of the device variable*/
__global__ void setvalue(int value){
    remind=value;
}




/* this function is used to calculate prefix sum of an array parallely. The maximum size of the array can be 1024 that is why 
we don't need a global barrier .Syncthreads will be sufficient.It takes as argument the array, size of the array over which
we want to calculate the prefix sum and also the startign index of the array
*/

__global__ void prefixSum(int *a,int a_size,int start){ 
    int tid=blockIdx.x*blockDim.x+threadIdx.x;
    int tempvar=0;
    if(tid==0)
        a[tid+start+1]+=remind;//Adding the last element of the previous block over which we calculated prefix summ to maintain consistency 
    __syncthreads();

    /*doing the obvious parallel prefix sum
    */
    for(int off=1;off<a_size;off*=2){
       
        if(tid>=off){
            tempvar=a[tid-off+start+1];
        }
        __syncthreads();
        if(tid>=off){
            a[tid+start+1]+=tempvar;
        }
        __syncthreads();
        
    }

    //again storing the last element into the device variable remind which acts as a remainder which is going to used when we 
    //calculate prefix sum of the next block
    if(tid==0){
        remind=a[start+a_size-1+1];
    } 
}






/*This is the main function and each instance of this kernel running on a thread represents a facility of a centre.It goes over all the
request of that facility and tries to accomodate that request according to maximum capacity of taht request.If request is granted
then increaments number of succesfull request for the corresponding centre
*/
__global__ void kernel1(int request,    int *d_req_id, int *d_req_cen, int *d_req_fac, int *d_req_start, int *d_req_slots,int startind,int lastind,int cen,int fac,int cap,int * d_succ_reqs  ){
    int tid=blockIdx.x*blockDim.x+threadIdx.x;
    if(tid<request){// check that it is valid request 
        int arr[24];//indicating capacity of each faciltiy at each timestamp initially it is original capacity of faciltiy
        for(int i=0;i<24;i++){
            arr[i]=cap;
        }
        
        //now accessing each request of that facility priority of smaller request id will be taken care of implicitly
        for(int j=startind;j<=lastind;j++){
            int start=d_req_start[j]-1;//start timestamp of that request
            int slot=d_req_slots[j];//number of slots for tha trequest
            int k;

            //for each request check availibilty of the request
            for(k=start;k<start+slot;k++){
                if(arr[k]==0){
                    //condition satisfies means request cant be granted
                    break;
                }
            }


            if(k==start+slot){//facility is available 
                atomicAdd(&d_succ_reqs[cen],1); //increamenting number of succesfull request for the corresponding centre                      
                for(int m=start;m<start+slot;m++){
                    arr[m]-=1;//decreamentign the capacity for that timestamp interval because is is occupies by the above successful request
                }
            }
        }
        
    }

}

//***********************************************


int main(int argc,char **argv)
{
	// variable declarations...
    int N, *centre ,*facility ,*capacity ,*fac_ids , *succ_reqs , *tot_reqs ;
    

    FILE *inputfilepointer;
    
    //File Opening for read
    char *inputfilename = argv[1];
    inputfilepointer    = fopen( inputfilename , "r");

    if ( inputfilepointer == NULL )  {
        printf( "input.txt file failed to open." );
        return 0; 
    }

    fscanf( inputfilepointer, "%d", &N ); // N is number of centres
	
    // Allocate memory on cpu
    centre=(int*)malloc(N * sizeof (int));  // Computer  centre numbers
    memset(centre, 0, N*sizeof(int));
    
    /*This one is important Anup*/
    facility=(int*)malloc(N * sizeof (int));  // Number of facilities in each computer centre
    memset(facility, 0, N*sizeof(int));

    fac_ids=(int*)malloc(max_P * N  * sizeof (int));  // Facility room numbers of each computer centre
    memset(fac_ids, 0, max_P * N*sizeof(int));

    capacity=(int*)malloc(max_P * N * sizeof (int));  // stores capacities of each facility for every computer centre 
    memset(capacity, 0, max_P * N*sizeof(int));


 
    int total=0;
    int success=0;  // total successful requests
    int fail = 0;   // total failed requests


    tot_reqs = (int *)malloc(N*sizeof(int));   // total requests for each centre
    memset(tot_reqs, 0, N*sizeof(int));

    succ_reqs = (int *)malloc(N*sizeof(int)); // total successful requests for each centre
    memset(succ_reqs, 0, N*sizeof(int));


    // Input the computer centres data
    int k1=0 , k2 = 0;
    for(int i=0;i<N;i++)
    {
      fscanf( inputfilepointer, "%d", &centre[i] );
      fscanf( inputfilepointer, "%d", &facility[i] );
      
      for(int j=0;j<facility[i];j++)
      {
        fscanf( inputfilepointer, "%d", &fac_ids[k1] );
        k1++;
      }
      for(int j=0;j<facility[i];j++)
      {
        fscanf( inputfilepointer, "%d", &capacity[k2]);
        k2++;     
      }
    }



    // variable declarations
    int *req_id, *req_cen, *req_fac, *req_start, *req_slots;   // Number of slots requested for every request

    // Allocate memory on CPU 
	int R;


  
	fscanf( inputfilepointer, "%d", &R); // Total requests
    req_id = (int *) malloc ( (R) * sizeof (int) );  // Request ids
    req_cen = (int *) malloc ( (R) * sizeof (int) );  // Requested computer centre
    req_fac = (int *) malloc ( (R) * sizeof (int) );  // Requested facility
    req_start = (int *) malloc ( (R) * sizeof (int) );  // Start slot of every request
    req_slots = (int *) malloc ( (R) * sizeof (int) );   // Number of slots requested for every request
    
    // Input the user request data
    for(int j = 0; j < R; j++)
    {
       fscanf( inputfilepointer, "%d", &req_id[j]);
       fscanf( inputfilepointer, "%d", &req_cen[j]);
       fscanf( inputfilepointer, "%d", &req_fac[j]);
       fscanf( inputfilepointer, "%d", &req_start[j]);
       fscanf( inputfilepointer, "%d", &req_slots[j]);
       tot_reqs[req_cen[j]]+=1;  
    }
        //code the kernels here
    //*******************************************************************************************************************************		


   ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
   //This is first stage of our overall solution that is finding prefix sum so that we can get the index of capacity array
   //when we know the centre number and facilty number so that we can get the corresponding maximum capacity in O(1) time
   /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// 
   
   
    setvalue<<<1,1>>>(0);//setting the value of our device variable remind to be initailly zero it represnets remainder of the
    // previous block over which prefix sum is calculated
    cudaDeviceSynchronize();


    /*d_newfacility is a device array of size N+1 .It is used to store the prefix sum of the array which stores the no. of facility
    of each centre. Our prefix sum will have its first element as zero. So that if we want to know the capacity of first facility
    ith centre we can get that by accessing index "d_newfacility[i]"" of capacity  array or if we want to access the second facility 
    of ith centre  we can get that by accessing index "d_newfacility[i]+1" of capacity  array ans so on
    */
    int *d_newfacility;
    cudaMalloc(&d_newfacility, (N+1)*sizeof(int));
    cudaMemset(d_newfacility, 0, (N+1)* sizeof(int));
    cudaMemcpy(d_newfacility+1,facility, N*sizeof(int), cudaMemcpyHostToDevice);

    /*we are deviding N that is no. of centres in blocks of size 1024 and then computing the prefix sum of each block one at a time
    after calculating for one block the last elelment of that block is stored in device  variable "remind" and then we are calculating
    the prefix sum of the next block and we are making sure that the remind variable is added appropriately so that we get the overall 
    prefix sum of the array and consistency is maintained.
    the reason I am calculating prefix sum over blocks of size 1024 because in that way I don't have to use Global barrier
    syncthreads will be sufficient
    */
    int p=((N)/1024);//number of blocks in our  each of size 1024
    int r=((N)%1024);//size of the last block 
    int starting=0;
    for(int i=0;i<p;i++){
         starting=i*1024;
         //calculating prefix sum of each block
         prefixSum<<<1,1024>>>(d_newfacility,1024,starting);//kernel takes device array,size of the block,and starting index of that block
         cudaDeviceSynchronize();
     }
     starting=p*1024;//starting index of the last block
     prefixSum<<<1,r>>>(d_newfacility,r,starting);//calculating prefix sum of last block
    
    
    //copying the prefix sum to host array cumfacility array
    int *cumfacility=(int*)malloc((N+1) * sizeof (int));
    cudaMemcpy(cumfacility,d_newfacility, (N+1)*sizeof(int), cudaMemcpyDeviceToHost);

    
    /*Now this nor array is very important and its similar to the capacity array that is first each element represents a facility of 
    a centre smilar to how facility array stores the capacity of each facility of each centre this nor array stores the number of
    request of each facility of each centre .And one more thing similar to how we were accessing the capacity of any facility of any 
    centre similar to that we can access the number of request corresponding to any facility of any centre .Indexing used will be
     similar  by takign help of the prefix sum we calculated
    */
    int *nor=(int*)malloc(max_P * N * sizeof (int));
    memset(nor, 0, max_P * N * sizeof(int));




/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This is the second stage of our overall solution that is sorting of the request array
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


/*In this step we are sorting our request array in such a way that a request requesting for a centre with smaller centre number
comes first and if centre number is same than a request requesting for a smaller facility comes first and at  last if a request 
is requestign for the same centre number and also the same facility number than the request with smaller request id should come
first ,according to question.

We are increamenting the number of request for each facility of each centre in the nor array also.
*/

    for(int i=0;i<R;i++){
      int id=req_id[i];
      int cen=req_cen[i];
      int fac=req_fac[i];
      int start=req_start[i];
      int slot=req_slots[i];


      nor[cumfacility[cen]+fac]+=1;//increamenting the number of request for that particukar facility

        int j;
      for( j=i-1;j>=0;j--){
        if( (req_cen[j]>cen) || (req_cen[j]==cen && req_fac[j]>fac )|| (req_cen[j]==cen && req_fac[j]==fac && req_id[j]>id ) ){
          req_id[j+1]=req_id[j];
          req_cen[j+1]=req_cen[j];
          req_fac[j+1]=req_fac[j];
          req_start[j+1]=req_start[j];
          req_slots[j+1]=req_slots[j];
        }
        else{
            break;
        }
      }
          req_id[j+1]=id;
          req_cen[j+1]=cen;
          req_fac[j+1]=fac;
          req_start[j+1]=start;
          req_slots[j+1]=slot;

    }

    int * d_succ_reqs;// decice array for number of successfull request for each centre
    cudaMalloc(&d_succ_reqs, N*sizeof(int));
    cudaMemset(d_succ_reqs, 0, N* sizeof(int));



/* 
 Here i am just transfering the request array info on to the device
*/
    int *d_req_id, *d_req_cen, *d_req_fac, *d_req_start, *d_req_slots;   
    cudaMalloc(&d_req_id, R*sizeof(int));
    cudaMalloc(&d_req_cen, R*sizeof(int));
    cudaMalloc(&d_req_fac, R*sizeof(int));
    cudaMalloc(&d_req_start, R*sizeof(int));
    cudaMalloc(&d_req_slots, R*sizeof(int));
    cudaMemcpy(d_req_id,req_id, R*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_req_cen,req_cen, R*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_req_fac,req_fac, R*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_req_start,req_start, R*sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_req_slots,req_slots, R*sizeof(int), cudaMemcpyHostToDevice);



    int totalfacility=cumfacility[N-1]+facility[N-1];//total no of facility in our problem
    int mycen=-1,myfac=-1,reqind=0;




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//This is the third stage of our overall solution that is finally serving each request 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




/*Now we are at the third stage of our solution. In this stage we are iterating over all the facility and considering
each request of that facility if there is any request for that facility.In each iteration of the loop variable represents the facility
variable mycen contains the corresponding centre number and myfac contains the corresponding request number and for 
each facility we have its number of request and we also have the starting index of the request and also the ending index of 
the request because we have sorted our request array in similar manner
*/
    for(int i=0;i<totalfacility;i++){
        int request=nor[i];// number of request for that facility
        if(fac_ids[i]==0){
            mycen+=1;
            myfac=0;
        }
        else{
            myfac+=1;
        }
        if(request==0)// if no request is for that facility than move on to the next facility
            continue;
       
        int mycapacity=capacity[cumfacility[mycen]+myfac];//maximum capacity of taht facility
        int startind=reqind;//in the request array starting request index of that facility
        int lastind=reqind+request-1;//in the request array last request index of tha facility
        /*launching a thread for each facility it takes as argumnet number of request for that facility all the request array info
        its center number so that if request is accepted it increaments the number of successful request in the d_succ_reqs array
        for that centre, it also takes facility number of that centre and also the maximum capacity of that facility 
        */
        kernel1<<<1,1>>>(request,   d_req_id,d_req_cen,d_req_fac,d_req_start,d_req_slots,startind,lastind,   mycen,myfac,mycapacity, d_succ_reqs);
        reqind+=request;
    }
    cudaDeviceSynchronize();

    //copying  back number of successfull request for each centre back to host 
    cudaMemcpy(succ_reqs,d_succ_reqs, N*sizeof(int), cudaMemcpyDeviceToHost);
    int total_succ=0;

    /*counting the total number of request for each centre and also the total number of successfull request for each centre
    */
    for(int i=0;i<N;i++){
        total+=tot_reqs[i];
        total_succ+=succ_reqs[i];
        
    }    
    success=total_succ;
    fail=total-total_succ;// total number of failed request

    //****************************************************************************************************************************





    // Output
    char *outputfilename = argv[2]; 
    FILE *outputfilepointer;
    outputfilepointer = fopen(outputfilename,"w");

    fprintf( outputfilepointer, "%d %d\n", success, fail);
    for(int j = 0; j < N; j++)
    {
        fprintf( outputfilepointer, "%d %d\n", succ_reqs[j], tot_reqs[j]-succ_reqs[j]);
    }
    fclose( inputfilepointer );
    fclose( outputfilepointer );
    cudaDeviceSynchronize();
	return 0;
}