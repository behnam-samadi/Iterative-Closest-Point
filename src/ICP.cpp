#include "ICP.h"
#include <iostream>
#include <omp.h>
#include <iomanip>
#include <fstream>
#include <vector>
#include <algorithm>
#include<cmath>
#include <pthread.h>
#include <queue>
#include <limits>
#include "numeric"
#include <limits>
#include <time.h>
#include <bits/stdc++.h>
int index_to_ckeck = 1237;
#define point_dim 3
using namespace std;
using namespace gs;
int num_calcs;
//double basis[3] = {26.19372829,12.88486764,7.69639231};
//double basis[3] = {4.41840962,2.42806299,1.49673601};
//double basis[3] = {6.89697076,3.52467804,2.04844848};
//double basis[3] = {1.87220823,0.76542574,0.41296428};
//double basis[3] = {1.09387599, -0.96979575 , 0.94387252};
//double basis[3] = {1,1,-1};
//double basis[3] = {1.51664444, 0.66704231, 0.51779803};
//double basis[3] = {4.07571969, 1.99038848, 1.32392432};
//double basis[3] = {1.12409238, -0.90933214,  0.94841502};
// double basis[3] = {0.11346113, -0.71007125,  1.34066972};
//double basis[3] = {0.29276996, -0.6530365 ,  1.33276279};
//double basis[3] = {-0.12370598, -0.10563619, -1.31032003};
//double basis[3]  = {
 //double basis[3] = {0.11346113, -0.71007125,  1.34066972}
//double basis[3] = {1,1,1};
//double basis[3] = {-0.51475309, -0.29207345,  0.58999927};
//double basis[3]  = {0.46084765, 0.18950231, 0.68975714};
//double basis[3] = {-0.06889851,  0.00933936, -0.71927585};
//double basis[3] = {0.02665797,  0.02011677, -0.11730251};
//double basis[3] = {-0.47216308, -0.06279136,  1.00655682};
//double basis[3] = {-0.10831018, -0.65576258,  2.53081975};
//double basis[3] = {-0.63842603,  0.02794823,  2.83228815};
//double basis[3] = {0.05957336, -0.35304188,  0.96837945};
//double basis[3] = {0,-1,-1};
//double basis[3] = {0.39195903, 0.57994389, 0.72309869};
//double basis[3] = {1,1,1};
//double basis[3] = {0.01105764, 0.50933852, 0.74604819};
double basis[3] = {1,1,1};
double basis_size = sqrt(basis[0] * basis[0] + basis[1] * basis[1] + basis[2] * basis[2]);
bool iter_save_points = true;
bool main_save_points = false;
bool random_init_points = false;
int num_initial_points =200;
int svd_size = 150;


vector<vector <double>> init_queries;
vector<vector <double>> init_references(num_initial_points);
vector <double> temp_ref(3);
bool use_proposed;
enum dist_metric
{
    Modified_Manhattan,
    Euclidean,
    Manhattan
};


void shuffle_array(int arr[], int n)
{
 
    // To obtain a time-based seed
    unsigned seed = 0;
 
    // Shuffling our array
    shuffle(arr, arr + n,
            default_random_engine(seed));
 
    // Printing our array
    //for (int i = 0; i < n; ++i)
        //cout << arr[i] << " ";
    //cout << endl;
}



double calc_distance (double *v1, vector<double>v2, dist_metric type)
{
    //cout<<"distance of ";
    //print_vector_double(v1);
    //print_vector_double(v2);
    if (type == Modified_Manhattan)
    {
        double sum1 = 0;
        double sum2 = 0;
        for(int i = 0; i<point_dim;i++)
            sum1+=v1[i];
        for(int i = 0; i<point_dim;i++)
            sum2+=v2[i];
        return (abs(sum2 - sum1));
    }
    else
    {
        double sum = 0;
        for(int i = 0; i<point_dim;i++)
        {
            if (type==Euclidean)
            {
            	sum+= pow(abs(v1[i] - v2[i]), 2);
            }
            if (type==Manhattan)
            sum+= abs(v1[i] - v2[i]);
        }
        double result = sum;
        if (type == Euclidean)
            result = sqrt(result);
        //cout<<result;
        
        return(result);
        }
}

class Frame{
    //later: change to private
public:
    vector<vector<double>> row_data;
    int num_points; 
    double ** data;
    Frame(string file_adress, int max_points=0)
    {
	    ifstream fin(file_adress, ios::binary);
	    fin.seekg(0, ios::end);
	    size_t num_elements = fin.tellg() / sizeof(double);
	    cout<<file_adress<<file_adress<< num_elements<<endl;
	    if (max_points!=0) num_elements = (max_points*point_dim);
	    num_points = num_elements/point_dim;
	    cout<<num_points<<endl;
	    fin.seekg(0, ios::beg);
	    row_data = vector<vector<double>> (num_points , vector<double> (point_dim, 0));
	    for (int c = 0 ; c<num_points; c++)
	    {
	        if (c%200 == 0) 
	            {cout<<c<<endl;}
	        fin.read(reinterpret_cast<char*>(&row_data[c][0]), point_dim*sizeof(double));
	    }
	    allocate_data();
	    }
	    void allocate_data()
	    {
	        //allocating 
	        double * temp = new double [num_points*(point_dim+1)];
	        data = new double*[num_points];
	        for (int i = 0 ; i < num_points;i++)
	        {
	            data[i] = (temp+i*(point_dim+1));
	        }
    }
    Frame(vector<Point*>* point_data)
    {
    	num_points = point_data->size();
    	row_data = vector<vector<double>> (point_data->size() , vector<double> (point_dim, 0));
    	for (int i = 0 ; i <point_data->size();i++)
    	{
    		for (int j = 0 ; j < point_dim; j++)
    		{
    			row_data[i][j] = (*point_data)[i]->pos[j];
    		}
    	}
    	allocate_data();
    }
    void create_reference_data()
    {
    vector<double> reference_projected(num_points);
    for (int i =0 ; i<num_points;i++)
    {
        reference_projected[i] = 0;
        for (int j = 0; j<point_dim;j++)
        {
        reference_projected[i] += ((row_data[i][j]) * basis[j]);
        }
    }
    vector<int> reference_order(num_points);
    iota(reference_order.begin(),reference_order.end(),0); //Initializing
    sort( reference_order.begin(),reference_order.end(), [&](int i,int j){return reference_projected[i]<reference_projected[j];} );
    double sum;
    for (int i = 0; i<num_points;i++)
    {
        sum = 0;
        for (int j = 0; j<point_dim;j++)
        {
            data[i][j] = row_data[reference_order[i]][j];
            sum += ((data[i][j])*basis[j]);
        }
        data[i][point_dim] = sum;
    }
    }

    void create_query_data(vector<double>* query_projected)
    {
        for (int i =0 ; i<num_points;i++)
    {
        (*query_projected)[i] = 0;
        for (int j = 0; j<point_dim;j++)
        {
        (*query_projected)[i] += row_data[i][j];
        }
    }
    }
};


int binary_search (double ** reference, double query, int begin, int end)
{
	//cout<<"query: "<<query<<" begin: "<<begin<<" end: "<<end<<endl;
    int length = end - begin+1;
    int end_orig = end;
    int middle_index = (begin + end) / 2;
    double middle = reference[(int)((begin + end) / 2)][point_dim];
    while (end >= begin)
    {
        middle_index = (begin + end) / 2;
        middle = reference[(int)((begin + end) / 2)][point_dim];

        if (query == middle) 
        {
        	//*(nearest_projected) = reference[middle_index][point_dim];
            return (middle_index);
        }
        else if (query > middle) 
        {
            begin = middle_index+1;
        }
        else if(query < middle) 
            {

                end = middle_index-1;
            }
        }
        double diff1 = abs(query - middle);
        double diff2;
        double diff3;
        if (middle_index < end_orig)
        {
            diff2 = abs(query - reference[(middle_index+1)][point_dim]);
        }
        else {
            diff2 =numeric_limits<double>::max() ;
        }
        if (middle_index > 0)
        {
            diff3 = abs(query - reference[middle_index-1][point_dim]);
        }
        else
        {
            diff3 = numeric_limits<double>::max();
        }
        if ((diff1 <= diff2) && (diff1 <= diff3))  {
        //*(nearest_projected) = reference[middle_index][point_dim];
        return(middle_index);
        }
        else if ((diff2 <= diff1) && (diff2 <= diff3))
        {
        	//*(nearest_projected) = reference[middle_index+1][point_dim];
            return(middle_index+1);
        }
        else if((diff3 <= diff2) && (diff3 <= diff1)) 
        {
        //*(nearest_projected) = reference[middle_index-1][point_dim];
        return(middle_index-1);
        }
}

int exact_knn_projected(const Frame* reference,vector<double>query, double query_projected, int nearest_index, int K, int row, int num_ref_points, double * NN_points, int result_index, vector<int>&search_range, int bucket_size = 0)
{
	if (iter_save_points) init_queries.push_back(query);
	search_range[0] = nearest_index;
	search_range[1] = nearest_index;
	//cout<<endl<<"init_references: "<<init_references.size();
	//cout<<endl<<"a new call for "<<nearest_index<<endl;
    double NN_dist = calc_distance(reference->data[nearest_index], query, Euclidean);
    
    /*if (dist_to_prev_NN != 0)
    {
    	if (dist_to_prev_NN < NN_dist)
    		NN_dist = dist_to_prev_NN;
    }*/
    int NN_index = nearest_index;    
    double dist;
    bool right_progress = true;
    bool left_progress = true;
    bool original_threshold_visited = false;
    bool bidirectional_cont = true;
    int right_candidate = nearest_index;
    int left_candidate = nearest_index;
    int next;
    bool search_cont = true;
    if (left_candidate == -1)
    {
        left_candidate++;
        bidirectional_cont = false;
        left_progress = false;
        search_range[0] = left_candidate;
    }
    if (right_candidate == num_ref_points)
    {
        right_candidate--;
        bidirectional_cont = false;
        right_progress = false;
        search_range[1] = right_candidate;
    }
    if (abs(reference->data[right_candidate][point_dim] - query_projected) <abs(reference->data[left_candidate][point_dim] - query_projected))
            next = right_candidate;
    else
        next = left_candidate;
    while(search_cont)
    {
    	//cout<<endl<<"in processing "<<nearest_index<<" "<< left_progress<<" , "<<right_progress<<" ----------: "<<next;
        dist = calc_distance(reference->data[next], query, Euclidean);
        num_calcs++;
        if ((bucket_size) && (num_calcs>bucket_size))
        {
        	search_cont = false;
        }
        if (dist < NN_dist)
        {
        	NN_index = next;
        	NN_dist = dist;
        }

        if  ((abs( reference->data[next][point_dim] - query_projected ) > (basis_size*NN_dist)))
        	{search_cont = false;}

        if (left_progress && right_progress)
        {
                        if (left_candidate == -1)
            {
                left_candidate++;
                bidirectional_cont = false;
                left_progress = false;
                search_range[0] = left_candidate;
            }
            if (right_candidate == num_ref_points)
            {
                right_candidate--;
                bidirectional_cont = false;
                right_progress = false;
                search_range[1] = right_candidate;
            }
             if (abs(reference->data[right_candidate][point_dim] - query_projected) <abs(reference->data[left_candidate][point_dim] - query_projected))
            {next = right_candidate;
                right_candidate++;
                search_range[1] =right_candidate;
            }
            else
            {
            next = left_candidate;
            left_candidate--;
            search_range[0] = left_candidate;
        }
        }
        else
        {
        
        if (!left_progress)
        {
            right_candidate++;
            next = right_candidate;
            search_range[1] = right_candidate;
            if (right_candidate == num_ref_points)
            {
                right_candidate--;
                search_range[1] = right_candidate;
                bidirectional_cont = false;
                right_progress = false;
            }

        }
        if (!right_progress)
        {
            left_candidate--;
            search_range[0] = left_candidate;
            next = left_candidate;
            if (left_candidate == -1)
            {
                left_candidate++;
                bidirectional_cont = false;
                left_progress = false;
                search_range[0] = left_candidate;
            }
        }
        if ((!left_progress) && (!right_progress))        	
        {
        	search_cont = 0;
        	break;
        }

    }
    }
    for (int i = 0 ; i < point_dim; i++)
    {
    	NN_points[result_index *point_dim+i] = reference->data[NN_index][i];
    	temp_ref[i]= reference->data[NN_index][i];
    }
    if (iter_save_points)
    {init_references[result_index] = temp_ref;}
    //cout<<endl<<num_calcs<<endl;
    return NN_index;
}





/*int exact_knn_projected_(const Frame* reference,vector<double>query, double query_projected, int nearest_index, int K, int row, int num_ref_points)
{
    
    double NN_dist = calc_distance(reference->data[nearest_index], query, Euclidean);
    int NN_index = nearest_index;
    double max_dist = dist;
    
    
    if (right_arrow<num_ref_points)
    {
    while( abs( reference->data[right_arrow][point_dim] - query_projected ) <= (sqrt(3)*max_dist)    )
    {
        dist = calc_distance(reference->data[right_arrow], query, Euclidean);
        if (dist < max_dist)
        {
        	NN_index = right_arrow;
        	NN_dist = dist;
        	max_dist = dist;
        }
        right_arrow++;
        if (right_arrow == num_ref_points)
            break;
    	}
	}
if (left_arrow>0)
{
        while(abs(reference->data[left_arrow][point_dim] - query_projected) <= (sqrt(3)*max_dist))
    {
        dist = calc_distance(reference->data[left_arrow], query, Euclidean);
        if (dist < max_dist)
        {
        	NN_index = left_arrow;
        	NN_dist = dist;
        	max_dist = dist;
        }
        left_arrow--;
        if (left_arrow<0) break;
    }
}
	int c = 0;

}
*/




void gs::icp(std::vector<Point*> &dynamicPointCloud, std::vector<Point*> &staticPointCloud)
{
	

	//fout.close();
	
	int test_sum_calcs = 0;
	float rotationMatrix[9];
	float translation[3];
	int number_of_in_ranges = 0;

	std::vector<Point*> staticPointCloudCopy;
	
	gs::Point dynamicMid(0.0,0.0,0.0);
	gs::Point staticMid(0.0,0.0,0.0);

	// copy the static point cloud
	for (int i = 0; i < staticPointCloud.size(); i++)
	{
		Point* pCopy = new Point(staticPointCloud[i]->pos[0], staticPointCloud[i]->pos[1], staticPointCloud[i]->pos[2]);
		staticPointCloudCopy.push_back(pCopy);
	}
      /*cout<<endl<<"staticPointCloudCopy "<< staticPointCloudCopy.size();
  cout<<endl<<"dynamicPointCloud "<< dynamicPointCloud.size();
  for (int i =0 ; i < 10; i++)
  {
    cout<<endl<<staticPointCloudCopy[i]->pos[0]<<" , "<<staticPointCloudCopy[i]->pos[1]<<" , " <<staticPointCloudCopy[i]->pos[2]<<endl;    
    cout<<endl<<dynamicPointCloud[i]->pos[0]<<" , "<<dynamicPointCloud[i]->pos[1]<<" , " <<dynamicPointCloud[i]->pos[2]<<endl;    

  }

  exit(0);*/

	// create the kd tree
	KdTree* tree = new KdTree(staticPointCloudCopy);
	Frame reference(&staticPointCloudCopy);
	double create_reference_data_time = -omp_get_wtime();
	reference.create_reference_data();
	create_reference_data_time += omp_get_wtime();
	cout<<endl<<create_reference_data_time<<endl;
	

	for (int i = 0 ; i < 5; i++)
	{
		for (int j =0 ; j < point_dim+1; j++)
		{
			//cout<<endl<<"i: "<<i<<" j: "<<j<<" "<<reference.data[i][j]<<endl;
		}
	}
	size_t numDynamicPoints = dynamicPointCloud.size();

	computeCloudMean(staticPointCloud, &staticMid);
	computeCloudMean(dynamicPointCloud, &dynamicMid);

	// initialize the translation vector
	clearTranslation(translation);

	// initialize the rotation matrix
	clearRotation(rotationMatrix);

	const int maxIterations = 50;
	//const int numRandomSamples = dynamicPointCloud.size();
	const int numRandomSamples = svd_size;
	const float eps = 1e-8;
	gs::Point p;
	gs::Point x;
	int j1 = 0;
	int point_cloud_size = dynamicPointCloud.size();

	//vector<int> random_indices1;
	int randSample;
	int * random_indices2 = new int[point_cloud_size];
	for (int i = 0 ; i < point_cloud_size; i++)
	{
		random_indices2[i] = i;
	}
	//shuffle_array(random_indices2, point_cloud_size);
	/*
	for (int i = 0; i < svd_size; i++)
		{
			randSample = std::rand() % dynamicPointCloud.size();
			//random_indices1.push_back(randSample);
			random_indices1.push_back(i);
			//for (int j2 = 0 ; j2<1000;j2++)
			//j1++;
			//cout<<"index: " <<random_indices1[random_indices1.size()-1]<<endl;
		}*/

		if (random_init_points)
		{


		//vector<vector <double>> init_first_pair(num_initial_points);
		//vector<vector <double>> init_second_pair(num_initial_points);
		//vector <double> temp_pair;
		vector<vector<double>> init_first_pair(num_initial_points , vector<double> (point_dim, 0));
		vector<vector<double>> init_second_pair(num_initial_points , vector<double> (point_dim, 0));
		vector<double> init_distances(num_initial_points);
		int randSample1;
		int randSample2;
		//vector<int> init_indices(num_initial_points*2) [2];
		for (int i = 0; i < num_initial_points; i++)
		{
			randSample1 = std::rand() % dynamicPointCloud.size();
			randSample2 = std::rand() % dynamicPointCloud.size();
			for (int j = 0 ; j < point_dim; j++)
			{
				init_first_pair[i][j] = reference.data[randSample1][j];
				init_second_pair[i][j] = reference.data[randSample2][j];
			}
			double sum = 0;
			for (int j = 0 ; j < point_dim; j++)
			{
				sum+= pow((reference.data[randSample1][j]) - (reference.data[randSample2][j]),2);
			}
			init_distances[i] = sqrt(sum);
		}
		ofstream random_points("/home/behnam/phd/Research/ICP/init_points/random_points");
		for(int i = 0; i < num_initial_points; i++)
			{
				cout<<endl<<"init first pair:  "<<i<<" 'th point: ["<<init_first_pair[i][0]<<" , "<<init_first_pair[i][1]<<" , "<<init_first_pair[i][2]<<"]"<<endl;
				cout<<endl<<"init second pair: "<<i<<" 'th point: ["<<init_second_pair[i][0]<<" , "<<init_second_pair[i][1]<<" , "<<init_second_pair[i][2]<<"]"<<endl;
				cout<<endl<<"distance        : "<<init_distances[i]<<endl;
			}


			for(int i = 0; i < num_initial_points; i++)
			{
				for (int j = 0 ; j <point_dim;j++)
				{
					random_points<<init_first_pair[i][j]<<endl;
					//cout<<endl<<"reference "<<i<<" "<<j<<" wrote"<<endl;
				}
				for (int j = 0 ; j <point_dim;j++)
				{
					random_points<<init_second_pair[i][j]<<endl;
				}
			}


			
exit(0);
		}



	float cost = 1.0;
	std::srand(std::time(0));

	gs::Point qd;
	gs::Point qs;

	float U[9];
	float w[9];
	float sigma[3];
	float V[9];

	float** uSvd = new float*[3];
	float** vSvd = new float*[3];
	uSvd[0] = new float[3];
	uSvd[1] = new float[3];
	uSvd[2] = new float[3];

	vSvd[0] = new float[3];
	vSvd[1] = new float[3];
	vSvd[2] = new float[3];
	int num_iterations = 0;
	int sum_calcs = 0;
	double iter_time = 0;
	double kd_tree_search_time;
	double sum_kd_search_time = 0;
	double sum_proposed_time = 0;
	//vector<double> NN_projected  (dynamicPointCloud.size());
	//to do : free this space
	double * NN_projected =  new double[dynamicPointCloud.size()];
	//double NN_projected[39000];
	double * NN_points = new double [dynamicPointCloud.size() * point_dim];
	//double NN_points[39000*3];
	//double * displacement = new double[dynamicPointCloud.size()];
	
	
	double x_prev,y_prev,z_prev;
	double itrative_threshold;
	double prev_NN_projected;
	double p_projected;
	int switch_iteration = -1;
	int num_rights = 0;
	int num_lefts = 0;
	use_proposed = false;
	vector<int> search_range;
	search_range.push_back(0);
	search_range.push_back(1);
	int * past_nearest_indices = new int[point_cloud_size];
	for (int iter = 0; iter < maxIterations && abs(cost) > eps; iter++)
	{
		number_of_in_ranges = 0;
		num_lefts = 0;
		num_rights = 0;
		if (iter > switch_iteration)
			use_proposed = true;
		iter_time = -omp_get_wtime();
		sum_calcs = 0;
		num_iterations++;
		cost = 0.0;
		//clearRotation(rotationMatrix);
		clearMatrix(U);
		clearMatrix(V);
		clearMatrix(w);
		computeCloudMean(dynamicPointCloud, &dynamicMid);
		int nearest_index = 0;
		if (!use_proposed)
		{
			for (int i = 0; i < numRandomSamples; i++)
			{
				int randSample = std::rand() % dynamicPointCloud.size();
				// sample the dynamic point cloud
				p = *dynamicPointCloud[randSample];

				// get the closest point in the static point cloud
				kd_tree_search_time = -omp_get_wtime();
				tree->search(&p, &x);
				kd_tree_search_time += omp_get_wtime();
				sum_kd_search_time += kd_tree_search_time;
				qd = p - dynamicMid;
				qs = x - staticMid;
				outerProduct(&qs, &qd, w);
				addMatrix(w, U, U);
				cost += error(&x, &p, rotationMatrix, translation);
			}
		}
		else
		{
			if (iter == maxIterations-1)
				svd_size = 150000;
			for (int i = 0; i < svd_size; i++)
			{
				if (iter == maxIterations-1)
				cout<<"i:"<<i<<endl;
				//if (iter < 2) cout<<endl<<"processing point: "<<i;
				int randSample = std::rand() % dynamicPointCloud.size();
				// sample the dynamic point cloud
				//p = *dynamicPointCloud[randSample];
				p = *dynamicPointCloud[random_indices2[i]];
				p_projected = (p.pos[0]*basis[0]) + (p.pos[1]*basis[1]) + (p.pos[2]*basis[2]);
				
				vector <double> p_vec(3);
				for (int d = 0 ; d < point_dim;d++)
				{
					p_vec[d] = p.pos[d];
				}
				double binary_search_time = -omp_get_wtime();
				int nearest_index = binary_search (reference.data,p_projected, 0, reference.num_points-1);
				binary_search_time += omp_get_wtime();
				/*if (i ==random_indices1[5])
				{
					cout<<endl<<"************************"<<endl;
					cout<<"in iteration: "<<iter<<" nearest_index is: "<<nearest_index;
					cout<<endl<<"************************"<<endl;
				}*/

				num_calcs = 0;
				if (i<num_initial_points) {iter_save_points = true;}
				else {iter_save_points = false;}
				double exact_search_time = -omp_get_wtime();
				int bucket_size = 0;
	        	exact_knn_projected(&reference,p_vec,p_projected, nearest_index, 1, 0,reference.num_points, NN_points, i, search_range, bucket_size);
	        	//cout<<endl<<num_calcs<<endl;
	        	//cout<<endl<<"iteration: "<<iter<<" in processing point: "<<i<<" "<<search_range[0]<<" to "<<search_range[1]<<" : " <<search_range[1] - search_range[0]<< " "<<num_calcs;
	        	if (iter>0)
	        	{
	        		//cout<<" "<<past_nearest_indices[i]<<" ";
	        		if ((past_nearest_indices[i] >= search_range[0]) && (past_nearest_indices[i] <= search_range[1]))
	        		{
	        			//cout<<"Yes"<<endl;
	        			number_of_in_ranges++;
	        		}
	        		else
	        		{
	        			//cout<<"No"<<endl;
	        		}

	        	}
	        	else
	        	{
	        		//cout<<endl;
	        	}
	        	past_nearest_indices[i] = nearest_index;
	        	exact_search_time += omp_get_wtime();
	        	//cout<<endl<<"point: "<<i<<" ratio: "<<binary_search_time / exact_search_time;
	        	//cout<<endl<<"in iteration: "<<iter<<" for point " <<i<<" : "<<nearest_index<<" , " <<output_index[i]<<" direction: "<< output_index[i] - nearest_index<<endl;
	        	sum_calcs += num_calcs;
	        	kd_tree_search_time = -omp_get_wtime();
				x.pos[0] = NN_points[i * point_dim];
				x.pos[1] = NN_points[i * point_dim + 1];
				x.pos[2] = NN_points[i * point_dim + 2];
			
				qd = p - dynamicMid;
				qs = x - staticMid;
				outerProduct(&qs, &qd, w);
				addMatrix(w, U, U);
				cost += error(&x, &p, rotationMatrix, translation);
				/*if (i<num_initial_points)
				{
					cout<<"point"<<i<<"added "<<init_queries.size()<<" "<<init_references.size()<<endl;;
					init_queries[i] = p_vec;
					temp_ref[0] = x.pos[0];
					temp_ref[1] = x.pos[1];
					temp_ref[2] = x.pos[2];
					init_references[i] = temp_ref;
				}*/

				//***********************************
				//saving the points when processing num_initial_points's point
				//***********************************

				if ((iter == 0) &&main_save_points && (i == num_initial_points-1))
				{
				ofstream fout("/home/behnam/phd/Research/ICP/init_points/points", ios::app);
				for (int i3 = 0 ;i3 < num_initial_points; i3++ )
				{
					cout<<endl<<i3<<":"<<endl;
					cout<<endl<<"inja: reference point: ";
					for (int i4 = 0 ; i4 < point_dim; i4++)
					{
						cout<<init_references[i3][i4]<<" , ";
						fout<<init_references[i3][i4]<<endl;
						cout<<"neveshte shod";

					}
					cout<<endl<<"inja: query point: ";
					for (int i4 = 0 ; i4 < point_dim; i4++)
					{
						cout<<init_queries[i3][i4]<<" , ";
						fout<<init_queries[i3][i4]<<endl;
					}
					cout<<endl;
				}
				fout.close();
				exit(0);
				}


		}
				//shuffle_array(random_indices2, point_cloud_size);
	}
		copyMatToUV(U, uSvd);
		dsvd(uSvd, 3, 3, sigma, vSvd);
		copyUVtoMat(uSvd, U);
		copyUVtoMat(vSvd, V);

		transpose(V);
		matrixMult(U, V, rotationMatrix);

		gs::Point t(0.0, 0.0, 0.0);
		rotate(&dynamicMid, rotationMatrix, &t);
		translation[0] = staticMid.pos[0] - t.pos[0];
		translation[1] = staticMid.pos[1] - t.pos[1];
		translation[2] = staticMid.pos[2] - t.pos[2];
		//update the point cloud
		for (int i = 0; i < point_cloud_size; i++)
		{
			x_prev = dynamicPointCloud[i]->pos[0];
			y_prev = dynamicPointCloud[i]->pos[1];
			z_prev = dynamicPointCloud[i]->pos[2];
			rotate(dynamicPointCloud[i], rotationMatrix, &p);
			translate(&p, translation, dynamicPointCloud[i]);			
			//dist_to_prev_NN[i] = sqrt(pow(NN_points[(i*point_dim)] - dynamicPointCloud[i]->pos[0],2) + pow(NN_points[(i*point_dim+1)]-dynamicPointCloud[i]->pos[1],2) +pow(NN_points[i*point_dim+2] - dynamicPointCloud[i]->pos[2],2));
			//displacement[i] = sqrt(pow(x_prev - dynamicPointCloud[i]->pos[0],2) + pow(y_prev-dynamicPointCloud[i]->pos[1],2) +pow(z_prev - dynamicPointCloud[i]->pos[2],2));
		}
		cout<<endl<<"in iteration "<<iter<<" average number of calcs: "<<((float)sum_calcs / svd_size)<<endl;
		cout<<endl<<"in iteration "<<iter<<" cost is: "<<cost<<endl;
		cout<<endl<<"in iteration "<<iter<<" number_of_in_ranges is: "<<(float)number_of_in_ranges/svd_size<<endl;
		test_sum_calcs += ((float)sum_calcs / svd_size);

		//cout<<endl<<"num_lefts: "<<num_lefts<<" num_rights: "<<num_rights<<endl;
		//cout<<endl<<"in iteration "<<iter<<" time: "<<iter_time+omp_get_wtime()<<endl;
		//cout<<endl<<"in iteration "<<iter<<" cost: "<<cost<<endl;
		//cout<<endl<<endl<<endl<<"first iteration is ok:"<<endl;
		//cout<<endl<<"and the results are:"<<endl;
		cout<<endl<<init_references.size()<<" "<<init_queries.size()<<endl;

		
			
		
		/*if (save_points)
		{
		for (int i3 = 0 ;i3 < num_initial_points; i3++ )
		{
			cout<<endl<<i3<<":"<<endl;
			cout<<endl<<"reference point: ";
			for (int i4 = 0 ; i4 < point_dim; i4++)
			{
				cout<<init_references[i3][i4]<<" , ";
				fout<<init_references[i3][i4]<<endl;

			}
			cout<<endl<<"query point: ";
			for (int i4 = 0 ; i4 < point_dim; i4++)
			{
				cout<<init_queries[i3][i4]<<" , ";
				fout<<init_queries[i3][i4]<<endl;
			}
			cout<<endl;

		}
		
		//exit(0);
		}*/
		//char temp;
		//cin>>temp;
	}

	//cout<<endl<<"sum_kd_search_time"<<sum_kd_search_time<<endl;

	cout<<endl<<num_iterations<<" maxIterations ";
	

	staticPointCloudCopy.clear();
	delete tree;
	
	delete[] uSvd[0]; 
	delete[] uSvd[1];
	delete[] uSvd[2];
	delete[] uSvd;

	delete[] vSvd[0];
	delete[] vSvd[1];
	delete[] vSvd[2];
	delete[] vSvd;
	cout<<endl<<"test_sum_calcs: "<<test_sum_calcs<<endl;
	cout<<endl<<endl<<"----deleting-----"<<endl<<endl;
	delete [] NN_projected;
	delete [] NN_points; 
	//delete [] displacement;
	//delete [] dist_to_prev_NN ;



}


/*

whole kd-tree search:

total_error0.00118645
whole_time: 0.277023
sh: 1: pause: not found



after 15 iterations:

Alignment Error: 0.00000 
total_error5.14011e-07
whole_time: 0.903092
sh: 1: pause: not found



after 17 iterations:
Alignment Error: 0.00002 
total_error0.000408422
whole_time: 0.690345
sh: 1: pause: not found


double basis[3] = {1,1,-1};
in iteration 0 average number of calcs: 15760.5

double basis[3] = {0.11346113, -0.71007125,  1.34066972};
in iteration 0 average number of calcs: 22121.8

double basis[3] = {0.29276996, -0.6530365 ,  1.33276279};
22216.2


double basis[3] = {-0.12370598, -0.10563619, -1.31032003};
10046.3

double basis[3] = {1,1,1};
12795.5

double basis[3] = {-0.51475309, -0.29207345,  0.58999927};
16699.8

double basis[3] = {-0.06889851,  0.00933936, -0.71927585};
14020.4

double basis[3] = {0.02665797,  0.02011677, -0.11730251};
21420.4

four times of expand:
in iteration 99 average number of calcs: 3014.09

in iteration 99 cost is: 206.711

200 15000

100 maxIterations 
test_sum_calcs: 379976


----deleting-----

Alignment Error: 103796.72656 

each_test_time: 550.876
total_error103797

whole_time: 550.895


when shufflig in each iteration: eval size= 50000
513412

519036

5639.71875


150000  not shuffling every time
in iteration 49 average number of calcs: 6029.03
in iteration 49 cost is: 1.57389e+06


150000, shuffling every time
in iteration 49 average number of calcs: 6084.35
in iteration 49 cost is: 1.54776e+06



49988
i:149989
i:149990
i:149991
i:149992
i:149993
i:149994
i:149995
i:149996
i:149997
i:149998
i:149999

in iteration 49 average number of calcs: -4393.93

in iteration 49 cost is: 6.71558e+08

in iteration 49 number_of_in_ranges is: 0.00348

200 10000

50 maxIterations 
test_sum_calcs: 259609


----deleting-----

Alignment Error: 6305.30127 

each_test_time: 7533.38
total_error6305.3

whole_time: 7533.38
sh: 1: pause: not found
*/
