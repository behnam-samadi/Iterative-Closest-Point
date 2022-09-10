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
int index_to_ckeck = 1237;
#define point_dim 3
using namespace std;
using namespace gs;
int num_calcs;
bool use_proposed;
enum dist_metric
{
    Modified_Manhattan,
    Euclidean,
    Manhattan
};

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
            sum+= pow(abs(v1[i] - v2[i]), 2);
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
        reference_projected[i] += row_data[i][j];
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
            sum += data[i][j];
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

void exact_knn_projected(const Frame* reference,vector<double>query, double query_projected, int nearest_index, int K, int row, int num_ref_points, int result_index, double * NN_points, int * output_index, double displacement = 0, double itrative_center = 0, double dist_to_prev_NN = 0)
{
    double NN_dist = calc_distance(reference->data[nearest_index], query, Euclidean);
    if (NN_dist < dist_to_prev_NN) dist_to_prev_NN = NN_dist;
    /*if (dist_to_prev_NN != 0)
    {
    	if (dist_to_prev_NN < NN_dist)
    		NN_dist = dist_to_prev_NN;
    }*/
    int NN_index = nearest_index;    
    double dist;
    bool itrative_threshold_visitd = false;
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
    }
    if (right_candidate == num_ref_points)
    {
        right_candidate--;
        bidirectional_cont = false;
        right_progress = false;
    }
    if (abs(reference->data[right_candidate][point_dim] - query_projected) <abs(reference->data[left_candidate][point_dim] - query_projected))
            next = right_candidate;
    else
        next = left_candidate;
    while(search_cont)
    {
        dist = calc_distance(reference->data[next], query, Euclidean);
        num_calcs++;
        if (dist < NN_dist)
        {
        	NN_index = next;
        	NN_dist = dist;
        	if (NN_dist < dist_to_prev_NN) dist_to_prev_NN = NN_dist;
        }
        if (displacement > 0)
        {
        	if (abs(reference->data[next][point_dim] - itrative_center)> (sqrt(3)*(displacement+dist_to_prev_NN)) )
        	{	
        		itrative_threshold_visitd = true;
        		cout<<"itrative threshold visitd for point "<<result_index<<" "<<reference->data[next][point_dim]<<endl;
        	}
        }

        if  ((abs( reference->data[next][point_dim] - query_projected ) > (sqrt(3)*NN_dist)))
        	{original_threshold_visited = true;}

            search_cont = !(original_threshold_visited || itrative_threshold_visitd);
            if (itrative_threshold_visitd)
            {
            	cout<<"original_threshold_visited: "<<original_threshold_visited<<" itrative_threshold_visitd: "<<itrative_threshold_visitd<<endl;
            }
        if (left_progress && right_progress)
        {
                        if (left_candidate == -1)
            {
                left_candidate++;
                bidirectional_cont = false;
                left_progress = false;
            }
            if (right_candidate == num_ref_points)
            {
                right_candidate--;
                bidirectional_cont = false;
                right_progress = false;
            }
             if (abs(reference->data[right_candidate][point_dim] - query_projected) <abs(reference->data[left_candidate][point_dim] - query_projected))
            {next = right_candidate;
                right_candidate++;
            }
            else
            {
            next = left_candidate;
            left_candidate--;
        }
        }
        else
        {
        if (!left_progress)
        {
            right_candidate++;
            next = right_candidate;
        }
        if (!right_progress)
        {
            left_candidate--;
            next = left_candidate;
        }
    }
    }
    for (int i = 0 ; i < point_dim; i++)
    {
    	NN_points[result_index *point_dim+i] = reference->data[NN_index][i];
    	output_index[result_index] = NN_index;
    }
    //cout<<endl<<num_calcs<<endl;
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
	float rotationMatrix[9];
	float translation[3];

	std::vector<Point*> staticPointCloudCopy;
	
	gs::Point dynamicMid(0.0,0.0,0.0);
	gs::Point staticMid(0.0,0.0,0.0);

	// copy the static point cloud
	for (int i = 0; i < staticPointCloud.size(); i++)
	{
		Point* pCopy = new Point(staticPointCloud[i]->pos[0], staticPointCloud[i]->pos[1], staticPointCloud[i]->pos[2]);
		staticPointCloudCopy.push_back(pCopy);
	}

	// create the kd tree
	KdTree* tree = new KdTree(staticPointCloudCopy);
	Frame reference(&staticPointCloudCopy);
	reference.create_reference_data();
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

	const int maxIterations = 100;
	//const int numRandomSamples = dynamicPointCloud.size();
	const int numRandomSamples = 400;
	const float eps = 1e-8;
	gs::Point p;
	gs::Point x;
	int j1 = 0;
	int point_cloud_size = dynamicPointCloud.size();

	vector<int> random_indices1;
	int randSample;
	for (int i = 0; i < point_cloud_size; i++)
		{
			randSample = std::rand() % dynamicPointCloud.size();
			//random_indices1.push_back(randSample);
			random_indices1.push_back(i);
			//for (int j2 = 0 ; j2<1000;j2++)
			//j1++;
			//cout<<"index: " <<random_indices1[random_indices1.size()-1]<<endl;
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
	//double NN_projected[2048];
	double * NN_points = new double [dynamicPointCloud.size() * point_dim];
	//double NN_points[2048*3];
	//double * displacement = new double[dynamicPointCloud.size()];
	double displacement[2048];
	double dist_to_prev_NN[2048];
	double p_projected_values[2048];
	int output_index[2048];
	double x_prev,y_prev,z_prev;
	double itrative_threshold;
	double prev_NN_projected;
	double p_projected;
	int switch_iteration = -1;
	int sum_needed_points = 0;
	use_proposed = false;
	int num_wows = 0;
	for (int iter = 0; iter < maxIterations && abs(cost) > eps; iter++)
	{
		num_wows = 0;
		sum_needed_points = 0;
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
		for (int i = 0; i < random_indices1.size(); i++)
		{
			p = *dynamicPointCloud[random_indices1[i]];
			p_projected = p.pos[0] + p.pos[1] + p.pos[2];
			p_projected_values[i] = p_projected;
			vector <double> p_vec(3);
			for (int d = 0 ; d < point_dim;d++)
			{
				p_vec[d] = p.pos[d];
			}
			if (iter == 0)
			{
				int nearest_index = binary_search (reference.data,p_projected, 0, reference.num_points-1);
				num_calcs = 0;
	        	exact_knn_projected(&reference,p_vec,p_projected, nearest_index, 1, 0,reference.num_points, i, NN_points,output_index);
	        	sum_calcs += num_calcs;
	        	kd_tree_search_time = -omp_get_wtime();
				x.pos[0] = NN_points[i * point_dim];
				x.pos[1] = NN_points[i * point_dim + 1];
				x.pos[2] = NN_points[i * point_dim + 2];
			}
			else
			{
				//itrative_threshold = displacement[i] + dist_to_prev_NN[i];
				int nearest_index = binary_search (reference.data,p_projected, 0, reference.num_points-1);
	        	num_calcs = 0;

	        	int start = output_index[i];
	        	int end = start;
	        	int needed_points = 0;
	        	while( abs(reference.data[start][point_dim] - reference.data[end][point_dim]) < sqrt(3) * (displacement[i] * 2 + dist_to_prev_NN[i]))
	        		{
	        		end++;
	        		//cout<<endl<<"end chenged to "<<end<<endl;
	        		if (end > reference.num_points-1)
	        		{
	        			break;
	        		}
	        		}
	        	needed_points += (end - start);
	        	end = start;
	        	while( abs(reference.data[start][point_dim] - reference.data[end][point_dim]) < sqrt(3) * (displacement[i] * 2 + dist_to_prev_NN[i]))
	        		{
	        		end--;
	        		if (end<0)
	        			break;
	        		}
	        		needed_points += (start - end);
	        	//cout<<endl<<"needed_points:"<< needed_points<<endl;
	        	sum_needed_points += needed_points;


	        	//exact_knn_projected(&reference,p_vec,p_projected, nearest_index, 1, 0,reference.num_points, i, NN_points, NN_index, displacement[i], p_projected_values[i], dist_to_prev_NN[i]);
	        	exact_knn_projected(&reference,p_vec,p_projected, nearest_index, 1, 0,reference.num_points, i, NN_points, output_index);
	        	sum_calcs += num_calcs;
	        	//cout<<"iteration" <<iter<<" point ["<<i<<"]: "<<needed_points<<" , "<<num_calcs<<endl;
	        	if (num_calcs > needed_points)
	        		{num_wows++;}
	        	kd_tree_search_time = -omp_get_wtime();
				x.pos[0] = NN_points[i * point_dim];
				x.pos[1] = NN_points[i * point_dim + 1];
				x.pos[2] = NN_points[i * point_dim + 2];
			}
			qd = p - dynamicMid;
			qs = x - staticMid;
			outerProduct(&qs, &qd, w);
			addMatrix(w, U, U);
			cost += error(&x, &p, rotationMatrix, translation);
		}
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
			dist_to_prev_NN[i] = sqrt(pow(NN_points[(i*point_dim)] - dynamicPointCloud[i]->pos[0],2) + pow(NN_points[(i*point_dim+1)]-dynamicPointCloud[i]->pos[1],2) +pow(NN_points[i*point_dim+2] - dynamicPointCloud[i]->pos[2],2));
			displacement[i] = sqrt(pow(x_prev - dynamicPointCloud[i]->pos[0],2) + pow(y_prev-dynamicPointCloud[i]->pos[1],2) +pow(z_prev - dynamicPointCloud[i]->pos[2],2));
		}
		cout<<endl<<"in iteration "<<iter<<" average number of calcs: "<<((float)sum_calcs / point_cloud_size)<<endl;
		cout<<endl<<"in iteration "<<iter<<" average number of needed_points: "<<((float)sum_needed_points / point_cloud_size)<<endl;
		cout<<endl<<"in iteration "<<iter<<" number of wows:"<<num_wows<<endl;
		//cout<<endl<<"in iteration "<<iter<<" time: "<<iter_time+omp_get_wtime()<<endl;
		//cout<<endl<<"in iteration "<<iter<<" cost: "<<cost<<endl;
	}
	cout<<"num_wows: "<<num_wows<<endl;

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



*/
