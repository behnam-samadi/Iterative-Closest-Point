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
#define point_dim 3
using namespace std;
using namespace gs;
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
        //exit(0);
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
        return(middle_index);
        }
        else if ((diff2 <= diff1) && (diff2 <= diff3))
        {
            return(middle_index+1);
        }
        else if((diff3 <= diff2) && (diff3 <= diff1)) 
        {
        return(middle_index-1);
        }
}

void exact_knn_projected(vector<vector<int>>* output,const Frame* reference,vector<double>query, double query_projected, int nearest_index, int K, int row, int num_ref_points)
{
    
    int start_knn = nearest_index;
    int end_knn = nearest_index;
    while((end_knn - start_knn + 1) < K)
    {
        if (start_knn ==0)    
        {
            end_knn += (K - (end_knn - start_knn + 1));
            break;
        }
        if (end_knn == num_ref_points-1)
        {
            start_knn -= (K - (end_knn - start_knn + 1));
            break;
        }
        if ((abs((reference->data)[start_knn-1][point_dim]-query_projected)) < (abs((reference->data)[end_knn+1][point_dim]-query_projected)))
        {
            start_knn--;
        }
        else
        {
            end_knn++;
        }
    }

    double max_dist = calc_distance(reference->data[start_knn], query, Euclidean);	
    double dist;
    int calculated_distances_num = 0;
    priority_queue<pair<double, int>> knn;
    for(int c = start_knn; c<= end_knn; c++)
    {
        dist = calc_distance(reference->data[c], query, Euclidean);
        calculated_distances_num ++;
        knn.push(make_pair(dist, c));
        if (dist > max_dist)
        {
            max_dist = dist;
        }
    }
    int right_arrow = end_knn+1;
    int left_arrow = start_knn-1;
    max_dist = knn.top().first;
    
    if (right_arrow<num_ref_points)
        {
    while( abs( reference->data[right_arrow][point_dim] - query_projected ) <= (sqrt(3)*max_dist)    )
    {
        dist = calc_distance(reference->data[right_arrow], query, Euclidean);
        calculated_distances_num++;
        if (dist < max_dist)
        {
            knn.pop();
            knn.push(make_pair(dist, right_arrow));
            max_dist = knn.top().first;
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
        calculated_distances_num++;
        if (dist < max_dist)
        {
            
            knn.pop();
            knn.push(make_pair(dist, left_arrow));
            max_dist = knn.top().first;
        }
        left_arrow--;
        if (left_arrow<0) break;
    }
}
int c = 0;
    while(knn.size())
    {
        (*output)[row][c++] = knn.top().second;
        knn.pop();
    }
}





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
			cout<<endl<<"i: "<<i<<" j: "<<j<<" "<<reference.data[i][j]<<endl;
		}
	}
	size_t numDynamicPoints = dynamicPointCloud.size();

	computeCloudMean(staticPointCloud, &staticMid);
	computeCloudMean(dynamicPointCloud, &dynamicMid);

	// initialize the translation vector
	clearTranslation(translation);

	// initialize the rotation matrix
	clearRotation(rotationMatrix);

	const int maxIterations = 15;
	//const int numRandomSamples = dynamicPointCloud.size();
	const int numRandomSamples = 400;
	const float eps = 1e-8;
	gs::Point p;
	gs::Point x;

	vector<int> random_indices;
	int randSample;

	for (int i = 0; i < numRandomSamples; i++)
		{
			randSample = std::rand() % dynamicPointCloud.size();
			random_indices.push_back(randSample);
			//cout<<"index: " <<random_indices[random_indices.size()-1]<<endl;
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
	double kd_tree_search_time;
	double sum_kd_search_time = 0;
	double sum_proposed_time = 0;
	vector<vector<int>> result_test  (1 , vector<int> (1, 0));
	for (int iter = 0; iter < maxIterations && abs(cost) > eps; iter++)
	{
		num_iterations++;
		cost = 0.0;
		
		//clearRotation(rotationMatrix);
		clearMatrix(U);
		clearMatrix(V);
		clearMatrix(w);
		computeCloudMean(dynamicPointCloud, &dynamicMid);
		int nearest_index = 0;
		for (int i = 0; i < numRandomSamples; i++)
		{
			int randSample = std::rand() % dynamicPointCloud.size();
			// sample the dynamic point cloud
			p = *dynamicPointCloud[randSample];

			// get the closest point in the static point cloud

			
			double p_projected = p.pos[0] + p.pos[1] + p.pos[2];
			vector <double> p_vec(3);
			for (int d = 0 ; d < point_dim;d++)
			{
				p_vec[d] = p.pos[d];
			}
			//cout<<endl<<i<<" "<<"point created"<<endl;
			double proposed_time = -omp_get_wtime();
			int nearest_index = binary_search (reference.data,p_projected, 0, reference.num_points-1);
			//cout<<i<<" "<<"binary search done "<<nearest_index<<endl;
        	exact_knn_projected(&result_test,&reference,p_vec,p_projected, nearest_index, 1, 0,reference.num_points);
        	//cout<<endl<<result_test.size()<<" and: ";
        	//cout<<result_test[0].size()<<endl;
        	//cout<<endl<<"result_test[0][0]: "<<result_test[0][0]<<endl;
        	//cout<<" X: "<<reference.data[result_test[0][0]][0]<<" Y: "<<reference.data[result_test[0][0]][1]<<" Z: "<<reference.data[result_test[0][0]][2]<<endl;
        	proposed_time += omp_get_wtime();
        	kd_tree_search_time = -omp_get_wtime();
			tree->search(&p, &x);
			//cout<<" X: "<<x.pos[0]<<" Y: "<<x.pos[1]<<" Z: "<<x.pos[2]<<endl;
			x.pos[0] = reference.data[result_test[0][0]][0];
			x.pos[1] = reference.data[result_test[0][0]][1];
			x.pos[2] = reference.data[result_test[0][0]][2];
			kd_tree_search_time += omp_get_wtime();
			cout<<endl<<"proposed_time: "<<sum_proposed_time<<" and kd_tree_search_time: " <<sum_kd_search_time<<endl;
			sum_kd_search_time += kd_tree_search_time;
			sum_proposed_time += proposed_time;

			qd = p - dynamicMid;
			qs = x - staticMid;

			outerProduct(&qs, &qd, w);
			addMatrix(w, U, U);

			cost += error(&x, &p, rotationMatrix, translation);
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
		for (int i = 0; i < dynamicPointCloud.size(); i++)
		{
			rotate(dynamicPointCloud[i], rotationMatrix, &p);
			translate(&p, translation, dynamicPointCloud[i]);
		}
	}
	cout<<endl<<"sum_kd_search_time"<<sum_kd_search_time<<endl;

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
}