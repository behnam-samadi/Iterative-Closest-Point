#include "ICP.h"
#include <iostream>
#include <fstream>
#include <omp.h>
using namespace gs;
using namespace std;
#define point_dim 3


float total_error = 0;
/*
Create a test box of points.
*/

void print_vector_2D (vector<vector<float>>input){
    for (int i = 0; i< input.size();i++)
    {
        for(int j = 0; j<input[0].size();j++)
        {
            cout<<input[i][j]<<" ";
        }
        cout<<endl;
    }

}
class Frame{
    //later: change to private
public:
    vector<vector<float>> data;
    Frame(string file_adress)
    {
    ifstream fin(file_adress);
    
    bool finished = false;
    string a1 = "start";
    int counter = 0;
    while(!finished)
    {
        //cout<<a1<<" has been read" ;
        getline(fin, a1, ',');      
        if (a1[0]!='e')
        {
            //cout<<stof(a1)<<" has been read"<<endl ;
            counter++;
            if ((counter%100) == 0)
            {
                //cout<<counter<<endl;
            }
            //cout<<stof(a1)<<" has been read"<<endl ;
            data.push_back(vector<float>(point_dim));       
            data[data.size()-1][0] = stof(a1);
            for (int c = 1 ;c<point_dim;c++)
            {
                getline(fin, a1, ',');      
                data[data.size()-1][c] = stof(a1);
            }
        }
        else finished = true;
        
    }
    }
};

void createPoints(std::vector<Point*>& points, Frame * reference , Frame * query)
{
	for (int i = 0 ; i < reference->data.size(); i++)
	{
		points.push_back(new Point(reference->data[i][0], reference->data[i][1], reference->data[i][2]));		
	}
}


void createPoints1(std::vector<Point*>& points)
{
	points.push_back(new Point(3.0f, -2.0f, -7.0f));
	points.push_back(new Point(1.0f, -1.0f, -1.0f));
	points.push_back(new Point(0.0f, 1.0f, -5.0f));
	points.push_back(new Point(2.6f, 1.7f, 9.8f));
}

/*
Apply an afine transofrm to a point cloud
*/
void applyAffineTransform(std::vector<Point*>& points, float* rotationMatrix, float* translation)
{
	Point pRot;
	for (int i = 0; i < points.size(); i++)
	{
		pRot.pos[0] = rotationMatrix[0] * points[i]->pos[0] + rotationMatrix[1] * points[i]->pos[1] + rotationMatrix[2] * points[i]->pos[2] + translation[0];
		pRot.pos[1] = rotationMatrix[3] * points[i]->pos[0] + rotationMatrix[4] * points[i]->pos[1] + rotationMatrix[5] * points[i]->pos[2] + translation[1];
		pRot.pos[2] = rotationMatrix[6] * points[i]->pos[0] + rotationMatrix[7] * points[i]->pos[1] + rotationMatrix[8] * points[i]->pos[2] + translation[2];

		*points[i] = pRot;
	}
}

/*
ICP is used to minimize the distance between 2 point clouds.

This example computes a point cloud of a unit box (Static point cloud)
and a second unit box with a linear transform applied (dynamic point cloud).
ICP is used to transform the dynamic point cloud to best match the static point cloud.
*/
void icpExample(Frame *reference, Frame* query)
{
	
	//create a static box point cloud used as a reference.
	std::vector<Point*> staticPointCloud;
	createPoints(staticPointCloud, reference, query);

	//create a dynamic box point cloud.
	//this point cloud is transformed to match the static point cloud.
	std::vector<Point*> dynamicPointCloud;
	createPoints(dynamicPointCloud, reference , query);

	//apply an artitrary rotation and translation to the dynamic point cloud to misalign the point cloud.
	float rotation[] = { 1.0f, 0.0f, 0.0f,	0.0f, 0.70710678f, -0.70710678f,	0.0f, 0.70710678f, 0.70710678f };
	float translation[] = { -0.75f, 0.5f, -0.5f };
	applyAffineTransform(dynamicPointCloud, rotation, translation);

	/*printf("Static point Cloud: \n");
	for (int i = 0; i < staticPointCloud.size(); i++)
	{
		printf("%0.2f, %0.2f, %0.2f \n", staticPointCloud[i]->pos[0], staticPointCloud[i]->pos[1], staticPointCloud[i]->pos[2]);
	}
	printf("\n");*/

	/*printf("Dynamic point Cloud: \n");
	for (int i = 0; i < dynamicPointCloud.size(); i++)
	{
		printf("%0.2f, %0.2f, %0.2f \n", dynamicPointCloud[i]->pos[0], dynamicPointCloud[i]->pos[1], dynamicPointCloud[i]->pos[2]);
	}
	printf("\n");*/

	//use iterative closest point to transform the dynamic point cloud to best align the static point cloud.
	icp(dynamicPointCloud, staticPointCloud);

	/*printf("Dynamic point Cloud Transformed: \n");
	for (int i = 0; i < dynamicPointCloud.size(); i++)
	{
		printf("%0.2f, %0.2f, %0.2f \n", dynamicPointCloud[i]->pos[0], dynamicPointCloud[i]->pos[1], dynamicPointCloud[i]->pos[2]);
	}
	printf("\n");*/

	float alignmentError = 0.0f;
	for (int i = 0; i < dynamicPointCloud.size(); i++)
	{
		alignmentError += pow(dynamicPointCloud[i]->pos[0] - staticPointCloud[i]->pos[0], 2.0f);
		alignmentError += pow(dynamicPointCloud[i]->pos[1] - staticPointCloud[i]->pos[1], 2.0f);
		alignmentError += pow(dynamicPointCloud[i]->pos[2] - staticPointCloud[i]->pos[2], 2.0f);
	}

	alignmentError /= (float)dynamicPointCloud.size();
	total_error += alignmentError;

	printf("Alignment Error: %0.5f \n", alignmentError);
}

int main()
{
	//Frame reference("reformed_dataset/PC1.txt");
	//Frame reference("reformed_dataset/1_gr.txt");
	Frame reference("reformed_dataset/rad_and_black_0.txt");
	Frame query("reformed_dataset/1_gr.txt");
	cout<<reference.data.size()<<endl;
	cout<<query.data.size();
	double whole_time = -omp_get_wtime();
	double each_test_time = 0;
	int num_tests = 20;
	for (int i = 0 ;i< num_tests; i++)
	{
		cout<<"test number "<<i<<endl;
		each_test_time = -omp_get_wtime();
		icpExample(&reference , &query);
		each_test_time += omp_get_wtime();
		cout<<endl<<"each_test_time: "<<each_test_time<<endl;
	}
	total_error /= num_tests;
	cout<<"total_error"<<total_error<<endl;
	whole_time += omp_get_wtime();
	cout<<endl<<"whole_time: "<<whole_time<<endl;
	system("pause");
	return 0;
}